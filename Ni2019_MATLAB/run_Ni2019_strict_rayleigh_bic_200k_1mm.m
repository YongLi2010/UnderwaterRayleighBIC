%% Strict two-groove Rayleigh-BIC with a 1 mm minimum groove width at 200 kHz
clear; close all; clc;

rootDir=fileparts(mfilename('fullpath'));
outputDir=fullfile(rootDir,'results');
if ~exist(outputDir,'dir'), mkdir(outputDir); end

fTarget=200e3;                 % Hz
cWater=1500;                  % m/s, nominal water sound speed
minimumWidth=1e-3;             % m
kappaRange=[.04 .13];
OmegaMinimum=1-kappaRange(2);
% A constant normalized lower bound evaluated at the smallest possible
% period guarantees w_phys>=1 mm throughout the complete kappa range.
widthMinimumOverA=fTarget*minimumWidth/(cWater*OmegaMinimum);

old=load(fullfile(outputDir,'StrictRayleighBIC_optimized_practical.mat'));
xOld=old.xFinal;
xSeed=xOld;
xSeed(5)=max(xSeed(5),1.02*widthMinimumOverA);

lb=[kappaRange(1) .10 .10 .20 widthMinimumOverA .08];
ub=[kappaRange(2) .80 .80 .68 .40 .28];
fillMax=.96;
if sum(xSeed(4:6))>fillMax
    xSeed(6)=max(lb(6),fillMax-xSeed(4)-xSeed(5)-.005);
end

fprintf('200 kHz manufacturable strict-BIC search\n');
fprintf('  c=%.1f m/s, lambda=%.6f mm\n',cWater,1e3*cWater/fTarget);
fprintf('  guaranteed normalized width bound w/a >= %.9f\n',widthMinimumOverA);

searchCfg=make_cfg(xSeed,61,9);
broad=ni2019_optimize_strict_rayleigh_bic(searchCfg, ...
    'Truncations',[41 7;61 9;81 11], ...
    'KappaRange',kappaRange,'DepthRange',[lb(2) ub(2)], ...
    'WidthRange',[widthMinimumOverA ub(4)],'GapRange',[lb(6) ub(6)], ...
    'FillMax',fillMax,'Starts',10,'UseGlobal',true, ...
    'GlobalSwarmSize',48,'GlobalMaxIterations',45, ...
    'GlobalMaxStallIterations',14,'MaxFunctionEvaluations',1400, ...
    'MaxIterations',180,'RandomSeed',260824,'Display','off');

% Respect the asymmetric bounds used for the manufacturable small groove.
x=broad.x;
x=min(max(x,lb),ub);
rootTruncations=[121 15;185 23;281 35;313 39];
nRoot=size(rootTruncations,1);
xSequence=zeros(nRoot,6); rootSigma=nan(nRoot,1); rootRaw=nan(nRoot,1);
grooveFraction=nan(nRoot,1); rootResult=cell(nRoot,1);
for j=1:nRoot
    NK=rootTruncations(j,:);
    cfg=make_cfg(x,NK(1),NK(2));
    rootResult{j}=ni2019_refine_strict_rayleigh_bic(cfg,x, ...
        'Truncation',NK,'LowerBounds',lb,'UpperBounds',ub,'FillMax',fillMax, ...
        'MaxFunctionEvaluations',1500,'MaxIterations',180,'Display','off');
    x=rootResult{j}.x;
    xSequence(j,:)=x;
    rootSigma(j)=rootResult{j}.sigma_ratio;
    rootRaw(j)=rootResult{j}.strict_operator.strict_residual;
    grooveFraction(j)=rootResult{j}.strict_operator.groove.pressure_proxy_fraction;
    fprintf('  root N=%3d K=%2d: sigma %.3e, raw %.3e, w2/a %.9f\n', ...
        NK(1),NK(2),rootSigma(j),rootRaw(j),x(5));
end

xFinal=xSequence(end,:);
OmegaFinal=1-xFinal(1); thetaFinal=asind(xFinal(1)/OmegaFinal);
aPhysical=cWater*OmegaFinal/fTarget;
lambdaPhysical=cWater/fTarget;
depthsMm=1e3*aPhysical*xFinal(2:3);
widthsMm=1e3*aPhysical*xFinal(4:5);
gapMm=1e3*aPhysical*xFinal(6);

crossTruncations=[81 10;121 15;161 20;185 23;241 30;281 35;313 39;345 43];
nCross=size(crossTruncations,1); crossSigma=nan(nCross,1); crossRaw=nan(nCross,1);
for j=1:nCross
    R=ni2019_strict_rayleigh_operator(make_cfg(xFinal,crossTruncations(j,1), ...
        crossTruncations(j,2)),'TargetOrder',-1);
    crossSigma(j)=R.sigma_ratio; crossRaw(j)=R.strict_residual;
end

constraintMarginMm=min(widthsMm)-1;
parameterStep=max(abs(xSequence(end,:)-xSequence(end-1,:)));
accepted=rootSigma(end)<1e-8 && crossSigma(end)<1e-4 && ...
    parameterStep<1e-3 && constraintMarginMm>=-1e-8;

%% Summary figure
navy=[.035 .10 .30]; orange=[1 .36 .08]; red=[.90 .05 .06]; gray=[.5 .5 .5];
fig=figure('Color','w','Position',[80 70 1350 790]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
plot(rootTruncations(:,2),xSequence(:,2),'-o','Color',navy,'LineWidth',1.5); hold on;
plot(rootTruncations(:,2),xSequence(:,3),'-s','Color',orange,'LineWidth',1.5);
plot(rootTruncations(:,2),xSequence(:,4),'-^','Color',[.1 .55 .38],'LineWidth',1.5);
plot(rootTruncations(:,2),xSequence(:,5),'-v','Color',red,'LineWidth',1.5);
plot(rootTruncations(:,2),xSequence(:,6),'-d','Color',gray,'LineWidth',1.4);
xlabel('groove truncation K'); ylabel('normalized geometry'); grid on; box on;
title('(a) Manufacturable-root sequence');
legend('d_1/a','d_2/a','w_1/a','w_2/a','g/a','Location','eastoutside');

nexttile;
semilogy(crossTruncations(:,2),crossSigma,'-s','Color',navy,'LineWidth',1.6); hold on;
semilogy(rootTruncations(:,2),rootSigma,'o','Color',red,'MarkerFaceColor',red);
yline(1e-8,'--','Color',gray);
xlabel('groove truncation K'); ylabel('\sigma_{min}/\sigma_{max}');
title('(b) Strict compatibility convergence'); grid on; box on;
legend('fixed final geometry','reoptimized roots','10^{-8} guide','Location','best');

nexttile;
bar([depthsMm(:),widthsMm(:)]); hold on; yline(1,'--','1 mm');
set(gca,'XTickLabel',{'groove 1','groove 2'});
ylabel('physical dimension (mm)'); title('(c) Depths and widths at 200 kHz');
legend('depth','width','Location','best'); grid on; box on;

nexttile;
axis off;
text(.03,.92,sprintf('f = %.3f kHz',fTarget/1e3),'FontSize',13,'Color','k');
text(.03,.80,sprintf('a = %.6f mm',1e3*aPhysical),'FontSize',13,'Color','k');
text(.03,.68,sprintf('lambda = %.6f mm',1e3*lambdaPhysical),'FontSize',13,'Color','k');
text(.03,.56,sprintf('theta = %.6f deg',thetaFinal),'FontSize',13,'Color','k');
text(.03,.44,sprintf('w = [%.6f, %.6f] mm',widthsMm),'FontSize',13,'Color','k');
text(.03,.32,sprintf('d = [%.6f, %.6f] mm',depthsMm),'FontSize',13,'Color','k');
text(.03,.20,sprintf('gap = %.6f mm',gapMm),'FontSize',13,'Color','k');
text(.03,.08,sprintf('strict accepted = %d',accepted),'FontSize',13,'Color',red);
title('(d) Physical design');

sgtitle('Strict single-Rayleigh BIC with w >= 1 mm at 200 kHz','Color','k');
style_figure(fig);
figureFile=fullfile(outputDir,'StrictRayleighBIC_200kHz_min1mm.png');
exportgraphics(fig,figureFile,'Resolution',220);

parameterTable=table(xFinal(1),OmegaFinal,thetaFinal,1e3*aPhysical, ...
    1e3*lambdaPhysical,depthsMm(1),depthsMm(2),widthsMm(1),widthsMm(2),gapMm, ...
    rootSigma(end),rootRaw(end),parameterStep,constraintMarginMm,accepted, ...
    'VariableNames',{'kappa','Omega','theta_deg','period_mm','lambda_mm', ...
    'depth1_mm','depth2_mm','width1_mm','width2_mm','gap_mm', ...
    'sigma_ratio','raw_residual','parameter_step','width_margin_mm','accepted'});
csvFile=fullfile(outputDir,'StrictRayleighBIC_200kHz_min1mm.csv');
writetable(parameterTable,csvFile);
matFile=fullfile(outputDir,'StrictRayleighBIC_200kHz_min1mm.mat');
save(matFile,'fTarget','cWater','minimumWidth','widthMinimumOverA','lb','ub', ...
    'broad','xSequence','rootTruncations','rootSigma','rootRaw','grooveFraction', ...
    'xFinal','OmegaFinal','thetaFinal','aPhysical','lambdaPhysical','depthsMm', ...
    'widthsMm','gapMm','crossTruncations','crossSigma','crossRaw', ...
    'parameterStep','constraintMarginMm','accepted','parameterTable');

fprintf('\n200 kHz / 1 mm constrained result\n');
disp(parameterTable);
fprintf('  %s\n  %s\n  %s\n',figureFile,csvFile,matFile);

function cfg=make_cfg(x,N,K)
Omega=1-x(1);
cfg=struct('a',1,'lambda',1/Omega,'theta_i_deg',asind(x(1)/Omega), ...
    'depths',x(2:3),'widths',x(4:5),'gaps',x(6), ...
    'N',N,'K',K,'solve_scattering',false);
end

function style_figure(fig)
axesList=findall(fig,'Type','axes');
for j=1:numel(axesList)
    ax=axesList(j); ax.Color='w'; ax.XColor='k'; ax.YColor='k';
    ax.GridColor=[.72 .72 .72]; ax.MinorGridColor=[.84 .84 .84];
    ax.Title.Color='k'; ax.XLabel.Color='k'; ax.YLabel.Color='k';
end
legends=findall(fig,'Type','legend');
for j=1:numel(legends)
    legends(j).Color='w'; legends(j).TextColor='k'; legends(j).EdgeColor=[.3 .3 .3];
end
end
