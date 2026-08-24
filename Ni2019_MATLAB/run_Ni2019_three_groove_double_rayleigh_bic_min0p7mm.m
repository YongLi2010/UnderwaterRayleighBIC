%% Three-groove Gamma Rayleigh BIC with all lateral features above 0.7 mm
% Water: c=1500 m/s, f=200 kHz, a=lambda=7.5 mm.  The practical numerical
% lower bound is 0.71 mm, giving 0.01 mm clearance above the requested
% strict >0.70 mm manufacturing constraint.  A_-1, A_0, and A_+1 are
% removed exactly in the odd-parity homogeneous eigenproblem.
clear; close all; clc;

outputDir=fullfile(fileparts(mfilename('fullpath')),'results');
if ~exist(outputDir,'dir'), mkdir(outputDir); end

cWater=1500;
frequency=200e3;
aMm=1e3*cWater/frequency;
requiredMinimumMm=0.70;
numericalMinimumMm=0.71;
h=numericalMinimumMm/aMm;

% The lateral geometry was obtained by continuing the unconstrained odd
% root to the 0.71-mm feasible boundary.  Only the two independent depths
% need to be repolished as the transverse truncation K is increased.
lateral=[0.0946671632665097 0.5266646285814554 ...
    0.0946671632665052 0.0946666666666925 0.0946666666666903];
rootNK=[185 23;249 31;313 39;345 43;377 47; ...
    441 55;505 63;633 79;761 95];
rootDepth=[ ...
    0.1933834164414804 0.7119321892652759; ...
    0.1933757446511810 0.7120272602563635; ...
    0.1933713184296437 0.7120782796169131; ...
    0.1933681837665580 0.7120945638015009; ...
    0.1933684881099692 0.7121095628538512; ...
    0.1933665473415647 0.7121304513964285; ...
    0.1933651495756524 0.7121452549311726; ...
    0.1933633044690238 0.7121646427069011; ...
    0.1933621622537939 0.7121766633679519];
depthSeed=[0.1933989397543590 0.7117070599937559];

% Set NI2019_MIN0P7_REFINE=1 to rerun the sequential depth polishing.
redoDepthRefinement=strcmp(getenv('NI2019_MIN0P7_REFINE'),'1');
if redoDepthRefinement
    u=depthSeed;
    for j=1:size(rootNK,1)
        objective=@(v)depth_objective(v,lateral,rootNK(j,1),rootNK(j,2));
        options=optimset('Display','off','MaxFunEvals',700, ...
            'MaxIter',260,'TolX',1e-11,'TolFun',1e-12);
        u=fminsearch(objective,u,options);
        rootDepth(j,:)=u;
    end
end

nRoot=size(rootNK,1);
roots=cell(nRoot,1);
sigmaRatio=zeros(nRoot,1);
fullRaw=zeros(nRoot,1);
grooveFraction=zeros(nRoot,1);
xSequence=zeros(nRoot,8);
for j=1:nRoot
    xSequence(j,:)=[rootDepth(j,1) rootDepth(j,2) rootDepth(j,1) lateral];
    roots{j}=ni2019_three_groove_parity_strict_operator( ...
        make_cfg(xSequence(j,:),rootNK(j,1),rootNK(j,2)),'odd');
    sigmaRatio(j)=roots{j}.sigma_ratio;
    fullRaw(j)=roots{j}.strict_residual;
    grooveFraction(j)=roots{j}.diagnostics.physical_groove_fraction;
end

xFinal=xSequence(end,:);
finalRoot=roots{end};
genericRoot=ni2019_strict_rayleigh_operator( ...
    make_cfg(xFinal,rootNK(end,1),rootNK(end,2)), ...
    'TargetOrder',-1,'EnforceOtherThreshold',true);

% Quadratic extrapolation in 1/K over the five highest roots defines a
% single estimated continuum geometry.  Its fixed-geometry residual is
% then evaluated without retuning at increasingly large K.
fitId=(nRoot-4:nRoot).';
pOuter=polyfit(1./rootNK(fitId,2),rootDepth(fitId,1),2);
pCenter=polyfit(1./rootNK(fitId,2),rootDepth(fitId,2),2);
limitDepth=[polyval(pOuter,0),polyval(pCenter,0)];
xLimit=[limitDepth(1) limitDepth(2) limitDepth(1) lateral];
fixedK=[39 47 55 63 79 95 111].';
fixedN=8*fixedK+1;
nFixed=numel(fixedK);
fixedSigma=zeros(nFixed,1); fixedRaw=zeros(nFixed,1);
fixedGroove=zeros(nFixed,1);
for j=1:nFixed
    R=ni2019_three_groove_parity_strict_operator( ...
        make_cfg(xLimit,fixedN(j),fixedK(j)),'odd');
    fixedSigma(j)=R.sigma_ratio;
    fixedRaw(j)=R.strict_residual;
    fixedGroove(j)=R.diagnostics.physical_groove_fraction;
end

boundaryGap=1-sum(xFinal(4:8));
featuresMm=[xFinal(4:6),xFinal(7:8),boundaryGap]*aMm;
minimumFeatureMm=min(featuresMm);
centerWidthMaxMm=aMm-5*numericalMinimumMm;
centerCutoffMm=aMm/2;
centralBeta1a=finalRoot.transverse.central_beta(2)*aMm/aMm;

% The K=15 finite root is adequate for outgoing-sheet pole continuation.
xPole=[depthSeed(1) depthSeed(2) depthSeed(1) lateral];
poleCfg=make_cfg(xPole,121,15);
kappa=logspace(log10(3e-4),-3,17);
poles=ni2019_continue_complex_pole(poleCfg,kappa, ...
    'TargetOrder',1,'OuterIterations',14, ...
    'InitialScale',-.05-1.41421356237i, ...
    'PredictorExponent',.5,'Display','off');
poleId=2:numel(poles.kappa);
validPole=poleId(imag(poles.Omega(poleId))>0 & ...
    poles.continued_sigma_ratio(poleId)<1e-11 & ...
    real(poles.qbar(poleId))<0 & imag(poles.qbar(poleId))<0);
if numel(validPole)>=3
    fitQ=polyfit(log(poles.kappa(validPole)),log(poles.Q(validPole)),1);
    qExponent=-fitQ(1);
else
    fitQ=[nan nan];
    qExponent=nan;
end

parameterStep=[nan;max(abs(diff(xSequence)),[],2)];
convergenceTable=table(rootNK(:,1),rootNK(:,2),rootDepth(:,1), ...
    rootDepth(:,2),sigmaRatio,fullRaw,grooveFraction,parameterStep, ...
    'VariableNames',{'N','K','d_outer_over_a','d_center_over_a', ...
    'sigma_ratio','full_raw_residual','physical_groove_fraction', ...
    'max_parameter_step'});
fixedLimitTable=table(fixedN,fixedK,fixedSigma,fixedRaw,fixedGroove, ...
    'VariableNames',{'N','K','sigma_ratio','full_raw_residual', ...
    'physical_groove_fraction'});
parameterTable=table(aMm,frequency/1e3,requiredMinimumMm, ...
    numericalMinimumMm,xFinal(1)*aMm,xFinal(2)*aMm,xFinal(3)*aMm, ...
    xLimit(1)*aMm,xLimit(2)*aMm,xLimit(3)*aMm, ...
    xFinal(4)*aMm,xFinal(5)*aMm,xFinal(6)*aMm, ...
    xFinal(7)*aMm,xFinal(8)*aMm,boundaryGap*aMm, ...
    minimumFeatureMm,real(centralBeta1a), ...
    'VariableNames',{'period_mm','frequency_kHz','required_minimum_mm', ...
    'numerical_minimum_mm','d1_K95_mm','d2_K95_mm','d3_K95_mm', ...
    'd1_limit_mm','d2_limit_mm','d3_limit_mm','w1_mm','w2_mm', ...
    'w3_mm','g1_mm','g2_mm','boundary_gap_mm', ...
    'minimum_lateral_feature_mm','central_beta1_a'});
poleTable=table(poles.kappa(:),real(poles.Omega(:)), ...
    imag(poles.Omega(:)),poles.Q(:),poles.continued_sigma_ratio(:), ...
    'VariableNames',{'kappa','real_Omega','imag_Omega','Q','sigma_ratio'});

% ------------------------------- Figure -------------------------------
colors=ni2019_viridis(9);
purple=colors(2,:); blue=colors(4,:); green=colors(6,:);
yellow=colors(9,:); red=[.88 .12 .10]; gray=[.55 .55 .57];
fig=figure('Color','w','Position',[35 35 1510 870]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

nexttile;
draw_geometry(xFinal,green,.68);
axis equal; xlim([0 1]); ylim([-0.80 .08]);
xlabel('x/a'); ylabel('y/a');
title('(a) 0.71-mm-minimum strict cell');
text(.5,-.76,sprintf('minimum feature = %.6f mm',minimumFeatureMm), ...
    'HorizontalAlignment','center');
style_axes(gca);

nexttile;
plot(rootNK(:,2),rootDepth(:,1)*aMm,'-o','Color',blue, ...
    'MarkerFaceColor',blue,'LineWidth',1.5); hold on;
plot(rootNK(:,2),rootDepth(:,2)*aMm,'-s','Color',purple, ...
    'MarkerFaceColor',purple,'LineWidth',1.5);
yline(limitDepth(1)*aMm,'--','Color',blue,'LineWidth',1.0);
yline(limitDepth(2)*aMm,'--','Color',purple,'LineWidth',1.0);
grid on; xlabel('transverse truncation K'); ylabel('depth (mm)');
title('(b) Truncation-by-truncation root');
legend('outer-depth roots','center-depth roots','outer extrapolation', ...
    'center extrapolation','Location','best');
style_axes(gca);

nexttile;
semilogy(rootNK(:,2),sigmaRatio,'-o','Color',green, ...
    'MarkerFaceColor',green,'LineWidth',1.5); hold on;
semilogy(rootNK(:,2),fullRaw,'-s','Color',red, ...
    'MarkerFaceColor',red,'LineWidth',1.5);
semilogy(fixedK,fixedSigma,'--^','Color',blue, ...
    'MarkerFaceColor',blue,'LineWidth',1.2);
semilogy(fixedK,fixedRaw,'--v','Color',purple, ...
    'MarkerFaceColor',purple,'LineWidth',1.2);
yline(1e-10,':','Color',gray);
grid on; xlabel('transverse truncation K'); ylabel('residual');
title('(c) Strict eigenvalue compatibility');
legend('retuned-root \sigma ratio','retuned-root raw', ...
    'fixed-limit \sigma ratio','fixed-limit raw','10^{-10} guide', ...
    'Location','east');
style_axes(gca);

nexttile;
bar(1:6,featuresMm,.68,'FaceColor',blue);
hold on; yline(requiredMinimumMm,'--','Color',red,'LineWidth',1.3);
yline(numericalMinimumMm,':','Color',yellow,'LineWidth',1.4);
xticks(1:6); xticklabels({'w_1','w_2','w_3','g_1','g_2','g_b'});
ylabel('lateral feature (mm)'); ylim([0 max(featuresMm)*1.16]); grid on;
title('(d) Manufacturing constraints');
legend('optimized geometry','required >0.70 mm', ...
    'numerical minimum 0.71 mm','Location','northeast');
style_axes(gca);

nexttile;
qShow=(0:5).';
qBars=bar(qShow,[ ...
    finalRoot.transverse.surface_pressure_fraction_by_q(1:6), ...
    finalRoot.transverse.surface_velocity_fraction_by_q(1:6)],.82);
qBars(1).FaceColor=purple; qBars(2).FaceColor=green;
grid on; xlabel('transverse order q'); ylabel('fraction');
title('(e) Aperture modal composition, K=95');
legend('surface pressure','surface velocity','Location','northeast');
style_axes(gca);

nexttile;
if ~isempty(validPole)
    loglog(poles.kappa(validPole),poles.Q(validPole),'-o', ...
        'Color',purple,'MarkerFaceColor',purple,'LineWidth',1.5); hold on;
    if all(isfinite(fitQ))
        qFit=exp(polyval(fitQ,log(poles.kappa(validPole))));
        loglog(poles.kappa(validPole),qFit,'--','Color',yellow, ...
            'LineWidth',1.3);
    end
end
grid on; xlabel('|\kappa|'); ylabel('Q');
title(sprintf('(f) Outgoing-sheet poles, Q\\sim|\\kappa|^{-%.2f}', ...
    qExponent));
legend('continued pole','power-law fit','Location','southwest');
style_axes(gca);

allLegend=findobj(fig,'Type','Legend');
set(allLegend,'Color','w','TextColor',[.10 .10 .12], ...
    'EdgeColor',[.55 .55 .58]);

figureFile=fullfile(outputDir, ...
    'ThreeGrooveDoubleRayleighBIC_min0p7mm.png');
dataFile=fullfile(outputDir, ...
    'ThreeGrooveDoubleRayleighBIC_min0p7mm.mat');
convergenceFile=fullfile(outputDir, ...
    'ThreeGrooveDoubleRayleighBIC_min0p7mm_convergence.csv');
fixedLimitFile=fullfile(outputDir, ...
    'ThreeGrooveDoubleRayleighBIC_min0p7mm_fixed_limit.csv');
parameterFile=fullfile(outputDir, ...
    'ThreeGrooveDoubleRayleighBIC_min0p7mm_parameters.csv');
poleFile=fullfile(outputDir, ...
    'ThreeGrooveDoubleRayleighBIC_min0p7mm_poles.csv');
exportgraphics(fig,figureFile,'Resolution',220);
writetable(convergenceTable,convergenceFile);
writetable(fixedLimitTable,fixedLimitFile);
writetable(parameterTable,parameterFile);
writetable(poleTable,poleFile);
save(dataFile,'aMm','frequency','requiredMinimumMm', ...
    'numericalMinimumMm','h','lateral','rootNK','rootDepth', ...
    'xSequence','xFinal','roots','finalRoot','genericRoot','featuresMm', ...
    'minimumFeatureMm','centerWidthMaxMm','centerCutoffMm', ...
    'centralBeta1a','limitDepth','xLimit','fixedN','fixedK', ...
    'fixedSigma','fixedRaw','fixedGroove','poles','validPole', ...
    'fitQ','qExponent','convergenceTable','fixedLimitTable', ...
    'parameterTable','poleTable');

fprintf('\nThree-groove >0.70-mm Rayleigh-BIC result\n');
fprintf('  period %.6f mm, frequency %.6f kHz\n',aMm,frequency/1e3);
fprintf('  depths [mm] = [%.9f %.9f %.9f]\n',xFinal(1:3)*aMm);
fprintf('  extrapolated depths [mm] = [%.9f %.9f %.9f]\n', ...
    xLimit(1:3)*aMm);
fprintf('  widths [mm] = [%.9f %.9f %.9f]\n',xFinal(4:6)*aMm);
fprintf('  internal gaps [mm] = [%.9f %.9f]\n',xFinal(7:8)*aMm);
fprintf('  boundary gap [mm] = %.9f; minimum = %.9f\n', ...
    boundaryGap*aMm,minimumFeatureMm);
fprintf('  central beta1*a = %.9f (propagating)\n',real(centralBeta1a));
fprintf('  K=95 parity sigma %.3e, raw %.3e, groove %.6f\n', ...
    finalRoot.sigma_ratio,finalRoot.strict_residual, ...
    finalRoot.diagnostics.physical_groove_fraction);
fprintf('  K=95 generic sigma %.3e, raw %.3e\n', ...
    genericRoot.sigma_ratio,genericRoot.strict_residual);
fprintf('  fixed extrapolated geometry at K=111: sigma %.3e, raw %.3e\n', ...
    fixedSigma(end),fixedRaw(end));
fprintf('  outgoing-sheet valid poles %d/%d, Qmax %.3e, exponent %.4f\n', ...
    numel(validPole),numel(poleId),max(poles.Q(validPole)),qExponent);
fprintf('  figure: %s\n  data:   %s\n',figureFile,dataFile);

assignin('base','ThreeGrooveMin0p7mm',struct( ...
    'parameters',parameterTable,'convergence',convergenceTable, ...
    'fixedLimit',fixedLimitTable,'poles',poles,'xFinal',xFinal, ...
    'xLimit',xLimit,'finalRoot',finalRoot, ...
    'genericRoot',genericRoot,'figure',figureFile,'data',dataFile));

function f=depth_objective(u,lateral,N,K)
if any(~isfinite(u))||any(u<.03)||any(u>3)
    f=10+sum(max(0,.03-u).^2)+sum(max(0,u-3).^2);
    return;
end
x=[u(1) u(2) u(1) lateral];
R=ni2019_three_groove_parity_strict_operator(make_cfg(x,N,K),'odd');
f=log10(max(R.sigma_ratio,1e-16));
end

function cfg=make_cfg(x,N,K)
cfg=struct('a',1,'lambda',1,'theta_i_deg',0, ...
    'depths',x(1:3),'widths',x(4:6),'gaps',x(7:8), ...
    'N',N,'K',K,'solve_scattering',false);
end

function draw_geometry(x,color,alphaValue)
occupied=sum(x(4:8)); x0=.5*(1-occupied);
xl=[x0,x0+x(4)+x(7),x0+x(4)+x(7)+x(5)+x(8)];
for ell=1:3
    patch([xl(ell) xl(ell)+x(3+ell) xl(ell)+x(3+ell) xl(ell)], ...
        [0 0 -x(ell) -x(ell)],color,'FaceAlpha',alphaValue, ...
        'EdgeColor',color,'LineWidth',1.1); hold on;
end
plot([0 1],[0 0],'-','Color',[.15 .15 .17],'LineWidth',1.4);
end

function style_axes(ax)
set(ax,'Color','w','Box','on','LineWidth',.8,'FontName','Arial', ...
    'FontSize',10.5,'XColor',[.10 .10 .12],'YColor',[.10 .10 .12], ...
    'GridColor',[.82 .82 .84],'GridAlpha',.35,'Layer','top');
end
