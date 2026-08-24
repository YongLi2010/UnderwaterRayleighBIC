%% Three-groove Gamma-point double-Rayleigh-BIC optimization
% This run is independent of the saved off-Gamma two-groove result.
% It fixes kappa=0 and Omega=1, constrains A_-1=A_0=A_+1=0, follows one
% mirror-symmetric odd root through increasing modal truncations, and then
% checks parity, exterior localization, and the outgoing-sheet pole trend.
clear; close all; clc;

runGlobalSearch=strcmp(getenv('NI2019_GAMMA3_GLOBAL_SEARCH'),'1');
runPoleTrack=~strcmp(getenv('NI2019_GAMMA3_SKIP_POLES'),'1');
outputDir=fullfile(fileparts(mfilename('fullpath')),'results');
if ~exist(outputDir,'dir'), mkdir(outputDir); end

navy=[.035 .10 .30]; orange=[1.00 .36 .08]; gray=[.55 .55 .55];
red=[.90 .05 .06]; green=[.10 .55 .38]; purple=[.43 .23 .64];

% Physics-guided basin seed: opposite q=0 amplitudes in the two narrow
% outer grooves and a q=1 transverse mode in the wide central groove.
xSeed=[.25 .65 .25 .10 .52 .10 .07 .07];
lb=[.08 .40 .08 .045 .38 .045 .025 .025];
ub=[.35 .90 .35 .18 .65 .18 .15 .15];

if runGlobalSearch
    searchCfg=make_cfg(xSeed,41,3);
    globalResult=ni2019_optimize_three_groove_double_rayleigh_bic( ...
        searchCfg,'Truncations',[31 3;41 5;57 7], ...
        'Parameterization','both','DepthRange',[.08 .90], ...
        'WidthRange',[.045 .65],'GapRange',[.025 .15], ...
        'FillMax',.94,'Starts',8,'UseGlobal',true, ...
        'GlobalSwarmSize',48,'GlobalMaxIterations',40, ...
        'MaxFunctionEvaluations',1600,'MaxIterations',240, ...
        'BoundaryPenaltyWeight',0,'RandomSeed',390831,'Display','off');
    xSeed=globalResult.x;
else
    globalResult=[];
end

% Each finite discretization is polished to its own strict compatibility
% zero. Parameter convergence, rather than one fixed-discretization zero,
% is the continuum diagnostic.
rootTruncations=[41 3;41 5;73 9;121 15;185 23;249 31;313 39];
nRoot=size(rootTruncations,1);
xSequence=zeros(nRoot,8); rootSigma=zeros(nRoot,1);
rootRaw=zeros(nRoot,1); grooveFraction=zeros(nRoot,1);
oddResidual=zeros(nRoot,1); rootResult=cell(nRoot,1);
x=xSeed;
for j=1:nRoot
    N=rootTruncations(j,1); K=rootTruncations(j,2);
    cfg=make_cfg(x,N,K);
    rootResult{j}=ni2019_refine_three_groove_double_rayleigh_bic( ...
        cfg,x,'Truncation',[N K],'LowerBounds',lb,'UpperBounds',ub, ...
        'FillMax',.94,'MaxFunctionEvaluations',1800, ...
        'MaxIterations',220,'Display','off');
    x=rootResult{j}.x;
    parity=ni2019_three_groove_parity_diagnostic( ...
        rootResult{j}.strict_operator);
    xSequence(j,:)=x;
    rootSigma(j)=rootResult{j}.sigma_ratio;
    rootRaw(j)=rootResult{j}.strict_operator.strict_residual;
    grooveFraction(j)=rootResult{j}.strict_operator.groove.pressure_proxy_fraction;
    oddResidual(j)=parity.coefficient_residuals.odd;
    fprintf('root N=%3d K=%2d: sigma %.3e, raw %.3e, odd %.3e\n', ...
        N,K,rootSigma(j),rootRaw(j),oddResidual(j));
end

xFinal=xSequence(end,:);
finalCfg=make_cfg(xFinal,rootTruncations(end,1),rootTruncations(end,2));
finalStrict=ni2019_strict_rayleigh_operator(finalCfg,'TargetOrder',-1, ...
    'EnforceOtherThreshold',true);
finalParity=ni2019_three_groove_parity_diagnostic(finalStrict);

% Cross-evaluate the final fixed geometry. The dip at K=39 is the polished
% finite root; the decreasing envelope and xSequence are interpreted
% together rather than mistaking that isolated dip for convergence alone.
crossTruncations=[41 3;41 5;57 7;73 9;89 11;121 15;153 19; ...
    185 23;217 27;249 31;281 35;313 39;345 43];
nCross=size(crossTruncations,1);
crossSigma=zeros(nCross,1); crossRaw=zeros(nCross,1);
for j=1:nCross
    cfg=make_cfg(xFinal,crossTruncations(j,1),crossTruncations(j,2));
    r=ni2019_strict_rayleigh_operator(cfg,'TargetOrder',-1, ...
        'EnforceOtherThreshold',true);
    crossSigma(j)=r.sigma_ratio;
    crossRaw(j)=r.strict_residual;
end

% Analytic exterior-window energy of the strict mode. All finite/open and
% threshold amplitudes are exactly zero; only evanescent harmonics remain.
heightList=[.05 .10 .20 .35 .50 .75 1 1.5 2 3 5];
energy=exterior_energy(finalStrict,heightList);

% Use the exact K=15 root for field reconstruction and pole continuation;
% this avoids unnecessary high-order hyperbolic ratios in a visualization.
fieldId=find(rootTruncations(:,2)==15,1);
fieldStrict=rootResult{fieldId}.strict_operator;
xGrid=linspace(0,1,401);
yGrid=linspace(-max(fieldStrict.full_operator.depths),.55,361);
surface=surface_result(fieldStrict);
field=ni2019_reconstruct_field(surface,xGrid,yGrid,'scattered');
pScale=max(abs(field.p),[],'all','omitnan');
pFinite=field.p(isfinite(field.p));
[~,phaseId]=max(abs(pFinite));
fieldPhase=angle(pFinite(phaseId));
pNormalized=field.p*exp(-1i*fieldPhase)/max(pScale,eps);

if runPoleTrack
    poleCfg=rootResult{fieldId}.strict_operator.cfg;
    kappa=logspace(log10(3e-4),-3,17);
    poles=ni2019_continue_complex_pole(poleCfg,kappa, ...
        'TargetOrder',1,'OuterIterations',14, ...
        'InitialScale',-.05-1.41421356237i, ...
        'PredictorExponent',.5,'Display','off');
    poleId=2:numel(poles.kappa);
    validPole=poleId(imag(poles.Omega(poleId))>0 & ...
        poles.continued_sigma_ratio(poleId)<1e-13 & ...
        real(poles.qbar(poleId))<0 & imag(poles.qbar(poleId))<0);
    if numel(validPole)<5
        warning('Only %d pole points passed the outgoing-sheet filter.', ...
            numel(validPole));
        fitIm=[nan nan]; fitQ=[nan nan];
    else
        fitIm=polyfit(log(poles.kappa(validPole)), ...
            log(imag(poles.Omega(validPole))),1);
        fitQ=polyfit(log(poles.kappa(validPole)), ...
            log(poles.Q(validPole)),1);
    end
else
    poleCfg=[]; poles=[]; poleId=[]; validPole=[];
    fitIm=[nan nan]; fitQ=[nan nan];
end

% ------------------------------- Figure -------------------------------
fig=figure('Color','w','Position',[40 40 1540 900]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

nexttile;
imagesc(xGrid,yGrid,real(pNormalized)); axis xy; hold on;
plot([0 1],[0 0],'-','Color',[.2 .2 .2],'LineWidth',.8);
axis equal tight; colormap(gca,ni2019_viridis(256)); cb=colorbar;
cb.Label.String='Re(p)/max|p|'; cb.Color=[.1 .1 .1];
xlabel('x/a'); ylabel('y/a'); title('(a) Strict odd eigenmode');
style_axes(gca,[]);

nexttile;
Kroot=rootTruncations(:,2);
plot(Kroot,.5*(xSequence(:,1)+xSequence(:,3)),'-o','Color',navy,'LineWidth',1.5); hold on;
plot(Kroot,xSequence(:,2),'-s','Color',orange,'LineWidth',1.5);
plot(Kroot,.5*(xSequence(:,4)+xSequence(:,6)),'-^','Color',green,'LineWidth',1.4);
plot(Kroot,xSequence(:,5),'-v','Color',red,'LineWidth',1.4);
plot(Kroot,.5*(xSequence(:,7)+xSequence(:,8)),'-d','Color',gray,'LineWidth',1.4);
grid on; box on; xlabel('Groove truncation K'); ylabel('normalized geometry');
title('(b) Root-parameter convergence');
lg=legend('d_{out}/a','d_c/a','w_{out}/a','w_c/a','g/a', ...
    'Location','eastoutside'); style_axes(gca,lg);

nexttile;
semilogy(crossTruncations(:,2),crossSigma,'-o','Color',navy,'LineWidth',1.6); hold on;
semilogy(Kroot,rootSigma,'o','Color',red,'MarkerFaceColor',red,'MarkerSize',5);
yline(1e-8,'--','Color',gray,'LineWidth',1.0);
grid on; box on; xlabel('Groove truncation K');
ylabel('\sigma_{min}/\sigma_{max}'); title('(c) Three-channel compatibility');
lg=legend('fixed final geometry','polished roots','10^{-8} guide', ...
    'Location','best'); style_axes(gca,lg);

nexttile;
C=finalStrict.groove.C_by_groove;
C=C/max(abs(C(:)));
b=bar(1:3,[abs(C(1,:)).',abs(C(2,:)).'],.78);
b(1).FaceColor=navy; b(2).FaceColor=orange;
xticks(1:3); xticklabels({'outer 1','center','outer 3'});
ylabel('|C_q|/max|C_q|'); ylim([0 1.08]); grid on; box on;
title('(d) Odd-mode internal composition');
lg=legend('q=0','q=1','Location','north'); style_axes(gca,lg);

nexttile;
if runPoleTrack
    loglog(poles.kappa(validPole),poles.Q(validPole),'-o', ...
        'Color',orange,'MarkerFaceColor',orange,'MarkerSize',4, ...
        'LineWidth',1.5); hold on;
    qFit=exp(polyval(fitQ,log(poles.kappa(validPole))));
    loglog(poles.kappa(validPole),qFit,'--','Color',gray,'LineWidth',1.1);
    lg=legend('outgoing-sheet pole', ...
        sprintf('fit: Q\\sim|\\kappa|^{-%0.2f}',-fitQ(1)), ...
        'Location','best');
else
    lg=[];
end
grid on; box on; xlabel('|\kappa|'); ylabel('Q');
title('(e) Pole trend away from \Gamma'); style_axes(gca,lg);

nexttile;
plot(heightList,energy.total/energy.limit,'-o','Color',purple, ...
    'MarkerFaceColor',purple,'LineWidth',1.5); hold on;
yline(1,'--','Color',gray); grid on; box on; ylim([0 1.08]);
xlabel('Exterior height H/a'); ylabel('E_{ext}(H)/E_{ext}(\infty)');
title('(f) Evanescent-energy saturation'); style_axes(gca,[]);

figureFile=fullfile(outputDir,'ThreeGrooveDoubleRayleighBIC_optimized.png');
exportgraphics(fig,figureFile,'Resolution',220);

% ---------------------------- Compact output ---------------------------
rootTable=array2table([rootTruncations xSequence rootSigma rootRaw ...
    grooveFraction oddResidual],'VariableNames', ...
    {'N','K','d1','d2','d3','w1','w2','w3','g1','g2', ...
    'sigma_ratio','raw_residual','groove_fraction','odd_residual'});
crossTable=table(crossTruncations(:,1),crossTruncations(:,2), ...
    crossSigma,crossRaw,'VariableNames', ...
    {'N','K','sigma_ratio','raw_residual'});

fDesign=200e3; cWater=1500; periodMm=1e3*cWater/fDesign;
normalizedNames={'d1','d2','d3','w1','w2','w3','g1','g2'};
parameterTable=array2table(xFinal,'VariableNames',normalizedNames);
parameterTable.period_mm=periodMm;
for j=1:numel(normalizedNames)
    parameterTable.([normalizedNames{j},'_mm'])=xFinal(j)*periodMm;
end
parameterTable.fill_fraction=sum(xFinal(4:8));
parameterTable.frequency_kHz=fDesign/1e3;
writetable(parameterTable,fullfile(outputDir, ...
    'ThreeGrooveDoubleRayleighBIC_parameters.csv'));

fieldPlot=struct('x',xGrid,'y',yGrid,'p_normalized',pNormalized);
globalSummary=[];
if ~isempty(globalResult)
    globalSummary=struct('x',globalResult.x,'score',globalResult.score, ...
        'strict_residual',globalResult.strict_residual, ...
        'parameterization',globalResult.parameterization);
end
matFile=fullfile(outputDir,'ThreeGrooveDoubleRayleighBIC_optimized.mat');
save(matFile,'xSeed','xFinal','xSequence','rootTruncations','rootSigma', ...
    'rootRaw','grooveFraction','oddResidual','crossTruncations', ...
    'crossSigma','crossRaw','rootTable','crossTable','finalParity', ...
    'energy','heightList','fieldPlot','poleCfg','poles','validPole', ...
    'fitIm','fitQ','parameterTable','globalSummary');

fprintf('\nThree-groove Gamma double-Rayleigh optimization\n');
fprintf('  x = ['); fprintf(' %.12g',xFinal); fprintf(' ]\n');
fprintf('  removed orders ='); fprintf(' %d',finalStrict.removed_orders); fprintf('\n');
fprintf('  final finite-root sigma ratio = %.3e, raw residual = %.3e\n', ...
    rootSigma(end),rootRaw(end));
fprintf('  max |x(K=39)-x(K=31)| = %.3e\n', ...
    max(abs(xSequence(end,:)-xSequence(end-1,:))));
fprintf('  parity = %s, odd residual = %.3e\n', ...
    finalParity.parity_label,finalParity.coefficient_residuals.odd);
fprintf('  groove fraction = %.6f, propagating power = %.3e\n', ...
    grooveFraction(end),finalStrict.radiation.total_propagating_power);
fprintf('  E_ext(5a)/E_ext(infinity) = %.9f\n',energy.total(end)/energy.limit);
if runPoleTrack && ~isempty(validPole)
    fprintf('  outgoing-sheet Q max = %.3e, fitted exponent = %.4f\n', ...
        max(poles.Q(validPole)),-fitQ(1));
end
fprintf('  at 200 kHz water: period a = %.6f mm\n',periodMm);
fprintf('  figure: %s\n  data:   %s\n',figureFile,matFile);

assignin('base','ThreeGrooveDoubleRayleighBIC',struct( ...
    'x',xFinal,'root_table',rootTable,'cross_table',crossTable, ...
    'parity',finalParity,'energy',energy,'poles',poles));

function cfg=make_cfg(x,N,K)
cfg=struct('a',1,'lambda',1,'theta_i_deg',0, ...
    'depths',x(1:3),'widths',x(4:6),'gaps',x(7:8), ...
    'N',N,'K',K,'solve_scattering',false);
end

function s=surface_result(r)
op=r.full_operator;
s=struct('A',r.mode.A,'kx',op.kx,'ky',op.ky, ...
    'orders',op.orders,'k0',op.k0,'a',op.a,'lambda',op.lambda, ...
    'theta_i_deg',0,'ky_incident',0,'widths',op.widths, ...
    'depths',op.depths,'gaps',op.gaps,'xleft',op.xleft, ...
    'N',op.N,'K',op.K,'groove_surface_coefficients', ...
    reshape(r.mode.surface_coefficients,op.K,op.L));
end

function e=exterior_energy(r,heightList)
op=r.full_operator; A=r.mode.A(:); ky=op.ky(:); kx=op.kx(:);
weight=1+(abs(kx).^2+abs(ky).^2)/op.k0^2;
total=zeros(size(heightList)); limit=0;
for n=1:numel(A)
    if abs(A(n))<=10*eps, continue; end
    if imag(ky(n))<0
        gamma=-imag(ky(n));
        total=total+op.a*abs(A(n))^2*weight(n)* ...
            (1-exp(-2*gamma*heightList))/(2*gamma);
        limit=limit+op.a*abs(A(n))^2*weight(n)/(2*gamma);
    else
        total=total+op.a*abs(A(n))^2*weight(n)*heightList;
        limit=Inf;
    end
end
e=struct('height',heightList,'total',total,'limit',limit, ...
    'finite_open_power',r.radiation.total_propagating_power, ...
    'threshold_amplitude_norm',r.radiation.grazing_amplitude_norm);
end

function style_axes(ax,lg)
ax.Color='w'; ax.XColor=[.1 .1 .1]; ax.YColor=[.1 .1 .1];
ax.GridColor=[.72 .72 .72]; ax.MinorGridColor=[.84 .84 .84];
ax.Title.Color=[.1 .1 .1]; ax.XLabel.Color=[.1 .1 .1];
ax.YLabel.Color=[.1 .1 .1];
if ~isempty(lg)
    lg.Color='w'; lg.TextColor=[.1 .1 .1]; lg.EdgeColor=[.55 .55 .55];
end
end
