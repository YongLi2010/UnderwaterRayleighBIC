%% Practical optimization of a strict off-Gamma n=-1 Rayleigh-BIC
% The program separates two numerical questions:
%   1. Does each finite modal truncation admit a nonzero homogeneous mode
%      after A_-1 and A_0 are constrained to zero?
%   2. Do the corresponding geometry parameters stabilize as N,K grow?
%
% The optional global search is deliberately disabled by default.  The
% supplied basin seed makes the reproducible refinement take tens of
% seconds rather than several minutes.  Set environment variable
% NI2019_GLOBAL_SEARCH=1 to repeat the particle-swarm/multistart discovery
% stage, or NI2019_SKIP_POLES=1 to skip pole continuation.
clear; close all; clc;

runGlobalSearch=strcmp(getenv('NI2019_GLOBAL_SEARCH'),'1');
runPoleTrack=~strcmp(getenv('NI2019_SKIP_POLES'),'1');
outputDir=fullfile(fileparts(mfilename('fullpath')),'results');
if ~exist(outputDir,'dir'), mkdir(outputDir); end

navy=[.035 .10 .30]; orange=[1.00 .36 .08]; gray=[.55 .55 .55];
red=[.90 .05 .06]; green=[.10 .55 .38];

% Design vector: [kappa,d1/a,d2/a,w1/a,w2/a,g/a].
% This seed is inside the repeatable high-Q basin found by the robust
% multi-truncation search; it is not inserted as a claimed final answer.
xSeed=[.0825310229112419 .656866725522583 .221818379955834 ...
    .575916727584499 .037605905465381 .193285483674839];
lb=[.04 .16 .16 .02 .02 .14];
ub=[.13 .72 .72 .62 .62 .24];

if runGlobalSearch
    searchCfg=make_cfg(xSeed,61,9);
    broad=ni2019_optimize_strict_rayleigh_bic(searchCfg, ...
        'Truncations',[41 7;61 9;81 11], ...
        'KappaRange',[lb(1) ub(1)],'DepthRange',[lb(2) ub(2)], ...
        'WidthRange',[lb(4) ub(4)],'GapRange',[lb(6) ub(6)], ...
        'Starts',6,'UseGlobal',true,'GlobalSwarmSize',36, ...
        'GlobalMaxIterations',30,'MaxFunctionEvaluations',1000, ...
        'RandomSeed',260823,'Display','off');
    xSeed=broad.x;
else
    broad=[];
end

% Solve one exact compatibility zero at each increasing truncation.  The
% previous root is the predictor for the next discretization.
rootTruncations=[121 15;185 23;281 35;313 39];
nRoot=size(rootTruncations,1);
rootResult=cell(nRoot,1);
xSequence=zeros(nRoot,6);
rootSigma=zeros(nRoot,1); rootRaw=zeros(nRoot,1);
grooveFraction=zeros(nRoot,1);
x=xSeed;
for j=1:nRoot
    N=rootTruncations(j,1); K=rootTruncations(j,2);
    cfg=make_cfg(x,N,K);
    rootResult{j}=ni2019_refine_strict_rayleigh_bic(cfg,x, ...
        'Truncation',[N K],'LowerBounds',lb,'UpperBounds',ub, ...
        'MaxFunctionEvaluations',1000,'MaxIterations',120,'Display','off');
    x=rootResult{j}.x;
    xSequence(j,:)=x;
    rootSigma(j)=rootResult{j}.sigma_ratio;
    rootRaw(j)=rootResult{j}.strict_operator.strict_residual;
    grooveFraction(j)=rootResult{j}.strict_operator.groove.pressure_proxy_fraction;
    fprintf('root N=%3d K=%2d: sigma ratio %.3e, raw residual %.3e\n', ...
        N,K,rootSigma(j),rootRaw(j));
end

xFinal=xSequence(end,:);
OmegaFinal=1-xFinal(1);
thetaFinal=asind(xFinal(1)/OmegaFinal);

% Cross-evaluate the final fixed geometry.  A finite-discretization root is
% exact only at the truncation where it was polished; convergence is judged
% from the decreasing envelope and the stability of xSequence.
crossTruncations=[81 10;121 15;161 20;185 23;241 30;281 35;313 39;345 43];
nCross=size(crossTruncations,1);
crossSigma=zeros(nCross,1); crossRaw=zeros(nCross,1);
crossGroove=zeros(nCross,1);
crossResult=cell(nCross,1);
for j=1:nCross
    cfg=make_cfg(xFinal,crossTruncations(j,1),crossTruncations(j,2));
    crossResult{j}=ni2019_strict_rayleigh_operator(cfg,'TargetOrder',-1);
    crossSigma(j)=crossResult{j}.sigma_ratio;
    crossRaw(j)=crossResult{j}.strict_residual;
    crossGroove(j)=crossResult{j}.groove.pressure_proxy_fraction;
end

% A practical pole calculation uses the moderate (121,15) root.  It tracks
% the same nontrivial pole on the chosen outgoing sheet without relying on
% q=0 as a starting seed.  The endpoint is stopped before machine precision
% dominates; the local exponent is diagnostic rather than imposed.
if runPoleTrack
    poleX=xSequence(1,:);
    poleCfg=make_cfg(poleX,121,15);
    poleTrack=ni2019_track_leaky_pole_to_rayleigh(poleCfg,poleX(1), ...
        'DeltaStart',-3e-3,'DeltaEnd',-1e-5,'NumSteps',19, ...
        'OuterIterations',9,'ScanOuterIterations',11, ...
        'Verbose',false,'Display','off');
    poleId=1:numel(poleTrack.delta_kappa)-1;
    fitId=poleId(abs(poleTrack.delta_kappa(poleId))>=3e-5 & ...
        abs(poleTrack.delta_kappa(poleId))<=1e-3 & ...
        imag(poleTrack.Omega(poleId))>1e-14);
    fitIm=polyfit(log(abs(poleTrack.delta_kappa(fitId))), ...
        log(imag(poleTrack.Omega(fitId))),1);
    fitQ=polyfit(log(abs(poleTrack.delta_kappa(fitId))), ...
        log(poleTrack.Q(fitId)),1);
    practicalQ=max(poleTrack.Q(poleId));
else
    poleTrack=[]; poleCfg=[]; fitId=[]; fitIm=[nan nan]; fitQ=[nan nan];
    practicalQ=nan;
end

% Paper-compatible navy/orange/red palette and viridis trajectory colorbar.
fig=figure('Color','w','Position',[80 80 1320 850]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
Kroot=rootTruncations(:,2);
plot(Kroot,xSequence(:,2),'-o','Color',navy,'LineWidth',1.6); hold on;
plot(Kroot,xSequence(:,3),'-s','Color',orange,'LineWidth',1.6);
plot(Kroot,xSequence(:,4),'-^','Color',green,'LineWidth',1.4);
plot(Kroot,xSequence(:,5),'-v','Color',red,'LineWidth',1.4);
plot(Kroot,xSequence(:,6),'-d','Color',gray,'LineWidth',1.4);
grid on; box on; xlabel('Groove truncation K'); ylabel('normalized geometry');
title('(a) Root-parameter convergence');
lg=legend('d_1/a','d_2/a','w_1/a','w_2/a','g/a','Location','eastoutside');
style_axes(gca,lg);

nexttile;
semilogy(crossTruncations(:,2),crossSigma,'-o','Color',navy,'LineWidth',1.7); hold on;
semilogy(Kroot,rootSigma,'o','Color',red,'MarkerFaceColor',red,'MarkerSize',6);
yline(1e-8,'--','Color',gray,'LineWidth',1.1);
grid on; box on; xlabel('Groove truncation K');
ylabel('\sigma_{min}/\sigma_{max}');
title('(b) Strict compatibility residual');
lg=legend('fixed final geometry','polished finite roots','10^{-8} guide','Location','best');
style_axes(gca,lg);

nexttile;
if runPoleTrack
    d=abs(poleTrack.delta_kappa(poleId)); q=poleTrack.Q(poleId);
    loglog(d,q,'-o','Color',orange,'MarkerFaceColor',orange, ...
        'MarkerSize',4,'LineWidth',1.6); hold on;
    qFit=exp(polyval(fitQ,log(d(fitId))));
    loglog(d(fitId),qFit,'--','Color',gray,'LineWidth',1.2);
    lg=legend('complex-pole Q',sprintf('local fit: Q\\sim|\\Delta\\kappa|^{-%0.2f}', ...
        -fitQ(1)),'Location','best');
else
    lg=[];
end
grid on; box on; xlabel('|\Delta\kappa|'); ylabel('Q');
title('(c) Practical high-Q pole trend');
style_axes(gca,lg);

nexttile;
if runPoleTrack
    qp=poleTrack.q(poleId); colorValue=log10(abs(poleTrack.delta_kappa(poleId)));
    surface([real(qp);real(qp)],[imag(qp);imag(qp)],zeros(2,numel(qp)), ...
        [colorValue;colorValue],'EdgeColor','interp','FaceColor','none','LineWidth',2);
    hold on; plot(0,0,'o','Color',red,'MarkerFaceColor',red,'MarkerSize',6);
    colormap(gca,ni2019_viridis(256)); cb=colorbar;
    cb.Label.String='log_{10}|\Delta\kappa|';
    cb.Color=[.1 .1 .1];
end
axis equal; grid on; box on; xlabel('Re(qa/2\pi)'); ylabel('Im(qa/2\pi)');
title('(d) Outgoing-sheet pole trajectory');
style_axes(gca,[]);

figureFile=fullfile(outputDir,'StrictRayleighBIC_optimized_practical.png');
exportgraphics(fig,figureFile,'Resolution',220);

parameterNames={'kappa','d1_over_a','d2_over_a','w1_over_a', ...
    'w2_over_a','gap_over_a'};
parameterTable=array2table(xFinal,'VariableNames',parameterNames);
parameterTable.Omega=OmegaFinal;
parameterTable.theta_deg=thetaFinal;
parameterTable.fill_fraction=sum(xFinal(4:6));
writetable(parameterTable, ...
    fullfile(outputDir,'StrictRayleighBIC_optimized_parameters.csv'));

rootTable=array2table([rootTruncations xSequence rootSigma rootRaw grooveFraction], ...
    'VariableNames',{'N','K','kappa','d1_over_a','d2_over_a', ...
    'w1_over_a','w2_over_a','gap_over_a','sigma_ratio','raw_residual', ...
    'groove_fraction'});
crossTable=table(crossTruncations(:,1),crossTruncations(:,2),crossSigma, ...
    crossRaw,crossGroove,'VariableNames', ...
    {'N','K','sigma_ratio','raw_residual','groove_fraction'});

matFile=fullfile(outputDir,'StrictRayleighBIC_optimized_practical.mat');
save(matFile,'xSeed','xFinal','OmegaFinal','thetaFinal','rootTruncations', ...
    'xSequence','rootSigma','rootRaw','grooveFraction','crossTruncations', ...
    'crossSigma','crossRaw','crossGroove','rootTable','crossTable', ...
    'poleTrack','poleCfg','fitId','fitIm','fitQ','practicalQ','broad');

fprintf('\nPractical strict Rayleigh-BIC optimization\n');
fprintf('  x = ['); fprintf(' %.12g',xFinal); fprintf(' ]\n');
fprintf('  Omega = %.12f, theta = %.9f deg\n',OmegaFinal,thetaFinal);
fprintf('  final finite-root sigma ratio = %.3e\n',rootSigma(end));
fprintf('  max |x(K=39)-x(K=35)| = %.3e\n', ...
    max(abs(xSequence(end,:)-xSequence(end-1,:))));
fprintf('  adjacent K=43 fixed-geometry ratio = %.3e\n',crossSigma(end));
fprintf('  groove participation = %.6f\n',grooveFraction(end));
if runPoleTrack
    fprintf('  tracked practical Q = %.3e; Im(Omega) minimum = %.3e\n', ...
        practicalQ,min(imag(poleTrack.Omega(poleId))));
    fprintf('  observed finite-model Q exponent = %.4f\n',-fitQ(1));
end
fprintf('  figure: %s\n',figureFile);
fprintf('  data:   %s\n',matFile);

assignin('base','StrictRayleighBIC_optimized',struct( ...
    'x',xFinal,'Omega',OmegaFinal,'theta_deg',thetaFinal, ...
    'root_table',rootTable,'cross_table',crossTable,'pole',poleTrack));

function cfg=make_cfg(x,N,K)
Omega=1-x(1);
cfg=struct('a',1,'lambda',1/Omega, ...
    'theta_i_deg',asind(x(1)/Omega), ...
    'depths',x(2:3),'widths',x(4:5),'gaps',x(6), ...
    'N',N,'K',K,'solve_scattering',false);
end

function style_axes(ax,lg)
ax.Color='w';
ax.XColor=[.1 .1 .1]; ax.YColor=[.1 .1 .1];
ax.GridColor=[.72 .72 .72]; ax.MinorGridColor=[.84 .84 .84];
ax.Title.Color=[.1 .1 .1]; ax.XLabel.Color=[.1 .1 .1];
ax.YLabel.Color=[.1 .1 .1];
if ~isempty(lg)
    lg.Color='w'; lg.TextColor=[.1 .1 .1]; lg.EdgeColor=[.55 .55 .55];
end
end
