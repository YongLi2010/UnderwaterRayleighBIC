%% Complex-frequency pole and Q scaling at the Gamma double-Rayleigh point
% Color conventions reproduce Fig. 2 of 2607.07228v1:
% navy = eigenfrequency/pole, orange = Rayleigh/Q, gray = reference,
% red = BIC, and viridis = parameter colorbar.
clear; close all; clc;

navy=[.035 .10 .30]; orange=[1.00 .36 .08]; gray=[.55 .55 .55];
red=[.90 .05 .06];

lambda0=1;
thresholdCfg=struct('lambda',lambda0,'a',lambda0,'theta_i_deg',0, ...
    'widths',[.085 .527]*lambda0,'depths',[.419 .070]*lambda0, ...
    'gaps',.101*lambda0,'N',61,'K',10,'solve_scattering',false);

% Machine-precision BIC endpoint before complex continuation.
depthScan=ni2019_find_rayleigh_bic(thresholdCfg,'GridSize',31, ...
    'DepthRange',[.01 .49],'TargetOrder',-1);
cfg=thresholdCfg;
cfg.depths=depthScan.best_depths*lambda0;

% For kappa>0 the +1 threshold is the upper Rayleigh branch. Its pole lies
% on the second sheet (Re q<0, Im q<0) while Im Omega>0 gives temporal decay
% under exp(+j*omega*t). Dense logarithmic stepping prevents branch jumps.
kappa=logspace(-9,-3,49);
poles=ni2019_continue_complex_pole(cfg,kappa,'TargetOrder',1, ...
    'OuterIterations',8,'InitialScale',.1+.1i);

valid=poles.kappa>=1e-8 & poles.kappa<=1e-5 & ...
    imag(poles.Omega)>0 & poles.continued_sigma_ratio<1e-8;
fitIm=polyfit(log(poles.kappa(valid)),log(imag(poles.Omega(valid))),1);
fitQ=polyfit(log(poles.kappa(valid)),log(poles.Q(valid)),1);
pIm=fitIm(1); pQ=-fitQ(1);

fprintf('Complex pole on the continued +1 Rayleigh sheet\n');
fprintf('  BIC depths/lambda0 = [%.12f %.12f]\n',depthScan.best_depths);
fprintf('  endpoint sigma_min/sigma_max = %.3e\n',depthScan.best_sigma_ratio);
fprintf('  Im(Omega) -> 0: %.3e at kappa=%.1e\n', ...
    imag(poles.Omega(2)),poles.kappa(2));
fprintf('  fitted Im(Omega) exponent = %.6f\n',pIm);
fprintf('  fitted Q divergence exponent = %.6f\n',pQ);
fprintf('  expected single-Rayleigh exponent = 3; Gamma result is non-cubic.\n');

outputDir=fullfile(pwd,'results');
if ~exist(outputDir,'dir'), mkdir(outputDir); end

fig=figure('Color','w','Position',[100 100 1320 430]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% (a) Real dispersion and the split Rayleigh branches.
nexttile;
id=poles.kappa<=1e-3;
plot(poles.kappa(id),real(poles.Omega(id)),'Color',navy,'LineWidth',1.8); hold on;
plot(poles.kappa(id),1+poles.kappa(id),'Color',orange,'LineWidth',1.5);
plot(poles.kappa(id),1-poles.kappa(id),'Color',gray,'LineWidth',1.3);
plot(0,1,'o','Color',red,'MarkerFaceColor',red,'MarkerSize',6);
grid on; box on; xlabel('\kappa=k_Ba/(2\pi)'); ylabel('Re \Omega');
title('(a) Pole dispersion');
legend('Pole','+1 Rayleigh','-1 Rayleigh','BIC','Location','best');

% (b) Paper-style Im(omega) and Q curves, plus cubic reference.
nexttile;
kv=poles.kappa(valid); imv=imag(poles.Omega(valid)); qv=poles.Q(valid);
yyaxis left;
loglog(kv,imv,'Color',navy,'LineWidth',1.8); ylabel('Im \Omega');
ax=gca; ax.YAxis(1).Color=navy; hold on;
yyaxis right;
loglog(kv,qv,'Color',orange,'LineWidth',1.8); hold on;
anchor=ceil(numel(kv)/2);
qCubic=qv(anchor)*(kv/kv(anchor)).^(-3);
loglog(kv,qCubic,'--','Color',gray,'LineWidth',1.2);
ylabel('Q-factor'); ax.YAxis(2).Color=orange;
xlabel('|\Delta\kappa|'); grid on; box on;
title(sprintf('(b) Q \\propto |\\Delta\\kappa|^{-%0.3f}',pQ));
legend('Im \Omega','Q','cubic reference','Location','best');
ax.YAxis(1).Scale='log'; ax.YAxis(2).Scale='log';

% (c) q-plane Riemann-sheet trajectory with the paper's viridis colorbar.
nexttile;
q=poles.qbar(2:end); colorValue=log10(poles.kappa(2:end));
surface([real(q);real(q)],[imag(q);imag(q)],zeros(2,numel(q)), ...
    [colorValue;colorValue],'EdgeColor','interp','FaceColor','none','LineWidth',2);
hold on; plot(0,0,'o','Color',red,'MarkerFaceColor',red,'MarkerSize',6);
axis equal; grid on; box on; xlabel('Re(qa/2\pi)'); ylabel('Im(qa/2\pi)');
title('(c) Second-sheet pole trajectory');
colormap(gca,ni2019_viridis(256)); cb=colorbar;
cb.Label.String='log_{10}|\Delta\kappa|';

exportgraphics(fig,fullfile(outputDir,'RayleighBIC_complex_pole_Q_scaling.png'), ...
    'Resolution',200);
save(fullfile(outputDir,'RayleighBIC_complex_poles.mat'), ...
    'poles','depthScan','cfg','fitIm','fitQ','pIm','pQ','valid');
assignin('base','RayleighBIC_complex_poles',poles);
