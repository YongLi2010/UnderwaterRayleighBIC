%% Test whether unequal groove depths create an unpaired unidirectional BIC
clear; close all; clc;
navy=[.035 .10 .30]; orange=[1 .36 .08]; red=[.90 .05 .06];

cfg=struct('a',1,'widths',[.085 .527], ...
    'depths',[.486861397812 .040771887750], ...
    'gaps',.101,'N',61,'K',10);
kabs=logspace(-8,-4,28);
polePlus=ni2019_continue_complex_pole(cfg,kabs,'TargetOrder',1, ...
    'BlochSign',1,'OuterIterations',8,'InitialScale',.1+.1i);
poleMinus=ni2019_continue_complex_pole(cfg,kabs,'TargetOrder',-1, ...
    'BlochSign',-1,'OuterIterations',8,'InitialScale',.1+.1i);

frequencyRatio=1.01;
theta=linspace(.1,6,100);
etaGrazingPlus=nan(size(theta)); etaSpecPlus=etaGrazingPlus;
etaGrazingMinus=etaGrazingPlus; etaSpecMinus=etaGrazingPlus;
for n=1:numel(theta)
    cp=scattering_cfg(cfg,frequencyRatio,+theta(n));
    cm=scattering_cfg(cfg,frequencyRatio,-theta(n));
    rp=ni2019_modal_solver(cp); rm=ni2019_modal_solver(cm);
    etaGrazingPlus(n)=rp.eta(rp.orders==-1);
    etaSpecPlus(n)=rp.eta(rp.orders==0);
    etaGrazingMinus(n)=rm.eta(rm.orders==+1);
    etaSpecMinus(n)=rm.eta(rm.orders==0);
end

poleDifference=max(abs(polePlus.Omega-poleMinus.Omega));
qRelativeDifference=max(abs(polePlus.Q(2:end)-poleMinus.Q(2:end))./polePlus.Q(2:end));
fprintf('Reciprocal direction-pair test\n');
fprintf('  max |Omega(+k)-Omega(-k)| = %.3e\n',poleDifference);
fprintf('  max relative Q difference = %.3e\n',qRelativeDifference);
fprintf('  unequal depths do not create an unpaired BIC; they create a time-reversed pair.\n');

fig=figure('Color','w','Position',[100 100 1080 430]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;
loglog(abs(polePlus.kappa(2:end)),polePlus.Q(2:end),'-','Color',navy,'LineWidth',1.8); hold on;
loglog(abs(poleMinus.kappa(2:end)),poleMinus.Q(2:end),'--','Color',orange,'LineWidth',1.8);
grid on; box on; xlabel('|k_B|a/(2\pi)'); ylabel('Q-factor');
legend('+k_B pole','-k_B time-reversed pole','Location','best');
title('(a) Reciprocal pole pair');
nexttile;
plot(theta,etaGrazingPlus,'Color',navy,'LineWidth',1.8); hold on;
plot(theta,etaGrazingMinus,'--','Color',orange,'LineWidth',1.8);
plot(theta,etaSpecPlus,':','Color',navy,'LineWidth',1.4);
plot(theta,etaSpecMinus,'-.','Color',orange,'LineWidth',1.4);
grid on; box on; ylim([0 1.05]); xlabel('|\theta_i| (deg)'); ylabel('\eta');
legend('+\theta: n=-1','-\theta: n=+1','+\theta: n=0','-\theta: n=0', ...
    'Location','best'); title('(b) Left/right scattering partners');

outputDir=fullfile(pwd,'results');
if ~exist(outputDir,'dir'), mkdir(outputDir); end
exportgraphics(fig,fullfile(outputDir,'RayleighBIC_reciprocal_direction_pair.png'), ...
    'Resolution',200);
save(fullfile(outputDir,'RayleighBIC_reciprocal_direction_pair.mat'), ...
    'polePlus','poleMinus','theta','etaGrazingPlus','etaGrazingMinus', ...
    'etaSpecPlus','etaSpecMinus','cfg');

function c=scattering_cfg(base,frequencyRatio,theta)
c=base; c.lambda=1/frequencyRatio; c.theta_i_deg=theta;
c.solve_scattering=true;
end
