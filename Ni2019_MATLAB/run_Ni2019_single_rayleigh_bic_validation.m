%% Validate the converged off-Gamma single-Rayleigh BIC
clear; close all; clc;
navy=[.035 .10 .30]; orange=[1 .36 .08]; gray=[.45 .45 .45]; red=[.90 .05 .06];

% N=101, K=14 optimized geometry. All lengths are normalized by a.
cfg=struct('a',1,'lambda',1/0.759964678141,'theta_i_deg',0, ...
    'depths',[0.682459207769 0.122072875277], ...
    'widths',[0.658186442030 0.052212626051], ...
    'gaps',0.144800546474,'N',101,'K',14,'solve_scattering',false);
kappa0=0.240035321859; Omega0=1-kappa0;
cfg.theta_i_deg=asind(kappa0/Omega0);

% The physical leaky branch is deltaKappa<0 for exp(+i*omega*t).
dk=-logspace(-6,-3,37);
poles=ni2019_continue_offgamma_pole(cfg,kappa0,dk,'TargetOrder',-1, ...
    'OuterIterations',10,'InitialScale',.1+.1i);
fitId=abs(poles.delta_kappa)>=1e-4 & abs(poles.delta_kappa)<=1e-3;
fitQ=polyfit(log(abs(poles.delta_kappa(fitId))),log(poles.Q(fitId)),1);
fitIm=polyfit(log(abs(poles.delta_kappa(fitId))), ...
    log(imag(poles.Omega(fitId))),1);
pQ=-fitQ(1); pIm=fitIm(1);

op0=ni2019_full_eigen_operator(cfg); ss=svd(op0.Fscaled);
[~,~,V]=svd(op0.Fscaled,'econ'); z=V(:,end)./op0.column_scale.'; z=z/norm(z);
A=z(1:op0.N); finiteOpen=abs(imag(op0.ky))<1e-9*op0.k0 & ...
    real(op0.ky)>1e-7*op0.k0;
rad=sum(abs(A(finiteOpen)).^2)/sum(abs(A).^2);

fprintf('Validated single-Rayleigh BIC\n');
fprintf('  kappa_BIC = %.12f, Omega_BIC = %.12f, theta = %.6f deg\n', ...
    kappa0,Omega0,asind(kappa0/Omega0));
fprintf('  [d1 d2]/a = [%.12f %.12f]\n',cfg.depths);
fprintf('  [w1 w2 g]/a = [%.12f %.12f %.12f]\n',cfg.widths,cfg.gaps);
fprintf('  endpoint sigma_min/sigma_max = %.3e\n',ss(end)/ss(1));
fprintf('  endpoint finite-channel radiation fraction = %.3e\n',rad);
fprintf('  Im(Omega) = %.3e at |delta kappa| = %.1e\n', ...
    imag(poles.Omega(2)),abs(poles.delta_kappa(2)));
fprintf('  Q exponent = %.6f; Im(Omega) exponent = %.6f\n',pQ,pIm);

fig=figure('Color','w','Position',[80 80 1130 440]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile;
plot(poles.kappa,real(poles.Omega),'-','Color',navy,'LineWidth',1.8); hold on;
plot(poles.kappa,1-poles.kappa,'-','Color',orange,'LineWidth',1.5);
plot(kappa0,Omega0,'o','MarkerFaceColor',red,'MarkerEdgeColor','w','MarkerSize',8);
grid on; box on; xlabel('\kappa=k_Ba/(2\pi)'); ylabel('Re \Omega');
lg=legend('complex pole','n=-1 Rayleigh anomaly','BIC','Location','best');
paper_legend(lg);
title('(a) Pole reaches Rayleigh branch point');
paper_axes(gca);

nexttile;
x=abs(poles.delta_kappa(fitId)); q=poles.Q(fitId);
loglog(x,q,'-','Color',orange,'LineWidth',1.8); hold on;
reference=q(end)*(x/x(end)).^-3;
loglog(x,reference,'--','Color',gray,'LineWidth',1.4);
grid on; box on; xlabel('|\Delta\kappa|'); ylabel('Q-factor');
lg=legend('computed pole','|\Delta\kappa|^{-3}','Location','best');
paper_legend(lg);
title(['(b) Q \propto |\Delta\kappa|^{-',num2str(pQ,'%.3f'),'}']);
paper_axes(gca);

nexttile;
qbar=1e3*poles.qbar(2:end); c=log10(abs(poles.delta_kappa(2:end)));
surface([real(qbar);real(qbar)],[imag(qbar);imag(qbar)],zeros(2,numel(qbar)), ...
    [c;c],'EdgeColor','interp','FaceColor','none','LineWidth',2); hold on;
plot(0,0,'o','MarkerFaceColor',red,'MarkerEdgeColor','w','MarkerSize',8);
axis equal; grid on; box on;
xlabel('10^3 Re(q_ya/2\pi)'); ylabel('10^3 Im(q_ya/2\pi)');
colormap(ni2019_viridis(256)); cb=colorbar; cb.Label.String='log_{10}|\Delta\kappa|';
title('(c) Physical-sheet pole trajectory');
paper_axes(gca);

outputDir=fullfile(pwd,'results');
if ~exist(outputDir,'dir'), mkdir(outputDir); end
exportgraphics(fig,fullfile(outputDir,'SingleRayleighBIC_validation.png'),'Resolution',220);
save(fullfile(outputDir,'SingleRayleighBIC_validation.mat'), ...
    'cfg','kappa0','Omega0','poles','fitQ','fitIm','pQ','pIm','rad');

function paper_axes(ax)
set(ax,'Color','w','XColor','k','YColor','k','GridColor',[.75 .75 .75], ...
    'MinorGridColor',[.86 .86 .86]);
ax.Title.Color='k'; ax.XLabel.Color='k'; ax.YLabel.Color='k';
end

function paper_legend(lg)
set(lg,'Color','w','TextColor','k','EdgeColor',[.25 .25 .25]);
end
