%% Two-depth eigenvalue search at the first-order Rayleigh threshold
clear; close all; clc;

lambda0=1;
cfg=struct('lambda',lambda0,'a',lambda0,'theta_i_deg',0, ...
    'widths',[.085 .527]*lambda0,'depths',[.419 .070]*lambda0, ...
    'gaps',.101*lambda0,'N',61,'K',10,'solve_scattering',false);

scan=ni2019_find_rayleigh_bic(cfg,'GridSize',41,'DepthRange',[.01 .49]);
fprintf('Rayleigh-threshold two-depth search (a=lambda0)\n');
fprintf('  best d/lambda = [%.9f %.9f]\n',scan.best_depths);
fprintf('  sigma_min/sigma_max = %.3e\n',scan.best_sigma_ratio);
fprintf('  zeroth-order modal fraction = %.3e\n',scan.best_radiation_fraction);
opBest=ni2019_full_eigen_operator(set_depths(cfg,scan.best_depths*lambda0));
fprintf('  grazing |A_-1|, |A_+1| = %.6f, %.6f\n', ...
    abs(scan.mode(opBest.orders==-1)),abs(scan.mode(opBest.orders==1)));

outputDir=fullfile(pwd,'results');
if ~exist(outputDir,'dir'), mkdir(outputDir); end

fig1=figure('Color','w','Position',[100 100 1120 440]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;
imagesc(scan.depth_values,scan.depth_values,log10(scan.sigma_ratio)); axis xy equal tight;
xlabel('d_2/\lambda_0'); ylabel('d_1/\lambda_0'); colorbar;
title('log_{10}(\sigma_{min}/\sigma_{max})'); hold on;
plot(scan.best_depths(2),scan.best_depths(1),'wp','MarkerFaceColor','r','MarkerSize',12);
nexttile;
imagesc(scan.depth_values,scan.depth_values,log10(scan.radiation_fraction)); axis xy equal tight;
xlabel('d_2/\lambda_0'); ylabel('d_1/\lambda_0'); colorbar;
title('log_{10}|A_0|^2 of null vector'); hold on;
plot(scan.best_depths(2),scan.best_depths(1),'wp','MarkerFaceColor','r','MarkerSize',12);
exportgraphics(fig1,fullfile(outputDir,'RayleighBIC_depth_maps.png'),'Resolution',180);

% Argand trajectory of a locally projected characteristic function. Fixed
% left/right null vectors avoid arbitrary swapping of ordinary matrix
% eigenvalues and guarantee that a true operator zero maps to the origin.
d2=scan.best_depths(2)+linspace(-.04,.04,161);
[U0,~,V0]=svd(scan.operator,'econ');
uRef=U0(:,end); vRef=V0(:,end);
mu=complex(nan(size(d2)));
for n=1:numel(d2)
    local=cfg; local.depths=[scan.best_depths(1),d2(n)]*lambda0;
    out=ni2019_full_eigen_operator(local);
    D=out.Fscaled;
    mu(n)=uRef'*D*vRef;
end
fig2=figure('Color','w','Position',[100 100 600 520]);
surface([real(mu);real(mu)],[imag(mu);imag(mu)],zeros(2,numel(mu)),[d2;d2], ...
    'EdgeColor','interp', ...
    'FaceColor','none','LineWidth',2); hold on;
plot(0,0,'kp','MarkerFaceColor','y','MarkerSize',12); axis equal; grid on; box on;
xlabel('Re \mu'); ylabel('Im \mu'); title('Argand trajectory of tracked characteristic value');
cb=colorbar; cb.Label.String='d_2/\lambda_0';
exportgraphics(fig2,fullfile(outputDir,'RayleighBIC_Argand.png'),'Resolution',180);

assignin('base','RayleighBIC_depth_scan',scan);

% Scattering and homogeneous-operator diagnostics across the threshold.
freqRatio=[linspace(.97,.999,100),1,linspace(1.001,1.05,160)];
etaM=nan(size(freqRatio)); eta0=etaM; etaP=etaM; sFreq=etaM;
muFreq=complex(nan(size(freqRatio)));
candidateCfg=set_depths(cfg,scan.best_depths*lambda0);
[Uc,~,Vc]=svd(scan.operator,'econ'); uc=Uc(:,end); vc=Vc(:,end);
for n=1:numel(freqRatio)
    local=candidateCfg;
    local.lambda=lambda0/freqRatio(n);
    if abs(freqRatio(n)-1)>1e-12
        local.solve_scattering=true;
        rs=ni2019_modal_solver(local);
        etaM(n)=rs.eta(rs.orders==-1); eta0(n)=rs.eta(rs.orders==0);
        etaP(n)=rs.eta(rs.orders==1);
    end
    oe=ni2019_full_eigen_operator(local);
    sval=svd(oe.Fscaled);
    sFreq(n)=sval(end)/sval(1);
    muFreq(n)=uc'*oe.Fscaled*vc;
end
fig3=figure('Color','w','Position',[100 100 1120 440]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;
plot(freqRatio,etaM,'LineWidth',1.5); hold on;
plot(freqRatio,eta0,'LineWidth',1.5); plot(freqRatio,etaP,'LineWidth',1.5);
xline(1,'k--','Rayleigh'); grid on; box on; ylim([0 1.05]);
xlabel('f/f_0'); ylabel('\eta_n'); legend('n=-1','n=0','n=+1','Location','best');
title('Scattering across first-order threshold');
nexttile;
semilogy(freqRatio,sFreq,'LineWidth',1.5); xline(1,'k--','Rayleigh');
grid on; box on; xlabel('f/f_0'); ylabel('\sigma_{min}/\sigma_{max}');
title('Homogeneous-operator singular value');
exportgraphics(fig3,fullfile(outputDir,'RayleighBIC_frequency_diagnostics.png'),'Resolution',180);

fig4=figure('Color','w','Position',[100 100 600 520]);
surface([real(muFreq);real(muFreq)],[imag(muFreq);imag(muFreq)], ...
    zeros(2,numel(muFreq)),[freqRatio;freqRatio],'EdgeColor','interp', ...
    'FaceColor','none','LineWidth',2); hold on;
plot(0,0,'kp','MarkerFaceColor','y','MarkerSize',12); axis equal; grid on; box on;
xlabel('Re g'); ylabel('Im g'); title('Frequency-swept Argand trajectory');
cb=colorbar; cb.Label.String='f/f_0';
exportgraphics(fig4,fullfile(outputDir,'RayleighBIC_frequency_Argand.png'),'Resolution',180);

save(fullfile(outputDir,'RayleighBIC_depth_study.mat'), ...
    'scan','d2','mu','cfg','freqRatio','etaM','eta0','etaP','sFreq','muFreq');

function cfg=set_depths(cfg,depths)
cfg.depths=depths;
end
