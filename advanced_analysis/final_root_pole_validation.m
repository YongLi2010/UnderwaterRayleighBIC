%% Same-geometry pole continuation for the final (N,K)=(313,39) strict root
% The earlier manuscript cache belongs to an intermediate (121,15) root.
% This script deliberately starts from the final retained geometry and uses
% the same truncation for both the endpoint and the continued leaky pole.
clear; close all; clc;

thisFile = mfilename('fullpath');
analysisDir = fileparts(thisFile);
repoDir = fileparts(analysisDir);
solverDir = fullfile(repoDir,'Ni2019_MATLAB');
outDir = fullfile(analysisDir,'final_root_pole');
if ~exist(outDir,'dir'), mkdir(outDir); end
addpath(solverDir);

D = load(fullfile(solverDir,'results', ...
    'StrictRayleighBIC_200kHz_min1mm.mat'));
x = D.xFinal(:).';
kappa0 = x(1);
Omega0 = 1-kappa0;
N = 313; K = 39;
cfg = struct('a',1,'lambda',1/Omega0, ...
    'theta_i_deg',asind(kappa0/Omega0), ...
    'depths',x(2:3),'widths',x(4:5),'gaps',x(6), ...
    'N',N,'K',K,'solve_scattering',false);

track = ni2019_track_leaky_pole_to_rayleigh(cfg,kappa0, ...
    'TargetOrder',-1,'DeltaStart',-1e-3,'DeltaEnd',-1e-5, ...
    'NumSteps',13,'ScanRealFactors',[-.15,0,.15], ...
    'ScanImagFactors',[-.55,-.40,-.25,-.10], ...
    'OuterIterations',7,'ScanOuterIterations',8, ...
    'MinOverlap',.90,'Verbose',true,'Display','off');

finite = track.delta_kappa~=0 & isfinite(track.Q) & track.Q>0 & ...
    track.sigma_ratio<1e-7;
fitBand = finite & abs(track.delta_kappa)>=3e-5 & ...
    abs(track.delta_kappa)<=1e-3;
if nnz(fitBand)<4
    error('Too few converged same-geometry pole points for a local fit.');
end
fitQ = polyfit(log(abs(track.delta_kappa(fitBand))), ...
    log(track.Q(fitBand)),1);

% Direct endpoint check with target ky imposed as exactly zero.
opEnd = ni2019_full_eigen_operator_complex(cfg,2*pi*Omega0, ...
    2*pi*kappa0,-1,0);
[~,S,V] = svd(opEnd.Fscaled,'econ');
sv = diag(S);
zEnd = V(:,end)./transpose(opEnd.column_scale);
zEnd = zEnd/max(norm(zEnd),eps);
fullEndpointSigma = sv(end)/max(sv(1),eps);
fullEndpointRaw = norm(opEnd.F*zEnd)/max(norm(zEnd),eps);

pointTable = table(track.delta_kappa(:),track.kappa(:), ...
    real(track.Omega(:)),imag(track.Omega(:)),track.Q(:), ...
    track.sigma_ratio(:),track.mode_overlap(:), ...
    track.radiation_fraction(:),track.groove_content(:), ...
    'VariableNames',{'delta_kappa','kappa','Re_Omega','Im_Omega','Q', ...
    'sigma_ratio','mode_overlap','finite_channel_radiation_fraction', ...
    'groove_content'});
writetable(pointTable,fullfile(outDir,'final_root_pole_track.csv'));

[~,nearIndex] = min(abs(abs(track.delta_kappa)-5e-5));
[~,ordinaryIndex] = min(abs(abs(track.delta_kappa)-1e-3));
save(fullfile(outDir,'final_root_pole_track.mat'),'D','cfg','track', ...
    'fitQ','fitBand','fullEndpointSigma','fullEndpointRaw','zEnd', ...
    'nearIndex','ordinaryIndex','pointTable');

fprintf('\nFINAL-ROOT POLE AUDIT\n');
fprintf('  same branch = %d; warning = %d; min overlap = %.9f\n', ...
    track.same_branch,track.numerical_warning, ...
    min(track.mode_overlap(2:end-1)));
fprintf('  endpoint full sigma = %.6e; raw = %.6e\n', ...
    fullEndpointSigma,fullEndpointRaw);
fprintf('  nearest finite Q = %.6e at |dk| = %.3e\n', ...
    max(track.Q(finite)),min(abs(track.delta_kappa(finite))));
fprintf('  local log-log Q slope = %.8f\n',fitQ(1));

%% PRL-style vector diagnostic
blue=[.10,.31,.55]; orange=[.86,.38,.08]; red=[.76,.15,.18];
teal=[.08,.48,.46]; gray=[.43,.47,.51];
fig=figure('Visible','off','Color','w','Units','inches', ...
    'Position',[.5,.5,7.0,5.25]);
set(fig,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica', ...
    'DefaultAxesFontSize',8,'DefaultTextFontSize',8, ...
    'DefaultAxesLineWidth',.65,'DefaultLineLineWidth',1.15);
tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

nexttile;
plot(track.kappa(finite),real(track.Omega(finite)),'o-', ...
    'Color',blue,'MarkerFaceColor',blue,'MarkerSize',3); hold on;
plot(track.kappa(finite),1-track.kappa(finite),'--','Color',orange);
plot(kappa0,Omega0,'o','Color',red,'MarkerFaceColor',red,'MarkerSize',6);
xlabel('\kappa'); ylabel('Re \Omega_p');
legend('same-geometry pole','n=-1 Rayleigh line','strict endpoint', ...
    'Location','best','Box','off'); title('(a) final-root continuation');

nexttile;
loglog(abs(track.delta_kappa(finite)),track.Q(finite),'o-', ...
    'Color',blue,'MarkerFaceColor',blue,'MarkerSize',3); hold on;
dkFit=sort(abs(track.delta_kappa(fitBand)));
loglog(dkFit,exp(fitQ(2))*dkFit.^fitQ(1),'--','Color',orange);
xlabel('|\Delta\kappa|'); ylabel('Q');
legend('continued pole',sprintf('local slope %.2f',fitQ(1)), ...
    'Location','southwest','Box','off'); title('(b) finite-model linewidth');

nexttile;
semilogy(abs(track.delta_kappa(finite)),track.sigma_ratio(finite),'o-', ...
    'Color',teal,'MarkerFaceColor',teal,'MarkerSize',3); hold on;
yline(1e-7,'--','Color',gray);
xlabel('|\Delta\kappa|'); ylabel('\sigma_{min}/\sigma_{max}');
title('(c) complex-root residual');

nexttile;
plot(abs(track.delta_kappa(finite)),track.mode_overlap(finite),'o-', ...
    'Color',blue,'MarkerFaceColor',blue,'MarkerSize',3); hold on;
yline(.90,'--','Color',gray);
xlabel('|\Delta\kappa|'); ylabel('consecutive mode overlap');
ylim([.88,1.005]); title('(d) branch continuity');

axs=findall(fig,'Type','axes');
for j=1:numel(axs)
    set(axs(j),'Box','on','Layer','top','TickDir','out', ...
        'XColor',[.12,.14,.18],'YColor',[.12,.14,.18], ...
        'GridColor',[.78,.80,.83],'GridAlpha',.34, ...
        'TitleFontWeight','normal');
    grid(axs(j),'on');
end
exportgraphics(fig,fullfile(outDir,'fig_final_root_pole.pdf'), ...
    'ContentType','vector','BackgroundColor','white');
close(fig);
