%% Continue the target Rayleigh-pole sheet across the positive half zone
clear; clc;
rootDir=fileparts(mfilename('fullpath'));
outDir=fullfile(rootDir,'results','full_complex_band_180k');
S=load(fullfile(rootDir,'results','StrictRayleighBIC_180kHz_7p10deg_final.mat'));
cfg=S.cfg; cfg.N=121; cfg.K=15; cfg.solve_scattering=false;
k0=S.x(1);
dLeft=-(0.005:0.005:k0);
dRight=0.005:0.005:(0.5-k0);
fprintf('Target-sheet continuation from kappa=%.8f: %d left, %d right.\n', ...
    k0,numel(dLeft),numel(dRight));
left=ni2019_continue_offgamma_pole(cfg,k0,dLeft,'TargetOrder',-1, ...
    'OuterIterations',10,'InitialScale',.08+.08i,'Display','off');
right=ni2019_continue_offgamma_pole(cfg,k0,dRight,'TargetOrder',-1, ...
    'OuterIterations',10,'InitialScale',.08+.08i,'Display','off');

kappa=[fliplr(left.kappa(2:end)),k0,right.kappa(2:end)];
Omega=[fliplr(left.Omega(2:end)),1-k0,right.Omega(2:end)];
Q=real(Omega)./(2*max(abs(imag(Omega)),realmin));
sigma=[fliplr(left.continued_sigma_ratio(2:end)), ...
    left.continued_sigma_ratio(1),right.continued_sigma_ratio(2:end)];
save(fullfile(outDir,'target_rayleigh_band.mat'),'cfg','kappa','Omega','Q', ...
    'sigma','left','right','k0');
fprintf('Target sheet: kappa %.3f--%.3f, ReOmega %.4f--%.4f, max sigma %.3e\n', ...
    min(kappa),max(kappa),min(real(Omega)),max(real(Omega)),max(sigma));
