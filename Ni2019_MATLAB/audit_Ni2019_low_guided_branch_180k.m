%% Audit the near-Gamma continuation of the lowest guided branch
% The independent pole census misses this branch for 0 < kappa < 0.04
% because it lies extremely close to the n=0 sound line.  This script
% continues the pole directly and records both the low- and high-truncation
% checks used by the corrected full-band figure.
clear; clc;

rootDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(rootDir,'results','full_complex_band_180k');
S = load(fullfile(rootDir,'results', ...
    'StrictRayleighBIC_180kHz_7p10deg_final.mat'));
D = load(fullfile(dataDir,'positive_pole_census.mat'), ...
    'rootsByK','kappa','tracks');

cfgLow = S.cfg;
cfgLow.N = 89;
cfgLow.K = 11;
cfgLow.solve_scattering = false;

kNear = [0.039; 0.030; 0.020; 0.010];
nearRows = nan(numel(kNear),7);
for j = 1:numel(kNear)
    kap = kNear(j);
    pole = ni2019_refine_outgoing_pole_kappa(cfgLow,kap,0.95*kap, ...
        'OuterIterations',12,'Display','off');
    nearRows(j,:) = [kap,real(pole.Omega),imag(pole.Omega),pole.Q, ...
        pole.sigma_ratio,pole.raw_residual,abs(kap-real(pole.Omega))];
end
nearGamma = array2table(nearRows,'VariableNames', ...
    {'kappa','Omega_real','Omega_imag','Q','sigma_ratio', ...
     'raw_residual','distance_below_n0_line'});
writetable(nearGamma,fullfile(dataDir, ...
    'low_guided_branch_near_gamma.csv'));

cfgHigh = S.cfg;
cfgHigh.N = 121;
cfgHigh.K = 15;
cfgHigh.solve_scattering = false;
kCheck = [0.04; 0.11; 0.30; 0.50];
convRows = nan(numel(kCheck),9);
for j = 1:numel(kCheck)
    kap = kCheck(j);
    [~,ik] = min(abs(D.kappa-kap));
    roots = D.rootsByK{ik};
    [~,it] = min(abs(D.tracks(4).kappa-D.kappa(ik)));
    seed = D.tracks(4).Omega(it);
    [~,ir] = min(abs([roots.Omega]-seed));
    poleLow = roots(ir);
    poleHigh = ni2019_refine_outgoing_pole_kappa( ...
        cfgHigh,D.kappa(ik),poleLow.Omega, ...
        'OuterIterations',10,'Display','off');
    convRows(j,:) = [D.kappa(ik),real(poleLow.Omega),imag(poleLow.Omega), ...
        real(poleHigh.Omega),imag(poleHigh.Omega), ...
        abs(poleHigh.Omega-poleLow.Omega),poleHigh.sigma_ratio, ...
        poleHigh.raw_residual,abs(D.kappa(ik)-real(poleHigh.Omega))];
end
convergence = array2table(convRows,'VariableNames', ...
    {'kappa','Omega_low_real','Omega_low_imag','Omega_high_real', ...
     'Omega_high_imag','absolute_shift','sigma_high','raw_high', ...
     'distance_below_n0_line_high'});
writetable(convergence,fullfile(dataDir, ...
    'low_guided_branch_convergence.csv'));

strictBIC = table(S.x(1),1-S.x(1),0,S.fullsig,S.fullraw, ...
    'VariableNames',{'kappa','Omega_real','Omega_imag', ...
    'full_operator_sigma_ratio','full_operator_raw_residual'});
writetable(strictBIC,fullfile(dataDir,'strict_bic_endpoint.csv'));

save(fullfile(dataDir,'low_guided_branch_audit.mat'), ...
    'nearGamma','convergence','strictBIC','cfgLow','cfgHigh');

disp(nearGamma);
disp(convergence);
disp(strictBIC);
