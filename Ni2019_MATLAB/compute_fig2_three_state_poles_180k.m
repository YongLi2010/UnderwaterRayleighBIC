%% High-truncation pole samples used by the three-state Fig. 2e field strip
% The existing full-band scan was intentionally performed at a lower
% truncation.  Recompute the two off-BIC samples with the final production
% truncation so that their eigenvectors and residuals are internally
% consistent with the strict endpoint.
clear; clc;

rootDir = fileparts(mfilename('fullpath'));
outDir = fullfile(rootDir,'results','fig2_three_state_fields_180k');
if ~exist(outDir,'dir'), mkdir(outDir); end
S = load(fullfile(rootDir,'results', ...
    'StrictRayleighBIC_180kHz_7p10deg_final.mat'));
cfg = S.cfg;
cfg.solve_scattering = false;
kappaBIC = S.x(1);

% Approach Gamma through small steps on the fixed n=-1 sheet.  The final
% sample is at kappa ~= 10^-7 because the optimized BIC kappa contains that
% offset from 0.1100000000.
dLeft = -(0.005:0.005:0.11);
left = ni2019_continue_offgamma_pole(cfg,kappaBIC,dLeft, ...
    'TargetOrder',-1,'OuterIterations',12, ...
    'InitialScale',0.08+0.08i,'Display','off');

% A clearly post-BIC but still nearby physical-sheet pole.
dRight = 0.01;
right = ni2019_continue_offgamma_pole(cfg,kappaBIC,dRight, ...
    'TargetOrder',-1,'OuterIterations',12, ...
    'InitialScale',0.08+0.08i,'Display','off');

gamma = sample_at(left,numel(left.kappa));
beyond = sample_at(right,numel(right.kappa));
save(fullfile(outDir,'high_truncation_pole_samples.mat'), ...
    'cfg','kappaBIC','left','right','gamma','beyond');

fprintf('High-truncation Fig. 2e pole samples (N=%d,K=%d):\n',cfg.N,cfg.K);
fprintf('  Gamma: kappa=%.12g, Omega=%.12g%+.5gi, q=%.12g%+.5gi, sigma=%.3e\n', ...
    gamma.kappa,real(gamma.Omega),imag(gamma.Omega), ...
    real(gamma.qbar),imag(gamma.qbar),gamma.sigma_ratio);
fprintf('  Beyond: kappa=%.12g, Omega=%.12g%+.5gi, q=%.12g%+.5gi, sigma=%.3e\n', ...
    beyond.kappa,real(beyond.Omega),imag(beyond.Omega), ...
    real(beyond.qbar),imag(beyond.qbar),beyond.sigma_ratio);

function s = sample_at(p,j)
s = struct('kappa',p.kappa(j),'delta_kappa',p.delta_kappa(j), ...
    'Omega',p.Omega(j),'qbar',p.qbar(j),'Q',p.Q(j), ...
    'sigma_ratio',p.continued_sigma_ratio(j), ...
    'radiation_fraction',p.radiation_fraction(j), ...
    'target_order',p.target_order,'convention',p.convention);
end
