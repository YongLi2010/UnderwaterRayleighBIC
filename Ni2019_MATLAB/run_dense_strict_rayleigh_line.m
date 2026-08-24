%% Full-truncation scan of the strict homogeneous residual on the Rayleigh line
% The scan keeps the audited final geometry fixed and enforces both pressure
% conditions A_0=0 and A_-1=0.  It is intended for PRL Fig. 2(c), where a
% vanishing normal flux alone would be an insufficient BIC diagnostic.
clear; clc;

rootDir = fileparts(mfilename('fullpath'));
resultDir = fullfile(rootDir,'results');
paperDir = fullfile(fileparts(rootDir),'arxiv_theory_paper');
dataDir = fullfile(paperDir,'figure_data');
if ~exist(dataDir,'dir'), mkdir(dataDir); end

D = load(fullfile(resultDir,'StrictRayleighBIC_200kHz_min1mm.mat'));
x = D.xFinal(:).';
kappaBIC = x(1);

% Symmetric 5e-6 spacing; the audited endpoint is included exactly.
deltaKappa = linspace(-4e-4,4e-4,161).';
kappa = kappaBIC + deltaKappa;
Omega = 1-kappa;
thetaDeg = asind(kappa./Omega);
sigmaRelative = nan(size(kappa));
rawResidual = nan(size(kappa));

N = 313;
K = 39;
fprintf('Strict Rayleigh-line scan: %d points, N=%d, K=%d\n',numel(kappa),N,K);
tic;
parfor j = 1:numel(kappa)
    cfg = struct('a',1,'lambda',1/Omega(j), ...
        'theta_i_deg',thetaDeg(j), ...
        'depths',x(2:3),'widths',x(4:5),'gaps',x(6), ...
        'N',N,'K',K,'solve_scattering',false);
    R = ni2019_strict_rayleigh_operator(cfg,'TargetOrder',-1, ...
        'KyTolerance',1e-6);
    sigmaRelative(j) = R.sigma_ratio;
    rawResidual(j) = R.residual.reduced_raw;
end
elapsedSeconds = toc;

T = table(deltaKappa,kappa,Omega,thetaDeg,sigmaRelative,rawResidual);
writetable(T,fullfile(dataDir,'fig2_strict_rayleigh_line.csv'));
save(fullfile(resultDir,'DenseStrictRayleighLine_161.mat'),'T','x','N','K', ...
    'kappaBIC','elapsedSeconds');

[minimumSigma,minimumIndex] = min(sigmaRelative);
fprintf('Completed in %.2f s. Minimum sigma_rel=%.16g at kappa=%.16g\n', ...
    elapsedSeconds,minimumSigma,kappa(minimumIndex));
