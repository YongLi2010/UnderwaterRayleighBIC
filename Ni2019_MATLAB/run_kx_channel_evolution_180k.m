%% Radiation-channel evolution along the 180-kHz Rayleigh-BIC pole branch
% Numerical calculation only. Plotting is performed by the companion Python
% script so that all manuscript rendering follows one graphics backend.
clear; clc;

rootDir = fileparts(mfilename('fullpath'));
resultDir = fullfile(rootDir,'results','kx_channel_evolution_180k');
if ~exist(resultDir,'dir'), mkdir(resultDir); end

S = load(fullfile(rootDir,'results', ...
    'StrictRayleighBIC_180kHz_7p10deg_final.mat'));
cfg = S.cfg;
cfg.solve_scattering = false;
kappaBIC = S.x(1);
OmegaBIC = 1-kappaBIC;

% Logarithmic spacing resolves the common radiation zero while retaining a
% sufficiently broad interval to expose the two sides of the Rayleigh point.
deltaAbs = logspace(-6,log10(2.5e-2),33);
poleNeg = ni2019_continue_offgamma_pole(cfg,kappaBIC,-deltaAbs, ...
    'TargetOrder',-1,'OuterIterations',10,'InitialScale',.1+.1i);
polePos = ni2019_continue_offgamma_pole(cfg,kappaBIC,+deltaAbs, ...
    'TargetOrder',-1,'OuterIterations',10,'InitialScale',.1+.1i);

neg = extract_channels(cfg,poleNeg);
pos = extract_channels(cfg,polePos);

% Insert the exact strict endpoint rather than a finite numerical plotting
% floor. Both pressure amplitudes are constrained to zero at this point.
endpoint = table(kappaBIC,0,OmegaBIC,0,0,0,0,0,0,0,0,0,0, ...
    'VariableNames',neg.Properties.VariableNames);
data = [flipud(neg); endpoint; pos];
data = sortrows(data,'kappa');

% Fit the near-BIC radiation-amplitude exponents on both sides together.
fitMask = abs(data.delta_kappa)>=1e-6 & abs(data.delta_kappa)<=1e-3;
p0 = polyfit(log10(abs(data.delta_kappa(fitMask))), ...
    log10(data.A0_over_C(fitMask)),1);
pm = polyfit(log10(abs(data.delta_kappa(fitMask))), ...
    log10(data.Am1_over_C(fitMask)),1);

writetable(data,fullfile(resultDir,'kx_channel_evolution.csv'));
save(fullfile(resultDir,'kx_channel_evolution.mat'), ...
    'data','poleNeg','polePos','p0','pm','kappaBIC','OmegaBIC','cfg');

fprintf('kx-channel evolution complete (N=%d, K=%d)\n',cfg.N,cfg.K);
fprintf('  kappa_BIC = %.15g, Omega_BIC = %.15g\n',kappaBIC,OmegaBIC);
fprintf('  |A0|/||C|| ~ |Delta kappa|^(%.6f)\n',p0(1));
fprintf('  |A-1|/||C|| ~ |Delta kappa|^(%.6f)\n',pm(1));
fprintf('  max continued sigma ratio = %.3e\n', ...
    max([data.sigma_ratio;0]));

function T = extract_channels(cfg,poles)
n = numel(poles.kappa)-1;
kappa = poles.kappa(2:end).';
delta_kappa = poles.delta_kappa(2:end).';
Omega_real = real(poles.Omega(2:end)).';
Omega_imag = imag(poles.Omega(2:end)).';
q_real = real(poles.qbar(2:end)).';
q_imag = imag(poles.qbar(2:end)).';
A0_over_C = nan(n,1);
Am1_over_C = nan(n,1);
sigma_ratio = nan(n,1);
A0_real = nan(n,1); A0_imag = nan(n,1);
Am1_real = nan(n,1); Am1_imag = nan(n,1);

reference = [];
for j = 1:n
    op = ni2019_full_eigen_operator_complex(cfg,2*pi*poles.Omega(j+1), ...
        2*pi*kappa(j),-1,2*pi*poles.qbar(j+1));
    [~,Sv,V] = svd(op.Fscaled,'econ');
    singular = diag(Sv);
    y = V(:,end);
    z = y./op.column_scale.';
    C = z(op.N+1:end);
    z = z/max(norm(C),eps);

    % Fix the arbitrary eigenvector phase continuously for export. Magnitudes
    % are phase invariant; this convention also makes the complex data usable.
    if isempty(reference)
        [~,idRef] = max(abs(z(op.N+1:end)));
        phase = angle(z(op.N+idRef));
    else
        phase = angle(reference'*z);
    end
    z = z*exp(-1i*phase);
    reference = z;

    A = z(1:op.N);
    id0 = find(op.orders==0,1);
    idm = find(op.orders==-1,1);
    A0_over_C(j) = abs(A(id0));
    Am1_over_C(j) = abs(A(idm));
    A0_real(j) = real(A(id0)); A0_imag(j) = imag(A(id0));
    Am1_real(j) = real(A(idm)); Am1_imag(j) = imag(A(idm));
    sigma_ratio(j) = singular(end)/singular(1);
end

T = table(kappa,delta_kappa,Omega_real,Omega_imag,q_real,q_imag, ...
    A0_over_C,Am1_over_C,sigma_ratio,A0_real,A0_imag,Am1_real,Am1_imag);
end
