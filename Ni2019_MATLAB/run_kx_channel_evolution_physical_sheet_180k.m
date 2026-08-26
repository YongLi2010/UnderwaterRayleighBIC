%% Physical-sheet channel evolution around the strict 180-kHz Rayleigh BIC
% Numerical calculation only.  The outgoing pole is continued away from the
% strict endpoint independently on the two sides.  All closed orders use the
% decaying physical-sheet sign Im(ky)<0; open orders use Re(ky)>0.
clear; clc;

rootDir = fileparts(mfilename('fullpath'));
resultDir = fullfile(rootDir,'results','kx_channel_evolution_physical_180k');
if ~exist(resultDir,'dir'), mkdir(resultDir); end

S = load(fullfile(rootDir,'results', ...
    'StrictRayleighBIC_180kHz_7p10deg_final.mat'));
cfg = S.cfg;
cfg.solve_scattering = false;
kappaBIC = S.x(1);
OmegaBIC = 1-kappaBIC;

% Ascending distance lets each converged pole seed the next calculation.
deltaAbs = logspace(-6,log10(2.5e-2),61).';
left = continue_side(cfg,kappaBIC,OmegaBIC,-deltaAbs);
right = continue_side(cfg,kappaBIC,OmegaBIC,+deltaAbs);

names = left.Properties.VariableNames;
endpoint = array2table(zeros(1,numel(names)),'VariableNames',names);
endpoint.kappa = kappaBIC;
endpoint.Omega_real = OmegaBIC;
endpoint.Q = Inf;
endpoint.sigma_ratio = S.fullsig;
endpoint.raw_residual = S.fullraw;
endpoint.is_m1_open = 0;
endpoint.is_m1_threshold = 1;

data = [flipud(left); endpoint; right];
data = sortrows(data,'kappa');
writetable(data,fullfile(resultDir,'physical_sheet_channel_evolution.csv'));
save(fullfile(resultDir,'physical_sheet_channel_evolution.mat'), ...
    'data','left','right','deltaAbs','kappaBIC','OmegaBIC','cfg');

assert(all(left.ky_m1_imag<0), ...
    'The left n=-1 order did not remain on the decaying physical sheet.');
assert(all(right.ky_m1_real>0), ...
    'The right n=-1 order did not remain on the outgoing open sheet.');
assert(max(data.sigma_ratio)<1e-9, ...
    'At least one continued pole failed the full-operator null test.');
assert(all(left.flux_m1==0), ...
    'A closed n=-1 order was assigned nonzero normal radiative flux.');

fprintf('Physical-sheet channel evolution complete (N=%d,K=%d)\n',cfg.N,cfg.K);
fprintf('  kappa_BIC=%.15g, Omega_BIC=%.15g\n',kappaBIC,OmegaBIC);
fprintf('  max sigma ratio %.3e, max raw residual %.3e\n', ...
    max(data.sigma_ratio),max(data.raw_residual));
fprintf('  left edge: ky_-1=%.5g%+.5gi, A0/C=%.3e, A-1/C=%.3e\n', ...
    left.ky_m1_real(end),left.ky_m1_imag(end), ...
    left.A0_over_C(end),left.Am1_over_C(end));
fprintf('  right edge: ky_-1=%.5g%+.5gi, A0/C=%.3e, A-1/C=%.3e\n', ...
    right.ky_m1_real(end),right.ky_m1_imag(end), ...
    right.A0_over_C(end),right.Am1_over_C(end));

function T = continue_side(cfg,kappaBIC,OmegaBIC,delta)
n = numel(delta);
kappa = kappaBIC+delta;
Omega_real = nan(n,1); Omega_imag = nan(n,1); Q = nan(n,1);
sigma_ratio = nan(n,1); raw_residual = nan(n,1);
A0_over_C = nan(n,1); Am1_over_C = nan(n,1);
flux_0 = nan(n,1); flux_m1 = nan(n,1);
ky_0_real = nan(n,1); ky_0_imag = nan(n,1);
ky_m1_real = nan(n,1); ky_m1_imag = nan(n,1);
is_m1_open = nan(n,1); is_m1_threshold = zeros(n,1);

% The real pole slope is finite at the endpoint and the linewidth vanishes
% faster than linearly.  This seed is close enough for the first Newton step.
seed = complex(OmegaBIC+0.062*delta(1),1e-13);
for j = 1:n
    p = ni2019_refine_outgoing_pole_kappa(cfg,kappa(j),seed, ...
        'OuterIterations',12,'Display','off', ...
        'FunctionTolerance',5e-12,'StepTolerance',5e-12);
    if p.sigma_ratio>1e-9 || p.raw_residual>1e-7
        error('Physical-sheet pole failed at kappa %.12g: sigma %.3e raw %.3e', ...
            kappa(j),p.sigma_ratio,p.raw_residual);
    end
    seed = p.Omega;
    op = p.operator;
    C = p.mode_vector(op.N+1:end);
    scale = max(norm(C),eps);
    id0 = find(op.orders==0,1);
    idm = find(op.orders==-1,1);
    a0 = abs(p.A(id0))/scale;
    am = abs(p.A(idm))/scale;
    ky0 = op.ky(id0);
    kym = op.ky(idm);
    k0 = op.k0;
    openM = real(kym)>0 && abs(imag(kym))<0.05*max(real(kym),eps);

    Omega_real(j) = real(p.Omega); Omega_imag(j) = imag(p.Omega);
    Q(j) = p.Q; sigma_ratio(j) = p.sigma_ratio;
    raw_residual(j) = p.raw_residual;
    A0_over_C(j) = a0; Am1_over_C(j) = am;
    flux_0(j) = max(real(ky0/k0),0)*a0^2;
    flux_m1(j) = double(openM)*max(real(kym/k0),0)*am^2;
    ky_0_real(j) = real(ky0); ky_0_imag(j) = imag(ky0);
    ky_m1_real(j) = real(kym); ky_m1_imag(j) = imag(kym);
    is_m1_open(j) = double(openM);
end

T = table(kappa,delta,Omega_real,Omega_imag,Q,sigma_ratio,raw_residual, ...
    A0_over_C,Am1_over_C,flux_0,flux_m1, ...
    ky_0_real,ky_0_imag,ky_m1_real,ky_m1_imag, ...
    is_m1_open,is_m1_threshold);
end
