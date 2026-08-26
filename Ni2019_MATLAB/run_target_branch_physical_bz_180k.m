%% Continue the target Rayleigh-BIC eigenvalue branch across the physical BZ
% The dense near-BIC segment is reused from the verified physical-boundary
% calculation.  Its two edges seed outgoing/decaying pole continuations to
% Gamma and to the Brillouin-zone boundary.  Plotting remains Python-only.
clear; clc;

rootDir = fileparts(mfilename('fullpath'));
nearDir = fullfile(rootDir,'results','kx_channel_evolution_physical_180k');
outDir = fullfile(rootDir,'results','target_branch_physical_bz_180k');
if ~exist(outDir,'dir'), mkdir(outDir); end
S = load(fullfile(nearDir,'physical_sheet_channel_evolution.mat'));
cfg = S.cfg;

near = S.data;
near = near(near.kappa>=0,:);
leftEdge = near(1,:);
rightEdge = near(end,:);

dk = 0.005;
kLow = (leftEdge.kappa-dk):-dk:0;
if isempty(kLow) || kLow(end)>1e-10, kLow = [kLow,0]; end
kHigh = (rightEdge.kappa+dk):dk:0.5;
if isempty(kHigh) || kHigh(end)<0.5-1e-10, kHigh = [kHigh,0.5]; end

low = continue_path(cfg,kLow, ...
    complex(leftEdge.Omega_real,leftEdge.Omega_imag));
high = continue_path(cfg,kHigh, ...
    complex(rightEdge.Omega_real,rightEdge.Omega_imag));

nearOut = near(:,{'kappa','Omega_real','Omega_imag','Q', ...
    'sigma_ratio','raw_residual'});
nearOut.mode_overlap = ones(height(nearOut),1);
nearOut.source_code = ones(height(nearOut),1); % dense near-BIC source

data = [flipud(low);nearOut;high];
data = sortrows(data,'kappa');
[~,uniqueId] = unique(round(data.kappa,12),'stable');
data = data(uniqueId,:);

assert(data.kappa(1)<1e-12 && abs(data.kappa(end)-0.5)<1e-12, ...
    'The physical target branch does not span Gamma to the BZ boundary.');
assert(max(data.sigma_ratio)<1e-9, ...
    'At least one physical-boundary eigenvalue failed the null test.');
assert(min(data.mode_overlap)>0.65, ...
    'Mode-overlap tracking indicates a branch jump.');

writetable(data,fullfile(outDir,'target_branch_physical_bz.csv'));
save(fullfile(outDir,'target_branch_physical_bz.mat'), ...
    'data','near','low','high','cfg','dk');

fprintf('Physical target branch completed: %d points, kappa %.3f to %.3f\n', ...
    height(data),data.kappa(1),data.kappa(end));
fprintf('  Omega(Gamma)=%.10f%+.3ei, Omega(BZ)=%.10f%+.3ei\n', ...
    data.Omega_real(1),data.Omega_imag(1), ...
    data.Omega_real(end),data.Omega_imag(end));
fprintf('  max sigma %.3e, max raw %.3e, min overlap %.6f\n', ...
    max(data.sigma_ratio),max(data.raw_residual),min(data.mode_overlap));

function T = continue_path(cfg,kappa,seed)
n = numel(kappa);
Omega_real = nan(n,1); Omega_imag = nan(n,1); Q = nan(n,1);
sigma_ratio = nan(n,1); raw_residual = nan(n,1);
mode_overlap = nan(n,1); source_code = 2*ones(n,1);
reference = [];
for j = 1:n
    p = ni2019_refine_outgoing_pole_kappa(cfg,kappa(j),seed, ...
        'OuterIterations',12,'Display','off', ...
        'FunctionTolerance',5e-12,'StepTolerance',5e-12);
    if p.sigma_ratio>1e-9 || p.raw_residual>1e-7
        error('Pole failed at kappa %.12g: sigma %.3e raw %.3e', ...
            kappa(j),p.sigma_ratio,p.raw_residual);
    end
    current = p.mode_vector/max(norm(p.mode_vector),eps);
    if isempty(reference)
        mode_overlap(j) = 1;
    else
        mode_overlap(j) = abs(reference'*current);
    end
    reference = current;
    seed = p.Omega;
    Omega_real(j) = real(p.Omega); Omega_imag(j) = imag(p.Omega);
    Q(j) = p.Q; sigma_ratio(j) = p.sigma_ratio;
    raw_residual(j) = p.raw_residual;
end
T = table(kappa(:),Omega_real,Omega_imag,Q,sigma_ratio,raw_residual, ...
    mode_overlap,source_code,'VariableNames', ...
    {'kappa','Omega_real','Omega_imag','Q','sigma_ratio', ...
    'raw_residual','mode_overlap','source_code'});
end
