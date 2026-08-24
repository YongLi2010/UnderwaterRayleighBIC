function result = ni2019_strict_rayleigh_operator(cfg,varargin)
%NI2019_STRICT_RAYLEIGH_OPERATOR Build a strict zero-radiation modal system.
%   R = NI2019_STRICT_RAYLEIGH_OPERATOR(CFG) constructs the pole-free
%   homogeneous operator used by NI2019_FULL_EIGEN_OPERATOR, removes the
%   Floquet-amplitude columns that must vanish at a strict Rayleigh BIC,
%   and returns the smallest right-singular-vector candidate of the resulting
%   tall system.
%
%   The full unknown is z = [A; Cscaled].  A(n) is the exterior pressure
%   amplitude of Floquet order n and Cscaled are the groove coefficients with
%   the vertical scaling used by NI2019_FULL_EIGEN_OPERATOR.  At a target
%   Rayleigh point the strict condition is imposed directly on pressure
%   amplitudes, not on normal flux:
%
%       A(target) = 0,
%       A(n) = 0 for every finite-flux open order n.
%
%   The reduced matrix is formed from the raw full matrix F.  Its row and
%   column scaling is recomputed after the amplitude columns are removed;
%   the pre-existing full Fscaled is deliberately not sliced.  If
%       FscaledReduced = diag(1./rowScale) * Fred * diag(1./colScale),
%   and y is the smallest right singular vector, the physical reduced mode
%   is u = y./colScale.  The full mode is reconstructed by inserting exact
%   zeros in the removed A entries.
%
%   Name/value options:
%       TargetOrder           Rayleigh order to enforce (default -1)
%       EnforceAllFiniteOpen  Enforce all finite-flux open orders as well as
%                             TargetOrder (default true)
%       EnforceFiniteOpen     Alias for EnforceAllFiniteOpen
%       EnforceOtherThreshold Also remove other |ky|~0 orders (default false)
%       KyTolerance           Absolute threshold tolerance (default
%                             1e-8*max(abs(k0),1))
%       NormalizeMode         Normalize full [A;Cscaled] mode (default true)
%       Verbose               Print a compact diagnostic (default false)
%
%   The result contains raw/reduced/scaled matrices, scaling vectors,
%   singular values, residual consistency checks, exact reconstructed
%   amplitudes, channel classification, and groove/exterior diagnostics.

p = inputParser;
addParameter(p,'TargetOrder',-1,@(x)isscalar(x) && isfinite(x));
addParameter(p,'EnforceAllFiniteOpen',true,@(x)islogical(x) && isscalar(x));
addParameter(p,'EnforceFiniteOpen',[],@(x)isempty(x) || (islogical(x) && isscalar(x)));
addParameter(p,'EnforceOtherThreshold',false,@(x)islogical(x) && isscalar(x));
addParameter(p,'KyTolerance',[],@(x)isempty(x) || (isscalar(x) && isfinite(x) && x>=0));
addParameter(p,'NormalizeMode',true,@(x)islogical(x) && isscalar(x));
addParameter(p,'Verbose',false,@(x)islogical(x) && isscalar(x));
parse(p,varargin{:});
opts = p.Results;

if ~isfield(cfg,'N'), cfg.N = 101; end
if ~isfield(cfg,'K'), cfg.K = 10; end
if mod(cfg.N,2) ~= 1
    error('cfg.N must be odd so that the retained Floquet orders are symmetric.');
end

% This call supplies the raw pole-free operator and the physical geometry.
% Only op.F is used below; op.Fscaled must not be sliced after columns are
% removed because its original column normalization is no longer valid.
fullOp = ni2019_full_eigen_operator(cfg);
target = round(opts.TargetOrder);
targetId = (fullOp.orders == target);
if ~any(targetId)
    error('TargetOrder %d is outside the retained Floquet orders.',target);
end

if isempty(opts.KyTolerance)
    kyTol = 1e-8*max(abs(fullOp.k0),1);
else
    kyTol = opts.KyTolerance;
end

% At the real threshold ky is real in the existing outgoing-square-root
% convention.  finiteOpen excludes the target threshold itself by requiring
% a strictly positive real ky.  threshold includes any order at |ky|~0 so
% that the caller can see an unremoved second threshold at Gamma.
finiteOpen = abs(imag(fullOp.ky)) <= kyTol & real(fullOp.ky) > kyTol;
threshold = abs(fullOp.ky) <= kyTol;
if ~isempty(opts.EnforceFiniteOpen)
    enforceFinite = opts.EnforceFiniteOpen;
else
    enforceFinite = opts.EnforceAllFiniteOpen;
end

removedA = targetId;
if enforceFinite
    removedA = removedA | finiteOpen;
end
if opts.EnforceOtherThreshold
    removedA = removedA | (threshold & ~targetId);
end
keepA = ~removedA;
nGroove = fullOp.n_groove;
keepUnknown = [keepA; true(nGroove,1)];

Ffull = fullOp.F;
Fred = Ffull(:,keepUnknown);

% Recompute the conditioning transform for the reduced matrix.  This is the
% mathematically correct transform because deleting columns changes every
% reduced column norm.
rowScale = max(vecnorm(Fred,2,2),sqrt(eps));
Frow = Fred./rowScale;
columnScale = max(vecnorm(Frow,2,1),sqrt(eps));
Fscaled = Frow./columnScale;

[~,S,V] = svd(Fscaled,'econ');
singularValues = diag(S);
if isempty(singularValues) || isempty(V)
    error('The reduced operator has no singular-vector solution.');
end
y = V(:,end);
u = y./columnScale(:);
if opts.NormalizeMode
    u = u/max(norm(u),eps);
end

z = complex(zeros(size(Ffull,2),1));
z(keepUnknown) = u;

% The scaling identity is Fscaled*y = diag(1./rowScale)*Fred*u with
% y = diag(columnScale)*u.  Check it explicitly in both scaled and raw
% coordinates to catch accidental use of the old full-matrix scaling.
yFromU = columnScale(:).*u;
scaledResidual = norm(Fscaled*yFromU)/max(norm(yFromU),eps);
rawResidual = norm(Fred*u)/max(norm(u),eps);
fullResidual = norm(Ffull*z)/max(norm(z),eps);
consistencyResidual = norm(Fred*u - rowScale.*(Fscaled*yFromU))/ ...
    max(norm(Fred*u),eps);

nKeptA = nnz(keepA);
Cscaled = u(nKeptA+1:end);
[Cphysical,verticalScale,beta] = physical_groove_coefficients(Cscaled,fullOp);
depthVector = groove_depth_vector(fullOp);
if isfield(fullOp,'cos_depth_normalized')
    surfaceCoeff = Cscaled.*fullOp.cos_depth_normalized;
else
    surfaceCoeff = Cphysical.*cos(beta.*depthVector);
end
A = z(1:fullOp.N);

% Channel diagnostics use pressure amplitudes and the corresponding normal
% power weight.  A threshold amplitude is reported separately because its
% normal flux is zero although a nonzero pressure field is not square
% integrable in the exterior.
propagatingPowerByOrder = max(real(fullOp.ky),0)./max(abs(fullOp.k0),eps).*abs(A).^2;
finiteOpenPower = sum(propagatingPowerByOrder(finiteOpen));
targetAmplitude = A(targetId);
removedAmplitude = A(removedA);
keptThreshold = threshold & keepA;
keptFiniteOpen = finiteOpen & keepA;
evanescent = ~(finiteOpen | threshold);

grooveNorm = norm(Cphysical);
exteriorNorm = norm(A);
grooveFraction = grooveNorm^2/max(grooveNorm^2 + exteriorNorm^2,eps);
exterior = struct();
exterior.amplitudes = A;
exterior.orders = fullOp.orders;
exterior.kx = fullOp.kx;
exterior.ky = fullOp.ky;
exterior.finite_open_orders = fullOp.orders(finiteOpen);
exterior.threshold_orders = fullOp.orders(threshold);
exterior.evanescent_orders = fullOp.orders(evanescent);
exterior.finite_open_amplitudes = A(finiteOpen);
exterior.threshold_amplitudes = A(threshold);
exterior.kept_finite_open_amplitudes = A(keptFiniteOpen);
exterior.kept_threshold_amplitudes = A(keptThreshold);
exterior.power_by_order = propagatingPowerByOrder;
exterior.finite_open_power = finiteOpenPower;
exterior.total_propagating_power = sum(propagatingPowerByOrder);
exterior.grazing_amplitude_norm = norm(A(threshold));
exterior.evanescent_amplitude_norm = norm(A(evanescent));
exterior.pressure_norm = exteriorNorm;

groove = struct();
groove.C_scaled = Cscaled;
groove.C_physical = Cphysical;
groove.C_by_groove = reshape(Cphysical,fullOp.K,fullOp.L);
groove.surface_coefficients = surfaceCoeff;
groove.beta = beta;
groove.vertical_scale = verticalScale;
groove.pressure_coefficient_norm = grooveNorm;
groove.pressure_proxy_fraction = grooveFraction;
groove.n_groove = nGroove;

result = struct();
result.cfg = cfg;
result.target_order = target;
result.target_index = find(targetId,1);
result.target_amplitude = targetAmplitude;
result.target_is_rayleigh = any(targetId & threshold);
result.ky_tolerance = kyTol;
result.enforce_all_finite_open = enforceFinite;
result.enforce_other_threshold = opts.EnforceOtherThreshold;
result.full_operator = fullOp;
result.F_full_raw = Ffull;
result.F_reduced_raw = Fred;
result.F_reduced_scaled = Fscaled;
result.matrices = struct('full_raw',Ffull,'reduced_raw',Fred, ...
    'reduced_scaled',Fscaled);
result.raw = Ffull;
result.reduced = Fred;
result.scaled = Fscaled;
result.row_scale = rowScale;
result.column_scale = columnScale;
result.singular_values = singularValues;
result.sigma_ratio = singularValues(end)/max(singularValues(1),eps);
result.removed_amplitude_indices = find(removedA);
result.removed_orders = fullOp.orders(removedA);
result.finite_open_orders = fullOp.orders(finiteOpen);
result.threshold_orders = fullOp.orders(threshold);
result.retained_orders = fullOp.orders(keepA);
result.keep_amplitude = keepA;
result.keep_unknown = keepUnknown;
result.dimensions = struct('full',size(Ffull), ...
    'reduced',size(Fred), ...
    'equations',size(Fred,1), ...
    'unknowns',size(Fred,2), ...
    'removed_channels',nnz(removedA), ...
    'is_tall',size(Fred,1)>size(Fred,2));
result.mode = struct('z_full',z,'u_reduced',u,'y_scaled',yFromU, ...
    'A',A,'C_scaled',Cscaled,'C_physical',Cphysical, ...
    'surface_coefficients',surfaceCoeff,'beta',beta, ...
    'vertical_scale',verticalScale,'coefficient_norm',norm(z));
result.residual = struct('scaled',scaledResidual,'reduced_raw',rawResidual, ...
    'full_raw',fullResidual,'scaling_consistency',consistencyResidual);
result.strict_residual = rawResidual;
result.radiation = exterior;
result.groove = groove;
result.exterior_diagnostic = exterior;
result.groove_diagnostic = groove;
result.diagnostics = struct('groove_pressure_proxy_fraction',grooveFraction, ...
    'finite_open_power',finiteOpenPower, ...
    'finite_open_amplitude_norm',norm(A(finiteOpen)), ...
    'threshold_amplitude_norm',norm(A(threshold)), ...
    'removed_amplitude_norm',norm(removedAmplitude), ...
    'removed_amplitude_is_exact_zero',all(removedAmplitude==0));

if opts.Verbose
    fprintf('strict Rayleigh operator: target n=%d, removed n=[%s], ', ...
        target,strtrim(sprintf('%d ',result.removed_orders)));
    fprintf('size %d x %d, sigma ratio %.3e, raw residual %.3e\n', ...
        result.dimensions.equations,result.dimensions.unknowns, ...
        result.sigma_ratio,result.residual.reduced_raw);
    fprintf('  finite-open power %.3e, threshold amplitude %.3e, groove proxy %.3f\n', ...
        finiteOpenPower,exterior.grazing_amplitude_norm,grooveFraction);
end
end

function [Cphysical,scale,beta] = physical_groove_coefficients(Cscaled,op)
if isfield(op,'vertical_scale') && isfield(op,'beta')
    scale=op.vertical_scale;
    beta=op.beta;
    Cphysical=Cscaled./scale;
    return;
end
scale = zeros(op.n_groove,1);
beta = zeros(op.n_groove,1);
row = 0;
for ell = 1:op.L
    d = op.depths(ell);
    t = op.widths(ell);
    for q = 0:op.K-1
        row = row+1;
        alpha = q*pi/t;
        beta(row) = groove_sqrt_local(op.k0^2-alpha^2);
        scale(row) = max([abs(cos(beta(row)*d)), ...
            abs((beta(row)/op.k0)*sin(beta(row)*d)),sqrt(eps)]);
    end
end
Cphysical = Cscaled./scale;
end

function d = groove_depth_vector(op)
d = zeros(op.n_groove,1);
row = 0;
for ell = 1:op.L
    for q = 1:op.K
        row = row+1;
        d(row) = op.depths(ell);
    end
end
end

function value = groove_sqrt_local(z)
if real(z)>=0
    value = sqrt(real(z));
else
    value = 1i*sqrt(-real(z));
end
end
