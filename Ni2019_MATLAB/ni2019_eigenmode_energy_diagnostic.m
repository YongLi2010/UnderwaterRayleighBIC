function diagnostic = ni2019_eigenmode_energy_diagnostic(cfg,varargin)
%NI2019_EIGENMODE_ENERGY_DIAGNOSTIC Diagnose localization of a homogeneous mode.
%   D = NI2019_EIGENMODE_ENERGY_DIAGNOSTIC(CFG) takes the smallest-singular-
%   vector mode of NI2019_FULL_EIGEN_OPERATOR(CFG), reconstructs its physical
%   pressure/velocity field, and separates the finite groove energy from the
%   exterior energy above the periodic surface.  The latter is evaluated in a
%   finite-height window because a nonzero grazing Floquet coefficient has
%   zero normal flux but non-integrable energy at a Rayleigh threshold.
%
%   The unknowns in NI2019_FULL_EIGEN_OPERATOR are [A; Cscaled].  A are
%   outgoing Floquet pressure amplitudes.  Cscaled are groove cosine-mode
%   amplitudes rescaled internally by that operator to avoid tan(beta*d)
%   poles.  This function converts Cscaled back to physical bottom
%   amplitudes before integrating the groove field.
%
%   The dimensionless energy used here is
%       E = integral ( |p|^2 + (|rho*omega*v_x|^2+
%             |rho*omega*v_y|^2)/k0^2 ) dA .
%   It differs from the time-averaged acoustic energy only by one common
%   positive factor, so all localization and window-divergence conclusions
%   are unchanged.  Integration over one period is exact: cosine/sine
%   groove modes and distinct Floquet harmonics are orthogonal in x.
%
%   Optional name/value pairs:
%       TargetOrder       Rayleigh order to diagnose (default -1)
%       HeightList        Exterior window heights (default [1 2 5 10 20])
%       NList             Floquet truncations for convergence
%       KList             Groove truncations for convergence
%       ReturnField       Also return sampled pressure/velocity arrays
%                         (default false)
%       X                 x samples when ReturnField is true
%       Y                 y samples when ReturnField is true
%       Verbose            Print a compact numerical verdict (default true)
%
%   The output contains base.mode (A, physical groove coefficients and a
%   reconstruction-ready result), base.energy, and separate NSeries and
%   KSeries convergence tables.  A strict square-integrability verdict is
%   deliberately conservative: any nonzero finite-flux open coefficient or
%   any nonzero threshold/grazing coefficient makes the exterior norm
%   divergent.

p = inputParser;
addParameter(p,'TargetOrder',-1,@(x)isscalar(x) && isfinite(x));
addParameter(p,'HeightList',[1 2 5 10 20],@(x)isvector(x) && all(x>0));
addParameter(p,'NList',[],@(x)isvector(x) && all(x>=3));
addParameter(p,'KList',[],@(x)isvector(x) && all(x>=1));
addParameter(p,'ReturnField',false,@(x)islogical(x) && isscalar(x));
addParameter(p,'X',[],@(x)isvector(x));
addParameter(p,'Y',[],@(x)isvector(x));
addParameter(p,'Verbose',true,@(x)islogical(x) && isscalar(x));
parse(p,varargin{:});
opts = p.Results;

if ~isfield(cfg,'N'), cfg.N = 101; end
if ~isfield(cfg,'K'), cfg.K = 10; end
if mod(cfg.N,2)~=1
    error('cfg.N must be odd.');
end
if any(opts.HeightList<=0) || any(~isfinite(opts.HeightList))
    error('HeightList must contain finite positive heights.');
end

if isempty(opts.NList)
    opts.NList = unique([max(21,cfg.N-40),max(21,cfg.N-20),cfg.N]);
end
if isempty(opts.KList)
    opts.KList = unique([max(4,cfg.K-6),max(4,cfg.K-3),cfg.K]);
end
opts.NList = round(opts.NList(:).');
opts.KList = round(opts.KList(:).');
if any(mod(opts.NList,2)==0)
    error('Every entry in NList must be odd.');
end

target = round(opts.TargetOrder);
base = evaluate_one(cfg,target,opts.HeightList,opts);

% Keep the two convergence sweeps one-dimensional.  This both exposes which
% truncation controls the result and avoids pretending that a rectangular
% N-by-K table is a physical uncertainty estimate.
NSeries = repmat(summary_template(),1,numel(opts.NList));
for j = 1:numel(opts.NList)
    cj = cfg;
    cj.N = opts.NList(j);
    e = evaluate_one(cj,target,opts.HeightList(1),opts);
    NSeries(j) = summarize_one(e,opts.HeightList(1));
end

KSeries = repmat(summary_template(),1,numel(opts.KList));
for j = 1:numel(opts.KList)
    cj = cfg;
    cj.K = opts.KList(j);
    e = evaluate_one(cj,target,opts.HeightList(1),opts);
    KSeries(j) = summarize_one(e,opts.HeightList(1));
end

diagnostic = struct();
diagnostic.cfg = cfg;
diagnostic.target_order = target;
diagnostic.height_list = opts.HeightList(:).';
diagnostic.base = base;
diagnostic.NSeries = NSeries;
diagnostic.KSeries = KSeries;
diagnostic.assumptions = struct( ...
    'energy_definition','integral(|p|^2+(|rho*omega*v_x|^2+|rho*omega*v_y|^2)/k0^2)dA', ...
    'periodic_x_window','one period; distinct Floquet orders are x-orthogonal', ...
    'groove_integration','analytic cosine/sine modal integrals', ...
    'exterior_integration','analytic finite-height integral for each Floquet order', ...
    'time_convention','exp(+j*omega*t)', ...
    'normalization','Euclidean norm of [A;Cscaled] is one');

diagnostic.verdict = base.verdict;
if opts.Verbose
    print_diagnostic(diagnostic);
end
end

function e = evaluate_one(cfg,target,heightList,opts)
op = ni2019_full_eigen_operator(cfg);
[~,S,V] = svd(op.Fscaled,'econ');
sv = diag(S);
z = V(:,end)./transpose(op.column_scale);
if norm(z)==0 || any(~isfinite(z))
    error('The smallest-singular-vector mode is invalid.');
end
z = z/norm(z);
A = z(1:op.N);
Cscaled = z(op.N+1:end);

% The full operator rescales each groove coefficient by the largest of its
% pressure and velocity factors. Recover the coefficient multiplying
% cos(beta*(y+d)) in the physical field.
[Cphysical,verticalScale,beta] = physical_groove_coefficients(Cscaled,op);
CbyGroove = reshape(Cphysical,op.K,op.L);
surfaceByGroove = zeros(op.K,op.L,'like',Cphysical);
row = 0;
for ell = 1:op.L
    d = op.depths(ell);
    for q = 0:op.K-1
        row = row+1;
        surfaceByGroove(q+1,ell) = Cphysical(row)*cos(beta(row)*d);
    end
end

targetId = op.orders==target;
if ~any(targetId)
    error('TargetOrder %d is outside retained Floquet orders.',target);
end
targetKy = op.ky(targetId);
kyTol = 1e-7*max(abs(op.k0),1);
% At a threshold, classify the requested order as the grazing order even if
% roundoff leaves a tiny positive ky.  Other finite-flux orders remain open.
grazingId = targetId & abs(targetKy)<=kyTol;
finiteOpenId = ~targetId & abs(imag(op.ky))<=1e-8*max(abs(op.k0),1) & ...
    real(op.ky)>kyTol;
otherGrazingId = ~targetId & abs(op.ky)<=kyTol;
evanescentId = ~(finiteOpenId | grazingId | otherGrazingId);

energy = modal_energy(op,A,CbyGroove,beta,grazingId,finiteOpenId, ...
    evanescentId,heightList);
residual = norm(op.F*z)/max(norm(z),eps);
surface = struct('A',A,'kx',op.kx,'ky',op.ky,'orders',op.orders, ...
    'k0',op.k0,'a',op.a,'lambda',op.lambda,'theta_i_deg',cfg.theta_i_deg, ...
    'ky_incident',0,'widths',op.widths,'depths',op.depths,'gaps',op.gaps, ...
    'xleft',op.xleft,'N',op.N,'K',op.K, ...
    'groove_surface_coefficients',surfaceByGroove);

mode = struct('A',A,'C_scaled',Cscaled,'C_physical',Cphysical, ...
    'C_by_groove',CbyGroove,'surface_coefficients',surfaceByGroove, ...
    'orders',op.orders,'kx',op.kx,'ky',op.ky,'beta',beta, ...
    'vertical_scale',verticalScale,'coefficient_norm',norm(z), ...
    'z',z,'reconstruction_result',surface);
if opts.ReturnField
    if isempty(opts.X), x = linspace(0,op.a,241); else, x = opts.X; end
    if isempty(opts.Y)
        y = linspace(-max(op.depths),max(op.a,max(op.depths)),301);
    else
        y = opts.Y;
    end
    mode.field = ni2019_reconstruct_field(surface,x,y,'scattered');
end

targetAmplitude = A(targetId);
coeffScale = max(abs(A));
targetRelative = abs(targetAmplitude)/max(coeffScale,eps);
coefficientFraction = abs(targetAmplitude)^2/max(sum(abs(A).^2),eps);
finiteCoefficientFraction = sum(abs(A(finiteOpenId)).^2)/max(sum(abs(A).^2),eps);
finitePower = energy.finite_open_power;
isThreshold = abs(targetKy)<=kyTol;
strictSquare = isThreshold && targetRelative<=1e-8 && ...
    all(abs(A(finiteOpenId))<=1e-8*max(coeffScale,eps));
if ~isThreshold
    verdict = 'requested order is not at a Rayleigh threshold at this cfg';
elseif targetRelative>1e-8
    verdict = 'not square-integrable: nonzero grazing Floquet amplitude gives E_ext(H)~H';
elseif finitePower>1e-12*max(abs(op.a)*max(abs(A)).^2*abs(op.k0),eps)
    verdict = 'not a BIC: finite-flux open-order radiation remains';
else
    verdict = 'no finite-height divergence detected at this numerical tolerance';
end

e = struct('operator',op,'mode',mode,'singular_values',sv, ...
    'sigma_ratio',sv(end)/max(sv(1),eps),'residual',residual, ...
    'target_ky',targetKy,'target_amplitude',targetAmplitude, ...
    'target_relative_amplitude',targetRelative, ...
    'target_coefficient_fraction',coefficientFraction, ...
    'finite_open_coefficient_fraction',finiteCoefficientFraction, ...
    'finite_open_orders',op.orders(finiteOpenId), ...
    'grazing_orders',op.orders(grazingId|otherGrazingId), ...
    'energy',energy,'strict_square_integrable',strictSquare, ...
    'verdict',verdict);
end

function [Cphysical,scale,beta] = physical_groove_coefficients(Cscaled,op)
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

function energy = modal_energy(op,A,CbyGroove,beta,grazingId,finiteOpenId, ...
    evanescentId,heightList)
% Exact modal integrals in one x period and y in [0,H] outside.
EgrooveByGroove = zeros(1,op.L);
EgrooveByMode = zeros(op.K,op.L);
for ell = 1:op.L
    t = op.widths(ell); d = op.depths(ell);
    for q = 0:op.K-1
        alpha = q*pi/t;
        b = beta(q+1 + (ell-1)*op.K);
        [Jc,Js] = abs_trig_integrals(b,d);
        xc = t*(q==0) + 0.5*t*(q>0);
        xs = 0.5*t*(q>0);
        c = CbyGroove(q+1,ell);
        E = abs(c)^2*(xc*Jc + ...
            (alpha^2/op.k0^2)*xs*Jc + ...
            (abs(b)^2/op.k0^2)*xc*Js);
        EgrooveByMode(q+1,ell) = real(E);
    end
    EgrooveByGroove(ell) = sum(EgrooveByMode(:,ell));
end

kx = op.kx(:); ky = op.ky(:); aa = op.a;
weight = 1 + (abs(kx).^2 + abs(ky).^2)/op.k0^2;
heightList = heightList(:).';
J = zeros(numel(ky),numel(heightList));
for n = 1:numel(ky)
    if abs(imag(ky(n)))<=1e-10*max(abs(op.k0),1) && real(ky(n))>=0
        J(n,:) = heightList;
    elseif imag(ky(n))<0
        gamma = -imag(ky(n));
        J(n,:) = (1-exp(-2*gamma*heightList))/(2*gamma);
    else
        % This branch is only defensive for a non-outgoing operator.
        gamma = abs(imag(ky(n)));
        J(n,:) = (exp(2*gamma*heightList)-1)/(2*gamma);
    end
end
orderEnergy = aa*(abs(A).^2.*weight).*J;
targetEnergy = sum(orderEnergy(grazingId,:),1);
finiteEnergy = sum(orderEnergy(finiteOpenId,:),1);
evanescentEnergy = sum(orderEnergy(evanescentId,:),1);
totalExterior = sum(orderEnergy,1);
powerByOrder = aa/2*real(ky).*abs(A).^2;
energy = struct('groove_total',sum(EgrooveByGroove), ...
    'groove_by_groove',EgrooveByGroove,'groove_by_mode',EgrooveByMode, ...
    'height_list',heightList,'exterior_total',totalExterior, ...
    'exterior_target_grazing',targetEnergy, ...
    'exterior_finite_open',finiteEnergy, ...
    'exterior_evanescent',evanescentEnergy, ...
    'target_energy_slope',sum(orderEnergy(grazingId,:),1)./heightList, ...
    'finite_open_energy_slope',sum(orderEnergy(finiteOpenId,:),1)./heightList, ...
    'target_power',sum(powerByOrder(grazingId)), ...
    'finite_open_power',sum(powerByOrder(finiteOpenId)), ...
    'power_by_order',powerByOrder,'order_energy_weight',weight);
end

function [Jc,Js] = abs_trig_integrals(b,d)
% Integrate |cos(b*s)|^2 and |sin(b*s)|^2 from 0 to d analytically.
bp = b+conj(b); bm = b-conj(b);
E = @(z) exp_integral(z,d);
Jc = 0.25*(E(1i*bp)+E(1i*bm)+E(-1i*bm)+E(-1i*bp));
Js = 0.25*(-E(1i*bp)+E(1i*bm)+E(-1i*bm)-E(-1i*bp));
Jc = real(Jc); Js = real(Js);
% Roundoff can make a theoretically nonnegative integral very slightly
% negative for a nearly zero beta.
Jc = max(Jc,0); Js = max(Js,0);
end

function value = exp_integral(z,d)
if abs(z*d)<1e-8
    value = d*(1 + z*d/2 + (z*d)^2/6 + (z*d)^3/24);
else
    value = expm1(z*d)/z;
end
end

function b = groove_sqrt_local(z)
if real(z)>=0
    b = sqrt(real(z));
else
    b = 1i*sqrt(-real(z));
end
end

function s = summary_template()
s = struct('N',nan,'K',nan,'sigma_ratio',nan,'residual',nan, ...
    'groove_total',nan,'target_relative_amplitude',nan, ...
    'target_coefficient_fraction',nan, ...
    'target_energy_at_height',nan,'finite_open_energy_at_height',nan, ...
    'evanescent_energy_at_height',nan,'finite_open_power',nan, ...
    'target_power',nan,'strict_square_integrable',false);
end

function s = summarize_one(e,H)
s = summary_template();
s.N = e.operator.N; s.K = e.operator.K;
s.sigma_ratio = e.sigma_ratio; s.residual = e.residual;
s.groove_total = e.energy.groove_total;
s.target_relative_amplitude = e.target_relative_amplitude;
s.target_coefficient_fraction = e.target_coefficient_fraction;
s.target_energy_at_height = value_at_height(e.energy.height_list, ...
    e.energy.exterior_target_grazing,H);
s.finite_open_energy_at_height = value_at_height(e.energy.height_list, ...
    e.energy.exterior_finite_open,H);
s.evanescent_energy_at_height = value_at_height(e.energy.height_list, ...
    e.energy.exterior_evanescent,H);
s.finite_open_power = e.energy.finite_open_power;
s.target_power = e.energy.target_power;
s.strict_square_integrable = e.strict_square_integrable;
end

function value = value_at_height(x,y,xq)
if isscalar(x)
    value = y(1);
else
    value = interp1(x,y,xq,'nearest','extrap');
end
end

function print_diagnostic(d)
b = d.base;
fprintf('Eigenmode energy diagnostic\n');
fprintf('  N=%d, K=%d, target order n=%d, ky_target/k0 = %.3e%+.3ei\n', ...
    b.operator.N,b.operator.K,d.target_order,real(b.target_ky/b.operator.k0), ...
    imag(b.target_ky/b.operator.k0));
fprintf('  sigma_min/sigma_max = %.3e; homogeneous residual = %.3e\n', ...
    b.sigma_ratio,b.residual);
fprintf('  target |A|/max|A| = %.6e; finite-open orders = %s\n', ...
    b.target_relative_amplitude,mat2str(b.finite_open_orders.'));
fprintf('  target coefficient fraction = %.6e; finite-open coefficient fraction = %.6e\n', ...
    b.target_coefficient_fraction,b.finite_open_coefficient_fraction);
fprintf('  groove energy = %.6e; target grazing power = %.6e; finite-open power = %.6e\n', ...
    b.energy.groove_total,b.energy.target_power,b.energy.finite_open_power);
fprintf('  H      E_target_grazing      E_finite_open       E_evanescent\n');
for j = 1:numel(b.energy.height_list)
    fprintf('  %-5.3g  %-18.8e  %-18.8e  %-18.8e\n', ...
        b.energy.height_list(j),b.energy.exterior_target_grazing(j), ...
        b.energy.exterior_finite_open(j),b.energy.exterior_evanescent(j));
end
fprintf('  verdict: %s\n',b.verdict);
fprintf('  N-series: ');
for j = 1:numel(d.NSeries)
    fprintf('%d:%0.3e ',d.NSeries(j).N,d.NSeries(j).target_relative_amplitude);
end
fprintf('\n  K-series: ');
for j = 1:numel(d.KSeries)
    fprintf('%d:%0.3e ',d.KSeries(j).K,d.KSeries(j).target_relative_amplitude);
end
fprintf('\n');
end
