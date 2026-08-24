function result = ni2019_radiation_zero_diagnostic(cfg,varargin)
%NI2019_RADIATION_ZERO_DIAGNOSTIC Test the numerator condition separately.
%
%   RESULT = NI2019_RADIATION_ZERO_DIAGNOSTIC(CFG) evaluates the real-axis
%   homogeneous two-groove operator near the supplied Rayleigh candidate.
%   It reports two different quantities:
%
%     pole_denominator
%       The normalized smallest singular value of the full pole-free
%       homogeneous operator.  A real-axis pole requires this quantity to
%       vanish.
%
%     radiation_numerator
%       The Fourier coupling of the normalized near-null groove aperture
%       velocity into every Floquet channel.  The finite-flux, non-target
%       component is the acoustic analogue of the numerator of a pole
%       residue.  It is computed before division by k_y, so the zero is not
%       manufactured by the Rayleigh factor 1/k_y.
%
%   The block equations used by NI2019_FULL_EIGEN_OPERATOR are
%
%       [ B       -C_d ] [ A ] = 0,
%       [ Y       i R  ] [ C ]
%
%   where Y=diag(k_y/k_0), A are outgoing Floquet amplitudes, and C are
%   groove cosine-mode amplitudes with the vertical scaling used by the
%   full operator.  The lower-right block gives
%
%       R C = (1/a) int_aperture v_y(x,0) exp(-i k_x,n x) dx,
%
%   up to the common rho*omega normalization.  For a finite open order,
%
%       A_n = -i k_0 (R C)_n/k_y,n.
%
%   Consequently the raw Fourier coefficient (R C)_n is the useful
%   numerator-like diagnostic.  At the selected Rayleigh order k_y=0, the
%   homogeneous velocity equation itself forces (R C)_target=0.  That
%   threshold identity is reported as target_raw_coupling but is explicitly
%   excluded from the anapole/numerator verdict.  Only finite-flux
%   non-target orders are used for finite_channel_numerator_residual.
%
%   Optional name/value arguments:
%     'Kappa'          Bloch number k_B*a/(2*pi); inferred from CFG.
%     'Omega'          k_0*a/(2*pi); inferred from CFG.
%     'TargetOrder'    Rayleigh order, default -1.
%     'KappaOffsets'   local scan offsets, default [-.03 ... .03].
%     'OmegaOffsets'   local frequency offsets, default [-.02 ... .02].
%     'NList'          external-mode convergence list, default CFG.N.
%     'KList'          groove-mode convergence list, default CFG.K.
%     'GridN/GridK'    truncation used for the local 2-D diagnostic grid.
%     'Verbose'        print a compact summary, default false.
%
%   This diagnostic deliberately does not call a total-reflection minimum
%   an anapole.  A lossless one-port scattering coefficient cannot vanish
%   on the real axis.  The returned radiation numerator is instead the
%   intrinsic aperture-source projection of the near-pole mode.  If the
%   caller needs an isolated-metaatom polarizability, that is a separate
%   model and cannot be extracted from the present periodic boundary-value
%   equations without an additional definition.

p = inputParser;
addParameter(p,'Kappa',[]);
addParameter(p,'Omega',[]);
addParameter(p,'TargetOrder',-1);
addParameter(p,'KappaOffsets',[-.03 -.02 -.01 0 .01 .02 .03]);
addParameter(p,'OmegaOffsets',[-.02 -.01 0 .01 .02]);
addParameter(p,'NList',[]);
addParameter(p,'KList',[]);
addParameter(p,'GridN',[]);
addParameter(p,'GridK',[]);
addParameter(p,'Verbose',false);
parse(p,varargin{:});
opt = p.Results;

if ~isfield(cfg,'a') || ~isscalar(cfg.a) || ~(cfg.a > 0)
    error('cfg.a must be a positive scalar.');
end
a = cfg.a;
if ~isfield(cfg,'N'), cfg.N = 101; end
if ~isfield(cfg,'K'), cfg.K = 10; end
if ~isfield(cfg,'gaps')
    cfg.gaps = zeros(1,max(0,numel(cfg.widths)-1));
end

if isempty(opt.Omega)
    if ~isfield(cfg,'lambda') || ~(cfg.lambda > 0)
        error('Supply Omega or cfg.lambda.');
    end
    Omega0 = a/cfg.lambda;
else
    Omega0 = opt.Omega;
end
if isempty(opt.Kappa)
    if ~isfield(cfg,'theta_i_deg')
        error('Supply Kappa or cfg.theta_i_deg.');
    end
    kappa0 = Omega0*sind(cfg.theta_i_deg);
else
    kappa0 = opt.Kappa;
end
if ~isscalar(Omega0) || ~isreal(Omega0) || ~(Omega0 > 0)
    error('Omega must be a positive real scalar.');
end
if ~isscalar(kappa0) || ~isreal(kappa0)
    error('Kappa must be a real scalar.');
end

target = opt.TargetOrder;
if ~isscalar(target) || target ~= round(target)
    error('TargetOrder must be an integer scalar.');
end

NList = opt.NList;
if isempty(NList), NList = cfg.N; end
KList = opt.KList;
if isempty(KList), KList = cfg.K; end
NList = unique(round(NList(:).'),'stable');
KList = unique(round(KList(:).'),'stable');
if any(NList < 3) || any(mod(NList,2) ~= 1)
    error('Every N in NList must be an odd integer >= 3.');
end
if any(KList < 1)
    error('Every K in KList must be a positive integer.');
end
if isempty(opt.GridN), gridN = NList(end); else, gridN = opt.GridN; end
if isempty(opt.GridK), gridK = KList(end); else, gridK = opt.GridK; end
if ~isscalar(gridN) || mod(gridN,2) ~= 1 || gridN < 3
    error('GridN must be an odd integer >= 3.');
end
if ~isscalar(gridK) || gridK < 1
    error('GridK must be a positive integer.');
end

kappaOffsets = unique(opt.KappaOffsets(:).','stable');
omegaOffsets = unique(opt.OmegaOffsets(:).','stable');
if any(~isreal(kappaOffsets)) || any(~isreal(omegaOffsets))
    error('Local scan offsets must be real.');
end

% The exact target Rayleigh line is Omega=abs(kappa+target).  The selected
% candidate is allowed a small mismatch; the mismatch is returned rather
% than silently projecting the point onto the line.
rayleighOmega0 = abs(kappa0 + target);

centerByN = cell(numel(NList),1);
for ii = 1:numel(NList)
    centerByN{ii} = evaluate_point(kappa0,Omega0,NList(ii),gridK);
end
centerByK = cell(numel(KList),1);
for ii = 1:numel(KList)
    centerByK{ii} = evaluate_point(kappa0,Omega0,gridN,KList(ii));
end
center = evaluate_point(kappa0,Omega0,gridN,gridK);

% Evaluate a rectangular neighborhood.  This is a local diagnostic of the
% two independent conditions; it is not a global pole continuation.
nk = numel(kappaOffsets); no = numel(omegaOffsets);
grid = struct();
grid.kappa = kappa0 + kappaOffsets(:);
grid.Omega = Omega0 + omegaOffsets(:).';
grid.rayleigh_residual = nan(nk,no);
grid.pole_denominator = nan(nk,no);
grid.finite_channel_numerator_residual = nan(nk,no);
grid.finite_channel_power_residual = nan(nk,no);
grid.target_raw_coupling = nan(nk,no);
grid.groove_fraction = nan(nk,no);
grid.external_finite_amplitude = nan(nk,no);
grid.valid = false(nk,no);
if opt.Verbose
    fprintf('Radiation-zero diagnostic grid: %d x %d at N=%d, K=%d\n', ...
        nk,no,gridN,gridK);
end
for ii = 1:nk
    for jj = 1:no
        Om = grid.Omega(jj);
        kap = grid.kappa(ii);
        if Om <= 0 || abs(kap/Om) > 1
            continue;
        end
        ev = evaluate_point(kap,Om,gridN,gridK);
        grid.rayleigh_residual(ii,jj) = ev.rayleigh_residual;
        grid.pole_denominator(ii,jj) = ev.pole_denominator;
        grid.finite_channel_numerator_residual(ii,jj) = ...
            ev.finite_channel_numerator_residual;
        grid.finite_channel_power_residual(ii,jj) = ...
            ev.finite_channel_power_residual;
        grid.target_raw_coupling(ii,jj) = ev.target_raw_coupling;
        grid.groove_fraction(ii,jj) = ev.groove_fraction;
        grid.external_finite_amplitude(ii,jj) = ev.external_finite_amplitude;
        grid.valid(ii,jj) = true;
    end
end

% Locate the minima independently.  If the two minima do not coincide, the
% candidate is only a pole or only a radiation cancellation.  A coarse grid
% cannot prove a zero, so the result also exposes the raw maps and their
% coordinates for refinement by the caller.
[poleMin,poleId] = finite_minimum(grid.pole_denominator,grid.valid);
[radMin,radId] = finite_minimum( ...
    grid.finite_channel_numerator_residual,grid.valid);
if isempty(poleId)
    poleMinCoord = [nan nan];
else
    [ip,jp] = ind2sub(size(grid.pole_denominator),poleId);
    poleMinCoord = [grid.kappa(ip),grid.Omega(jp)];
end
if isempty(radId)
    radMinCoord = [nan nan];
else
    [ir,jr] = ind2sub(size(grid.finite_channel_numerator_residual),radId);
    radMinCoord = [grid.kappa(ir),grid.Omega(jr)];
end
if numel(kappaOffsets) > 1
    dkGrid = min(diff(sort(grid.kappa)));
else
    dkGrid = inf;
end
if numel(omegaOffsets) > 1
    dOmegaGrid = min(diff(sort(grid.Omega)));
else
    dOmegaGrid = inf;
end
if any(isfinite([poleMinCoord radMinCoord]))
    minimaDistance = hypot(poleMinCoord(1)-radMinCoord(1), ...
        poleMinCoord(2)-radMinCoord(2));
else
    minimaDistance = nan;
end
minimaOneCell = minimaDistance <= 1.01*hypot(dkGrid,dOmegaGrid);

% An apparent small residual at one truncation is not enough.  Compare the
% final two entries of each requested convergence sequence.  The loose
% factor-of-ten gate is intentional: it flags accidental cancellation while
% not pretending to provide a rigorous error bound.
convN = convergence_summary(centerByN,'N');
convK = convergence_summary(centerByK,'K');
converged = convN.stable && convK.stable;
centerSmall = center.pole_denominator < 1e-5 && ...
    center.finite_channel_numerator_residual < 1e-5;
candidateWorkingPoint = centerSmall && minimaOneCell;
if candidateWorkingPoint && converged
    verdict = ['At the requested truncation, the real-axis pole and the ', ...
        'finite-channel aperture-radiation zero coincide within the local ', ...
        'grid.  This supports, but does not by itself prove, a Rayleigh BIC.'];
elseif candidateWorkingPoint && ~converged
    verdict = ['The pole and finite-channel radiation numerator are both ', ...
        'small at the requested truncation, but N/K convergence is not ', ...
        'stable; a distinct converged anapole/radiation zero is not ', ...
        'established.'];
elseif ~candidateWorkingPoint
    verdict = ['The current candidate does not show a coincident pole and ', ...
        'finite-channel radiation zero on the supplied local grid.'];
else
    verdict = 'Inconclusive: no finite non-target channel was available for a numerator test.';
end

result = struct();
result.kappa = kappa0;
result.Omega = Omega0;
result.rayleigh_order = target;
result.rayleigh_Omega = rayleighOmega0;
result.rayleigh_residual = center.rayleigh_residual;
result.center = center;
result.center_by_N = centerByN;
result.center_by_K = centerByK;
result.convergence = struct('N',convN,'K',convK,'stable',converged);
result.grid = grid;
result.pole_minimum = struct('value',poleMin,'coordinate',poleMinCoord);
result.radiation_minimum = struct('value',radMin,'coordinate',radMinCoord);
result.minima_distance = minimaDistance;
result.minima_within_one_grid_cell = minimaOneCell;
result.candidate_working_point = candidateWorkingPoint;
result.converged_intersection = candidateWorkingPoint && converged;
result.verdict = verdict;
result.definition = struct( ...
    'pole_denominator','sigma_min(Fscaled)/sigma_max(Fscaled)', ...
    'radiation_map','R = lower_right_block(F)/1i', ...
    'radiation_numerator','R*C for finite-flux non-target Floquet orders', ...
    'threshold_warning','R(target)=0 at ky(target)=0 is a boundary-equation identity, not an anapole test');
result.settings = struct('NList',NList,'KList',KList,'GridN',gridN, ...
    'GridK',gridK,'KappaOffsets',kappaOffsets,'OmegaOffsets',omegaOffsets);

if opt.Verbose
    print_summary(result);
end

    function ev = evaluate_point(kap,Om,N,K)
        if Om <= 0 || abs(kap/Om) > 1+1e-12
            ev = invalid_point(kap,Om,N,K);
            return;
        end
        local = cfg;
        local.N = N;
        local.K = K;
        local.lambda = a/Om;
        local.theta_i_deg = asind(kap/Om);
        local.solve_scattering = false;
        op = ni2019_full_eigen_operator(local);
        [~,S,V] = svd(op.Fscaled,'econ');
        sv = diag(S);
        v = V(:,end)./transpose(op.column_scale);
        v = v/max(norm(v),eps);
        A = v(1:N);
        C = v(N+1:end);
        % In the full operator, the bottom block has columns [A,C], so its
        % C block is exactly i*R.  This extraction avoids duplicating the
        % projection convention and its oblique-k correction here.
        R = op.F(op.n_groove+1:end,N+1:end)/1i;
        coupling = R*C;
        Cnorm = max(norm(C),eps);
        k0 = op.k0;
        ky = op.ky;
        orders = op.orders;
        finite = abs(imag(ky)) <= 1e-8*max(abs(k0),1) & ...
            real(ky) > 1e-8*max(abs(k0),1);
        targetId = orders == target;
        finiteNonTarget = finite & ~targetId;
        % The selected target is excluded even on the side where it becomes
        % propagating: its threshold behavior is diagnosed separately.
        openOrders = orders(finiteNonTarget);
        if any(finiteNonTarget)
            rr = coupling(finiteNonTarget);
            kyOpen = real(ky(finiteNonTarget));
            finiteNumerator = norm(rr)/Cnorm;
            finitePowerNumerator = norm(sqrt(k0./kyOpen).*rr)/Cnorm;
            externalFiniteAmplitude = norm(k0./kyOpen.*rr)/Cnorm;
            extFromR = -1i*k0./kyOpen.*rr;
            consistency = max(abs((kyOpen/k0).*extFromR + 1i*rr));
            consistency = consistency/max(norm(rr),eps);
        else
            finiteNumerator = nan;
            finitePowerNumerator = nan;
            externalFiniteAmplitude = nan;
            consistency = nan;
        end
        if any(targetId)
            targetRaw = abs(coupling(targetId))/Cnorm;
            targetAmplitude = abs(A(targetId));
            targetKy = ky(targetId);
            if finite(targetId)
                targetPower = sqrt(real(targetKy)/k0)*targetAmplitude;
            else
                targetPower = 0;
            end
        else
            targetRaw = nan;
            targetAmplitude = nan;
            targetPower = nan;
        end
        ev = struct();
        ev.kappa = kap;
        ev.Omega = Om;
        ev.rayleigh_residual = Om-abs(kap+target);
        ev.N = N;
        ev.K = K;
        ev.orders = orders;
        ev.kx = op.kx;
        ev.ky = ky;
        ev.pole_denominator = sv(end)/max(sv(1),eps);
        ev.sigma_min = sv(end);
        ev.sigma_max = sv(1);
        ev.mode = v;
        ev.floquet_amplitudes = A;
        ev.groove_amplitudes = C;
        ev.radiation_map = R;
        ev.radiation_coupling = coupling;
        ev.radiation_coupling_normalized = coupling/Cnorm;
        ev.open_non_target_orders = openOrders;
        ev.finite_channel_numerator = coupling(finiteNonTarget);
        ev.finite_channel_numerator_residual = finiteNumerator;
        ev.finite_channel_power_residual = finitePowerNumerator;
        ev.external_finite_amplitude = externalFiniteAmplitude;
        ev.channel_equation_consistency = consistency;
        ev.target_raw_coupling = targetRaw;
        ev.target_external_amplitude = targetAmplitude;
        ev.target_power_coupling = targetPower;
        ev.groove_fraction = norm(C)^2/max(norm(v)^2,eps);
        ev.finite_external_fraction = sum(abs(A(finiteNonTarget)).^2)/ ...
            max(norm(v)^2,eps);
        ev.valid = true;
    end

    function ev = invalid_point(kap,Om,N,K)
        ev = struct('kappa',kap,'Omega',Om,'rayleigh_residual',nan, ...
            'N',N,'K',K,'orders',[],'kx',[],'ky',[], ...
            'pole_denominator',nan,'sigma_min',nan,'sigma_max',nan, ...
            'mode',[],'floquet_amplitudes',[],'groove_amplitudes',[], ...
            'radiation_map',[],'radiation_coupling',[], ...
            'radiation_coupling_normalized',[],'open_non_target_orders',[], ...
            'finite_channel_numerator',[],'finite_channel_numerator_residual',nan, ...
            'finite_channel_power_residual',nan,'external_finite_amplitude',nan, ...
            'channel_equation_consistency',nan,'target_raw_coupling',nan, ...
            'target_external_amplitude',nan,'target_power_coupling',nan, ...
            'groove_fraction',nan,'finite_external_fraction',nan,'valid',false);
    end

    function [value,id] = finite_minimum(x,valid)
        mask = valid & isfinite(x);
        if ~any(mask(:))
            value = nan; id = [];
            return;
        end
        y = x; y(~mask) = inf;
        [value,id] = min(y(:));
    end

    function summary = convergence_summary(entries,label)
        summary = struct('label',label,'values',nan(numel(entries),1), ...
            'numerator_values',nan(numel(entries),1),'stable',false, ...
            'last_ratio_pole',nan,'last_ratio_numerator',nan);
        for mm = 1:numel(entries)
            summary.values(mm) = entries{mm}.pole_denominator;
            summary.numerator_values(mm) = ...
                entries{mm}.finite_channel_numerator_residual;
        end
        if numel(entries) >= 2
            summary.last_ratio_pole = safe_ratio(summary.values(end), ...
                summary.values(end-1));
            summary.last_ratio_numerator = safe_ratio( ...
                summary.numerator_values(end),summary.numerator_values(end-1));
            summary.stable = ratio_is_reasonable(summary.last_ratio_pole) && ...
                ratio_is_reasonable(summary.last_ratio_numerator);
        end
        summary.stable = summary.stable && all(isfinite(summary.values)) && ...
            all(isfinite(summary.numerator_values));
    end

    function r = safe_ratio(x,y)
        if y == 0
            if x == 0, r = 1; else, r = inf; end
        else
            r = x/y;
        end
    end

    function tf = ratio_is_reasonable(r)
        tf = isfinite(r) && r >= .1 && r <= 10;
    end

    function print_summary(res)
        fprintf('Radiation-zero diagnostic\n');
        fprintf('  kappa = %.12f, Omega = %.12f, target order = %d\n', ...
            res.kappa,res.Omega,res.rayleigh_order);
        fprintf('  Rayleigh residual Omega-|kappa+n| = %.3e\n', ...
            res.rayleigh_residual);
        fprintf('  pole denominator = %.3e\n',res.center.pole_denominator);
        fprintf('  finite-channel numerator residual = %.3e\n', ...
            res.center.finite_channel_numerator_residual);
        fprintf('  target raw coupling (diagnostic only) = %.3e\n', ...
            res.center.target_raw_coupling);
        fprintf('  pole minimum at [kappa Omega] = [%.8f %.8f]\n', ...
            res.pole_minimum.coordinate);
        fprintf('  radiation minimum at [kappa Omega] = [%.8f %.8f]\n', ...
            res.radiation_minimum.coordinate);
        fprintf('  N/K convergence stable = %d\n',res.convergence.stable);
        fprintf('  verdict: %s\n',res.verdict);
    end
end
