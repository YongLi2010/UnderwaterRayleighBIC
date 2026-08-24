function diagnostic = ni2019_two_groove_parity_diagnostic(inputData,varargin)
%NI2019_TWO_GROOVE_PARITY_DIAGNOSTIC Diagnose two-groove Gamma parity.
%   D = NI2019_TWO_GROOVE_PARITY_DIAGNOSTIC(R) accepts the output R of
%   NI2019_STRICT_RAYLEIGH_OPERATOR. It also accepts a struct containing
%   MODE and OPERATOR fields, or a compatible operator struct alone. The
%   function is read-only and does not write files or alter its input.
%
%   For two identical grooves mirrored by x -> a-x, the local cosine basis
%   gives the exact reflection law
%
%       C(q,2) = p*(-1)^q C(q,1),       p=+1 (even), p=-1 (odd),
%       A(-n)  = p A(n)                  at Gamma.
%
%   Consequently an odd mode has A(0)=0 automatically, but its two first
%   Rayleigh channels are related by A(+1)=-A(-1), leaving one independent
%   threshold radiation amplitude that must still be zero for a strict
%   double-Rayleigh BIC. The diagnostic also evaluates q-wise aperture
%   source projections. A q=0-only odd mode has one scalar amplitude; it
%   cannot cancel a nonzero q=0 projection against itself. It can be
%   non-radiating only at a special geometric projection zero or with its
%   q=0 coefficient equal to zero. Additional q components can interfere
%   with it and cancel the threshold projection.
%
%   Name/value options:
%       SymmetryTolerance  Relative geometry tolerance (default 1e-8)
%       ParityTolerance    Normalized modal residual tolerance (default 1e-6)
%       Verbose            Print a compact diagnostic (default false)

p = inputParser;
addParameter(p,'SymmetryTolerance',1e-8, ...
    @(x)isscalar(x) && isfinite(x) && x>=0);
addParameter(p,'ParityTolerance',1e-6, ...
    @(x)isscalar(x) && isfinite(x) && x>=0);
addParameter(p,'Verbose',false,@(x)islogical(x) && isscalar(x));
parse(p,varargin{:});
opts = p.Results;

[mode,op] = unpack_input(inputData);
if isempty(op)
    error('The input must contain a compatible operator or full_operator struct.');
end
if ~isfield(op,'L') || op.L~=2
    error('This diagnostic requires exactly two grooves (op.L==2).');
end
required = {'N','K','orders','widths','depths','gaps','xleft','a','k0'};
for j = 1:numel(required)
    if ~isfield(op,required{j})
        error('The operator is missing metadata field %s.',required{j});
    end
end

[widths,depths,gaps,xleft,a] = geometry_from_operator(op);
geometry = geometry_symmetry(widths,depths,gaps,xleft,a,opts.SymmetryTolerance);

if isempty(mode)
    mode = smallest_mode_from_operator(op);
end
[A,Cphysical] = coefficients_from_mode(mode,op);
if numel(A)~=op.N || numel(Cphysical)~=op.L*op.K
    error('Mode coefficient dimensions do not match the operator metadata.');
end
C = reshape(Cphysical,op.K,op.L);

% Reflection of cos(q*pi*xi/w) under xi -> w-xi contributes (-1)^q.
localReflectionSign = reshape((-1).^(0:op.K-1),[],1);
Cref = localReflectionSign.*C(:,[2 1]);
coeffScale = max(norm(C(:)),eps);
evenCoefficientResidual = norm(C(:)-Cref(:))/coeffScale;
oddCoefficientResidual = norm(C(:)+Cref(:))/coeffScale;
qEvenResidual = abs(C(:,2)-localReflectionSign.*C(:,1))/coeffScale;
qOddResidual = abs(C(:,2)+localReflectionSign.*C(:,1))/coeffScale;

orders = op.orders(:);
minusId = orders==-1;
zeroId = orders==0;
plusId = orders==1;
if ~any(minusId) || ~any(zeroId) || ~any(plusId)
    error('The retained Floquet orders must include -1, 0, and +1.');
end
Aminus = A(minusId);
Azero = A(zeroId);
Aplus = A(plusId);
channelScale = max(norm(A),eps);
evenChannelResidual = abs(Aplus-Aminus)/channelScale;
oddChannelResidual = abs(Aplus+Aminus)/channelScale;
zeroChannelAmplitude = abs(Azero)/channelScale;
evenCombination = (Aplus+Aminus)/sqrt(2);
oddCombination = (Aplus-Aminus)/sqrt(2);

% q-wise source projections at n=-1 and n=+1. The source factor is the
% physical groove normal-velocity coefficient beta/k0*sin(beta*d), and the
% projection matches the Cmap convention of the modal operator.
[qSourceCoupling,sourceByQ] = threshold_source_couplings(C,op);
qSourceScale = max(norm(qSourceCoupling(:)),eps);
qSourceNormByQ = sqrt(sum(abs(qSourceCoupling).^2,2));
qSourceNormByQNormalized = qSourceNormByQ/qSourceScale;
qEnergyProxy = zeros(op.K,1);
for q = 0:op.K-1
    xNorm = widths(1)*(q==0) + 0.5*widths(1)*(q>0);
    qEnergyProxy(q+1) = xNorm*sum(abs(C(q+1,:)).^2);
end

% A q=0 odd pair contains one scalar coefficient. Its projected threshold
% coupling is a geometric scalar; a nonzero value cannot be self-cancelled.
q0Coupling = qSourceCoupling(1,:).';
q0CouplingNorm = norm(q0Coupling);
q0BasisCoupling = q0_odd_basis_coupling(op);
q0BasisCouplingNorm = norm(q0BasisCoupling);
q0GeometryZero = q0BasisCouplingNorm <= opts.ParityTolerance;

if geometry.is_symmetric
    if evenCoefficientResidual<=opts.ParityTolerance && ...
            evenChannelResidual<=opts.ParityTolerance
        parityLabel = 'even';
        paritySign = 1;
    elseif oddCoefficientResidual<=opts.ParityTolerance && ...
            oddChannelResidual<=opts.ParityTolerance && ...
            zeroChannelAmplitude<=opts.ParityTolerance
        parityLabel = 'odd';
        paritySign = -1;
    elseif evenCoefficientResidual<oddCoefficientResidual && ...
            evenChannelResidual<=oddChannelResidual
        parityLabel = 'even-like / residual above tolerance';
        paritySign = 1;
    elseif oddCoefficientResidual<evenCoefficientResidual && ...
            oddChannelResidual<=evenChannelResidual
        parityLabel = 'odd-like / residual above tolerance';
        paritySign = -1;
    else
        parityLabel = 'mixed or undetermined';
        paritySign = NaN;
    end
else
    parityLabel = 'asymmetric geometry: parity not assigned';
    paritySign = NaN;
end

grooveNorms = sqrt(sum(abs(C).^2,1));
grooveNormsNormalized = grooveNorms/max(norm(grooveNorms),eps);

diagnostic = struct();
diagnostic.geometry = geometry;
diagnostic.geometry_mirror_error = geometry.mirror_error;
diagnostic.is_symmetric = geometry.is_symmetric;
diagnostic.parity_label = parityLabel;
diagnostic.parity_sign = paritySign;
diagnostic.parity_tolerance = opts.ParityTolerance;
diagnostic.coefficient_residuals = struct('even',evenCoefficientResidual, ...
    'odd',oddCoefficientResidual,'normalization',coeffScale, ...
    'q_even',qEvenResidual,'q_odd',qOddResidual);
diagnostic.channel_residuals = struct( ...
    'even_Aplus_minus_Aminus',evenChannelResidual, ...
    'odd_Aplus_plus_Aminus',oddChannelResidual, ...
    'A0',zeroChannelAmplitude,'normalization',channelScale);
diagnostic.channel_combination_amplitudes = struct( ...
    'A_minus1',Aminus,'A_0',Azero,'A_plus1',Aplus, ...
    'even_Aplus_plus_Aminus',evenCombination, ...
    'odd_Aplus_minus_Aminus',oddCombination, ...
    'odd_independent_threshold_amplitude',Aminus);
diagnostic.coefficients = struct('A',A,'C_physical',Cphysical, ...
    'C_by_groove',C,'C_reflected',Cref, ...
    'local_reflection_sign',localReflectionSign);
diagnostic.q_wise = struct('coefficient_norm',sqrt(sum(abs(C).^2,2)), ...
    'coefficient_energy_proxy',qEnergyProxy, ...
    'source_coupling_n_minus1_plus1',qSourceCoupling, ...
    'source_coupling_norm',qSourceNormByQ, ...
    'source_coupling_norm_normalized',qSourceNormByQNormalized, ...
    'source_by_q',sourceByQ);
diagnostic.q0_only = struct('threshold_coupling_minus1_plus1',q0Coupling, ...
    'coupling_norm',q0CouplingNorm, ...
    'unit_odd_q0_threshold_coupling_minus1_plus1',q0BasisCoupling, ...
    'unit_odd_q0_coupling_norm',q0BasisCouplingNorm, ...
    'geometric_projection_zero',q0GeometryZero, ...
    'nontrivial_q0_can_be_nonradiating',q0GeometryZero, ...
    'interpretation',q0_interpretation(q0GeometryZero));
diagnostic.groove_by_groove_norms = grooveNorms;
diagnostic.groove_by_groove_norms_normalized = grooveNormsNormalized;
diagnostic.orders = orders;
diagnostic.operator_dimensions = struct('N',op.N,'K',op.K,'L',op.L);
diagnostic.theory = struct( ...
    'mirror_law','C(q,2)=p*(-1)^q*C(q,1)', ...
    'channel_law','A(-n)=p*A(n) at Gamma', ...
    'odd_A0','p=-1 gives A(0)=-A(0), hence A(0)=0', ...
    'odd_threshold','p=-1 gives A(+1)=-A(-1), leaving one independent threshold amplitude', ...
    'q0_conclusion',q0_interpretation(q0GeometryZero), ...
    'q_ge1_role','q>=1 supplies additional internal coefficients whose threshold projections can interfere with q=0');

if opts.Verbose
    fprintf('two-groove Gamma parity: %s, mirror error %.3e\n', ...
        parityLabel,geometry.mirror_error);
    fprintf('  coefficient residuals: even %.3e, odd %.3e\n', ...
        evenCoefficientResidual,oddCoefficientResidual);
    fprintf('  channel residuals: even %.3e, odd %.3e, A0 %.3e\n', ...
        evenChannelResidual,oddChannelResidual,zeroChannelAmplitude);
    fprintf('  odd independent |A_-1|=%.3e, unit-q0 coupling norm=%.3e (%s)\n', ...
        abs(Aminus),q0BasisCouplingNorm,q0_interpretation(q0GeometryZero));
    fprintf('  groove norms: [%s]\n',strtrim(sprintf('%.3e ',grooveNorms)));
end
end

function [mode,op] = unpack_input(inputData)
mode = [];
op = [];
if ~isstruct(inputData)
    error('Input must be a struct returned by the strict operator or contain mode/operator fields.');
end
if isfield(inputData,'full_operator')
    op = inputData.full_operator;
elseif isfield(inputData,'operator')
    op = inputData.operator;
elseif isfield(inputData,'F') && isfield(inputData,'N')
    op = inputData;
end
if isfield(inputData,'mode')
    mode = inputData.mode;
elseif isfield(inputData,'A')
    mode = inputData;
end
end

function [widths,depths,gaps,xleft,a] = geometry_from_operator(op)
widths = op.widths(:).';
depths = op.depths(:).';
gaps = op.gaps(:).';
xleft = op.xleft(:).';
a = op.a;
if numel(widths)~=2 || numel(depths)~=2 || numel(xleft)~=2 || numel(gaps)~=1
    error('Two-groove geometry requires two widths/depths/xleft values and one gap.');
end
end

function g = geometry_symmetry(widths,depths,gaps,xleft,a,tol)
scale = max([abs(a),abs(widths),abs(depths),abs(gaps),1]);
widthError = abs(widths(1)-widths(2))/scale;
depthError = abs(depths(1)-depths(2))/scale;
marginError = abs(xleft(1)-(a-(xleft(2)+widths(2))))/scale;
centerError = abs((xleft(1)+widths(1)+xleft(2)) - ...
    (a))/scale;
% The center relation above is equivalent to xleft(1)+widths(1)=a-xleft(2)
% for equal widths; retaining both terms makes the reported error explicit.
centerError = max(centerError, ...
    abs((xleft(1)+widths(1))-(a-xleft(2)))/scale);
err = max([widthError,depthError,marginError,centerError]);
g = struct('mirror_error',err,'symmetry_error',err, ...
    'is_symmetric',err<=tol,'tolerance',tol, ...
    'width_error',widthError,'depth_error',depthError, ...
    'margin_error',marginError,'center_error',centerError);
end

function mode = smallest_mode_from_operator(op)
if isfield(op,'Fscaled') && isfield(op,'column_scale')
    [~,~,V] = svd(op.Fscaled,'econ');
    z = V(:,end)./op.column_scale(:);
elseif isfield(op,'F')
    [~,~,V] = svd(op.F,'econ');
    z = V(:,end);
else
    error('No mode coefficients are available and the operator has no F matrix.');
end
mode = struct('A',z(1:op.N),'C_scaled',z(op.N+1:end));
end

function [A,Cphysical] = coefficients_from_mode(mode,op)
if isfield(mode,'A')
    A = mode.A(:);
else
    error('The mode is missing Floquet amplitudes A.');
end
if isfield(mode,'C_physical')
    Cphysical = mode.C_physical(:);
elseif isfield(mode,'C_by_groove')
    Cphysical = mode.C_by_groove(:);
elseif isfield(mode,'C_scaled')
    Cphysical = physical_from_scaled(mode.C_scaled(:),op);
else
    error('The mode is missing C_physical, C_by_groove, or C_scaled.');
end
end

function Cphysical = physical_from_scaled(Cscaled,op)
scale = zeros(op.n_groove,1);
row = 0;
for ell = 1:op.L
    d = op.depths(ell);
    t = op.widths(ell);
    for q = 0:op.K-1
        row = row+1;
        alpha = q*pi/t;
        beta = groove_sqrt_local(op.k0^2-alpha^2);
        scale(row) = max([abs(cos(beta*d)), ...
            abs((beta/op.k0)*sin(beta*d)),sqrt(eps)]);
    end
end
Cphysical = Cscaled./scale;
end

function [coupling,sourceByQ] = threshold_source_couplings(C,op)
coupling = complex(zeros(op.K,2));
sourceByQ = complex(zeros(op.K,2,2));
orders = [-1,1];
for iq = 1:op.K
    q = iq-1;
    for ell = 1:op.L
        t = op.widths(ell);
        d = op.depths(ell);
        alpha = q*pi/t;
        beta = groove_sqrt_local(op.k0^2-alpha^2);
        source = C(iq,ell)*(beta/op.k0)*sin(beta*d);
        for in = 1:2
            n = orders(in);
            kx = op.kx(op.orders==n);
            proj = projection_local(-kx,q,op.xleft(ell),t);
            contribution = (t/op.a)*proj*source;
            sourceByQ(iq,in,ell) = contribution;
            coupling(iq,in) = coupling(iq,in)+contribution;
        end
    end
end
end

function coupling = q0_odd_basis_coupling(op)
% Coupling of a unit physical q=0 odd pair [C_1,C_2]=[1,-1].
coupling = complex(zeros(2,1));
for in = 1:2
    n = [-1,1];
    kx = op.kx(op.orders==n(in));
    for ell = 1:2
        t = op.widths(ell);
        d = op.depths(ell);
        beta = groove_sqrt_local(op.k0^2);
        sourceSign = 1-2*(ell==2);
        source = sourceSign*(beta/op.k0)*sin(beta*d);
        coupling(in) = coupling(in)+(t/op.a)* ...
            projection_local(-kx,0,op.xleft(ell),t)*source;
    end
end
end

function value = projection_local(kx,q,xl,t)
alpha = q*pi/t;
value = exp(-1i*kx*xl)*0.5*(expint_local(kx-alpha,t)+ ...
    expint_local(kx+alpha,t));
end

function value = expint_local(kappa,t)
z = kappa*t;
if abs(z)<1e-8
    value = 1-1i*z/2-z^2/6;
else
    value = exp(-1i*z/2)*sin(z/2)/(z/2);
end
end

function text = q0_interpretation(isZero)
if isZero
    text = 'q=0 projection is geometrically zero; a nontrivial q=0 null is possible in principle';
else
    text = 'q=0 projection is nonzero; q=0 alone cannot self-cancel the threshold channel';
end
end

function value = groove_sqrt_local(z)
if real(z)>=0
    value = sqrt(real(z));
else
    value = 1i*sqrt(-real(z));
end
end
