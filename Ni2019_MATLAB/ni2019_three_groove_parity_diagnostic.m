function diagnostic = ni2019_three_groove_parity_diagnostic(inputData,varargin)
%NI2019_THREE_GROOVE_PARITY_DIAGNOSTIC Diagnose Gamma-point modal parity.
%   D = NI2019_THREE_GROOVE_PARITY_DIAGNOSTIC(R) accepts the output R of
%   NI2019_STRICT_RAYLEIGH_OPERATOR. It also accepts a struct containing
%   fields MODE and OPERATOR, or an operator struct with a compatible mode.
%   The function is read-only and does not write files or modify its input.
%
%   For a three-groove structure symmetric under x -> a-x, grooves 1 and 3
%   are exchanged and groove 2 maps to itself. With the local cosine basis
%       cos(q*pi*xi/w),
%   reflection contributes (-1)^q. Thus a parity-p mode obeys
%       C(q,mirror(ell)) = p*(-1)^q*C(q,ell),
%       A(-n) = p*A(n),
%   where p=+1 is even and p=-1 is odd. At Gamma the channel diagnostics
%   therefore test A(+1)-A(-1) (even), A(+1)+A(-1) (odd), and A(0).
%
%   If the geometry is not symmetric within SymmetryTolerance, no parity
%   label is assigned. This avoids interpreting an asymmetric mode as an
%   even/odd mode merely because one residual happens to be smaller.
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

if ~isfield(op,'L') || op.L~=3
    error('This diagnostic requires exactly three grooves (op.L==3).');
end
if ~isfield(op,'N') || ~isfield(op,'K') || ~isfield(op,'orders')
    error('The operator is missing N, K, or orders metadata.');
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

% Reflection of the local cosine basis: q=0,1,...,K-1 has sign (-1)^q.
localReflectionSign = reshape((-1).^(0:op.K-1),[],1);
Cref = localReflectionSign.*C(:,[3 2 1]);
coeffNorm = norm(C(:));
coeffScale = max(coeffNorm,eps);
evenCoefficientResidual = norm(C(:)-Cref(:))/coeffScale;
oddCoefficientResidual = norm(C(:)+Cref(:))/coeffScale;

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

if geometry.is_symmetric
    if evenCoefficientResidual <= opts.ParityTolerance && ...
            evenChannelResidual <= opts.ParityTolerance
        parityLabel = 'even';
        paritySign = 1;
    elseif oddCoefficientResidual <= opts.ParityTolerance && ...
            oddChannelResidual <= opts.ParityTolerance && ...
            zeroChannelAmplitude <= opts.ParityTolerance
        parityLabel = 'odd';
        paritySign = -1;
    elseif evenCoefficientResidual < oddCoefficientResidual && ...
            evenChannelResidual <= oddChannelResidual
        parityLabel = 'even-like / residual above tolerance';
        paritySign = 1;
    elseif oddCoefficientResidual < evenCoefficientResidual && ...
            oddChannelResidual <= evenChannelResidual
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
diagnostic.geometry_symmetry_error = geometry.symmetry_error;
diagnostic.is_symmetric = geometry.is_symmetric;
diagnostic.parity_label = parityLabel;
diagnostic.parity_sign = paritySign;
diagnostic.parity_tolerance = opts.ParityTolerance;
diagnostic.coefficient_residuals = struct( ...
    'even',evenCoefficientResidual,'odd',oddCoefficientResidual, ...
    'normalization',coeffScale);
diagnostic.channel_residuals = struct( ...
    'even_Aplus_minus_Aminus',evenChannelResidual, ...
    'odd_Aplus_plus_Aminus',oddChannelResidual, ...
    'A0',zeroChannelAmplitude,'normalization',channelScale);
diagnostic.channel_combination_amplitudes = struct( ...
    'A_minus1',Aminus,'A_0',Azero,'A_plus1',Aplus, ...
    'even_Aplus_plus_Aminus',evenCombination, ...
    'odd_Aplus_minus_Aminus',oddCombination);
diagnostic.groove_by_groove_norms = grooveNorms;
diagnostic.groove_by_groove_norms_normalized = grooveNormsNormalized;
diagnostic.coefficients = struct('A',A,'C_physical',Cphysical, ...
    'C_by_groove',C,'C_reflected',Cref, ...
    'local_reflection_sign',localReflectionSign);
diagnostic.orders = orders;
diagnostic.operator_dimensions = struct('N',op.N,'K',op.K,'L',op.L);

if opts.Verbose
    fprintf('three-groove parity: %s, geometry error %.3e\n', ...
        parityLabel,geometry.symmetry_error);
    fprintf('  coefficient residuals: even %.3e, odd %.3e\n', ...
        evenCoefficientResidual,oddCoefficientResidual);
    fprintf('  channel residuals: even %.3e, odd %.3e, A0 %.3e\n', ...
        evenChannelResidual,oddChannelResidual,zeroChannelAmplitude);
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
required = {'widths','depths','gaps','xleft','a'};
for j = 1:numel(required)
    if ~isfield(op,required{j})
        error('The operator is missing geometry field %s.',required{j});
    end
end
widths = op.widths(:).';
depths = op.depths(:).';
gaps = op.gaps(:).';
xleft = op.xleft(:).';
a = op.a;
if numel(widths)~=3 || numel(depths)~=3 || numel(xleft)~=3 || numel(gaps)~=2
    error('Three-groove geometry requires three widths/depths/xleft values and two gaps.');
end
end

function g = geometry_symmetry(widths,depths,gaps,xleft,a,tol)
scale = max([abs(a),abs(widths),abs(depths),abs(gaps),1]);
pairWidth = abs(widths(1)-widths(3))/scale;
pairDepth = abs(depths(1)-depths(3))/scale;
gapSymmetry = abs(gaps(1)-gaps(2))/scale;
leftMargin = xleft(1);
rightMargin = a-(xleft(3)+widths(3));
marginSymmetry = abs(leftMargin-rightMargin)/scale;
outerReflection = abs((xleft(1)+widths(1))-(a-xleft(3)))/scale;
centerReflection = abs((2*xleft(2)+widths(2))-a)/scale;
err = max([pairWidth,pairDepth,gapSymmetry,marginSymmetry, ...
    outerReflection,centerReflection]);
g = struct('symmetry_error',err,'is_symmetric',err<=tol, ...
    'tolerance',tol,'width_pair_error',pairWidth, ...
    'depth_pair_error',pairDepth,'gap_error',gapSymmetry, ...
    'margin_error',marginSymmetry,'outer_reflection_error',outerReflection, ...
    'center_reflection_error',centerReflection);
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
if isfield(op,'N')
    mode = struct('A',z(1:op.N),'C_scaled',z(op.N+1:end));
else
    error('The operator is missing N metadata.');
end
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

function value = groove_sqrt_local(z)
if real(z)>=0
    value = sqrt(real(z));
else
    value = 1i*sqrt(-real(z));
end
end
