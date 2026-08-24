function result=ni2019_two_groove_odd_strict_operator(cfg,varargin)
%NI2019_TWO_GROOVE_ODD_STRICT_OPERATOR Strict odd double-Rayleigh system.
%   R=NI2019_TWO_GROOVE_ODD_STRICT_OPERATOR(CFG) restricts the homogeneous
%   two-groove modal-matching operator to the mirror-odd subspace at
%   kappa=0 and Omega=a/lambda=1.  The two grooves must be identical and
%   centered as a mirror pair.  The basis imposes
%
%       A(-n)=-A(n),                  n>=2,
%       A(-1)=A(0)=A(+1)=0,
%       C(q,2)=-(-1)^q C(q,1).
%
%   Hence a vanishing singular value is precisely the odd-sector
%   compatibility condition for a Gamma-point double-Rayleigh BIC.
%
%   Name/value options:
%       GeometryTolerance  normalized mirror tolerance, default 1e-10
%       NormalizeMode      normalize the reconstructed full mode, true
%       Verbose            print diagnostics, false

p=inputParser;
addParameter(p,'GeometryTolerance',1e-10, ...
    @(x)isscalar(x)&&isfinite(x)&&x>=0);
addParameter(p,'NormalizeMode',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'Verbose',false,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});
opt=p.Results;

if ~isfield(cfg,'N'), cfg.N=101; end
if ~isfield(cfg,'K'), cfg.K=10; end
if ~isfield(cfg,'a')||~isfield(cfg,'lambda')||~isfield(cfg,'theta_i_deg')
    error('cfg must provide a, lambda, and theta_i_deg.');
end
if mod(cfg.N,2)~=1||cfg.N<5||cfg.K<1
    error('The odd strict operator requires odd N>=5 and K>=1.');
end
if numel(cfg.widths)~=2||numel(cfg.depths)~=2||numel(cfg.gaps)~=1
    error('cfg must describe exactly two grooves and one internal gap.');
end

a=cfg.a;
tol=opt.GeometryTolerance;
normalizedError=max([abs(cfg.widths(1)-cfg.widths(2)), ...
    abs(cfg.depths(1)-cfg.depths(2))])/a;
occupied=sum(cfg.widths)+cfg.gaps;
centeredX0=.5*(a-occupied);
if isfield(cfg,'x0')
    normalizedError=max(normalizedError,abs(cfg.x0-centeredX0)/a);
end
normalizedError=max(normalizedError,abs(cfg.lambda-a)/a);
normalizedError=max(normalizedError,abs(sind(cfg.theta_i_deg)));
if normalizedError>tol
    error(['Geometry/frequency is not the required identical centered pair ', ...
        'at kappa=0, Omega=1 (normalized error %.3e).'],normalizedError);
end

op=ni2019_full_eigen_operator(cfg);
kyTolerance=1e-8*max(abs(op.k0),1);
finiteOpen=abs(imag(op.ky))<=kyTolerance&real(op.ky)>kyTolerance;
threshold=abs(op.ky)<=kyTolerance;
removedOrders=op.orders(finiteOpen|threshold);
if ~isequal(removedOrders(:),[-1;0;1])
    error('Expected removed orders [-1 0 1], obtained [%s].', ...
        strtrim(sprintf('%d ',removedOrders)));
end

positiveOrders=(2:(op.N-1)/2).';
nOddA=numel(positiveOrders);
T=complex(zeros(op.N+op.n_groove,nOddA+op.K));
for jj=1:nOddA
    plusId=find(op.orders==positiveOrders(jj),1);
    minusId=find(op.orders==-positiveOrders(jj),1);
    T(plusId,jj)=1/sqrt(2);
    T(minusId,jj)=-1/sqrt(2);
end
for q=0:op.K-1
    column=nOddA+q+1;
    T(op.N+q+1,column)=1/sqrt(2);
    T(op.N+op.K+q+1,column)=-(-1)^q/sqrt(2);
end

% Keep one equation from each mirror-related row pair.  The left-groove
% pressure equations and the positive-order velocity equations form the
% complete independent odd sector.  In particular, the n=+1 row is the
% one remaining Rayleigh compatibility equation.
FoddFull=op.F*T;
pressureRows=(1:op.K).';
velocityRows=2*op.K+find(op.orders>=1);
independentRows=[pressureRows;velocityRows];
Fodd=FoddFull(independentRows,:);
rowScale=max(vecnorm(Fodd,2,2),sqrt(eps));
Frow=Fodd./rowScale;
columnScale=max(vecnorm(Frow,2,1),sqrt(eps));
Fscaled=Frow./columnScale;
if any(~isfinite(Fscaled(:)))
    error('Odd reduced operator contains NaN or Inf.');
end
[~,S,V]=svd(Fscaled,'econ');
singularValues=diag(S);
v=V(:,end);
uOdd=v./columnScale(:);
z=T*uOdd;
if opt.NormalizeMode
    z=z/max(norm(z),eps);
    uOdd=uOdd/max(norm(T*uOdd),eps);
end

A=z(1:op.N);
Cscaled=z(op.N+1:end);
Cphysical=Cscaled./op.vertical_scale;
CbyGroove=reshape(Cphysical,op.K,2);
surfacePressure=Cscaled.*op.cos_depth_normalized;
surfaceVelocity=Cscaled.*op.beta_sin_normalized;
pressureByGroove=reshape(surfacePressure,op.K,2);
velocityByGroove=reshape(surfaceVelocity,op.K,2);

qNorm=sqrt(sum(abs(CbyGroove).^2,2));
pressureQNorm=sqrt(sum(abs(pressureByGroove).^2,2));
velocityQNorm=sqrt(sum(abs(velocityByGroove).^2,2));
qFraction=normalized_square(qNorm);
pressureQFraction=normalized_square(pressureQNorm);
velocityQFraction=normalized_square(velocityQNorm);

reflectedC=(-1).^(0:op.K-1).'.*CbyGroove(:,[2 1]);
oddCoefficientResidual=norm(CbyGroove+reflectedC,'fro')/ ...
    max(norm(CbyGroove,'fro'),eps);
evenCoefficientResidual=norm(CbyGroove-reflectedC,'fro')/ ...
    max(norm(CbyGroove,'fro'),eps);

rawResidual=norm(op.F*z)/max(norm(z),eps);
scaledResidual=norm(Fscaled*v)/max(norm(v),eps);
sigmaRatio=singularValues(end)/max(singularValues(1),eps);
removedAmplitude=A(finiteOpen|threshold);
propagatingPower=max(real(op.ky),0)./max(abs(op.k0),eps).*abs(A).^2;

result=struct();
result.cfg=cfg;
result.full_operator=op;
result.T_odd=T;
result.F_odd_full_raw=FoddFull;
result.F_odd_raw=Fodd;
result.F_odd_scaled=Fscaled;
result.row_scale=rowScale;
result.column_scale=columnScale;
result.singular_values=singularValues;
result.sigma_ratio=sigmaRatio;
result.strict_residual=rawResidual;
result.residual=struct('scaled',scaledResidual,'raw',rawResidual);
result.positive_orders=positiveOrders;
result.independent_rows=independentRows;
result.removed_orders=removedOrders;
result.mode=struct('z_full',z,'u_odd',uOdd,'v_scaled',v,'A',A, ...
    'C_scaled',Cscaled,'C_physical',Cphysical, ...
    'C_by_groove',CbyGroove, ...
    'surface_pressure_coefficients',surfacePressure, ...
    'surface_velocity_coefficients',surfaceVelocity);
result.transverse=struct('q',(0:op.K-1).','norm_by_q',qNorm, ...
    'fraction_by_q',qFraction, ...
    'surface_pressure_norm_by_q',pressureQNorm, ...
    'surface_pressure_fraction_by_q',pressureQFraction, ...
    'surface_velocity_norm_by_q',velocityQNorm, ...
    'surface_velocity_fraction_by_q',velocityQFraction, ...
    'bottom_higher_order_fraction',sum(qFraction(2:end)), ...
    'surface_pressure_higher_order_fraction', ...
    sum(pressureQFraction(2:end)), ...
    'surface_velocity_higher_order_fraction', ...
    sum(velocityQFraction(2:end)));
result.parity=struct('label','odd','geometry_error',normalizedError, ...
    'odd_coefficient_residual',oddCoefficientResidual, ...
    'even_coefficient_residual',evenCoefficientResidual);
result.radiation=struct('amplitudes',A,'orders',op.orders, ...
    'removed_amplitudes',removedAmplitude, ...
    'removed_amplitudes_exact_zero',all(removedAmplitude==0), ...
    'power_by_order',propagatingPower, ...
    'total_propagating_power',sum(propagatingPower), ...
    'evanescent_amplitude_norm',norm(A(~(finiteOpen|threshold))));
result.dimensions=struct('full',size(op.F), ...
    'odd_full_rows',size(FoddFull),'odd',size(Fodd), ...
    'odd_scaled',size(Fscaled));

if opt.Verbose
    fprintf(['two-groove odd strict operator: size %d x %d, ', ...
        'sigma ratio %.3e, raw residual %.3e\n'], ...
        size(Fodd,1),size(Fodd,2),sigmaRatio,rawResidual);
    fprintf(['  q>=1 aperture-pressure %.3f, aperture-velocity %.3f, ', ...
        'odd residual %.3e\n'], ...
        sum(pressureQFraction(2:end)),sum(velocityQFraction(2:end)), ...
        oddCoefficientResidual);
end
end

function fraction=normalized_square(value)
fraction=abs(value).^2/max(sum(abs(value).^2),eps);
end
