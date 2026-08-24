function result=ni2019_three_groove_parity_strict_operator(cfg,parity,varargin)
%NI2019_THREE_GROOVE_PARITY_STRICT_OPERATOR Gamma strict parity sector.
%   R=NI2019_THREE_GROOVE_PARITY_STRICT_OPERATOR(CFG,PARITY) constructs the
%   independent mirror-even or mirror-odd strict system for three centered
%   mirror-symmetric grooves at kappa=0, Omega=1. PARITY is 'odd', 'even',
%   -1, or +1. The basis imposes
%
%       A(-n)=p*A(n),
%       C(q,3)=p*(-1)^q*C(q,1),
%
%   and retains a central-groove coefficient only when
%
%       p*(-1)^q=1.
%
%   A(-1), A(0), and A(+1) are identically zero. The independent equation
%   set contains the left-groove pressure rows, the parity-allowed central
%   pressure rows, and positive-order velocity rows; the even sector also
%   retains the n=0 velocity compatibility row.

pIn=inputParser;
addParameter(pIn,'GeometryTolerance',1e-10, ...
    @(x)isscalar(x)&&isfinite(x)&&x>=0);
addParameter(pIn,'NormalizeMode',true,@(x)islogical(x)&&isscalar(x));
addParameter(pIn,'Verbose',false,@(x)islogical(x)&&isscalar(x));
parse(pIn,varargin{:});
opt=pIn.Results;

paritySign=parse_parity(parity);
parityName='even';
if paritySign<0, parityName='odd'; end
if ~isfield(cfg,'N'), cfg.N=101; end
if ~isfield(cfg,'K'), cfg.K=10; end
if ~isfield(cfg,'a')||~isfield(cfg,'lambda')||~isfield(cfg,'theta_i_deg')
    error('cfg must provide a, lambda, and theta_i_deg.');
end
if mod(cfg.N,2)~=1||cfg.N<5||cfg.K<1
    error('The parity strict operator requires odd N>=5 and K>=1.');
end
if numel(cfg.widths)~=3||numel(cfg.depths)~=3||numel(cfg.gaps)~=2
    error('cfg must describe exactly three grooves and two internal gaps.');
end

a=cfg.a;
occupied=sum(cfg.widths)+sum(cfg.gaps);
centeredX0=.5*(a-occupied);
geometryError=max([abs(cfg.widths(1)-cfg.widths(3)), ...
    abs(cfg.depths(1)-cfg.depths(3)),abs(cfg.gaps(1)-cfg.gaps(2))])/a;
if isfield(cfg,'x0')
    geometryError=max(geometryError,abs(cfg.x0-centeredX0)/a);
end
geometryError=max(geometryError,abs(cfg.lambda-a)/a);
geometryError=max(geometryError,abs(sind(cfg.theta_i_deg)));
if geometryError>opt.GeometryTolerance
    error(['Geometry/frequency is not a centered mirror-symmetric ', ...
        'three-groove cell at kappa=0, Omega=1 (error %.3e).'], ...
        geometryError);
end

op=ni2019_full_eigen_operator(cfg);
kyTolerance=1e-8*max(abs(op.k0),1);
finiteOpen=abs(imag(op.ky))<=kyTolerance&real(op.ky)>kyTolerance;
threshold=abs(op.ky)<=kyTolerance;
removedMask=finiteOpen|threshold;
removedOrders=op.orders(removedMask);
if ~isequal(removedOrders(:),[-1;0;1])
    error('Expected removed orders [-1 0 1], obtained [%s].', ...
        strtrim(sprintf('%d ',removedOrders)));
end

P=(op.N-1)/2;
positiveOrders=(2:P).';
q=(0:op.K-1).';
centralAllowed=paritySign*(-1).^q==1;
nUnknown=numel(positiveOrders)+op.K+nnz(centralAllowed);
T=complex(zeros(op.N+op.n_groove,nUnknown));
column=0;
for jj=1:numel(positiveOrders)
    column=column+1;
    plusId=find(op.orders==positiveOrders(jj),1);
    minusId=find(op.orders==-positiveOrders(jj),1);
    T(plusId,column)=1/sqrt(2);
    T(minusId,column)=paritySign/sqrt(2);
end
for iq=1:op.K
    column=column+1;
    T(op.N+iq,column)=1/sqrt(2);
    T(op.N+2*op.K+iq,column)=paritySign*(-1)^(iq-1)/sqrt(2);
end
for iq=find(centralAllowed).'
    column=column+1;
    T(op.N+op.K+iq,column)=1;
end
if column~=nUnknown
    error('Internal three-groove parity-map dimension mismatch.');
end

FparityFull=op.F*T;
outerPressureRows=(1:op.K).';
centerPressureRows=op.K+find(centralAllowed);
if paritySign<0
    velocityRows=3*op.K+find(op.orders>=1);
else
    velocityRows=3*op.K+find(op.orders>=0);
end
independentRows=[outerPressureRows;centerPressureRows;velocityRows];
Fparity=FparityFull(independentRows,:);
rowScale=max(vecnorm(Fparity,2,2),sqrt(eps));
Frow=Fparity./rowScale;
columnScale=max(vecnorm(Frow,2,1),sqrt(eps));
Fscaled=Frow./columnScale;
if any(~isfinite(Fscaled(:)))
    error('Parity-reduced operator contains NaN or Inf.');
end
[~,S,V]=svd(Fscaled,'econ');
singularValues=diag(S);
v=V(:,end);
uParity=v./columnScale(:);
z=T*uParity;
if opt.NormalizeMode
    normFactor=max(norm(z),eps);
    z=z/normFactor;
    uParity=uParity/normFactor;
end

A=z(1:op.N);
Cscaled=z(op.N+1:end);
Cphysical=Cscaled./op.vertical_scale;
CbyGroove=reshape(Cphysical,op.K,3);
surfacePressure=Cscaled.*op.cos_depth_normalized;
surfaceVelocity=Cscaled.*op.beta_sin_normalized;
pressureByGroove=reshape(surfacePressure,op.K,3);
velocityByGroove=reshape(surfaceVelocity,op.K,3);
pressureFraction=normalized_square(sum(abs(pressureByGroove).^2,2));
velocityFraction=normalized_square(sum(abs(velocityByGroove).^2,2));
bottomFraction=normalized_square(sum(abs(CbyGroove).^2,2));

reflectedC=(-1).^q.*CbyGroove(:,[3 2 1]);
parityResidual=norm(CbyGroove-paritySign*reflectedC,'fro')/ ...
    max(norm(CbyGroove,'fro'),eps);
oppositeResidual=norm(CbyGroove+paritySign*reflectedC,'fro')/ ...
    max(norm(CbyGroove,'fro'),eps);
rawResidual=norm(op.F*z)/max(norm(z),eps);
scaledResidual=norm(Fscaled*v)/max(norm(v),eps);
sigmaRatio=singularValues(end)/max(singularValues(1),eps);
grooveFraction=norm(Cphysical)^2/ ...
    max(norm(Cphysical)^2+norm(A)^2,eps);
propagatingPower=max(real(op.ky),0)./max(abs(op.k0),eps).*abs(A).^2;

betaByGroove=reshape(op.beta,op.K,3);
result=struct();
result.cfg=cfg;
result.parity=parityName;
result.parity_sign=paritySign;
result.full_operator=op;
result.T_parity=T;
result.F_parity_full_raw=FparityFull;
result.F_parity_raw=Fparity;
result.F_parity_scaled=Fscaled;
result.independent_rows=independentRows;
result.row_scale=rowScale;
result.column_scale=columnScale;
result.singular_values=singularValues;
result.sigma_ratio=sigmaRatio;
result.strict_residual=rawResidual;
result.residual=struct('scaled',scaledResidual,'raw',rawResidual);
result.removed_orders=removedOrders;
result.central_allowed_q=q(centralAllowed);
result.mode=struct('z_full',z,'u_parity',uParity,'v_scaled',v, ...
    'A',A,'C_scaled',Cscaled,'C_physical',Cphysical, ...
    'C_by_groove',CbyGroove, ...
    'surface_pressure_coefficients',surfacePressure, ...
    'surface_velocity_coefficients',surfaceVelocity);
result.transverse=struct('q',q,'bottom_fraction_by_q',bottomFraction, ...
    'surface_pressure_fraction_by_q',pressureFraction, ...
    'surface_velocity_fraction_by_q',velocityFraction, ...
    'beta_by_groove',betaByGroove, ...
    'central_beta',betaByGroove(:,2), ...
    'central_propagating',abs(imag(betaByGroove(:,2)))<=kyTolerance);
result.diagnostics=struct('geometry_error',geometryError, ...
    'parity_residual',parityResidual, ...
    'opposite_parity_residual',oppositeResidual, ...
    'physical_groove_fraction',grooveFraction);
result.radiation=struct('orders',op.orders,'amplitudes',A, ...
    'removed_amplitudes',A(removedMask), ...
    'removed_amplitudes_exact_zero',all(A(removedMask)==0), ...
    'power_by_order',propagatingPower, ...
    'total_propagating_power',sum(propagatingPower));
result.dimensions=struct('full',size(op.F), ...
    'parity_full_rows',size(FparityFull),'parity',size(Fparity), ...
    'extra_rows',size(Fparity,1)-size(Fparity,2));

if opt.Verbose
    fprintf(['three-groove %s strict operator: size %d x %d, ', ...
        'sigma %.3e, raw %.3e\n'],parityName,size(Fparity,1), ...
        size(Fparity,2),sigmaRatio,rawResidual);
    fprintf('  physical groove %.3f, parity residual %.3e, central q=[%s]\n', ...
        grooveFraction,parityResidual, ...
        strtrim(sprintf('%d ',result.central_allowed_q)));
end
end

function p=parse_parity(parity)
if isnumeric(parity)&&isscalar(parity)&&ismember(parity,[-1 1])
    p=parity;
    return;
end
name=lower(char(parity));
if strcmp(name,'odd')
    p=-1;
elseif strcmp(name,'even')
    p=1;
else
    error('parity must be odd/even or -1/+1.');
end
end

function fraction=normalized_square(value)
fraction=abs(value)/max(sum(abs(value)),eps);
end
