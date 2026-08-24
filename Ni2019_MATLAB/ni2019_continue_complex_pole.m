function poles=ni2019_continue_complex_pole(cfg,kappa,varargin)
%NI2019_CONTINUE_COMPLEX_POLE Continue the Rayleigh pole for kappa>0.
%   kappa=kBloch*a/(2*pi). The n=-1 grazing-channel normal wavenumber q is
%   the nonlinear spectral variable, which regularizes the branch point.
%   PredictorExponent controls the initial/continuation q scale.  Keep the
%   historical default 1 for regular branches; use 1/2 when the split
%   double-Rayleigh endpoint has q=O(sqrt(kappa)).

p=inputParser;
addParameter(p,'TargetOrder',-1,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'OuterIterations',5,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'Display','off',@(x)ischar(x)||isstring(x));
addParameter(p,'InitialScale',.35+.35i,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'PredictorExponent',1,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'BlochSign',1,@(x)isnumeric(x)&&isscalar(x)&&abs(x)==1);
parse(p,varargin{:});
if any(kappa<=0)
    error('Use strictly positive kappa; the BIC endpoint kappa=0 is appended internally.');
end
kappa=sort(kappa(:).'); a=cfg.a; target=p.Results.TargetOrder;
blochSign=p.Results.BlochSign;
predictorExponent=p.Results.PredictorExponent;

% BIC endpoint: normalized Omega=1 and q=0.
op0=ni2019_full_eigen_operator_complex(cfg,2*pi/a,0,target,0);
[U0,S0,V0]=svd(op0.Fscaled,'econ');
u=U0(:,end); v=V0(:,end);
s0=diag(S0);

n=numel(kappa); qbar=complex(nan(1,n)); Omega=qbar; residual=nan(1,n);
sigmaRatio=residual; radiationFraction=residual; exitflag=nan(1,n);
xPrevious=[];
opts=optimoptions('fsolve','Display',p.Results.Display, ...
    'FunctionTolerance',1e-11,'StepTolerance',1e-11, ...
    'OptimalityTolerance',1e-11,'MaxFunctionEvaluations',800);

for m=1:n
    kap=blochSign*kappa(m);
    if isempty(xPrevious)
        % qbar with equal real/imaginary parts gives Im(Omega)>0 under the
        % exp(+j*omega*t) convention used by the modal formulation.
        seed=p.Results.InitialScale*abs(kap)^predictorExponent;
        x=[real(seed),imag(seed)];
    else
        x=xPrevious*(abs(kap)/kappa(m-1))^predictorExponent;
    end
    flag=0;
    for outer=1:p.Results.OuterIterations
        [x,~,flag]=fsolve(@projected_equation,x,opts);
        [op,~]=operator_from_x(x,kap);
        [Un,Sn,Vn]=svd(op.Fscaled,'econ');
        u=phase_align(Un(:,end),u); v=phase_align(Vn(:,end),v);
        singular=diag(Sn);
        if singular(end)/singular(1)<1e-9, break; end
    end
    [op,om]=operator_from_x(x,kap);
    [~,Ss,Vs]=svd(op.Fscaled,'econ'); singular=diag(Ss);
    z=Vs(:,end)./transpose(op.column_scale); z=z/norm(z);
    id0=find(op.orders==0,1);
    qbar(m)=complex(x(1),x(2)); Omega(m)=om;
    sigmaRatio(m)=singular(end)/singular(1);
    residual(m)=norm(op.Fscaled*Vs(:,end));
    radiationFraction(m)=abs(z(id0))^2; exitflag(m)=flag;
    xPrevious=x;
end

Q=real(Omega)./(2*abs(imag(Omega)));
poles=struct('kappa',[0,blochSign*kappa],'qbar',[0,qbar],'Omega',[1,Omega], ...
    'Q',[Inf,Q],'continued_sigma_ratio',[s0(end)/s0(1),sigmaRatio], ...
    'residual',[s0(end),residual], ...
    'radiation_fraction',[0,radiationFraction], ...
    'exitflag',[1,exitflag],'target_order',target,'convention','exp(+j*omega*t)');

    function y=projected_equation(x)
        [op,~]=operator_from_x(x,kap);
        g=u'*op.Fscaled*v;
        y=[real(g);imag(g)];
    end
    function [op,om]=operator_from_x(x,kap)
        qb=complex(x(1),x(2)); kxTarget=kap+target;
        om=sqrt(kxTarget^2+qb^2);
        if real(om)<0, om=-om; end
        op=ni2019_full_eigen_operator_complex(cfg,2*pi*om/a, ...
            2*pi*kap/a,target,2*pi*qb/a);
    end
end

function x=phase_align(x,reference)
overlap=reference'*x;
if overlap~=0, x=x*exp(-1i*angle(overlap)); end
end
