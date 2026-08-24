function poles=ni2019_continue_offgamma_pole(cfg,kappa0,deltaKappa,varargin)
%NI2019_CONTINUE_OFFGAMMA_POLE Continue a pole away from a Rayleigh endpoint.
%   deltaKappa is the signed displacement from kappa0. The endpoint obeys
%   Omega0=abs(kappa0+targetOrder), and qbar regularizes its square root.

p=inputParser;
addParameter(p,'TargetOrder',-1);
addParameter(p,'OuterIterations',8);
addParameter(p,'Display','off');
addParameter(p,'InitialScale',.2+.2i);
parse(p,varargin{:});
deltaKappa=deltaKappa(:).';
if any(deltaKappa==0), error('deltaKappa must exclude zero; it is appended internally.'); end
[~,sequence]=sort(abs(deltaKappa)); deltaKappa=deltaKappa(sequence);
a=cfg.a; target=p.Results.TargetOrder; Omega0=abs(kappa0+target);
op0=ni2019_full_eigen_operator_complex(cfg,2*pi*Omega0/a, ...
    2*pi*kappa0/a,target,0);
[U0,S0,V0]=svd(op0.Fscaled,'econ'); u=U0(:,end); v=V0(:,end);
s0=diag(S0); n=numel(deltaKappa);
qbar=complex(nan(1,n)); Omega=qbar; sigmaRatio=nan(1,n);
radiationFraction=nan(1,n); exitflag=nan(1,n); xPrevious=[];
opts=optimoptions('fsolve','Display',p.Results.Display, ...
    'FunctionTolerance',1e-12,'StepTolerance',1e-12, ...
    'OptimalityTolerance',1e-12,'MaxFunctionEvaluations',1000);
for m=1:n
    kap=kappa0+deltaKappa(m);
    if isempty(xPrevious)
        seed=p.Results.InitialScale*sqrt(abs(deltaKappa(m)));
        x=[real(seed),imag(seed)];
    else
        x=xPrevious*sqrt(abs(deltaKappa(m)/deltaKappa(m-1)));
    end
    flag=0;
    for outer=1:p.Results.OuterIterations
        [x,~,flag]=fsolve(@projected_equation,x,opts);
        [op,~]=operator_from_x(x,kap);
        [Un,Sn,Vn]=svd(op.Fscaled,'econ');
        u=phase_align(Un(:,end),u); v=phase_align(Vn(:,end),v);
        ss=diag(Sn);
        if ss(end)/ss(1)<1e-10, break; end
    end
    [op,om]=operator_from_x(x,kap);
    [~,Ss,Vs]=svd(op.Fscaled,'econ'); ss=diag(Ss);
    z=Vs(:,end)./transpose(op.column_scale); z=z/norm(z);
    finiteOpen=abs(imag(op.ky))<1e-7*abs(op.k0) & real(op.ky)>1e-7*abs(op.k0);
    qbar(m)=complex(x(1),x(2)); Omega(m)=om;
    sigmaRatio(m)=ss(end)/ss(1);
    radiationFraction(m)=sum(abs(z(finiteOpen)).^2);
    exitflag(m)=flag; xPrevious=x;
end
Q=real(Omega)./(2*abs(imag(Omega)));
poles=struct('kappa',[kappa0,kappa0+deltaKappa], ...
    'delta_kappa',[0,deltaKappa],'qbar',[0,qbar], ...
    'Omega',[Omega0,Omega],'Q',[Inf,Q], ...
    'continued_sigma_ratio',[s0(end)/s0(1),sigmaRatio], ...
    'radiation_fraction',[0,radiationFraction], ...
    'exitflag',[1,exitflag],'target_order',target, ...
    'convention','exp(+j*omega*t)');

    function y=projected_equation(x)
        [op,~]=operator_from_x(x,kap);
        g=u'*op.Fscaled*v; y=[real(g);imag(g)];
    end
    function [op,om]=operator_from_x(x,kap)
        qb=complex(x(1),x(2)); longitudinal=kap+target;
        om=sqrt(longitudinal^2+qb^2);
        if real(om)<0, om=-om; end
        op=ni2019_full_eigen_operator_complex(cfg,2*pi*om/a, ...
            2*pi*kap/a,target,2*pi*qb/a);
    end
end

function x=phase_align(x,reference)
overlap=reference'*x;
if overlap~=0, x=x*exp(-1i*angle(overlap)); end
end
