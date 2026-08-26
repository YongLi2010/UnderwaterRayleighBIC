function pole = ni2019_refine_outgoing_pole_kappa(cfg,kappa,OmegaSeed,varargin)
%NI2019_REFINE_OUTGOING_POLE_KAPPA Refine one complex Bloch eigenfrequency.
%   The normalized variables are Omega=k0*a/(2*pi) and
%   kappa=kBloch*a/(2*pi).  With exp(+i*omega*t), passive outgoing poles have
%   Im(Omega)>=0.  The nonlinear solve freezes the current equilibrating
%   transforms and singular vectors within each projected-Newton step, then
%   refreshes them in an outer iteration.

p = inputParser;
addParameter(p,'OuterIterations',8,@(x)isscalar(x)&&x>=1);
addParameter(p,'Display','off',@(x)ischar(x)||isstring(x));
addParameter(p,'FunctionTolerance',2e-11,@(x)isscalar(x)&&x>0);
addParameter(p,'StepTolerance',2e-11,@(x)isscalar(x)&&x>0);
parse(p,varargin{:}); opt=p.Results;

Omega = complex(OmegaSeed);
history = complex(nan(opt.OuterIterations+1,1)); history(1)=Omega;
sigmaHistory = nan(opt.OuterIterations+1,1); exitflag=nan;
for outer = 1:opt.OuterIterations
    op0 = operator_at(Omega);
    rowScale = op0.row_scale(:);
    columnScale = op0.column_scale(:).';
    F0 = op0.F./rowScale./columnScale;
    [U,S,V] = svd(F0,'econ'); sv=diag(S);
    u=U(:,end); v=V(:,end);
    sigmaHistory(outer)=sv(end)/max(sv(1),eps);
    options=optimoptions('fsolve','Display',char(opt.Display), ...
        'FunctionTolerance',opt.FunctionTolerance, ...
        'StepTolerance',opt.StepTolerance, ...
        'OptimalityTolerance',opt.FunctionTolerance, ...
        'MaxFunctionEvaluations',120,'MaxIterations',60);
    [x,~,exitflag]=fsolve(@projected_equation,[real(Omega),imag(Omega)],options);
    OmegaNew=complex(x(1),x(2)); history(outer+1)=OmegaNew;
    if abs(OmegaNew-Omega)<=2e-10*max(abs(OmegaNew),1)
        Omega=OmegaNew; break;
    end
    Omega=OmegaNew;
end

op=operator_at(Omega);
[~,S,V]=svd(op.Fscaled,'econ'); sv=diag(S);
vScaled=V(:,end); z=vScaled./transpose(op.column_scale);
z=z/max(norm(z),eps);
A=z(1:op.N); Cscaled=z(op.N+1:end);
Cphysical=Cscaled./op.vertical_scale;
surface=Cscaled.*op.cos_depth_normalized;
grooveFraction=abs(Cphysical).^2/max(sum(abs(Cphysical).^2),eps);
surfaceFraction=abs(surface).^2/max(sum(abs(surface).^2),eps);

pole=struct('kappa',kappa,'Omega',Omega, ...
    'Q',real(Omega)/(2*max(abs(imag(Omega)),realmin)), ...
    'sigma_ratio',sv(end)/max(sv(1),eps), ...
    'raw_residual',norm(op.F*z)/max(norm(z),eps), ...
    'exitflag',exitflag,'history',history(isfinite(history)), ...
    'sigma_history',sigmaHistory(isfinite(sigmaHistory)), ...
    'mode_vector',z,'A',A,'C_physical',Cphysical, ...
    'groove_fraction',grooveFraction, ...
    'surface_pressure_fraction',surfaceFraction,'operator',op);

    function opLocal=operator_at(om)
        opLocal=ni2019_outgoing_eigen_operator_complex(cfg, ...
            2*pi*om/cfg.a,2*pi*kappa/cfg.a);
    end

    function residual=projected_equation(x)
        opLocal=operator_at(complex(x(1),x(2)));
        value=u'*(opLocal.F./rowScale./columnScale)*v;
        residual=[real(value);imag(value)];
    end
end
