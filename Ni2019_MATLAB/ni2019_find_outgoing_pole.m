function pole=ni2019_find_outgoing_pole(cfg,cWater,fSeed,varargin)
%NI2019_FIND_OUTGOING_POLE Refine a complex-frequency outgoing eigenpole.
%   fSeed is complex and is expressed in Hz.  With exp(+i*omega*t), a
%   passive leaky mode has Im(f)>0.

p=inputParser;
addParameter(p,'OuterIterations',10);
addParameter(p,'Display','off');
addParameter(p,'FunctionTolerance',1e-11);
addParameter(p,'StepTolerance',1e-11);
parse(p,varargin{:}); opt=p.Results;

f=fSeed; history=complex(nan(opt.OuterIterations+1,1));
sigmaHistory=nan(opt.OuterIterations+1,1); history(1)=f;
for outer=1:opt.OuterIterations
    op0=operator_at(f);
    rowScale=op0.row_scale(:); columnScale=op0.column_scale(:).';
    F0=op0.F./rowScale./columnScale;
    [U,S,V]=svd(F0,'econ'); sv=diag(S);
    u=U(:,end); v=V(:,end);
    sigmaHistory(outer)=sv(end)/max(sv(1),eps);
    x0=[real(f),imag(f)]/1e3;
    options=optimoptions('fsolve','Display',opt.Display, ...
        'FunctionTolerance',opt.FunctionTolerance, ...
        'StepTolerance',opt.StepTolerance, ...
        'OptimalityTolerance',opt.FunctionTolerance, ...
        'MaxFunctionEvaluations',100,'MaxIterations',50);
    [x,~,exitflag]=fsolve(@projected_equation,x0,options);
    fNew=(x(1)+1i*x(2))*1e3;
    history(outer+1)=fNew;
    if abs(fNew-f)<=1e-9*max(abs(fNew),1)
        f=fNew; break;
    end
    f=fNew;
end

op=operator_at(f);
[~,S,V]=svd(op.Fscaled,'econ'); sv=diag(S);
vScaled=V(:,end); z=vScaled./transpose(op.column_scale);
z=z/max(norm(z),eps);
Cscaled=z(op.N+1:end); Cphysical=Cscaled./op.vertical_scale;
surfacePressure=Cscaled.*op.cos_depth_normalized;
fraction=abs(Cphysical).^2/max(sum(abs(Cphysical).^2),eps);
surfaceFraction=abs(surfacePressure).^2/ ...
    max(sum(abs(surfacePressure).^2),eps);

pole=struct('frequency_hz',f,'Q',real(f)/(2*imag(f)), ...
    'sigma_ratio',sv(end)/max(sv(1),eps), ...
    'raw_residual',norm(op.F*z)/max(norm(z),eps), ...
    'exitflag',exitflag,'history_hz',history(isfinite(history)), ...
    'sigma_history',sigmaHistory(isfinite(sigmaHistory)), ...
    'operator',op,'mode_vector',z,'C_physical',Cphysical, ...
    'transverse_fraction',fraction,'surface_pressure_fraction',surfaceFraction);

    function opLocal=operator_at(frequency)
        k0=2*pi*frequency/cWater;
        opLocal=ni2019_outgoing_eigen_operator_complex(cfg,k0,0);
    end

    function residual=projected_equation(x)
        frequency=(x(1)+1i*x(2))*1e3;
        opLocal=operator_at(frequency);
        F=opLocal.F./rowScale./columnScale;
        value=u'*F*v;
        residual=[real(value);imag(value)];
    end
end
