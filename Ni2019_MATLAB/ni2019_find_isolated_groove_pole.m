function pole=ni2019_find_isolated_groove_pole(width,depth,K,cWater,fSeed,varargin)
%NI2019_FIND_ISOLATED_GROOVE_POLE Complex pole of one baffled open groove.
p=inputParser;
addParameter(p,'OuterIterations',10);
addParameter(p,'Display','off');
addParameter(p,'SMax',max(120,(K+12)*pi));
addParameter(p,'PanelsPerPi',8);
addParameter(p,'QParity','all');
parse(p,varargin{:}); opt=p.Results;
switch lower(char(opt.QParity))
    case 'all', modeIndices=(1:K).';
    case 'even', modeIndices=(1:2:K).';       % q=0,2,...
    case 'odd', modeIndices=(2:2:K).';        % q=1,3,...
    otherwise, error('QParity must be all, even, or odd.');
end

f=fSeed; history=complex(nan(opt.OuterIterations+1,1));
sigmaHistory=nan(opt.OuterIterations+1,1); history(1)=f; exitflag=nan;
for outer=1:opt.OuterIterations
    op0=operator_at(f); Fraw=op0.F(modeIndices,modeIndices);
    rowScale=max(vecnorm(Fraw,2,2),sqrt(eps));
    Frow=Fraw./rowScale;
    columnScale=max(vecnorm(Frow,2,1),sqrt(eps));
    F0=Frow./columnScale;
    [U,S,V]=svd(F0,'econ'); sv=diag(S); u=U(:,end); v=V(:,end);
    sigmaHistory(outer)=sv(end)/max(sv(1),eps);
    x0=[real(f),imag(f)]/1e3;
    options=optimoptions('fsolve','Display',opt.Display, ...
        'FunctionTolerance',1e-12,'StepTolerance',1e-12, ...
        'OptimalityTolerance',1e-12,'MaxFunctionEvaluations',100, ...
        'MaxIterations',50);
    [x,~,exitflag]=fsolve(@projected_equation,x0,options);
    fNew=(x(1)+1i*x(2))*1e3; history(outer+1)=fNew;
    if abs(fNew-f)<=1e-10*max(abs(fNew),1), f=fNew; break; end
    f=fNew;
end
op=operator_at(f); Fraw=op.F(modeIndices,modeIndices);
rowScale=max(vecnorm(Fraw,2,2),sqrt(eps)); Frow=Fraw./rowScale;
columnScale=max(vecnorm(Frow,2,1),sqrt(eps)); Fscaled=Frow./columnScale;
[~,S,V]=svd(Fscaled,'econ'); sv=diag(S);
vScaled=V(:,end); Cselected=vScaled./transpose(columnScale);
Cscaled=complex(zeros(K,1)); Cscaled(modeIndices)=Cselected;
Cphysical=Cscaled./op.vertical_scale;
fraction=abs(Cphysical).^2/max(sum(abs(Cphysical).^2),eps);
surface=Cscaled.*op.cos_depth_normalized;
surfaceFraction=abs(surface).^2/max(sum(abs(surface).^2),eps);
pole=struct('frequency_hz',f,'Q',real(f)/(2*imag(f)), ...
    'sigma_ratio',sv(end)/max(sv(1),eps), ...
    'raw_residual',norm(op.F*Cscaled)/max(norm(Cscaled),eps), ...
    'exitflag',exitflag,'history_hz',history(isfinite(history)), ...
    'sigma_history',sigmaHistory(isfinite(sigmaHistory)), ...
    'operator',op,'mode_indices',modeIndices,'Q_parity',opt.QParity, ...
    'C_scaled',Cscaled,'C_physical',Cphysical, ...
    'transverse_fraction',fraction,'surface_pressure_fraction',surfaceFraction);

    function opLocal=operator_at(frequency)
        opLocal=ni2019_isolated_groove_operator_complex(width,depth,K, ...
            2*pi*frequency/cWater,'SMax',opt.SMax, ...
            'PanelsPerPi',opt.PanelsPerPi);
    end
    function residual=projected_equation(x)
        opLocal=operator_at((x(1)+1i*x(2))*1e3);
        F=opLocal.F(modeIndices,modeIndices)./rowScale./columnScale;
        value=u'*F*v;
        residual=[real(value);imag(value)];
    end
end
