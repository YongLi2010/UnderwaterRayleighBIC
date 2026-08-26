function track = ni2019_continue_pole_fixed_sheet(cfg,kappaPath,OmegaSeed,targetOrder,varargin)
%NI2019_CONTINUE_POLE_FIXED_SHEET Continue a pole through one Floquet cut.
%   The selected target-order transverse wavenumber q is used as the
%   nonlinear variable, so its sign is continued rather than reselected
%   from Re(Omega) at every kappa.  This is intended for diagnosing apparent
%   jumps caused by a physical-sheet branch cut.

p = inputParser;
addParameter(p,'OuterIterations',10,@(x)isscalar(x)&&x>=1);
addParameter(p,'Display','off',@(x)ischar(x)||isstring(x));
addParameter(p,'FunctionTolerance',1e-12,@(x)isscalar(x)&&x>0);
addParameter(p,'StepTolerance',1e-12,@(x)isscalar(x)&&x>0);
parse(p,varargin{:}); opt=p.Results;

kappaPath = kappaPath(:).';
if isempty(kappaPath), error('kappaPath must not be empty.'); end
if ~isfield(cfg,'a'), error('cfg.a is required.'); end
a = cfg.a;

% Match the selected order to the outgoing operator at the supplied seed.
opSeed = ni2019_outgoing_eigen_operator_complex(cfg, ...
    2*pi*OmegaSeed/a,2*pi*kappaPath(1)/a);
targetId = find(opSeed.orders==targetOrder,1);
if isempty(targetId), error('Target order %d is outside the basis.',targetOrder); end
qPrevious = opSeed.ky(targetId)*a/(2*pi);

n = numel(kappaPath);
q = complex(nan(1,n)); Omega=q;
sigmaRatio=nan(1,n); residual=nan(1,n); overlap=nan(1,n);
exitflag=nan(1,n); mode=cell(1,n);

[op0,~]=operator_at(qPrevious,kappaPath(1));
[U,S,V]=svd(op0.Fscaled,'econ');
uPrevious=U(:,end); vPrevious=V(:,end);
vModePrevious=[];

options=optimoptions('fsolve','Display',char(opt.Display), ...
    'FunctionTolerance',opt.FunctionTolerance, ...
    'StepTolerance',opt.StepTolerance, ...
    'OptimalityTolerance',opt.FunctionTolerance, ...
    'MaxFunctionEvaluations',300,'MaxIterations',100);

for m=1:n
    kap=kappaPath(m);
    x=[real(qPrevious),imag(qPrevious)]; flag=0;
    uRef=uPrevious; vRef=vPrevious;
    for outer=1:opt.OuterIterations
        [x,~,flag]=fsolve(@(xx)projected_equation(xx,kap,uRef,vRef),x,options);
        qNow=complex(x(1),x(2));
        [opNow,~]=operator_at(qNow,kap);
        [Un,Sn,Vn]=svd(opNow.Fscaled,'econ'); s=diag(Sn);
        uNew=phase_align(Un(:,end),uRef);
        vNew=phase_align(Vn(:,end),vRef);
        uRef=uNew; vRef=vNew;
        if s(end)/max(s(1),eps)<1e-11, break; end
    end
    qNow=complex(x(1),x(2));
    [opNow,omNow]=operator_at(qNow,kap);
    [Un,Sn,Vn]=svd(opNow.Fscaled,'econ'); s=diag(Sn);
    vScaled=Vn(:,end);
    z=vScaled./transpose(opNow.column_scale); z=z/max(norm(z),eps);
    if ~isempty(vModePrevious), overlap(m)=abs(vModePrevious'*vScaled); end
    q(m)=qNow; Omega(m)=omNow;
    sigmaRatio(m)=s(end)/max(s(1),eps);
    residual(m)=norm(opNow.Fscaled*vScaled);
    exitflag(m)=flag; mode{m}=z;
    qPrevious=qNow; uPrevious=Un(:,end); vPrevious=vScaled;
    vModePrevious=vScaled;
end

track=struct('kappa',kappaPath,'Omega',Omega,'q',q, ...
    'Q',real(Omega)./(2*max(abs(imag(Omega)),realmin)), ...
    'sigma_ratio',sigmaRatio,'residual',residual,'mode_overlap',overlap, ...
    'exitflag',exitflag,'mode',{mode},'target_order',targetOrder, ...
    'convention','fixed target-order q sign; exp(+j*omega*t)');

    function value=projected_equation(x,kap,uRef,vRef)
        qTrial=complex(x(1),x(2));
        [opTrial,~]=operator_at(qTrial,kap);
        g=uRef'*opTrial.Fscaled*vRef;
        value=[real(g);imag(g)];
    end

    function [opLocal,om]=operator_at(qbar,kap)
        om=sqrt((kap+targetOrder)^2+qbar^2);
        if real(om)<0, om=-om; end
        opLocal=ni2019_full_eigen_operator_complex(cfg,2*pi*om/a, ...
            2*pi*kap/a,targetOrder,2*pi*qbar/a);
    end
end

function value=phase_align(value,reference)
overlap=reference'*value;
if overlap~=0, value=value*exp(-1i*angle(overlap)); end
end
