function scan = ni2019_find_rayleigh_bic(cfg,varargin)
%NI2019_FIND_RAYLEIGH_BIC Search two groove depths at a first-order threshold.
%   The period and wavelength must satisfy a=lambda at normal incidence.
%   Widths and gap remain fixed. The objective combines the normalized
%   smallest singular value of the homogeneous modal operator and the
%   zeroth-order radiation amplitude of its approximate null vector.

p=inputParser;
addParameter(p,'DepthRange',[.01 .49],@(x)isnumeric(x)&&numel(x)==2);
addParameter(p,'GridSize',41,@(x)isnumeric(x)&&isscalar(x)&&x>=5);
addParameter(p,'Display','off',@(x)ischar(x)||isstring(x));
addParameter(p,'TargetOrder',-1,@(x)isnumeric(x)&&isscalar(x));
parse(p,varargin{:});
if numel(cfg.widths)~=2
    error('This depth scan is specifically for two-groove systems.');
end
k0check=2*pi/cfg.lambda;
kBlochCheck=k0check*sind(cfg.theta_i_deg);
kxTargetCheck=kBlochCheck+2*pi*p.Results.TargetOrder/cfg.a;
if abs(abs(kxTargetCheck)-k0check)>1e-9*k0check
    error('The supplied frequency and Bloch wavenumber do not put targetOrder at Rayleigh threshold.');
end

cfg.solve_scattering=false;
dvals=linspace(p.Results.DepthRange(1),p.Results.DepthRange(2),p.Results.GridSize);
score=nan(numel(dvals)); sigmaRatio=score; radiationFraction=score;
for i=1:numel(dvals)
    for j=1:numel(dvals)
        [score(i,j),sigmaRatio(i,j),radiationFraction(i,j)]= ...
            evaluate([dvals(i),dvals(j)]);
    end
end
[~,idx]=min(score,[],'all','linear');
[i0,j0]=ind2sub(size(score),idx);
x0=[dvals(i0),dvals(j0)];

if license('test','optimization_toolbox')
    opts=optimoptions('fmincon','Algorithm','sqp','Display',p.Results.Display, ...
        'MaxFunctionEvaluations',1500,'MaxIterations',250, ...
        'OptimalityTolerance',1e-12,'StepTolerance',1e-12);
    lb=p.Results.DepthRange(1)*ones(1,2);
    ub=p.Results.DepthRange(2)*ones(1,2);
    xbest=fmincon(@(x)evaluate(x),x0,[],[],[],[],lb,ub,[],opts);
    % Convert the minimum-singular-value estimate into a genuine complex
    % characteristic zero using alternating projected solves and SVD null
    % vector updates. Two real depths solve Re(g)=Im(g)=0.
    local=cfg; local.depths=xbest*local.lambda;
    opDepth=ni2019_full_eigen_operator(local);
    [uDepth,~,vDepth]=smallest_vectors(opDepth.Fscaled);
    rootOpts=optimoptions('fsolve','Display',p.Results.Display, ...
        'FunctionTolerance',1e-13,'StepTolerance',1e-13, ...
        'OptimalityTolerance',1e-13,'MaxFunctionEvaluations',1000);
    for refine=1:8
        xbest=fsolve(@depth_characteristic,xbest,rootOpts);
        xbest=min(max(xbest,lb),ub);
        local.depths=xbest*local.lambda;
        opDepth=ni2019_full_eigen_operator(local);
        [uNew,Snew,vNew]=smallest_vectors(opDepth.Fscaled);
        uDepth=align(uNew,uDepth); vDepth=align(vNew,vDepth);
        ss=diag(Snew);
        if ss(end)/ss(1)<1e-13, break; end
    end
else
    xbest=fminsearch(@penalized,x0,optimset('Display',p.Results.Display));
end
[bestScore,bestSigma,bestRadiation,mode,op]=evaluate(xbest);

scan=struct('cfg',cfg,'depth_values',dvals,'score',score, ...
    'sigma_ratio',sigmaRatio,'radiation_fraction',radiationFraction, ...
    'coarse_best_depths',x0,'best_depths',xbest, ...
    'best_score',bestScore,'best_sigma_ratio',bestSigma, ...
    'best_radiation_fraction',bestRadiation,'mode',mode,'operator',op);

    function f=penalized(x)
        lo=p.Results.DepthRange(1); hi=p.Results.DepthRange(2);
        penalty=1e3*sum(max(lo-x,0).^2+max(x-hi,0).^2);
        f=evaluate(min(max(x,lo),hi))+penalty;
    end
    function [f,srel,rad,v,Dscaled]=evaluate(depths)
        local=cfg; local.depths=depths(:).'*local.lambda;
        out=ni2019_full_eigen_operator(local);
        Dscaled=out.Fscaled;
        [~,S,V]=svd(Dscaled,'econ');
        singular=diag(S);
        srel=singular(end)/max(singular(1),eps);
        v=V(:,end)./transpose(out.column_scale);
        v=v/norm(v);
        open=abs(imag(out.ky))<1e-10*out.k0 & real(out.ky)>1e-8*out.k0;
        openIds=find(open);
        rad=sum(abs(v(openIds)).^2);
        % A true threshold BIC needs both a null operator and no radiation
        % into the only finite-flux open channel n=0.
        f=sqrt(srel^2+rad^2);
    end
    function y=depth_characteristic(depths)
        localDepth=cfg; localDepth.depths=depths(:).'*localDepth.lambda;
        depthOp=ni2019_full_eigen_operator(localDepth);
        g=uDepth'*depthOp.Fscaled*vDepth;
        y=[real(g);imag(g)];
    end
end

function [u,S,v]=smallest_vectors(F)
[U,S,V]=svd(F,'econ'); u=U(:,end); v=V(:,end);
end
function x=align(x,reference)
z=reference'*x;
if z~=0, x=x*exp(-1i*angle(z)); end
end
