function track=ni2019_track_leaky_pole_to_rayleigh(cfg,kappa0,varargin)
%NI2019_TRACK_LEAKY_POLE_TO_RAYLEIGH Track a pole without using q=0 data.
%   TRACK=NI2019_TRACK_LEAKY_POLE_TO_RAYLEIGH(CFG,KAPPA0) starts at a
%   finite distance from the n=-1 Rayleigh endpoint, locates a complex pole
%   by a local q-plane scan, and then continues that same pole toward the
%   endpoint.  KAPPA0 is the normalized endpoint Bloch number k_B*a/(2*pi).
%
%   This routine is deliberately independent of the endpoint singular
%   vector.  At the first point it scans nonzero q seeds, constructs the
%   local smallest-singular-vector pair, and refines the projected
%   characteristic equation.  Subsequent points use only the preceding
%   pole/mode as a predictor.  q is normalized as q_y*a/(2*pi); the
%   continuation convention is exp(+i*omega*t).
%
%   The returned structure contains the pole path (Omega and q), the
%   smallest-singular-value residual, the relative singular value, modal
%   overlaps, and a basis-dependent groove-content diagnostic.  The final
%   row is the directly evaluated Rayleigh endpoint, not a root obtained by
%   forcing q=0 during continuation.

p=inputParser;
addParameter(p,'TargetOrder',-1,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'DeltaStart',-1e-2,@(x)isnumeric(x)&&isscalar(x)&&x~=0);
addParameter(p,'DeltaEnd',-1e-6,@(x)isnumeric(x)&&isscalar(x)&&x~=0);
addParameter(p,'NumSteps',31,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
addParameter(p,'ScanRealFactors',-.25:.125:.25,@isnumeric);
addParameter(p,'ScanImagFactors',-.60:.10:-.10,@isnumeric);
addParameter(p,'OuterIterations',8,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'ScanOuterIterations',10,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'Display','off',@(x)ischar(x)||isstring(x));
addParameter(p,'Verbose',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'MinOverlap',.80,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<=1);
parse(p,varargin{:});

target=p.Results.TargetOrder;
deltaStart=p.Results.DeltaStart;
deltaEnd=p.Results.DeltaEnd;
if sign(deltaStart)~=sign(deltaEnd)
    error('DeltaStart and DeltaEnd must have the same sign for one-sheet continuation.');
end
if abs(deltaStart)<=abs(deltaEnd)
    error('|DeltaStart| must be larger than |DeltaEnd|.');
end
if ~isfield(cfg,'a'), error('cfg.a is required.'); end
a=cfg.a;
if ~isfield(cfg,'N'), cfg.N=101; end
if ~isfield(cfg,'K'), cfg.K=10; end

% The endpoint obeys Omega0=|kappa0+target|.  The default target=-1
% and 0<kappa0<1/2 therefore give Omega0=1-kappa0.
Omega0=abs(kappa0+target);
if Omega0==0, error('The chosen endpoint has zero normalized frequency.'); end

% Logarithmic path, ordered from the independently located leaky pole to
% the endpoint.  Do not include delta=0 in the nonlinear solve.
deltaPath=sign(deltaStart)*logspace(log10(abs(deltaStart)), ...
    log10(abs(deltaEnd)),p.Results.NumSteps);

% --- Independent starting-pole search ---------------------------------
% q=0 is intentionally not part of this scan.  The scale sqrt(|delta|) is
% only a dimensional predictor for a local q-plane search; no endpoint
% singular vector or endpoint mode is used.
qScale=sqrt(abs(deltaStart));
realFactors=p.Results.ScanRealFactors(:).';
imagFactors=p.Results.ScanImagFactors(:).';
seedList=zeros(numel(realFactors)*numel(imagFactors),1);
seedId=0;
scanSigma=nan(size(seedList));
scanMode=cell(size(seedList));
kapStart=kappa0+deltaStart;
for ir=1:numel(realFactors)
    for ii=1:numel(imagFactors)
        seedId=seedId+1;
        qseed=qScale*(realFactors(ir)+1i*imagFactors(ii));
        seedList(seedId)=qseed;
        [opSeed,~,~,~,~,srel]=mode_at_q(qseed,kapStart);
        scanSigma(seedId)=srel;
        scanMode{seedId}=opSeed;
    end
end

% Refine every few best seeds locally.  This makes the selected root an
% independently located pole rather than a continuation from the endpoint.
[~,order]=sort(scanSigma,'ascend');
nRefine=min(5,numel(order));
candidate=cell(nRefine,1);
for n=1:nRefine
    id=order(n);
    qseed=seedList(id);
    opSeed=scanMode{id};
    [U,S,V]=svd(opSeed.Fscaled,'econ');
    uRef=U(:,end); vRef=V(:,end);
    candidate{n}=refine_pole(kapStart,qseed,uRef,vRef, ...
        p.Results.ScanOuterIterations);
end

% Prefer a converged pole on the physical exp(+i omega t) branch (Re q<0,
% Im q<0 gives Im Omega>0 for the present delta<0 path).  If several roots
% meet that condition, choose the one closest to the best scan seed.
ok=false(nRefine,1);
for n=1:nRefine
    c=candidate{n};
    ok(n)=c.sigma_ratio<1e-7 && isfinite(c.sigma_ratio) && ...
        imag(c.q)<=1e-8 && imag(c.Omega)>=-1e-8;
end
if any(ok)
    ids=find(ok);
    [~,j]=min(cellfun(@(c)c.sigma_ratio,candidate(ids)));
    c0=candidate{ids(j)};
else
    [~,j]=min(cellfun(@(c)c.sigma_ratio,candidate));
    c0=candidate{j};
end

% --- Continuation ------------------------------------------------------
nPath=numel(deltaPath);
qPath=complex(nan(1,nPath)); OmegaPath=qPath;
sigmaMin=nan(1,nPath); sigmaRatio=nan(1,nPath);
residual=nan(1,nPath); grooveContent=nan(1,nPath);
radiationFraction=nan(1,nPath);
modeOverlap=nan(1,nPath); leftOverlap=nan(1,nPath);
exitflag=nan(1,nPath); jump=false(1,nPath);
sheetLabel=cell(1,nPath);
modeRight=cell(1,nPath); modeLeft=cell(1,nPath);

qPath(1)=c0.q; OmegaPath(1)=c0.Omega;
sigmaMin(1)=c0.sigma_min; sigmaRatio(1)=c0.sigma_ratio;
residual(1)=c0.residual; grooveContent(1)=c0.groove_content;
radiationFraction(1)=c0.radiation_fraction;
exitflag(1)=c0.exitflag; sheetLabel{1}=sheet_label(c0.q,c0.Omega);
modeRight{1}=c0.z; modeLeft{1}=c0.u;

qPrevious=c0.q; zPrevious=c0.z; uPrevious=c0.u; vPrevious=c0.v;
for m=2:nPath
    kap=kappa0+deltaPath(m);
    % Secant predictor in q is safer than reusing q=0.  The first point is
    % already independently located at |deltaStart|>=O(1e-2).
    if m==2
        qPredict=qPrevious;
    else
        d1=deltaPath(m-1)-deltaPath(m-2);
        d2=deltaPath(m)-deltaPath(m-1);
        qPredict=qPrevious+(qPrevious-qPath(m-2))*d2/d1;
    end
    c=refine_pole(kap,qPredict,uPrevious,vPrevious, ...
        p.Results.OuterIterations);
    qPath(m)=c.q; OmegaPath(m)=c.Omega;
    sigmaMin(m)=c.sigma_min; sigmaRatio(m)=c.sigma_ratio;
    residual(m)=c.residual; grooveContent(m)=c.groove_content;
    radiationFraction(m)=c.radiation_fraction;
    exitflag(m)=c.exitflag; sheetLabel{m}=sheet_label(c.q,c.Omega);
    zNow=c.z; uNow=c.u;
    modeOverlap(m)=abs(zPrevious'*zNow);
    leftOverlap(m)=abs(uPrevious'*uNow);
    % The magnitude of q changes as sqrt(|delta kappa|), so a large q
    % displacement is expected on the logarithmic path.  A branch jump is
    % therefore flagged from the normalized modal overlap, not from an
    % absolute q step.
    jump(m)=modeOverlap(m)<p.Results.MinOverlap;
    modeRight{m}=zNow; modeLeft{m}=uNow;
    qPrevious=c.q; zPrevious=zNow; uPrevious=uNow; vPrevious=c.v;
end

% Direct endpoint diagnostic.  This row is only an evaluation of the
% supplied candidate geometry at (kappa0,Omega0,q=0), so its residual is a
% meaningful test of the candidate and is not silently forced to zero.
opEnd=ni2019_full_eigen_operator_complex(cfg,2*pi*Omega0/a, ...
    2*pi*kappa0/a,target,0);
[Ue,Se,Ve]=svd(opEnd.Fscaled,'econ');
ze=Ve(:,end)./transpose(opEnd.column_scale); ze=ze/norm(ze);
ue=Ue(:,end);
endpointSigma=diag(Se);
endpointFinite=abs(imag(opEnd.ky))<1e-8*abs(opEnd.k0) & ...
    real(opEnd.ky)>1e-8*abs(opEnd.k0);
Aend=ze(1:opEnd.N); endpointRad=sum(abs(Aend(endpointFinite)).^2)/ ...
    max(sum(abs(Aend).^2),eps);

% Append endpoint as a diagnostic point.  Its mode overlap is with the last
% continued pole and is useful for detecting a branch jump near q=0.
deltaOut=[deltaPath,0]; qOut=[qPath,0]; OmegaOut=[OmegaPath,Omega0];
sigMinOut=[sigmaMin,endpointSigma(end)];
sigRatioOut=[sigmaRatio,endpointSigma(end)/endpointSigma(1)];
residualOut=[residual,endpointSigma(end)];
grooveOut=[grooveContent,groove_content_from_mode(ze,opEnd)];
radiationOut=[radiationFraction,endpointRad];
overlapOut=[modeOverlap,abs(zPrevious'*ze)];
leftOverlapOut=[leftOverlap,abs(uPrevious'*ue)];
exitOut=[exitflag,1]; jumpOut=[jump,overlapOut(end)<p.Results.MinOverlap];
sheetOut=[sheetLabel,{'Rayleigh endpoint (q=0)'}];

track=struct();
track.kappa=kappa0+deltaOut;
track.delta_kappa=deltaOut;
track.q=qOut;
track.Omega=OmegaOut;
track.Q=real(OmegaOut)./(2*abs(imag(OmegaOut)));
track.sigma_min=sigMinOut;
track.sigma_ratio=sigRatioOut;
track.singular_residual=residualOut;
track.mode_overlap=overlapOut;
track.left_mode_overlap=leftOverlapOut;
track.groove_content=grooveOut;
track.radiation_fraction=radiationOut;
track.exitflag=exitOut;
track.jump_detected=jumpOut;
track.sheet=sheetOut;
track.target_order=target;
track.kappa0=kappa0;
track.Omega0=Omega0;
track.convention='exp(+j*omega*t)';
track.endpoint=struct('sigma_min',endpointSigma(end), ...
    'sigma_ratio',endpointSigma(end)/endpointSigma(1), ...
    'finite_channel_radiation_fraction',endpointRad, ...
    'groove_content',grooveOut(end), ...
    'mode',ze);
track.independent_start=struct('delta_kappa',deltaStart, ...
    'kappa',kapStart,'q_seed',seedList(order(1)), ...
    'scan_seed',seedList,'scan_sigma_ratio',scanSigma, ...
    'candidate_sigma_ratio',cellfun(@(c)c.sigma_ratio,candidate), ...
    'q',c0.q,'Omega',c0.Omega,'sigma_min',c0.sigma_min, ...
    'sigma_ratio',c0.sigma_ratio,'residual',c0.residual, ...
    'groove_content',c0.groove_content,'radiation_fraction', ...
    c0.radiation_fraction,'mode',c0.z);
track.same_branch=~any(jumpOut) && all(isfinite(sigRatioOut(1:end-1))) && ...
    abs(qPath(end))<max(20*abs(qPath(1))*sqrt(abs(deltaEnd/deltaStart)),1e-5);
track.q_to_zero=abs(qPath(end));
track.Omega_to_endpoint=abs(OmegaPath(end)-Omega0);
track.numerical_warning=~all(sigRatioOut(1:end-1)<1e-7) || any(jumpOut);

if p.Results.Verbose
    fprintf('Independent leaky-pole to Rayleigh tracking\n');
    fprintf('  start delta_kappa = %+.6e, kappa = %.12f\n', ...
        deltaStart,kapStart);
    fprintf('  selected nonzero seed q = %+.6e%+.6ei\n', ...
        real(track.independent_start.q_seed),imag(track.independent_start.q_seed));
    fprintf('  start pole q = %+.12e%+.12ei, Omega = %+.12e%+.12ei\n', ...
        real(c0.q),imag(c0.q),real(c0.Omega),imag(c0.Omega));
    fprintf('  start sigma_min/sigma_max = %.3e; groove content = %.6f\n', ...
        c0.sigma_ratio,c0.groove_content);
    fprintf('  start finite-channel radiation fraction = %.3e\n', ...
        c0.radiation_fraction);
    fprintf('  endpoint Omega0 = %.12f; endpoint sigma_min/sigma_max = %.3e\n', ...
        Omega0,track.endpoint.sigma_ratio);
    fprintf('  endpoint finite-channel radiation fraction = %.3e\n',endpointRad);
    fprintf('  final continued q = %+.12e%+.12ei; |q| = %.3e\n', ...
        real(qPath(end)),imag(qPath(end)),abs(qPath(end)));
    fprintf('  final continued Omega = %+.12e%+.12ei; |Omega-Omega0| = %.3e\n', ...
        real(OmegaPath(end)),imag(OmegaPath(end)),abs(OmegaPath(end)-Omega0));
    fprintf('  min consecutive right-mode overlap = %.6f\n', ...
        min(overlapOut(2:end-1)));
    fprintf('  max consecutive right-mode jump flag = %d\n',any(jumpOut));
    fprintf('  same nontrivial branch reaches endpoint = %d\n',track.same_branch);
end

    function [op,om,z,u,v,srel]=mode_at_q(q,kap)
        om=omega_from_q(q,kap);
        op=ni2019_full_eigen_operator_complex(cfg,2*pi*om/a, ...
            2*pi*kap/a,target,2*pi*q/a);
        [U,S,V]=svd(op.Fscaled,'econ');
        s=diag(S); srel=s(end)/max(s(1),eps);
        v=V(:,end); u=U(:,end);
        z=v./transpose(op.column_scale); z=z/norm(z);
    end

    function c=refine_pole(kap,qSeed,uRef,vRef,nOuter)
        x=[real(qSeed),imag(qSeed)]; flag=0;
        % Rebuild a local projected characteristic equation after each SVD
        % update.  This is a two-dimensional complex-root solve in q and
        % does not use the endpoint singular vectors.
        for no=1:nOuter
            fun=@(x) projected_equation(x,kap,uRef,vRef);
            opts=optimoptions('fsolve','Display',p.Results.Display, ...
                'FunctionTolerance',1e-12,'StepTolerance',1e-12, ...
                'OptimalityTolerance',1e-12,'MaxFunctionEvaluations',1500);
            [x,~,flag]=fsolve(fun,x,opts);
            q=complex(x(1),x(2));
            [~,~,~,uNew,vNew,srel]=mode_at_q(q,kap);
            if no>1 && srel<=1e-10, break; end
            uRef=phase_align(uNew,uRef);
            vRef=phase_align(vNew,vRef);
        end
        q=complex(x(1),x(2));
        [op,om,z,u,v,srel]=mode_at_q(q,kap);
        ss=svd(op.Fscaled);
        c=struct('q',q,'Omega',om,'z',z,'u',u,'v',v, ...
            'sigma_min',ss(end),'sigma_ratio',srel, ...
            'residual',norm(op.Fscaled*v),'groove_content', ...
            groove_content_from_mode(z,op),'radiation_fraction', ...
            radiation_fraction_from_mode(z,op),'exitflag',flag);
    end

    function y=projected_equation(x,kap,uRef,vRef)
        q=complex(x(1),x(2)); om=omega_from_q(q,kap);
        op=ni2019_full_eigen_operator_complex(cfg,2*pi*om/a, ...
            2*pi*kap/a,target,2*pi*q/a);
        g=uRef'*op.Fscaled*vRef;
        y=[real(g);imag(g)];
    end

    function om=omega_from_q(q,kap)
        om=sqrt((kap+target)^2+q^2);
        if real(om)<0, om=-om; end
    end

    function z=phase_align(z,ref)
        ov=ref'*z;
        if ov~=0, z=z*exp(-1i*angle(ov)); end
    end

    function g=groove_content_from_mode(z,op)
        A=z(1:op.N); C=z(op.N+1:end);
        g=sum(abs(C).^2)/max(sum(abs(A).^2)+sum(abs(C).^2),eps);
    end

    function r=radiation_fraction_from_mode(z,op)
        A=z(1:op.N);
        % At a complex pole use the real-axis channel classification at the
        % same kappa/Omega point.  This retains the finite-flux n=0 channel
        % while excluding the analytically continued Rayleigh channel.
        finiteOpen=real(op.k0)^2-op.kx.^2>1e-9*abs(op.k0)^2;
        r=sum(abs(A(finiteOpen)).^2)/max(sum(abs(A).^2),eps);
    end

    function label=sheet_label(q,om)
        if imag(q)<-1e-9 && imag(om)>1e-12
            label='physical outgoing continuation (exp(+j*omega*t))';
        elseif imag(q)>1e-9 && imag(om)<-1e-12
            label='time-reversed / opposite pole sheet';
        else
            label='near Rayleigh branch point';
        end
    end

end
