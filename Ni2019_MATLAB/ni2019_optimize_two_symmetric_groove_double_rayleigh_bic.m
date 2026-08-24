function result = ni2019_optimize_two_symmetric_groove_double_rayleigh_bic(cfg,varargin)
%NI2019_OPTIMIZE_TWO_SYMMETRIC_GROOVE_DOUBLE_RAYLEIGH_BIC Odd Gamma BIC.
%
%   The structure is two identical rectangular grooves in one centered
%   period.  The design vector is
%       y = [d/a,w/a,g/a],
%   with expanded geometry
%       depths=[d,d], widths=[w,w], gaps=g.
%   The period and frequency are fixed at the double-Rayleigh Gamma point:
%       kappa=0, Omega=k0*a/(2*pi)=1, lambda=a.
%
%   The optimizer uses an explicit mirror-odd parity map.  In the Fourier
%   basis of NI2019_FULL_EIGEN_OPERATOR, A_-n=-A_n and A_0=0.  For the
%   groove cosine modes, C_(right,q)=-(-1)^q C_(left,q), because reflection
%   reverses the local coordinate of the second aperture.  The A_-1,
%   A_0 and A_+1 columns are absent from this parity map, so all three
%   Gamma-point Rayleigh channels are constrained to zero before SVD.
%
%   The objective is the worst equilibrated sigma_min/sigma_max over the
%   requested (N,K) truncations, plus a mean term and barriers against
%   weak groove participation or lower-bound geometry collapse.  K>=2 and
%   N>=5 are required so that the retained n=+/-2 evanescent orders are
%   represented.  No files are written by this function.
%
%   Required cfg fields:
%       cfg.a, cfg.depths (2), cfg.widths (2), cfg.gaps (scalar)
%
%   Example:
%       cfg=struct('a',1,'depths',[.4 .4],'widths',[.16 .16],'gaps',.18);
%       r=ni2019_optimize_two_symmetric_groove_double_rayleigh_bic(cfg, ...
%           'Truncations',[15 3;25 5],'Starts',3, ...
%           'GlobalSwarmSize',12,'GlobalMaxIterations',6, ...
%           'MaxFunctionEvaluations',250);

p=inputParser;
addParameter(p,'Truncations',[31 5;41 7;61 9]);
addParameter(p,'DepthRange',[.03 1.20]);
addParameter(p,'WidthRange',[.05 .48]);
addParameter(p,'GapRange',[.02 .75]);
addParameter(p,'FillMax',.94);
addParameter(p,'Starts',6);
addParameter(p,'UseGlobal',true);
addParameter(p,'GlobalMethod','particleswarm');
addParameter(p,'GlobalSwarmSize',24);
addParameter(p,'GlobalMaxIterations',18);
addParameter(p,'GlobalMaxStallIterations',8);
addParameter(p,'GlobalFunctionTolerance',1e-8);
addParameter(p,'MaxFunctionEvaluations',1200);
addParameter(p,'MaxIterations',300);
addParameter(p,'RobustMeanWeight',.15);
addParameter(p,'MinGrooveFraction',.05);
addParameter(p,'GroovePenaltyWeight',1.0);
addParameter(p,'BoundaryMargin',.05);
addParameter(p,'BoundaryPenaltyWeight',.20);
addParameter(p,'RandomSeed',53);
addParameter(p,'Display','off');
addParameter(p,'UseParallel',false);
parse(p,varargin{:});
opt=p.Results;
validate_inputs(cfg,opt);

a=cfg.a;
lb=[opt.DepthRange(1),opt.WidthRange(1),opt.GapRange(1)];
ub=[opt.DepthRange(2),opt.WidthRange(2),opt.GapRange(2)];
occIdx=2:3;
occWeights=[2,1];
Aineq=zeros(1,3);
Aineq(occIdx)=occWeights;
bineq=opt.FillMax;

yBase=[mean(cfg.depths(:))/a,mean(cfg.widths(:))/a,mean(cfg.gaps(:))/a];
yBase=min(max(yBase,lb),ub);
yBase=make_feasible(yBase,lb,ub,occIdx,occWeights,opt.FillMax);

rng(opt.RandomSeed,'twister');
startPool=zeros(opt.Starts,3);
startPool(1,:)=yBase;
for ii=2:opt.Starts
    startPool(ii,:)=random_feasible(lb,ub,occIdx,occWeights,opt.FillMax);
end

historyY=zeros(0,3);
historyScore=zeros(0,1);
historyRaw=zeros(0,1);
historyMinGroove=zeros(0,1);

[baselineScore,baselineDiag]=objective(yBase);
globalInfo=struct('used',false,'method','','y',[],'score',[], ...
    'exitflag',[],'output',struct());
if opt.UseGlobal
    if ~strcmpi(opt.GlobalMethod,'particleswarm')
        error('Unsupported GlobalMethod "%s". Use particleswarm or false.',opt.GlobalMethod);
    end
    if exist('particleswarm','file')~=2
        error('particleswarm is unavailable; set UseGlobal=false or install Global Optimization Toolbox.');
    end
    gopts=optimoptions('particleswarm','Display',opt.Display, ...
        'SwarmSize',max(4,round(opt.GlobalSwarmSize)), ...
        'MaxIterations',max(1,round(opt.GlobalMaxIterations)), ...
        'MaxStallIterations',max(1,round(opt.GlobalMaxStallIterations)), ...
        'FunctionTolerance',opt.GlobalFunctionTolerance, ...
        'UseParallel',logical(opt.UseParallel));
    try
        gopts.InitialPoints=startPool(1:min(size(startPool,1), ...
            max(1,round(opt.GlobalSwarmSize))),:);
    catch
        % InitialPoints is a convenience only; older releases can omit it.
    end
    [yg,fg,efg,og]=particleswarm(@objective,3,lb,ub,gopts);
    yg=make_feasible(yg,lb,ub,occIdx,occWeights,opt.FillMax);
    globalInfo=struct('used',true,'method',opt.GlobalMethod,'y',yg, ...
        'score',fg,'exitflag',efg,'output',og);
    if opt.Starts>=2
        startPool(2,:)=yg;
    else
        startPool(1,:)=yg;
    end
end

if exist('fmincon','file')~=2
    error('fmincon is unavailable; Optimization Toolbox is required.');
end
fopts=optimoptions('fmincon','Algorithm','sqp','Display',opt.Display, ...
    'MaxFunctionEvaluations',max(50,round(opt.MaxFunctionEvaluations)), ...
    'MaxIterations',max(20,round(opt.MaxIterations)), ...
    'OptimalityTolerance',1e-10,'StepTolerance',1e-10, ...
    'ConstraintTolerance',1e-10);
localSolutions=nan(opt.Starts,3);
localScores=inf(opt.Starts,1);
localRaw=inf(opt.Starts,1);
localExitflags=nan(opt.Starts,1);
localOutputs=cell(opt.Starts,1);
localDiagnostics=cell(opt.Starts,1);
for ii=1:opt.Starts
    y0=make_feasible(startPool(ii,:),lb,ub,occIdx,occWeights,opt.FillMax);
    [localSolutions(ii,:),localScores(ii),localExitflags(ii),localOutputs{ii}]= ...
        fmincon(@objective,y0,Aineq,bineq,[],[],lb,ub,[],fopts);
    localSolutions(ii,:)=make_feasible(localSolutions(ii,:),lb,ub, ...
        occIdx,occWeights,opt.FillMax);
    [~,localDiagnostics{ii}]=objective(localSolutions(ii,:));
    localRaw(ii)=localDiagnostics{ii}.max_residual;
end

[~,bestId]=min(localScores);
y=localSolutions(bestId,:);
[score,diagnostics,modes]=objective(y);
xFull=expand_geometry(y,a);
expandedCfg=cfg;
expandedCfg.a=a;
expandedCfg.lambda=a;
expandedCfg.theta_i_deg=0;
expandedCfg.depths=xFull.depths;
expandedCfg.widths=xFull.widths;
expandedCfg.gaps=xFull.gaps;
expandedCfg.solve_scattering=false;

result=struct();
result.y=y;
result.x=[y(1),y(1),y(2),y(2),y(3)];
result.score=score;
result.strict_residual=diagnostics.max_residual;
result.aggregate_residual=diagnostics.aggregate_residual;
result.kappa=0;
result.Omega=1;
result.lambda_over_a=1;
result.depth_over_a=y(1);
result.width_over_a=y(2);
result.gap_over_a=y(3);
result.fill_fraction=2*y(2)+y(3);
result.parity='odd';
result.removed_orders=[-1 0 1];
result.expanded_geometry=xFull;
result.expanded_cfg=expandedCfg;
result.baseline=struct('y',yBase,'x', ...
    [yBase(1),yBase(1),yBase(2),yBase(2),yBase(3)], ...
    'score',baselineScore,'strict_residual',baselineDiag.max_residual, ...
    'diagnostics',baselineDiag);
result.improvement_factor=baselineDiag.max_residual/max(result.strict_residual,eps);
result.truncations=opt.Truncations;
result.diagnostics=diagnostics;
result.modes=modes;
result.global=globalInfo;
result.local=struct('solutions',localSolutions,'scores',localScores, ...
    'strict_residuals',localRaw,'exitflags',localExitflags, ...
    'outputs',{localOutputs},'diagnostics',{localDiagnostics});
result.options=opt;
result.bounds=struct('lb',lb,'ub',ub,'FillMax',opt.FillMax, ...
    'occupancy_indices',occIdx,'occupancy_weights',occWeights);
result.history=struct('y',historyY,'score',historyScore, ...
    'strict_residual',historyRaw,'min_groove_fraction',historyMinGroove);

    function [score,dg,modesOut]=objective(yIn)
        if numel(yIn)~=3 || any(~isfinite(yIn)) || ...
                any(yIn<lb) || any(yIn>ub)
            score=1e3;
            dg=invalid_diag(yIn,'bounds');
            modesOut={};
            record(score,dg);
            return;
        end
        occupancy=sum(occWeights.*yIn(occIdx));
        if occupancy>opt.FillMax
            score=10+100*(occupancy-opt.FillMax)^2;
            dg=invalid_diag(yIn,'occupancy');
            dg.fill_fraction=occupancy;
            modesOut={};
            record(score,dg);
            return;
        end
        x=expand_geometry(yIn,a);
        nT=size(opt.Truncations,1);
        residuals=zeros(nT,1);
        rawResiduals=zeros(nT,1);
        grooveFractions=zeros(nT,1);
        highOrderFractions=zeros(nT,1);
        targetAmplitudeMax=zeros(nT,1);
        parityResiduals=zeros(nT,1);
        sigmaMin=zeros(nT,1);
        sigmaMax=zeros(nT,1);
        finiteOrders=cell(nT,1);
        operators=cell(nT,1);
        modesOut=cell(nT,1);
        for it=1:nT
            local=cfg;
            local.a=a;
            local.lambda=a;
            local.theta_i_deg=0;
            local.depths=x.depths;
            local.widths=x.widths;
            local.gaps=x.gaps;
            local.N=opt.Truncations(it,1);
            local.K=opt.Truncations(it,2);
            local.solve_scattering=false;
            out=ni2019_full_eigen_operator(local);
            T=odd_parity_map(local.N,local.K);
            Fmap=out.F*T;
            [Fscaled,rowScale,colScale]=equilibrate(Fmap);
            [~,S,V]=svd(Fscaled,'econ');
            sv=diag(S);
            if isempty(sv) || any(~isfinite(sv))
                residuals(it)=1e3;
                rawResiduals(it)=1e3;
                sigmaMin(it)=1e3;
                sigmaMax(it)=1;
                grooveFractions(it)=0;
                highOrderFractions(it)=0;
                parityResiduals(it)=1;
                finiteOrders{it}=out.orders([]);
                operators{it}=out;
                modesOut{it}=[];
                continue;
            end
            sigmaMin(it)=sv(end);
            sigmaMax(it)=max(sv(1),eps);
            residuals(it)=sigmaMin(it)/sigmaMax(it);
            vred=V(:,end)./colScale(:);
            vfull=T*vred;
            vfull=vfull/max(norm(vfull),eps);
            rawResiduals(it)=norm(out.F*vfull)/ ...
                max(norm(out.F,'fro')*norm(vfull),eps);
            Acoef=vfull(1:out.n_floquet);
            Ccoef=vfull(out.n_floquet+1:end);
            target=(abs(out.orders)<=1);
            targetAmplitudeMax(it)=max(abs(Acoef(target)));
            grooveFractions(it)=sum(abs(Ccoef).^2)/ ...
                max(sum(abs(vfull).^2),eps);
            retained=(abs(out.orders)>=2);
            highOrderFractions(it)=sum(abs(Acoef(retained)).^2)/ ...
                max(sum(abs(vfull).^2),eps);
            parityResiduals(it)=parity_error(vfull,out.N,out.K);
            finite=abs(imag(out.ky))<1e-8*out.k0 & ...
                real(out.ky)>1e-7*out.k0;
            finiteOrders{it}=out.orders(finite);
            modesOut{it}=struct('v',vfull,'A',Acoef,'C',Ccoef, ...
                'independent',vred,'T',T,'row_scale',rowScale, ...
                'column_scale_reduced',colScale,'raw_residual',rawResiduals(it), ...
                'sigma_ratio',residuals(it),'parity_error',parityResiduals(it), ...
                'high_order_fraction',highOrderFractions(it), ...
                'groove_fraction',grooveFractions(it), ...
                'target_amplitude_max',targetAmplitudeMax(it), ...
                'removed_orders',[-1 0 1]);
            operators{it}=out;
        end
        robust=max(residuals);
        aggregate=robust+opt.RobustMeanWeight*mean(residuals);
        minGroove=min(grooveFractions);
        groovePenalty=opt.GroovePenaltyWeight* ...
            (max(opt.MinGrooveFraction-minGroove,0)/ ...
            max(opt.MinGrooveFraction,eps))^2;
        if opt.BoundaryMargin>0 && opt.BoundaryPenaltyWeight>0
            lowerFraction=(yIn(2:3)-lb(2:3))./ ...
                max(ub(2:3)-lb(2:3),eps);
            boundaryPenalty=opt.BoundaryPenaltyWeight*mean( ...
                max((opt.BoundaryMargin-lowerFraction)/ ...
                opt.BoundaryMargin,0).^2);
        else
            boundaryPenalty=0;
        end
        score=aggregate+groovePenalty+boundaryPenalty;
        dg=struct('y',yIn,'x',x,'residuals',residuals, ...
            'raw_residuals',rawResiduals,'max_residual',robust, ...
            'aggregate_residual',aggregate,'sigma_min',sigmaMin, ...
            'sigma_max',sigmaMax,'groove_fractions',grooveFractions, ...
            'high_order_fractions',highOrderFractions, ...
            'target_amplitude_max',targetAmplitudeMax, ...
            'parity_residuals',parityResiduals, ...
            'min_groove_fraction',minGroove, ...
            'groove_penalty',groovePenalty, ...
            'boundary_penalty',boundaryPenalty, ...
            'finite_open_orders',{finiteOrders}, ...
            'removed_orders',[-1 0 1],'operators',{operators}, ...
            'fill_fraction',2*yIn(2)+yIn(3),'kappa',0,'Omega',1);
        record(score,dg);
    end

    function record(score,dg)
        historyY(end+1,:)=dg.y(:).';
        historyScore(end+1,1)=score;
        if isfield(dg,'max_residual')
            historyRaw(end+1,1)=dg.max_residual;
        else
            historyRaw(end+1,1)=score;
        end
        if isfield(dg,'min_groove_fraction')
            historyMinGroove(end+1,1)=dg.min_groove_fraction;
        else
            historyMinGroove(end+1,1)=0;
        end
    end
end

function validate_inputs(cfg,opt)
if ~isfield(cfg,'a') || ~isscalar(cfg.a) || ~isfinite(cfg.a) || cfg.a<=0
    error('cfg.a must be a positive finite scalar.');
end
if numel(cfg.depths)~=2 || numel(cfg.widths)~=2 || numel(cfg.gaps)~=1
    error('cfg must contain two depths, two widths, and one scalar gap.');
end
if ~isnumeric(opt.Truncations) || size(opt.Truncations,2)~=2 || ...
        any(~isfinite(opt.Truncations(:))) || any(opt.Truncations(:)<=0)
    error('Truncations must be a finite positive M-by-2 array.');
end
if any(mod(opt.Truncations(:,1),2)~=1) || any(opt.Truncations(:,1)<5)
    error('Every N must be odd and at least 5.');
end
if any(mod(opt.Truncations(:,2),1)~=0) || any(opt.Truncations(:,2)<2)
    error('Every K must be an integer at least 2.');
end
check_range(opt.DepthRange,'DepthRange',true);
check_range(opt.WidthRange,'WidthRange',true);
check_range(opt.GapRange,'GapRange',true);
if opt.FillMax<=0 || ~isfinite(opt.FillMax)
    error('FillMax must be positive and finite.');
end
if opt.FillMax<=2*opt.WidthRange(1)+opt.GapRange(1)
    error('FillMax is too small for the width/gap lower bounds.');
end
if opt.Starts<1 || opt.Starts~=round(opt.Starts)
    error('Starts must be a positive integer.');
end
if opt.MinGrooveFraction<0 || opt.MinGrooveFraction>=1
    error('MinGrooveFraction must lie in [0,1).');
end
if opt.BoundaryMargin<0 || opt.BoundaryMargin>=1
    error('BoundaryMargin must lie in [0,1).');
end
if opt.BoundaryPenaltyWeight<0 || ~isfinite(opt.BoundaryPenaltyWeight)
    error('BoundaryPenaltyWeight must be finite and nonnegative.');
end
end

function check_range(r,name,positive)
if ~isnumeric(r) || numel(r)~=2 || any(~isfinite(r)) || r(1)>=r(2)
    error('%s must be an ordered finite two-element range.',name);
end
if positive && r(1)<=0
    error('%s must be positive.',name);
end
end

function x=expand_geometry(y,a)
x=struct('depths',[y(1),y(1)]*a,'widths',[y(2),y(2)]*a, ...
    'gaps',y(3)*a);
end

function T=odd_parity_map(N,K)
P=(N-1)/2;
nIndependent=P-1+K;
T=complex(zeros(N+2*K,nIndependent));
column=0;
for n=2:P
    column=column+1;
    T(n+P+1,column)=1;
    T(-n+P+1,column)=-1;
end
for q=0:K-1
    column=column+1;
    T(N+q+1,column)=1;
    T(N+K+q+1,column)=-(-1)^q;
end
if column~=nIndependent
    error('Internal odd-parity map dimension mismatch.');
end
end

function err=parity_error(v,N,K)
A=v(1:N);
C=v(N+1:end);
Aref=flipud(A);
C1=C(1:K);
C2=C(K+1:2*K);
targetC=-((-1).^(0:K-1)).'.*C1;
errA=norm(A-Aref*(-1))/max(norm(A),eps);
errC=norm(C2-targetC)/max(norm(C),eps);
err=max(errA,errC);
end

function y=make_feasible(y,lb,ub,occIdx,occWeights,fillMax)
y=min(max(y(:).',lb),ub);
occ=sum(occWeights.*y(occIdx));
if occ>fillMax
    lower=lb(occIdx);
    free=sum(occWeights.*(y(occIdx)-lower));
    room=fillMax-sum(occWeights.*lower);
    if free>0 && room>0
        y(occIdx)=lower+(y(occIdx)-lower)*(room/free);
    else
        y(occIdx)=lower;
    end
end
if sum(occWeights.*y(occIdx))>fillMax
    y(occIdx)=y(occIdx)*(fillMax/max(sum(occWeights.*y(occIdx)),eps))*(1-10*eps);
end
end

function y=random_feasible(lb,ub,occIdx,occWeights,fillMax)
for attempt=1:1000
    y=lb+(ub-lb).*rand(1,numel(lb));
    if sum(occWeights.*y(occIdx))<=fillMax
        return;
    end
end
y=lb;
end

function [Fs,rowScale,colScale]=equilibrate(F)
rowScale=max(vecnorm(F,2,2),sqrt(eps));
Fr=F./rowScale;
colScale=max(vecnorm(Fr,2,1),sqrt(eps));
Fs=Fr./colScale;
end

function dg=invalid_diag(y,reason)
dg=struct('y',y,'x',y,'reason',reason,'residuals',inf, ...
    'raw_residuals',inf,'max_residual',inf,'aggregate_residual',inf, ...
    'sigma_min',inf,'sigma_max',1,'groove_fractions',0, ...
    'high_order_fractions',0,'parity_residuals',1, ...
    'target_amplitude_max',inf, ...
    'min_groove_fraction',0,'groove_penalty',inf, ...
    'boundary_penalty',inf,'finite_open_orders',{{}}, ...
    'removed_orders',[-1 0 1],'operators',{{}});
end
