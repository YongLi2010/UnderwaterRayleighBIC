function result = ni2019_optimize_three_groove_double_rayleigh_bic(cfg,varargin)
%NI2019_OPTIMIZE_THREE_GROOVE_DOUBLE_RAYLEIGH_BIC Gamma double-RA search.
%
%   This routine fixes the normalized frequency and Bloch wave number at
%       Omega = k0*a/(2*pi) = 1,   kappa = kBloch*a/(2*pi) = 0.
%   Thus n=-1, 0, +1 are the three Gamma-point Rayleigh/open channels.
%   The homogeneous modal-matching matrix is formed by
%   NI2019_FULL_EIGEN_OPERATOR, then the A_-1, A_0 and A_+1 columns are
%   deleted before row/column equilibration and SVD.  A nonzero null vector
%   of this reduced matrix is therefore a strict three-channel candidate.
%
%   Full asymmetric parameterization (default):
%       x = [d1/a,d2/a,d3/a,w1/a,w2/a,w3/a,g1/a,g2/a].
%   Optional mirror-symmetric parameterization uses
%       y = [dOuter/a,dCenter/a,wOuter/a,wCenter/a,gOuter/a],
%   expanded internally to
%       x = [dOuter,dCenter,dOuter,wOuter,wCenter,wOuter,gOuter,gOuter].
%   Use 'Parameterization','both' to run both searches and return the
%   better candidate while retaining both result structures.
%
%   The objective is the worst strict sigma_min/sigma_max over the
%   requested (N,K) truncations plus a small mean term.  A penalty prevents
%   a singular vector with negligible groove-mode participation.  Geometry
%   is constrained by sum(widths)+sum(gaps)<=FillMax.  No files are written.
%
%   Example:
%       cfg=struct('a',1,'depths',[.25 .45 .25], ...
%           'widths',[.12 .20 .12],'gaps',[.12 .12]);
%       r=ni2019_optimize_three_groove_double_rayleigh_bic(cfg, ...
%           'Truncations',[15 3;25 5],'Starts',3, ...
%           'GlobalSwarmSize',12,'GlobalMaxIterations',6, ...
%           'MaxFunctionEvaluations',250);

p=inputParser;
addParameter(p,'Truncations',[31 5;41 7;61 9]);
addParameter(p,'Parameterization','full');
addParameter(p,'DepthRange',[.01 .82]);
addParameter(p,'WidthRange',[.04 .65]);
addParameter(p,'GapRange',[.02 .30]);
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
addParameter(p,'RandomSeed',37);
addParameter(p,'Display','off');
addParameter(p,'UseParallel',false);
parse(p,varargin{:});
opt=p.Results;
validate_inputs(cfg,opt);

family=lower(char(opt.Parameterization));
if strcmp(family,'both')
    optS=opt;
    optS.Parameterization='symmetric';
    optS.RandomSeed=opt.RandomSeed;
    symmetric=search_one(cfg,optS);
    optF=opt;
    optF.Parameterization='full';
    optF.RandomSeed=opt.RandomSeed+1009;
    full=search_one(cfg,optF);
    if full.score<=symmetric.score
        result=full;
        result.best_family='full';
    else
        result=symmetric;
        result.best_family='symmetric';
    end
    result.parameterization='both';
    result.symmetric_result=symmetric;
    result.full_result=full;
else
    if ~ismember(family,{'full','symmetric'})
        error('Parameterization must be full, symmetric, or both.');
    end
    opt.Parameterization=family;
    result=search_one(cfg,opt);
end
end

function result=search_one(cfg,opt)
family=lower(char(opt.Parameterization));
if strcmp(family,'full')
    nvar=8;
    lb=[repmat(opt.DepthRange(1),1,3),repmat(opt.WidthRange(1),1,3), ...
        repmat(opt.GapRange(1),1,2)];
    ub=[repmat(opt.DepthRange(2),1,3),repmat(opt.WidthRange(2),1,3), ...
        repmat(opt.GapRange(2),1,2)];
    occIdx=4:8;
    occWeights=ones(1,5);
    geometryIdx=4:8;
else
    nvar=5;
    lb=[opt.DepthRange(1),opt.DepthRange(1), ...
        opt.WidthRange(1),opt.WidthRange(1),opt.GapRange(1)];
    ub=[opt.DepthRange(2),opt.DepthRange(2), ...
        opt.WidthRange(2),opt.WidthRange(2),opt.GapRange(2)];
    occIdx=3:5;
    occWeights=[2,1,2];
    geometryIdx=3:5;
end
Aineq=zeros(1,nvar);
Aineq(occIdx)=occWeights;
bineq=opt.FillMax;
xBaseFull=base_full_design(cfg,cfg.a,family);
yBase=xBaseFull;
yBase=min(max(yBase,lb),ub);
yBase=make_feasible(yBase,lb,ub,occIdx,occWeights,opt.FillMax);

rng(opt.RandomSeed,'twister');
startPool=zeros(opt.Starts,nvar);
startPool(1,:)=yBase;
for ii=2:opt.Starts
    startPool(ii,:)=random_feasible(lb,ub,occIdx,occWeights,opt.FillMax);
end

historyY=zeros(0,nvar);
historyScore=zeros(0,1);
historyRaw=zeros(0,1);
historyMinGroove=zeros(0,1);
[baselineScore,baselineDiag]=objective(yBase);

globalInfo=struct('used',false,'method','','x',[],'x_full',[], ...
    'score',[],'exitflag',[],'output',struct());
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
    [yg,fg,efg,og]=particleswarm(@objective,nvar,lb,ub,gopts);
    yg=make_feasible(yg,lb,ub,occIdx,occWeights,opt.FillMax);
    globalInfo=struct('used',true,'method',opt.GlobalMethod,'x',yg, ...
        'x_full',expand_design(yg,family),'score',fg, ...
        'exitflag',efg,'output',og);
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
localSolutions=nan(opt.Starts,nvar);
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
xFull=expand_design(y,family);
result=struct();
result.x=xFull;
result.search_vector=y;
result.parameterization=family;
result.score=score;
result.strict_residual=diagnostics.max_residual;
result.aggregate_residual=diagnostics.aggregate_residual;
result.kappa=0;
result.Omega=1;
result.depths_over_a=xFull(1:3);
result.widths_over_a=xFull(4:6);
result.gaps_over_a=xFull(7:8);
result.fill_fraction=sum(xFull(4:8));
result.removed_orders=[-1 0 1];
result.baseline=struct('x',expand_design(yBase,family), ...
    'search_vector',yBase,'score',baselineScore, ...
    'strict_residual',baselineDiag.max_residual, ...
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
result.history=struct('search_vector',historyY,'score',historyScore, ...
    'strict_residual',historyRaw,'min_groove_fraction',historyMinGroove);

    function [score,dg,modesOut]=objective(yIn)
        if numel(yIn)~=nvar || any(~isfinite(yIn)) || ...
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
        x=expand_design(yIn,family);
        nT=size(opt.Truncations,1);
        residuals=zeros(nT,1);
        grooveFractions=zeros(nT,1);
        sigmaMin=zeros(nT,1);
        sigmaMax=zeros(nT,1);
        finiteOrders=cell(nT,1);
        operators=cell(nT,1);
        modesOut=cell(nT,1);
        for it=1:nT
            local=cfg;
            local.a=cfg.a;
            local.lambda=cfg.a;
            local.theta_i_deg=0;
            local.depths=x(1:3)*cfg.a;
            local.widths=x(4:6)*cfg.a;
            local.gaps=x(7:8)*cfg.a;
            local.N=opt.Truncations(it,1);
            local.K=opt.Truncations(it,2);
            local.solve_scattering=false;
            out=ni2019_full_eigen_operator(local);
            target=(out.orders==-1 | out.orders==0 | out.orders==1);
            keepA=~target;
            keep=[find(keepA);(out.n_floquet+1: ...
                out.n_floquet+out.n_groove).'];
            Fred=out.F(:,keep);
            [FredScaled,rowScale,colScale]=equilibrate(Fred);
            [~,S,V]=svd(FredScaled,'econ');
            sv=diag(S);
            if isempty(sv) || any(~isfinite(sv))
                residuals(it)=1e3;
                sigmaMin(it)=1e3;
                sigmaMax(it)=1;
                grooveFractions(it)=0;
                finiteOrders{it}=out.orders([]);
                operators{it}=out;
                modesOut{it}=[];
                continue;
            end
            sigmaMin(it)=sv(end);
            sigmaMax(it)=max(sv(1),eps);
            residuals(it)=sigmaMin(it)/sigmaMax(it);
            vred=V(:,end)./colScale(:);
            vfull=zeros(size(out.F,2),1);
            vfull(keep)=vred;
            vfull=vfull/max(norm(vfull),eps);
            Ccoef=vfull(out.n_floquet+1:end);
            grooveFractions(it)=sum(abs(Ccoef).^2)/ ...
                max(sum(abs(vfull).^2),eps);
            finite=abs(imag(out.ky))<1e-8*out.k0 & ...
                real(out.ky)>1e-7*out.k0;
            finiteOrders{it}=out.orders(finite);
            modesOut{it}=struct('v',vfull,'A',vfull(1:out.n_floquet), ...
                'C',Ccoef,'row_scale',rowScale, ...
                'column_scale_reduced',colScale,'kept_columns',keep, ...
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
            lowerFraction=(yIn(geometryIdx)-lb(geometryIdx))./ ...
                max(ub(geometryIdx)-lb(geometryIdx),eps);
            boundaryPenalty=opt.BoundaryPenaltyWeight*mean( ...
                max((opt.BoundaryMargin-lowerFraction)/ ...
                opt.BoundaryMargin,0).^2);
        else
            boundaryPenalty=0;
        end
        score=aggregate+groovePenalty+boundaryPenalty;
        dg=struct('x',x,'search_vector',yIn,'residuals',residuals, ...
            'max_residual',robust,'aggregate_residual',aggregate, ...
            'sigma_min',sigmaMin,'sigma_max',sigmaMax, ...
            'groove_fractions',grooveFractions, ...
            'min_groove_fraction',minGroove, ...
            'groove_penalty',groovePenalty, ...
            'boundary_penalty',boundaryPenalty, ...
            'finite_open_orders',{finiteOrders}, ...
            'removed_orders',[-1 0 1],'operators',{operators}, ...
            'fill_fraction',sum(x(4:8)),'kappa',0,'Omega',1);
        record(score,dg);
    end

    function record(score,dg)
        historyY(end+1,:)=dg.search_vector(:).';
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
if numel(cfg.depths)~=3 || numel(cfg.widths)~=3 || numel(cfg.gaps)~=2
    error('cfg must contain depths(3), widths(3), and gaps(2).');
end
if ~isnumeric(opt.Truncations) || size(opt.Truncations,2)~=2 || ...
        any(~isfinite(opt.Truncations(:))) || any(opt.Truncations(:)<=0)
    error('Truncations must be a finite positive M-by-2 array.');
end
if any(mod(opt.Truncations(:,1),2)~=1) || any(opt.Truncations(:,1)<3)
    error('Every N must be odd and at least 3.');
end
if any(mod(opt.Truncations(:,2),1)~=0) || any(opt.Truncations(:,2)<1)
    error('Every K must be a positive integer.');
end
if ~ismember(lower(char(opt.Parameterization)),{'full','symmetric','both'})
    error('Parameterization must be full, symmetric, or both.');
end
check_range(opt.DepthRange,'DepthRange',true);
check_range(opt.WidthRange,'WidthRange',true);
check_range(opt.GapRange,'GapRange',true);
if opt.FillMax<=0 || ~isfinite(opt.FillMax)
    error('FillMax must be positive and finite.');
end
if opt.FillMax<=3*opt.WidthRange(1)+2*opt.GapRange(1)
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

function x=base_full_design(cfg,a,family)
x=[cfg.depths(:).'/a,cfg.widths(:).'/a,cfg.gaps(:).'/a];
if strcmp(family,'symmetric')
    x=[mean([x(1),x(3)]),x(2),mean([x(4),x(6)]),x(5),mean(x(7:8))];
else
    return;
end
end

function x=expand_design(y,family)
if strcmp(family,'full')
    x=y(:).';
else
    x=[y(1),y(2),y(1),y(3),y(4),y(3),y(5),y(5)];
end
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
dg=struct('search_vector',y,'x',y,'reason',reason,'residuals',inf, ...
    'max_residual',inf,'aggregate_residual',inf,'sigma_min',inf, ...
    'sigma_max',1,'groove_fractions',0,'min_groove_fraction',0, ...
    'groove_penalty',inf,'finite_open_orders',{{}}, ...
    'removed_orders',[-1 0 1],'operators',{{}});
end
