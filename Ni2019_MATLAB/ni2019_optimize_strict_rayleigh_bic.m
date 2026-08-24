function result = ni2019_optimize_strict_rayleigh_bic(cfg,varargin)
%NI2019_OPTIMIZE_STRICT_RAYLEIGH_BIC Robust off-Gamma Rayleigh-BIC search.
%
%   The design vector is
%       x = [kappa,d1/a,d2/a,w1/a,w2/a,g/a],
%   and the n=-1 Rayleigh condition is imposed exactly through
%       Omega = k0*a/(2*pi) = 1-kappa.
%
%   This routine searches for a non-trivial null vector of the homogeneous
%   modal-matching operator after the A_-1 and A_0 columns have been
%   removed. Those two amplitudes are therefore constrained to zero before
%   the SVD is taken. The reduced matrix is row/column equilibrated after
%   the column deletion; using the scaling of the unreduced matrix can hide
%   the very constraint this optimizer is intended to enforce.
%
%   The objective is robust over the requested (N,K) truncations. It is
%   the worst strict singular-value ratio plus a small mean-ratio term,
%   with a barrier against modes having negligible groove content. The
%   evanescent A_+1 coefficient is retained and is not penalized.
%
%   The returned result contains the best parameters, per-truncation
%   residuals, modal diagnostics, global/local exit information, and the
%   sampled objective history. No files are written by this function.

p = inputParser;
addParameter(p,'Truncations',[31 5;41 7;61 9]);
addParameter(p,'KappaRange',[.08 .48]);
addParameter(p,'DepthRange',[.01 .82]);
addParameter(p,'WidthRange',[.025 .72]);
addParameter(p,'GapRange',[.01 .36]);
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
addParameter(p,'RandomSeed',19);
addParameter(p,'Display','off');
addParameter(p,'UseParallel',false);
parse(p,varargin{:});
opt = p.Results;

validateattributes(opt.Truncations,{'numeric'},{'2d','ncols',2,'finite','positive'},mfilename,'Truncations');
if any(mod(opt.Truncations(:,1),2)~=1)
    error('Every N in Truncations must be odd.');
end
if any(opt.Truncations(:,1)<3)
    error('Every N in Truncations must include at least the n=-1,0,+1 orders (N>=3).');
end
if any(opt.Truncations(:,2)<1 | mod(opt.Truncations(:,2),1)~=0)
    error('Every K in Truncations must be a positive integer.');
end
if numel(cfg.widths)~=2 || numel(cfg.depths)~=2 || numel(cfg.gaps)~=1
    error('The strict optimizer requires two groove widths, two depths, and one separating gap.');
end
if ~isfield(cfg,'a') || ~isscalar(cfg.a) || cfg.a<=0
    error('cfg.a must be a positive scalar.');
end
if opt.FillMax <= sum([opt.WidthRange(1),opt.WidthRange(1),opt.GapRange(1)])
    error('FillMax is too small for the specified lower width/gap bounds.');
end
if opt.KappaRange(1)<=0 || opt.KappaRange(2)>=.5 || opt.KappaRange(1)>=opt.KappaRange(2)
    error('KappaRange must satisfy 0<kappa_min<kappa_max<0.5 for a single n=-1 Rayleigh channel.');
end
if any(opt.DepthRange<=0) || opt.DepthRange(1)>=opt.DepthRange(2)
    error('DepthRange must be positive and ordered.');
end
if any(opt.WidthRange<=0) || opt.WidthRange(1)>=opt.WidthRange(2)
    error('WidthRange must be positive and ordered.');
end
if any(opt.GapRange<=0) || opt.GapRange(1)>=opt.GapRange(2)
    error('GapRange must be positive and ordered.');
end
if opt.Starts<1 || opt.Starts~=round(opt.Starts)
    error('Starts must be a positive integer.');
end

% Bounds and deterministic feasible starting points.
lb = [opt.KappaRange(1),repmat(opt.DepthRange(1),1,2), ...
    repmat(opt.WidthRange(1),1,2),opt.GapRange(1)];
ub = [opt.KappaRange(2),repmat(opt.DepthRange(2),1,2), ...
    repmat(opt.WidthRange(2),1,2),opt.GapRange(2)];
xBase = base_design(cfg,cfg.a,opt,lb,ub);
Aineq = zeros(1,6); Aineq(4:6)=1; bineq=opt.FillMax;

rng(opt.RandomSeed,'twister');
startPool = zeros(opt.Starts,6);
startPool(1,:) = xBase;
for i=2:opt.Starts
    startPool(i,:) = random_feasible(lb,ub,opt.FillMax);
end

% Nested objective history is deliberately compact: full mode vectors are
% retained only for the returned best candidate, not for every PS iteration.
historyX = zeros(0,6);
historyScore = zeros(0,1);
historyRaw = zeros(0,1);
historyMinGroove = zeros(0,1);

% Establish a comparable strict baseline before any optimization.
[baselineScore,baselineDiag] = objective(xBase);

globalInfo = struct('used',false,'method','','x',[],'score',[], ...
    'exitflag',[],'output',struct());
if opt.UseGlobal
    if ~strcmpi(opt.GlobalMethod,'particleswarm')
        error('Unsupported GlobalMethod "%s". Use "particleswarm" or set UseGlobal=false.',opt.GlobalMethod);
    end
    if exist('particleswarm','file')~=2
        error('particleswarm is unavailable; set UseGlobal=false or install Global Optimization Toolbox.');
    end
    gopts = optimoptions('particleswarm','Display',opt.Display, ...
        'SwarmSize',max(4,round(opt.GlobalSwarmSize)), ...
        'MaxIterations',max(1,round(opt.GlobalMaxIterations)), ...
        'MaxStallIterations',max(1,round(opt.GlobalMaxStallIterations)), ...
        'FunctionTolerance',opt.GlobalFunctionTolerance, ...
        'UseParallel',logical(opt.UseParallel));
    % InitialPoints is supported by current MATLAB releases. It is not
    % required for correctness; use it only when dimensions are compatible.
    try
        gopts.InitialPoints = startPool(1:min(size(startPool,1),max(1,round(opt.GlobalSwarmSize))),:);
    catch
        % Older releases silently omit this convenience initialization.
    end
    [xg,fg,efg,og] = particleswarm(@objective,6,lb,ub,gopts);
    globalInfo = struct('used',true,'method',opt.GlobalMethod,'x',xg, ...
        'score',fg,'exitflag',efg,'output',og);
    if opt.Starts>=2
        startPool(2,:) = make_feasible(xg,lb,ub,opt.FillMax);
    else
        startPool(1,:) = make_feasible(xg,lb,ub,opt.FillMax);
    end
end

% Local SQP polishing from deterministic base/global/random starts.
localSolutions = nan(opt.Starts,6);
localScores = inf(opt.Starts,1);
localRaw = inf(opt.Starts,1);
localExitflags = nan(opt.Starts,1);
localOutputs = cell(opt.Starts,1);
localDiagnostics = cell(opt.Starts,1);
if exist('fmincon','file')~=2
    error('fmincon is unavailable; Optimization Toolbox is required.');
end
fopts = optimoptions('fmincon','Algorithm','sqp','Display',opt.Display, ...
    'MaxFunctionEvaluations',max(50,round(opt.MaxFunctionEvaluations)), ...
    'MaxIterations',max(20,round(opt.MaxIterations)), ...
    'OptimalityTolerance',1e-10,'StepTolerance',1e-10, ...
    'ConstraintTolerance',1e-10);
for i=1:opt.Starts
    x0 = make_feasible(startPool(i,:),lb,ub,opt.FillMax);
    [localSolutions(i,:),localScores(i),localExitflags(i),localOutputs{i}] = ...
        fmincon(@objective,x0,Aineq,bineq,[],[],lb,ub,[],fopts);
    localSolutions(i,:) = make_feasible(localSolutions(i,:),lb,ub,opt.FillMax);
    [~,localDiagnostics{i}] = objective(localSolutions(i,:));
    localRaw(i) = localDiagnostics{i}.max_residual;
end

[~,bestId] = min(localScores);
x = localSolutions(bestId,:);
[score,diagnostics,modes] = objective(x);

result = struct();
result.x = x;
result.score = score;
result.strict_residual = diagnostics.max_residual;
result.aggregate_residual = diagnostics.aggregate_residual;
result.kappa = x(1);
result.Omega = 1-x(1);
result.theta_deg = asind(x(1)/(1-x(1)));
result.depths_over_a = x(2:3);
result.widths_over_a = x(4:5);
result.gap_over_a = x(6);
result.fill_fraction = sum(x(4:6));
result.baseline = struct('x',xBase,'score',baselineScore, ...
    'strict_residual',baselineDiag.max_residual,'diagnostics',baselineDiag);
result.improvement_factor = baselineDiag.max_residual/max(result.strict_residual,eps);
result.truncations = opt.Truncations;
result.diagnostics = diagnostics;
result.modes = modes;
result.global = globalInfo;
result.local = struct('solutions',localSolutions,'scores',localScores, ...
    'strict_residuals',localRaw,'exitflags',localExitflags, ...
    'outputs',{localOutputs},'diagnostics',{localDiagnostics});
result.options = opt;
result.bounds = struct('lb',lb,'ub',ub,'FillMax',opt.FillMax);
result.history = struct('x',historyX,'score',historyScore, ...
    'strict_residual',historyRaw,'min_groove_fraction',historyMinGroove);

    function [score,dg,modesOut] = objective(x)
        % Return a smooth finite penalty for PS points outside occupancy.
        if numel(x)~=6 || any(~isfinite(x)) || any(x<lb) || any(x>ub)
            score = 1e3; dg=invalid_diag(x,'bounds'); modesOut={}; record(score,dg); return;
        end
        occupancy = sum(x(4:6));
        if occupancy>opt.FillMax
            score = 10 + 100*(occupancy-opt.FillMax)^2;
            dg=invalid_diag(x,'occupancy'); dg.fill_fraction=occupancy;
            modesOut={}; record(score,dg); return;
        end
        nT=size(opt.Truncations,1);
        residuals=zeros(nT,1); grooveFractions=zeros(nT,1);
        sigmaMin=zeros(nT,1); sigmaMax=zeros(nT,1);
        finiteOrders=cell(nT,1); operators=cell(nT,1); modesOut=cell(nT,1);
        for it=1:nT
            local=cfg;
            Omega=1-x(1);
            local.lambda=cfg.a/Omega;
            local.theta_i_deg=asind(x(1)/Omega);
            local.depths=x(2:3)*cfg.a;
            local.widths=x(4:5)*cfg.a;
            local.gaps=x(6)*cfg.a;
            local.N=opt.Truncations(it,1);
            local.K=opt.Truncations(it,2);
            local.solve_scattering=false;
            out=ni2019_full_eigen_operator(local);
            target=(out.orders==-1 | out.orders==0);
            keepA=~target;
            keep=[find(keepA); (out.n_floquet+1:out.n_floquet+out.n_groove).'];
            Fred=out.F(:,keep);
            [FredScaled,rowScale,colScale]=equilibrate(Fred);
            [~,S,V]=svd(FredScaled,'econ');
            sv=diag(S);
            if isempty(sv) || any(~isfinite(sv))
                residuals(it)=1e3; sigmaMin(it)=1e3; sigmaMax(it)=1;
                grooveFractions(it)=0; modesOut{it}=[]; operators{it}=out;
                finiteOrders{it}=out.orders([]);
                continue;
            end
            sigmaMin(it)=sv(end); sigmaMax(it)=max(sv(1),eps);
            residuals(it)=sigmaMin(it)/sigmaMax(it);
            vred=V(:,end)./colScale(:);
            vfull=zeros(size(out.F,2),1); vfull(keep)=vred;
            vfull=vfull/max(norm(vfull),eps);
            Ccoef=vfull(out.n_floquet+1:end);
            grooveFractions(it)=sum(abs(Ccoef).^2)/max(sum(abs(vfull).^2),eps);
            finite=(abs(imag(out.ky))<1e-8*out.k0 & real(out.ky)>1e-7*out.k0);
            finiteOrders{it}=out.orders(finite);
            modesOut{it}=struct('v',vfull,'A',vfull(1:out.n_floquet), ...
                'C',Ccoef,'row_scale',rowScale,'column_scale_reduced',colScale, ...
                'kept_columns',keep,'target_orders',[-1 0]);
            operators{it}=out;
        end
        robust=max(residuals);
        aggregate=robust+opt.RobustMeanWeight*mean(residuals);
        minGroove=min(grooveFractions);
        groovePenalty=opt.GroovePenaltyWeight* ...
            (max(opt.MinGrooveFraction-minGroove,0)/max(opt.MinGrooveFraction,eps))^2;
        score=aggregate+groovePenalty;
        dg=struct('residuals',residuals,'max_residual',robust, ...
            'aggregate_residual',aggregate,'sigma_min',sigmaMin, ...
            'sigma_max',sigmaMax,'groove_fractions',grooveFractions, ...
            'min_groove_fraction',minGroove,'groove_penalty',groovePenalty, ...
            'finite_open_orders',{finiteOrders},'operators',{operators}, ...
            'fill_fraction',occupancy,'kappa',x(1),'Omega',1-x(1), ...
            'theta_deg',asind(x(1)/(1-x(1))), 'x',x);
        record(score,dg);
    end

    function record(score,dg)
        historyX(end+1,:) = dg.x(:).';
        historyScore(end+1,1) = score;
        if isfield(dg,'max_residual'), historyRaw(end+1,1)=dg.max_residual;
        else, historyRaw(end+1,1)=score; end
        if isfield(dg,'min_groove_fraction'), historyMinGroove(end+1,1)=dg.min_groove_fraction;
        else, historyMinGroove(end+1,1)=0; end
    end
end

function [Fs,rowScale,colScale] = equilibrate(F)
% Equilibrate only after strict columns have been deleted.
rowScale=max(vecnorm(F,2,2),sqrt(eps));
Fr=F./rowScale;
colScale=max(vecnorm(Fr,2,1),sqrt(eps));
Fs=Fr./colScale;
end

function x=base_design(cfg,a,opt,lb,ub)
if isfield(cfg,'kappa')
    kappa=cfg.kappa;
elseif isfield(cfg,'lambda') && isfield(cfg,'theta_i_deg')
    kappa=(a/cfg.lambda)*sind(cfg.theta_i_deg);
else
    kappa=.24;
end
x=[kappa,cfg.depths(:).'/a,cfg.widths(:).'/a,cfg.gaps(:).'/a];
x=min(max(x,lb),ub);
x=make_feasible(x,lb,ub,opt.FillMax);
end

function x=make_feasible(x,lb,ub,fillMax)
x=min(max(x(:).',lb),ub);
if sum(x(4:6))>fillMax
    lower=lb(4:6);
    free=sum(x(4:6))-sum(lower);
    room=fillMax-sum(lower);
    if free>0 && room>0
        x(4:6)=lower+(x(4:6)-lower)*(room/free);
    else
        x(4:6)=lower;
    end
end
% Roundoff protection for fmincon's linear inequality.
if sum(x(4:6))>fillMax
    x(4:6)=x(4:6)*(fillMax/sum(x(4:6)))*(1-10*eps);
end
end

function x=random_feasible(lb,ub,fillMax)
for attempt=1:1000
    x=lb+(ub-lb).*rand(1,numel(lb));
    if sum(x(4:6))<=fillMax, return; end
end
x=lb;
end

function dg=invalid_diag(x,reason)
dg=struct('x',x,'reason',reason,'residuals',inf, ...
    'max_residual',inf,'aggregate_residual',inf, ...
    'sigma_min',inf,'sigma_max',1,'groove_fractions',0, ...
    'min_groove_fraction',0,'groove_penalty',inf, ...
    'finite_open_orders',{{}},'operators',{{}});
end
