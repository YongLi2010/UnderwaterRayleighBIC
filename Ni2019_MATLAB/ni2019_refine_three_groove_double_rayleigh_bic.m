function result=ni2019_refine_three_groove_double_rayleigh_bic(cfg,x0,varargin)
%NI2019_REFINE_THREE_GROOVE_DOUBLE_RAYLEIGH_BIC Polish a Gamma-point root.
%   R=NI2019_REFINE_THREE_GROOVE_DOUBLE_RAYLEIGH_BIC(CFG,X0) refines
%
%       X=[d1/a,d2/a,d3/a,w1/a,w2/a,w3/a,g1/a,g2/a]
%
%   at the exact double-Rayleigh endpoint kappa=0, Omega=a/lambda=1.
%   Before each SVD, the exterior columns A_-1, A_0, and A_+1 are removed
%   from the raw homogeneous modal-matching operator.  The remaining tall
%   matrix is re-equilibrated, its smallest right singular vector is
%   eliminated analytically, and LSQNONLIN minimizes the complex
%   compatibility residual with respect to the eight real geometry
%   variables.
%
%   This function polishes one finite (N,K) discretization.  It does not by
%   itself prove continuum convergence; call it along an increasing
%   truncation sequence and cross-evaluate every returned geometry.
%
%   Name/value options:
%       Truncation              [N K] (default [cfg.N cfg.K])
%       LowerBounds             eight normalized lower bounds
%       UpperBounds             eight normalized upper bounds
%       FillMax                 max sum(widths)+sum(gaps), default .94
%       FillPenaltyWeight       default 1e4
%       MaxFunctionEvaluations  default 2200
%       MaxIterations           default 300
%       FunctionTolerance       default 1e-14
%       StepTolerance           default 1e-13
%       OptimalityTolerance     default 1e-13
%       FiniteDifferenceType    default 'central'
%       Display                 default 'off'

p=inputParser;
addParameter(p,'Truncation',[],@(x)isnumeric(x)&&isequal(size(x),[1 2]));
addParameter(p,'LowerBounds',[.01 .01 .01 .015 .015 .015 .01 .01], ...
    @(x)isnumeric(x)&&numel(x)==8&&all(isfinite(x)));
addParameter(p,'UpperBounds',[.95 .95 .95 .75 .75 .75 .40 .40], ...
    @(x)isnumeric(x)&&numel(x)==8&&all(isfinite(x)));
addParameter(p,'FillMax',.94,@(x)isscalar(x)&&isfinite(x)&&x>0);
addParameter(p,'FillPenaltyWeight',1e4,@(x)isscalar(x)&&isfinite(x)&&x>0);
addParameter(p,'MaxFunctionEvaluations',2200,@(x)isscalar(x)&&x>=10);
addParameter(p,'MaxIterations',300,@(x)isscalar(x)&&x>=1);
addParameter(p,'FunctionTolerance',1e-14,@(x)isscalar(x)&&x>0);
addParameter(p,'StepTolerance',1e-13,@(x)isscalar(x)&&x>0);
addParameter(p,'OptimalityTolerance',1e-13,@(x)isscalar(x)&&x>0);
addParameter(p,'FiniteDifferenceType','central',@(x)ischar(x)||isstring(x));
addParameter(p,'Display','off',@(x)ischar(x)||isstring(x));
parse(p,varargin{:});
opt=p.Results;

if ~isfield(cfg,'a')||~isscalar(cfg.a)||cfg.a<=0
    error('cfg.a must be a positive scalar.');
end
if numel(cfg.widths)~=3||numel(cfg.depths)~=3||numel(cfg.gaps)~=2
    error('cfg must describe exactly three grooves and two internal gaps.');
end
if isempty(opt.Truncation)
    if ~isfield(cfg,'N')||~isfield(cfg,'K')
        error('Supply Truncation or provide cfg.N and cfg.K.');
    end
    truncation=[cfg.N cfg.K];
else
    truncation=round(opt.Truncation);
end
if mod(truncation(1),2)~=1||truncation(2)<1
    error('Truncation must contain an odd N and a positive K.');
end
if exist('lsqnonlin','file')~=2
    error('lsqnonlin from Optimization Toolbox is required.');
end

lb=opt.LowerBounds(:).'; ub=opt.UpperBounds(:).'; x0=x0(:).';
if any(lb>=ub)||any(x0<lb)||any(x0>ub)
    error('Bounds must be ordered and contain x0.');
end
if sum(x0(4:8))>opt.FillMax
    error('Initial widths and gaps exceed FillMax.');
end

% Select a stable groove entry to fix the arbitrary complex SVD phase.
[~,initialState]=projected_residual(x0,[]);
grooveStart=initialState.n_kept_amplitudes+1;
[~,relativeAnchor]=max(abs(initialState.y(grooveStart:end)));
anchor=grooveStart+relativeAnchor-1;

historyX=zeros(0,8); historyNorm=zeros(0,1);
lsqopt=optimoptions('lsqnonlin','Display',opt.Display, ...
    'MaxFunctionEvaluations',round(opt.MaxFunctionEvaluations), ...
    'MaxIterations',round(opt.MaxIterations), ...
    'FunctionTolerance',opt.FunctionTolerance, ...
    'StepTolerance',opt.StepTolerance, ...
    'OptimalityTolerance',opt.OptimalityTolerance, ...
    'FiniteDifferenceType',char(opt.FiniteDifferenceType));

[x,resnorm,residual,exitflag,output]=lsqnonlin(@objective,x0,lb,ub,lsqopt);
[compatibility,state]=projected_residual(x,anchor);
strict=ni2019_strict_rayleigh_operator(state.local_cfg, ...
    'TargetOrder',-1,'EnforceOtherThreshold',true);

result=struct();
result.x=x;
result.depths_over_a=x(1:3);
result.widths_over_a=x(4:6);
result.gaps_over_a=x(7:8);
result.fill_fraction=sum(x(4:8));
result.kappa=0;
result.Omega=1;
result.theta_deg=0;
result.truncation=truncation;
result.compatibility_residual=compatibility;
result.compatibility_norm=norm(compatibility);
result.sigma_ratio=state.sigma_ratio;
result.sigma_min=state.sigma_min;
result.strict_operator=strict;
result.removed_orders=state.removed_orders;
result.resnorm=resnorm;
result.residual=residual;
result.exitflag=exitflag;
result.output=output;
result.anchor_index_reduced=anchor;
result.history=struct('x',historyX,'compatibility_norm',historyNorm);
result.options=opt;

    function r=objective(x)
        [rComplex,~]=projected_residual(x,anchor);
        fillViolation=max(sum(x(4:8))-opt.FillMax,0);
        r=[real(rComplex);imag(rComplex); ...
            sqrt(opt.FillPenaltyWeight)*fillViolation];
        historyX(end+1,:)=x(:).';
        historyNorm(end+1,1)=norm(r);
    end

    function [r,state]=projected_residual(x,anchorIndex)
        local=cfg;
        local.lambda=cfg.a;
        local.theta_i_deg=0;
        local.depths=x(1:3)*cfg.a;
        local.widths=x(4:6)*cfg.a;
        local.gaps=x(7:8)*cfg.a;
        local.N=truncation(1);
        local.K=truncation(2);
        local.solve_scattering=false;
        op=ni2019_full_eigen_operator(local);
        kyTolerance=1e-8*max(abs(op.k0),1);
        finiteOpen=abs(imag(op.ky))<=kyTolerance&real(op.ky)>kyTolerance;
        threshold=abs(op.ky)<=kyTolerance;
        remove=finiteOpen|threshold;
        removedOrders=op.orders(remove);
        if ~isequal(removedOrders(:),[-1;0;1])
            error('Expected removed orders [-1 0 1], obtained [%s].', ...
                strtrim(sprintf('%d ',removedOrders)));
        end
        keep=[find(~remove);(op.N+1:op.N+op.n_groove).'];
        Fred=op.F(:,keep);
        rowScale=max(vecnorm(Fred,2,2),sqrt(eps));
        Frow=Fred./rowScale;
        columnScale=max(vecnorm(Frow,2,1),sqrt(eps));
        Fscaled=Frow./columnScale;
        if any(~isfinite(Fscaled(:)))
            error('Reduced operator contains NaN or Inf.');
        end
        [~,S,V]=svd(Fscaled,'econ');
        singularValues=diag(S); y=V(:,end);
        if ~isempty(anchorIndex)
            if anchorIndex>numel(y)||abs(y(anchorIndex))<100*eps
                error('The fixed phase anchor became singular during refinement.');
            end
            y=y*exp(-1i*angle(y(anchorIndex)));
        end
        sigmaMax=max(singularValues(1),eps);
        r=(Fscaled*y)/sigmaMax;
        state=struct('local_cfg',local,'operator',op,'Fred',Fred, ...
            'Fscaled',Fscaled,'y',y,'singular_values',singularValues, ...
            'sigma_min',singularValues(end), ...
            'sigma_ratio',singularValues(end)/sigmaMax, ...
            'n_kept_amplitudes',nnz(~remove),'removed_orders',removedOrders);
    end
end
