function result = ni2019_refine_strict_rayleigh_bic(cfg,x0,varargin)
%NI2019_REFINE_STRICT_RAYLEIGH_BIC Variable-projection strict-BIC refinement.
%   RESULT = NI2019_REFINE_STRICT_RAYLEIGH_BIC(CFG,X0) refines
%
%       X = [kappa,d1/a,d2/a,w1/a,w2/a,g/a]
%
%   at the exact n=-1 Rayleigh condition Omega=1-kappa.  The exterior
%   amplitudes A_-1 and every finite-flux open A_n (A_0 in the off-Gamma
%   range used here) are removed before the solve.  At each geometry, the
%   smallest right singular vector of the equilibrated reduced operator is
%   eliminated analytically.  LSQNONLIN then minimizes the real and
%   imaginary parts of the remaining compatibility residual.
%
%   This variable-projection step is intended to polish one finite
%   truncation after a global/multistart search.  A sequence of increasing
%   (N,K) truncations must still be used to establish convergence.
%
%   Name/value options:
%       Truncation             [N K] (default [cfg.N cfg.K])
%       LowerBounds            six design lower bounds
%       UpperBounds            six design upper bounds
%       FillMax                maximum w1/a+w2/a+g/a (default .94)
%       FillPenaltyWeight      residual penalty (default 1e4)
%       MaxFunctionEvaluations default 1800
%       MaxIterations          default 250
%       FunctionTolerance      default 1e-14
%       StepTolerance          default 1e-13
%       OptimalityTolerance    default 1e-13
%       FiniteDifferenceType   default 'central'
%       Display                default 'off'

p = inputParser;
addParameter(p,'Truncation',[],@(x)isnumeric(x) && isequal(size(x),[1 2]));
addParameter(p,'LowerBounds',[.01 .01 .01 .01 .01 .005], ...
    @(x)isnumeric(x) && numel(x)==6 && all(isfinite(x)));
addParameter(p,'UpperBounds',[.49 .95 .95 .90 .90 .50], ...
    @(x)isnumeric(x) && numel(x)==6 && all(isfinite(x)));
addParameter(p,'FillMax',.94,@(x)isscalar(x) && isfinite(x) && x>0);
addParameter(p,'FillPenaltyWeight',1e4,@(x)isscalar(x) && isfinite(x) && x>0);
addParameter(p,'MaxFunctionEvaluations',1800,@(x)isscalar(x) && x>=10);
addParameter(p,'MaxIterations',250,@(x)isscalar(x) && x>=1);
addParameter(p,'FunctionTolerance',1e-14,@(x)isscalar(x) && x>0);
addParameter(p,'StepTolerance',1e-13,@(x)isscalar(x) && x>0);
addParameter(p,'OptimalityTolerance',1e-13,@(x)isscalar(x) && x>0);
addParameter(p,'FiniteDifferenceType','central',@(x)ischar(x) || isstring(x));
addParameter(p,'Display','off',@(x)ischar(x) || isstring(x));
parse(p,varargin{:});
opt = p.Results;

if ~isfield(cfg,'a') || ~isscalar(cfg.a) || cfg.a<=0
    error('cfg.a must be a positive scalar.');
end
if numel(cfg.widths)~=2 || numel(cfg.depths)~=2 || numel(cfg.gaps)~=1
    error('This refiner requires two widths, two depths, and one gap.');
end
if isempty(opt.Truncation)
    if ~isfield(cfg,'N') || ~isfield(cfg,'K')
        error('Supply Truncation or provide cfg.N and cfg.K.');
    end
    truncation = [cfg.N cfg.K];
else
    truncation = round(opt.Truncation);
end
if mod(truncation(1),2)~=1 || truncation(2)<1
    error('Truncation must contain an odd N and a positive K.');
end

lb = opt.LowerBounds(:).';
ub = opt.UpperBounds(:).';
x0 = x0(:).';
if any(lb>=ub) || any(x0<lb) || any(x0>ub)
    error('Bounds must be ordered and contain x0.');
end
if x0(1)<=0 || x0(1)>=.5
    error('Use 0<kappa<0.5 for the single n=-1 Rayleigh threshold.');
end
if exist('lsqnonlin','file')~=2
    error('lsqnonlin from Optimization Toolbox is required.');
end

% Fix the singular-vector phase with a stable groove-dominated anchor from
% the initial design.  This removes the arbitrary complex phase from the
% vector residual seen by the finite-difference Jacobian.
[~,initialState] = projected_residual(x0,[]);
grooveStart = initialState.n_kept_amplitudes+1;
[~,relativeAnchor] = max(abs(initialState.y(grooveStart:end)));
anchor = grooveStart+relativeAnchor-1;

historyX = zeros(0,6);
historyNorm = zeros(0,1);
lsqopt = optimoptions('lsqnonlin','Display',opt.Display, ...
    'MaxFunctionEvaluations',round(opt.MaxFunctionEvaluations), ...
    'MaxIterations',round(opt.MaxIterations), ...
    'FunctionTolerance',opt.FunctionTolerance, ...
    'StepTolerance',opt.StepTolerance, ...
    'OptimalityTolerance',opt.OptimalityTolerance, ...
    'FiniteDifferenceType',char(opt.FiniteDifferenceType));

[x,resnorm,residual,exitflag,output] = lsqnonlin(@objective,x0,lb,ub,lsqopt);
[compatibility,state] = projected_residual(x,anchor);
strict = ni2019_strict_rayleigh_operator(state.local_cfg,'TargetOrder',-1);

result = struct();
result.x = x;
result.kappa = x(1);
result.Omega = 1-x(1);
result.theta_deg = asind(x(1)/(1-x(1)));
result.depths_over_a = x(2:3);
result.widths_over_a = x(4:5);
result.gap_over_a = x(6);
result.fill_fraction = sum(x(4:6));
result.truncation = truncation;
result.compatibility_residual = compatibility;
result.compatibility_norm = norm(compatibility);
result.sigma_ratio = state.sigma_ratio;
result.sigma_min = state.sigma_min;
result.strict_operator = strict;
result.resnorm = resnorm;
result.residual = residual;
result.exitflag = exitflag;
result.output = output;
result.anchor_index_reduced = anchor;
result.history = struct('x',historyX,'compatibility_norm',historyNorm);
result.options = opt;

    function r = objective(x)
        [rComplex,~] = projected_residual(x,anchor);
        phaseFixed = [real(rComplex);imag(rComplex)];
        fillViolation = max(sum(x(4:6))-opt.FillMax,0);
        r = [phaseFixed;sqrt(opt.FillPenaltyWeight)*fillViolation];
        historyX(end+1,:) = x(:).';
        historyNorm(end+1,1) = norm(r);
    end

    function [r,state] = projected_residual(x,anchorIndex)
        local = cfg;
        Omega = 1-x(1);
        local.lambda = cfg.a/Omega;
        local.theta_i_deg = asind(x(1)/Omega);
        local.depths = x(2:3)*cfg.a;
        local.widths = x(4:5)*cfg.a;
        local.gaps = x(6)*cfg.a;
        local.N = truncation(1);
        local.K = truncation(2);
        local.solve_scattering = false;
        op = ni2019_full_eigen_operator(local);
        kyTolerance = 1e-8*max(abs(op.k0),1);
        finiteOpen = abs(imag(op.ky))<=kyTolerance & real(op.ky)>kyTolerance;
        target = op.orders==-1;
        remove = finiteOpen | target;
        keep = [find(~remove);(op.N+1:op.N+op.n_groove).'];
        Fred = op.F(:,keep);
        rowScale = max(vecnorm(Fred,2,2),sqrt(eps));
        Frow = Fred./rowScale;
        columnScale = max(vecnorm(Frow,2,1),sqrt(eps));
        Fscaled = Frow./columnScale;
        [~,S,V] = svd(Fscaled,'econ');
        singularValues = diag(S);
        y = V(:,end);
        if ~isempty(anchorIndex)
            if anchorIndex>numel(y) || abs(y(anchorIndex))<100*eps
                error('The fixed phase anchor became singular during refinement.');
            end
            y = y*exp(-1i*angle(y(anchorIndex)));
        end
        sigmaMax = max(singularValues(1),eps);
        r = (Fscaled*y)/sigmaMax;
        state = struct('local_cfg',local,'operator',op,'Fred',Fred, ...
            'Fscaled',Fscaled,'y',y,'singular_values',singularValues, ...
            'sigma_min',singularValues(end), ...
            'sigma_ratio',singularValues(end)/sigmaMax, ...
            'n_kept_amplitudes',nnz(~remove),'removed_orders',op.orders(remove));
    end
end
