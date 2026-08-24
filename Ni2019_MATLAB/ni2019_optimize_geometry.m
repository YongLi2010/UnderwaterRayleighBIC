function best = ni2019_optimize_geometry(template,targetOrders,targetEta,varargin)
%NI2019_OPTIMIZE_GEOMETRY Re-design a groove grating for target efficiencies.
%   BEST = NI2019_OPTIMIZE_GEOMETRY(TEMPLATE,ORDERS,ETA) optimizes normalized
%   groove widths, depths, and gaps. TEMPLATE is a solver cfg structure and
%   defines a, lambda, incidence, number of grooves, N, and K.
%
%   This utility requires Optimization Toolbox. It uses bounded fmincon
%   searches and preserves the paper's fabrication bounds: width/a>=0.05,
%   gap/a>=0.05, and 0<=depth/lambda<=0.5.

p=inputParser;
addParameter(p,'NumStarts',8,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'Display','off',@(x)ischar(x)||isstring(x));
parse(p,varargin{:});
if ~license('test','optimization_toolbox')
    error('ni2019:OptimizationToolboxRequired', ...
        'Optimization Toolbox is required for geometry re-optimization.');
end

L=numel(template.widths);
xPrinted=[template.widths/template.a, template.depths/template.lambda, ...
          template.gaps/template.a];
lb=[.05*ones(1,L), zeros(1,L), .05*ones(1,L-1)];
ub=[.70*ones(1,L), .50*ones(1,L), .70*ones(1,L-1)];
A=[ones(1,L),zeros(1,L),ones(1,L-1)]; b=.95;
opts=optimoptions('fmincon','Algorithm','sqp','Display',p.Results.Display, ...
    'MaxFunctionEvaluations',5000,'MaxIterations',500, ...
    'OptimalityTolerance',1e-10,'StepTolerance',1e-10);

rng(2019,'twister');
starts=zeros(p.Results.NumStarts,numel(xPrinted));
starts(1,:)=min(max(xPrinted,lb),ub);
for s=2:p.Results.NumStarts
    valid=false;
    while ~valid
        z=lb+(ub-lb).*rand(size(lb));
        valid=A*z.'<=b;
    end
    starts(s,:)=z;
end

bestF=inf; bestX=[]; bestExit=[];
for s=1:size(starts,1)
    [x,fval,exitflag]=fmincon(@objective,starts(s,:),A,b,[],[],lb,ub,[],opts);
    if fval<bestF
        bestF=fval; bestX=x; bestExit=exitflag;
    end
end
bestCfg=unpack(bestX);
bestResult=ni2019_modal_solver(bestCfg);
best=struct('cfg',bestCfg,'result',bestResult,'objective',bestF, ...
    'exitflag',bestExit,'normalized_variables',bestX, ...
    'target_orders',targetOrders,'target_eta',targetEta);

    function f=objective(x)
        try
            r=ni2019_modal_solver(unpack(x));
            calc=zeros(size(targetEta));
            for ii=1:numel(targetOrders)
                calc(ii)=r.eta(r.orders==targetOrders(ii));
            end
            f=sum((calc-targetEta).^2) + 1e-4*r.energy_error^2;
            if ~isfinite(f), f=1e6; end
        catch
            f=1e6;
        end
    end
    function cfg=unpack(x)
        x=x(:).';
        cfg=template;
        cfg.widths=x(1:L)*cfg.a;
        cfg.depths=x(L+1:2*L)*cfg.lambda;
        cfg.gaps=x(2*L+1:end)*cfg.a;
    end
end
