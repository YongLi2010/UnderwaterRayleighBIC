function result=ni2019_optimize_single_rayleigh_bic(cfg,varargin)
%NI2019_OPTIMIZE_SINGLE_RAYLEIGH_BIC Optimize an off-Gamma threshold BIC.
%   x=[kappa,d1/a,d2/a,w1/a,w2/a,g/a], where kappa=kBloch*a/(2*pi).
%   The n=-1 Rayleigh condition is imposed exactly: Omega=k0*a/(2*pi)=1-kappa.
%   At 0<kappa<1/2, n=0 is the only finite-flux open diffraction order.

p=inputParser;
addParameter(p,'KappaRange',[.08 .48]);
addParameter(p,'DepthRange',[.01 .75]);
addParameter(p,'WidthRange',[.04 .65]);
addParameter(p,'GapRange',[.015 .35]);
addParameter(p,'FillMax',.94);
addParameter(p,'Starts',12);
addParameter(p,'Display','off');
addParameter(p,'RandomSeed',19);
parse(p,varargin{:});

if numel(cfg.widths)~=2 || numel(cfg.gaps)~=1
    error('The optimizer requires two grooves and one separating gap.');
end
if ~license('test','optimization_toolbox')
    error('Optimization Toolbox is required for the constrained multistart search.');
end

a=cfg.a;
if isfield(cfg,'lambda') && isfield(cfg,'theta_i_deg')
    kappaBase=(a/cfg.lambda)*sind(cfg.theta_i_deg);
else
    kappaBase=.25;
end
xBase=[kappaBase,cfg.depths(:).'/a,cfg.widths(:).'/a,cfg.gaps/a];
lb=[p.Results.KappaRange(1),repmat(p.Results.DepthRange(1),1,2), ...
    repmat(p.Results.WidthRange(1),1,2),p.Results.GapRange(1)];
ub=[p.Results.KappaRange(2),repmat(p.Results.DepthRange(2),1,2), ...
    repmat(p.Results.WidthRange(2),1,2),p.Results.GapRange(2)];
xBase=min(max(xBase,lb),ub);
A=zeros(1,6); A(4:6)=1; b=p.Results.FillMax;

rng(p.Results.RandomSeed,'twister');
starts=zeros(p.Results.Starts,6); starts(1,:)=xBase;
for n=2:p.Results.Starts
    starts(n,:)=lb+(ub-lb).*rand(1,6);
    if sum(starts(n,4:6))>p.Results.FillMax
        starts(n,4:6)=starts(n,4:6)*(.9*p.Results.FillMax/sum(starts(n,4:6)));
        starts(n,4:6)=max(starts(n,4:6),lb(4:6));
    end
end

opts=optimoptions('fmincon','Algorithm','sqp','Display',p.Results.Display, ...
    'MaxFunctionEvaluations',3500,'MaxIterations',450, ...
    'OptimalityTolerance',1e-12,'StepTolerance',1e-12, ...
    'ConstraintTolerance',1e-12);
solutions=nan(p.Results.Starts,6); scores=inf(p.Results.Starts,1);
diagnostics=cell(p.Results.Starts,1); flags=nan(p.Results.Starts,1);
for n=1:p.Results.Starts
    [solutions(n,:),scores(n),flags(n)]=fmincon(@objective,starts(n,:), ...
        A,b,[],[],lb,ub,[],opts);
    [~,diagnostics{n}]=objective(solutions(n,:));
end
[~,bestId]=min(scores); x=solutions(bestId,:);
[score,diagBest,opBest,modeBest,cfgBest]=objective(x);

result=struct('x',x,'score',score,'kappa',x(1),'Omega',1-x(1), ...
    'theta_deg',asind(x(1)/(1-x(1))),'depths_over_a',x(2:3), ...
    'widths_over_a',x(4:5),'gap_over_a',x(6), ...
    'diagnostics',diagBest,'operator',opBest,'mode',modeBest, ...
    'cfg',cfgBest,'solutions',solutions,'scores',scores, ...
    'exitflags',flags,'all_diagnostics',{diagnostics});

    function [f,dg,out,v,local]=objective(x)
        Omega=1-x(1);
        local=cfg;
        local.lambda=a/Omega;
        local.theta_i_deg=asind(x(1)/Omega);
        local.depths=x(2:3)*a;
        local.widths=x(4:5)*a;
        local.gaps=x(6)*a;
        local.solve_scattering=false;
        out=ni2019_full_eigen_operator(local);
        [~,S,V]=svd(out.Fscaled,'econ');
        sv=diag(S); srel=sv(end)/max(sv(1),eps);
        v=V(:,end)./transpose(out.column_scale); v=v/norm(v);
        Acoef=v(1:out.n_floquet);
        finiteOpen=abs(imag(out.ky))<1e-9*out.k0 & real(out.ky)>1e-7*out.k0;
        rad=sum(abs(Acoef(finiteOpen)).^2)/max(sum(abs(Acoef).^2),eps);
        target=out.orders==-1;
        grazing=abs(Acoef(target))^2/max(sum(abs(Acoef).^2),eps);
        grazingPenalty=max(.02-grazing,0)^2;
        f=sqrt(srel^2+rad)+1e-2*grazingPenalty;
        dg=struct('sigma_ratio',srel,'radiation_fraction',rad, ...
            'grazing_fraction',grazing,'finite_open_orders',out.orders(finiteOpen));
    end
end
