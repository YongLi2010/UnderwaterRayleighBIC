%% Exact-physical-feature search: two-groove strict Rayleigh BIC at 200 kHz
% z = [kappa,d1/a,d2/a,w1/lambda,w2/lambda,g/lambda].
% At 200 kHz in c=1500 m/s water, 1.5 mm/lambda = 0.2 exactly.
function run_Ni2019_strict_rayleigh_bic_200k_1p5mm_search(branch,minimumMm)
close all; clc;
if nargin<1, branch='general'; end
if nargin<2, minimumMm=1.5; end
rootDir=fileparts(mfilename('fullpath')); outputDir=fullfile(rootDir,'results');
addpath(rootDir);
fTarget=200e3; cWater=1500; lambda=cWater/fTarget;
minimumFeature=minimumMm*1e-3; h=minimumFeature/lambda;
truncations=[31 5;41 7;61 9]; fillMax=.999;
if strcmpi(branch,'practical08')
    w2Minimum=0.8e-3/lambda; gapMinimum=1.0e-3/lambda;
    lb=[.025 .02 .02 .42 w2Minimum gapMinimum];
    ub=[.18 1.45 1.45 .78 .45 .45];
    fileTag='w2_0p8_gap1'; trialSeeds=[1801 1802 1803];
elseif strcmpi(branch,'q1')
    % Preserve the demonstrated mechanism: the wide groove must be at or
    % above the q=1 cutoff w1/lambda=0.5.  The three minimum lateral
    % fractions already occupy 0.9 lambda, so kappa must remain below 0.1.
    kappaMax=max(.01,min(.22,(fillMax-(.5+2*h))/fillMax));
    lb=[.002 .02 .02 .5 h h]; ub=[kappaMax 1.45 1.45 .78 .40 .40];
    fileTag='q1'; trialSeeds=[1561 1562 1563]+round(100*minimumMm);
else
    lb=[.002 .02 .02 h h h]; ub=[.22 1.45 1.45 .78 .55 .55];
    fileTag='general'; trialSeeds=[1531 1532 1533]+round(100*minimumMm);
end
Aineq=[fillMax 0 0 1 1 1]; bineq=fillMax;
D=load(fullfile(outputDir,'StrictRayleighBIC_200kHz_min1mm.mat'));
x0=D.xFinal(:).'; Omega0=1-x0(1);
zBase=[x0(1),x0(2:3),max(Omega0*x0(4),h),h,h];
if strcmpi(branch,'practical08')
    zBase=[x0(1),x0(2:3),Omega0*x0(4),max(Omega0*x0(5),w2Minimum), ...
        max(Omega0*x0(6),gapMinimum)];
end
zBase=make_feasible(zBase,lb,ub,Aineq,bineq);

trials=cell(numel(trialSeeds),1);
bestScore=inf; bestZ=[]; bestTrial=0;
for trial=1:numel(trialSeeds)
    rng(trialSeeds(trial),'twister');
    initial=zeros(12,6); initial(1,:)=zBase;
    for j=2:size(initial,1)
        initial(j,:)=random_feasible(lb,ub,Aineq,bineq);
    end
    psopt=optimoptions('particleswarm','Display','off','SwarmSize',96, ...
        'MaxIterations',110,'MaxStallIterations',28, ...
        'FunctionTolerance',1e-10,'InitialPoints',initial);
    [zg,fg,efg,og]=particleswarm(@objective,6,lb,ub,psopt);
    fopt=optimoptions('fmincon','Algorithm','sqp','Display','off', ...
        'MaxFunctionEvaluations',3200,'MaxIterations',320, ...
        'OptimalityTolerance',1e-12,'StepTolerance',1e-12, ...
        'ConstraintTolerance',1e-12);
    starts=[zg;initial(1:7,:)]; localZ=zeros(size(starts));
    localScore=inf(size(starts,1),1); localExit=zeros(size(starts,1),1);
    for j=1:size(starts,1)
        [localZ(j,:),localScore(j),localExit(j)]=fmincon(@objective, ...
            starts(j,:),Aineq,bineq,[],[],lb,ub,[],fopt);
    end
    [score,id]=min(localScore); z=localZ(id,:);
    trials{trial}=struct('seed',trialSeeds(trial),'global_z',zg, ...
        'global_score',fg,'global_exitflag',efg,'global_output',og, ...
        'local_z',localZ,'local_score',localScore, ...
        'local_exitflag',localExit,'best_z',z,'best_score',score);
    [x,periodMm,featuresMm]=convert(z);
    fprintf('trial %d seed %d: score %.6e, kappa %.9f, ', ...
        trial,trialSeeds(trial),score,z(1));
    fprintf('features [%.6f %.6f %.6f] mm, x=[',featuresMm);
    fprintf(' %.10g',x); fprintf(' ], a %.6f mm\n',periodMm);
    if score<bestScore, bestScore=score; bestZ=z; bestTrial=trial; end
end

[bestX,periodMm,featuresMm]=convert(bestZ);
diagnostic=cell(size(truncations,1),1); sigma=zeros(size(truncations,1),1);
raw=sigma; groove=sigma;
for j=1:size(truncations,1)
    diagnostic{j}=strict_result(bestX,truncations(j,1),truncations(j,2));
    sigma(j)=diagnostic{j}.sigma_ratio;
    raw(j)=diagnostic{j}.strict_residual;
    groove(j)=diagnostic{j}.groove.pressure_proxy_fraction;
end
summary=table(bestZ(1),1-bestZ(1),asind(bestZ(1)/(1-bestZ(1))), ...
    periodMm,bestX(2),bestX(3),featuresMm(1),featuresMm(2), ...
    featuresMm(3),bestScore,max(sigma),max(raw),bestTrial, ...
    'VariableNames',{'kappa','Omega','theta_deg','period_mm','d1_over_a', ...
    'd2_over_a','w1_mm','w2_mm','gap_mm','objective','max_sigma', ...
    'max_raw','best_trial'});
minimumTag=strrep(sprintf('%.2g',minimumMm),'.','p');
writetable(summary,fullfile(outputDir, ...
    ['StrictRayleighBIC_200kHz_min' minimumTag 'mm_search_' fileTag '.csv']));
save(fullfile(outputDir, ...
    ['StrictRayleighBIC_200kHz_min' minimumTag 'mm_search_' fileTag '.mat']), ...
    'fTarget','cWater','minimumFeature','h','truncations','fillMax', ...
    'lb','ub','trials','bestScore','bestZ','bestX','periodMm', ...
    'featuresMm','diagnostic','sigma','raw','groove','summary');
disp(summary);

    function score=objective(z)
        if any(~isfinite(z)) || any(z<lb) || any(z>ub) || Aineq*z(:)>bineq
            violation=max(Aineq*z(:)-bineq,0);
            score=10+1e3*violation^2; return;
        end
        x=convert(z); ratios=zeros(size(truncations,1),1);
        fractions=ratios;
        for it=1:size(truncations,1)
            try
                R=strict_result(x,truncations(it,1),truncations(it,2));
                ratios(it)=R.sigma_ratio;
                fractions(it)=R.groove.pressure_proxy_fraction;
            catch
                score=100; return;
            end
        end
        penalty=(max(.04-min(fractions),0)/.04)^2;
        score=max(ratios)+.15*mean(ratios)+penalty;
    end

    function R=strict_result(x,N,K)
        Omega=1-x(1);
        cfg=struct('a',1,'lambda',1/Omega, ...
            'theta_i_deg',asind(x(1)/Omega),'depths',x(2:3), ...
            'widths',x(4:5),'gaps',x(6),'N',N,'K',K, ...
            'solve_scattering',false);
        R=ni2019_strict_rayleigh_operator(cfg,'TargetOrder',-1);
    end

    function [x,periodMm,featuresMm]=convert(z)
        Omega=1-z(1); x=[z(1),z(2:3),z(4:6)/Omega];
        periodMm=1e3*lambda*Omega; featuresMm=1e3*lambda*z(4:6);
    end
end

function z=make_feasible(z,lb,ub,A,b)
z=min(max(z,lb),ub);
if A*z(:)>b
    available=b-A(1)*z(1); lateral=z(4:6);
    lateral=lb(4:6)+(lateral-lb(4:6))* ...
        max((available-sum(lb(4:6)))/max(sum(lateral-lb(4:6)),eps),0);
    z(4:6)=lateral;
end
end

function z=random_feasible(lb,ub,A,b)
for attempt=1:1000
    z=lb+(ub-lb).*rand(size(lb));
    if A*z(:)<=b, return; end
end
z=make_feasible(.5*(lb+ub),lb,ub,A,b);
end
