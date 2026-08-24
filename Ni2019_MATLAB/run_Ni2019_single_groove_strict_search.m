%% Global search for a strict single-groove, single-Rayleigh BIC
% The single rectangular groove retains K transverse cosine modes.  The
% search varies x=[kappa,d/a,w/a] at the exact n=-1 Rayleigh threshold
% Omega=1-kappa and minimizes the strict homogeneous compatibility residual
% after A_-1 and the finite-flux A_0 channel are removed.
clear; close all; clc;

rootDir=fileparts(mfilename('fullpath'));
outputDir=fullfile(rootDir,'results');
if ~exist(outputDir,'dir'), mkdir(outputDir); end
rng(260823,'twister');

% A genuinely large aperture is enforced so the q=1 transverse component
% can approach or cross its groove cutoff.  Kappa is kept away from 0.5,
% where the nominal n=0 radiation channel would itself become grazing and
% turn the optimization into a different double-threshold problem.
lb=[.04 .04 .50];
ub=[.25 2.00 .95];
globalTruncations=[41 8;61 12];

% Deterministic broad random survey supplies diverse basin seeds.  Logarithmic
% depth sampling gives comparable coverage of shallow and multi-quarter-wave
% cavities.
nSurvey=1800;
survey=zeros(nSurvey,3);
survey(:,1)=lb(1)+(ub(1)-lb(1))*rand(nSurvey,1);
survey(:,2)=exp(log(lb(2))+(log(ub(2))-log(lb(2)))*rand(nSurvey,1));
survey(:,3)=lb(3)+(ub(3)-lb(3))*rand(nSurvey,1);
surveyScore=nan(nSurvey,1);
fprintf('Single-groove survey: %d samples, K up to %d transverse modes\n', ...
    nSurvey,max(globalTruncations(:,2)));
for j=1:nSurvey
    surveyScore(j)=multi_score(survey(j,:),globalTruncations);
end
[surveyScore,order]=sort(surveyScore);
survey=survey(order,:);

% Local bounded refinement from separated low-residual seeds.  The objective
% is the worst normalized strict singular value across two discretizations.
nStarts=14;
starts=select_separated_starts(survey,surveyScore,nStarts,lb,ub);
opts=optimoptions('fmincon','Algorithm','sqp','Display','off', ...
    'MaxIterations',100,'MaxFunctionEvaluations',900, ...
    'StepTolerance',1e-11,'OptimalityTolerance',1e-9);
localX=zeros(nStarts,3); localScore=nan(nStarts,1); localExit=zeros(nStarts,1);
for j=1:nStarts
    [localX(j,:),localScore(j),localExit(j)]=fmincon( ...
        @(z)multi_score(z,globalTruncations),starts(j,:),[],[],[],[],lb,ub,[],opts);
    fprintf('  start %2d: score=%+.6f, x=[%.9f %.9f %.9f]\n', ...
        j,localScore(j),localX(j,:));
end
[localScore,order]=sort(localScore); localX=localX(order,:); localExit=localExit(order);
x=localX(1,:);

% At each increasing truncation, minimize the single finite-discretization
% compatibility residual.  Unlike the two-groove case, this is not assumed
% to possess an exact root; the residual floor is part of the result.
rootTruncations=[61 12;101 20;141 28;181 36;221 44];
nRoot=size(rootTruncations,1);
xSequence=zeros(nRoot,3); rootSigma=nan(nRoot,1); rootRaw=nan(nRoot,1);
transverseFraction=nan(nRoot,1); rootResult=cell(nRoot,1);
for j=1:nRoot
    NK=rootTruncations(j,:);
    obj=@(z)single_score(z,NK(1),NK(2));
    [x,~,~]=fmincon(obj,x,[],[],[],[],lb,ub,[],opts);
    R=strict_result(x,NK(1),NK(2));
    xSequence(j,:)=x; rootSigma(j)=R.sigma_ratio; rootRaw(j)=R.strict_residual;
    C=R.groove.C_by_groove(:,1);
    transverseFraction(j)=sum(abs(C(2:end)).^2)/max(sum(abs(C).^2),eps);
    rootResult{j}=R;
    fprintf('polish N=%3d K=%2d: sigma=%.6e raw=%.6e, transverse fraction=%.6f\n', ...
        NK(1),NK(2),rootSigma(j),rootRaw(j),transverseFraction(j));
end
xFinal=xSequence(end,:);
OmegaFinal=1-xFinal(1); thetaFinal=asind(xFinal(1)/OmegaFinal);

% Fixed-geometry cross-truncation check distinguishes a true convergent zero
% from a minimum that moves with K.
crossTruncations=[61 12;81 16;101 20;121 24;141 28;161 32;181 36;201 40;221 44;241 48];
nCross=size(crossTruncations,1); crossSigma=nan(nCross,1); crossRaw=nan(nCross,1);
crossTransverse=nan(nCross,1);
for j=1:nCross
    R=strict_result(xFinal,crossTruncations(j,1),crossTruncations(j,2));
    crossSigma(j)=R.sigma_ratio; crossRaw(j)=R.strict_residual;
    C=R.groove.C_by_groove(:,1);
    crossTransverse(j)=sum(abs(C(2:end)).^2)/max(sum(abs(C).^2),eps);
end

% Compare directly against the documented two-groove finite root.
two=load(fullfile(outputDir,'StrictRayleighBIC_optimized_practical.mat'));
twoSigma=two.rootSigma(1); twoRaw=two.rootRaw(1);

% Report the transverse modal spectrum of the best high-order single groove.
best=rootResult{end};
C=best.groove.C_by_groove(:,1);
P=reshape(best.groove.surface_coefficients,xFinal(1)*0+rootTruncations(end,2),1);
q=(0:numel(C)-1).';
modalTable=table(q,abs(C),abs(P),'VariableNames', ...
    {'q','C_physical_abs','surface_pressure_abs'});

% Conservative numerical classification.  Passing requires a small residual,
% stability under the last refinement, and a small adjacent cross value.
parameterStep=max(abs(xSequence(end,:)-xSequence(end-1,:)));
foundStrict=rootSigma(end)<1e-8 && crossSigma(end)<1e-7 && parameterStep<1e-4;
if foundStrict
    verdict='strict single-groove candidate found; independent pole and energy checks required';
else
    verdict=['no strict single-groove BIC found in the searched domain; ' ...
        'the optimized residual remains finite or nonconvergent'];
end

%% Summary figure
navy=[.035 .10 .30]; orange=[1 .36 .08]; red=[.90 .05 .06]; gray=[.5 .5 .5];
fig=figure('Color','w','Position',[70 70 1350 790]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
scatter(survey(:,1),survey(:,2),13,surveyScore,'filled'); hold on;
plot(xSequence(:,1),xSequence(:,2),'w-o','LineWidth',1.5,'MarkerFaceColor',red);
xlabel('\kappa'); ylabel('d/a'); title('(a) Global single-groove search');
cb=colorbar; cb.Label.String='log_{10} worst \sigma ratio'; grid on; box on;
colormap(gca,ni2019_viridis(256));

nexttile;
semilogy(rootTruncations(:,2),rootSigma,'o-','Color',red,'LineWidth',1.6, ...
    'MarkerFaceColor',red); hold on;
semilogy(crossTruncations(:,2),crossSigma,'s-','Color',navy,'LineWidth',1.5, ...
    'MarkerFaceColor',navy);
yline(twoSigma,'--','Color',orange,'LineWidth',1.4);
yline(1e-8,':','Color',gray,'LineWidth',1.2);
xlabel('transverse truncation K'); ylabel('\sigma_{min}/\sigma_{max}');
title('(b) Residual convergence'); grid on; box on;
legend('reoptimized single groove','fixed final single geometry', ...
    'two-groove finite root','10^{-8} criterion','Location','best');

nexttile;
plot(rootTruncations(:,2),xSequence(:,1),'-o','Color',navy,'LineWidth',1.5); hold on;
plot(rootTruncations(:,2),xSequence(:,2),'-s','Color',orange,'LineWidth',1.5);
plot(rootTruncations(:,2),xSequence(:,3),'-^','Color',red,'LineWidth',1.5);
xlabel('transverse truncation K'); ylabel('optimized parameter');
title('(c) Minimum location versus truncation'); grid on; box on;
legend('\kappa','d/a','w/a','Location','best');

nexttile;
stem(q,abs(C)/max(abs(C)),'.','Color',navy,'LineWidth',1.3); hold on;
stem(q,abs(P)/max(abs(P)),'.','Color',orange,'LineWidth',1.2);
set(gca,'YScale','log'); ylim([1e-8 1.2]);
xlabel('single-groove transverse order q'); ylabel('normalized coefficient');
title(sprintf('(d) Best mode: transverse fraction %.3f',transverseFraction(end)));
grid on; box on; legend('|C_q|','surface |P_q|','Location','best');

sgtitle(sprintf('Single-groove strict search: %s',verdict), ...
    'Interpreter','none','Color','k');
style_figure(fig);
figureFile=fullfile(outputDir,'SingleGroove_strict_Rayleigh_search.png');
exportgraphics(fig,figureFile,'Resolution',220);

summaryTable=table(rootTruncations(:,1),rootTruncations(:,2),xSequence(:,1), ...
    xSequence(:,2),xSequence(:,3),rootSigma,rootRaw,transverseFraction, ...
    'VariableNames',{'N','K','kappa','d_over_a','w_over_a','sigma_ratio', ...
    'raw_residual','transverse_fraction'});
crossTable=table(crossTruncations(:,1),crossTruncations(:,2),crossSigma,crossRaw, ...
    crossTransverse,'VariableNames',{'N','K','sigma_ratio','raw_residual', ...
    'transverse_fraction'});
csvFile=fullfile(outputDir,'SingleGroove_strict_Rayleigh_search.csv');
writetable(summaryTable,csvFile);
matFile=fullfile(outputDir,'SingleGroove_strict_Rayleigh_search.mat');
save(matFile,'lb','ub','globalTruncations','survey','surveyScore','starts', ...
    'localX','localScore','localExit','rootTruncations','xSequence','rootSigma', ...
    'rootRaw','transverseFraction','crossTruncations','crossSigma','crossRaw', ...
    'crossTransverse','xFinal','OmegaFinal','thetaFinal','modalTable', ...
    'summaryTable','crossTable','twoSigma','twoRaw','parameterStep','foundStrict','verdict');

fprintf('\nSINGLE-GROOVE STRICT SEARCH VERDICT\n');
fprintf('  %s\n',verdict);
fprintf('  searched: kappa=[%.3f,%.3f], d/a=[%.3f,%.3f], w/a=[%.3f,%.3f]\n', ...
    lb(1),ub(1),lb(2),ub(2),lb(3),ub(3));
fprintf('  best x=[%.12g %.12g %.12g], Omega=%.12g, theta=%.9f deg\n', ...
    xFinal,OmegaFinal,thetaFinal);
fprintf('  final optimized sigma=%.6e, raw=%.6e\n',rootSigma(end),rootRaw(end));
fprintf('  adjacent fixed-geometry sigma=%.6e\n',crossSigma(end));
fprintf('  parameter step=%.6e, transverse fraction=%.6f\n', ...
    parameterStep,transverseFraction(end));
fprintf('  two-groove reference sigma=%.6e, raw=%.6e\n',twoSigma,twoRaw);
fprintf('  %s\n  %s\n  %s\n',figureFile,csvFile,matFile);

function value=multi_score(x,truncations)
values=zeros(size(truncations,1),1);
for m=1:size(truncations,1)
    R=strict_result(x,truncations(m,1),truncations(m,2));
    values(m)=R.sigma_ratio;
end
value=log10(max(max(values),1e-16));
if ~isfinite(value), value=4; end
end

function value=single_score(x,N,K)
R=strict_result(x,N,K);
value=log10(max(R.sigma_ratio,1e-16));
if ~isfinite(value), value=4; end
end

function R=strict_result(x,N,K)
Omega=1-x(1);
cfg=struct('a',1,'lambda',1/Omega,'theta_i_deg',asind(x(1)/Omega), ...
    'depths',x(2),'widths',x(3),'gaps',[],'N',N,'K',K, ...
    'solve_scattering',false);
R=ni2019_strict_rayleigh_operator(cfg,'TargetOrder',-1, ...
    'EnforceAllFiniteOpen',true,'Verbose',false);
end

function starts=select_separated_starts(survey,scores,nStarts,lb,ub)
starts=zeros(nStarts,3); n=0;
scale=ub-lb;
for j=1:size(survey,1)
    candidate=survey(j,:);
    if n==0 || all(vecnorm((starts(1:n,:)-candidate)./scale,2,2)>.08)
        n=n+1; starts(n,:)=candidate;
        if n==nStarts, break; end
    end
end
if n<nStarts
    starts(n+1:end,:)=survey(1:nStarts-n,:);
end
end

function style_figure(fig)
axesList=findall(fig,'Type','axes');
for j=1:numel(axesList)
    ax=axesList(j); ax.Color='w'; ax.XColor='k'; ax.YColor='k';
    ax.GridColor=[.72 .72 .72]; ax.MinorGridColor=[.84 .84 .84];
    ax.Title.Color='k'; ax.XLabel.Color='k'; ax.YLabel.Color='k';
end
legends=findall(fig,'Type','legend');
for j=1:numel(legends)
    legends(j).Color='w'; legends(j).TextColor='k'; legends(j).EdgeColor=[.3 .3 .3];
end
end
