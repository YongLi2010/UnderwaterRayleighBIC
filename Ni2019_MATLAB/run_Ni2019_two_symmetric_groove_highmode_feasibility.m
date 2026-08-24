%% Two identical side grooves with transverse-mode completeness
% This is an ablation of the converged three-groove Gamma design.  The
% central wide groove is removed, while both remaining side grooves retain
% all transverse modes.  Odd and even mirror sectors are tested separately
% with A_-1=A_0=A_+1=0 imposed exactly.
clear; close all; clc;

outputDir=fullfile(fileparts(mfilename('fullpath')),'results');
if ~exist(outputDir,'dir'), mkdir(outputDir); end

% The first design literally retains the two outer grooves of the converged
% three-groove result and replaces the removed center by rigid material.
xThree=[0.192681441082727 0.676102882473460 0.192681441082525 ...
    0.0917192490036675 0.529780736125398 0.0917192490043602 ...
    0.0736588091900755 0.0736588091898476];
yOuterClone=[mean(xThree([1 3])),mean(xThree([4 6])), ...
    xThree(7)+xThree(5)+xThree(8)];

% Best representative non-vanishing-aperture points from independent broad
% searches.  They are comparison points, not claimed BIC roots.
yOddSearch=[1.001999 0.250000 0.2499999];
yEvenSearch=[1.649290228 0.100000 0.410701595];

Klist=[1 3 5 9 15 23 31].';
Nlist=[51 51 51 51 81 121 161].';
candidateNames={'outer-clone','outer-clone','odd-search','even-search'};
candidateY=[yOuterClone;yOuterClone;yOddSearch;yEvenSearch];
candidateParity={'odd','even','odd','even'};

records=cell(0,1);
cross=cell(size(candidateY,1),1);
for ic=1:size(candidateY,1)
    y=candidateY(ic,:);
    nK=numel(Klist);
    sigma=zeros(nK,1); raw=zeros(nK,1);
    pHigher=zeros(nK,1); vHigher=zeros(nK,1);
    modes=cell(nK,1);
    for ik=1:nK
        cfg=make_cfg(y,Nlist(ik),Klist(ik));
        R=evaluate_parity(cfg,candidateParity{ic});
        sigma(ik)=R.sigma_ratio;
        raw(ik)=R.strict_residual;
        pHigher(ik)=R.transverse.surface_pressure_higher_order_fraction;
        vHigher(ik)=R.transverse.surface_velocity_higher_order_fraction;
        modes{ik}=R;
        records{end+1,1}=table(string(candidateNames{ic}), ...
            string(candidateParity{ic}),Nlist(ik),Klist(ik), ...
            y(1),y(2),y(3),sigma(ik),raw(ik),pHigher(ik),vHigher(ik), ...
            'VariableNames',{'candidate','parity','N','K','d_over_a', ...
            'w_over_a','g_over_a','sigma_ratio','full_raw_residual', ...
            'pressure_q_ge_1_fraction','velocity_q_ge_1_fraction'}); %#ok<SAGROW>
    end
    cross{ic}=struct('name',candidateNames{ic},'parity',candidateParity{ic}, ...
        'y',y,'N',Nlist,'K',Klist,'sigma_ratio',sigma, ...
        'full_raw_residual',raw,'pressure_higher_fraction',pHigher, ...
        'velocity_higher_fraction',vHigher,'modes',{modes});
end
comparisonTable=vertcat(records{:});

% Evaluate both parity sectors for the literal outer-groove clone.
cloneOdd=cross{1};
cloneEven=cross{2};
cloneEvenSigma=cloneEven.sigma_ratio;
cloneEvenRaw=cloneEven.full_raw_residual;

% Modal composition at K=15.  Surface coefficients are the meaningful
% comparison for evanescent transverse modes; bottom coefficients alone can
% be exponentially small even when an aperture boundary layer is large.
modeId=find(Klist==15,1);
outerMode=cloneOdd.modes{modeId};
outerPressure=outerMode.transverse.surface_pressure_fraction_by_q;
outerVelocity=outerMode.transverse.surface_velocity_fraction_by_q;

threeCfg=struct('a',1,'lambda',1,'theta_i_deg',0, ...
    'depths',xThree(1:3),'widths',xThree(4:6),'gaps',xThree(7:8), ...
    'N',121,'K',15,'solve_scattering',false);
threeStrict=ni2019_strict_rayleigh_operator(threeCfg,'TargetOrder',-1, ...
    'EnforceOtherThreshold',true);
op3=threeStrict.full_operator;
p3=reshape(threeStrict.mode.C_scaled.*op3.cos_depth_normalized,op3.K,op3.L);
v3=reshape(threeStrict.mode.C_scaled.*op3.beta_sin_normalized,op3.K,op3.L);
p3=abs(p3).^2./max(sum(abs(p3).^2,1),eps);
v3=abs(v3).^2./max(sum(abs(v3).^2,1),eps);

% At Omega=1, beta_1^2 a^2=(2*pi)^2-(pi/(w/a))^2.  Two identical
% non-overlapping grooves require w/a<1/2, so their q=1 mode is always
% vertically evanescent.  The central groove of the three-groove root lies
% on the propagating side of this cutoff.
wCurve=linspace(.05,.60,500);
beta1Sq=(2*pi)^2-(pi./wCurve).^2;
outerGamma=sqrt((pi/yOuterClone(2))^2-(2*pi)^2);
outerDecay=exp(-outerGamma*yOuterClone(1));
centralBeta=sqrt((2*pi)^2-(pi/xThree(5))^2);

% ------------------------------- Figure -------------------------------
colors=ni2019_viridis(8);
navy=colors(2,:); teal=colors(4,:); green=colors(6,:);
yellow=colors(8,:); gray=[.52 .52 .52]; red=[.88 .13 .12];
fig=figure('Color','w','Position',[35 35 1500 850]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

nexttile;
hold on;
patch([0 1 1 0],[-.015 -.015 .015 .015],[.18 .20 .24], ...
    'EdgeColor','none');
s=.5*(1-2*yOuterClone(2)-yOuterClone(3));
for xl=[s,s+yOuterClone(2)+yOuterClone(3)]
    patch([xl xl+yOuterClone(2) xl+yOuterClone(2) xl], ...
        [0 0 -yOuterClone(1) -yOuterClone(1)],colors(5,:), ...
        'FaceAlpha',.75,'EdgeColor',navy,'LineWidth',1.2);
end
plot([0 1],[0 0],'-','Color',[.08 .08 .10],'LineWidth',1.4);
axis equal; xlim([0 1]); ylim([-max(.28,1.18*yOuterClone(1)) .08]);
xlabel('x/a'); ylabel('y/a'); title('(a) Two outer grooves only');
text(.5,-.55*yOuterClone(1),'rigid center (wide groove removed)', ...
    'HorizontalAlignment','center','Color',[.15 .15 .18]);
style_axes(gca);

nexttile;
plot(wCurve,beta1Sq,'Color',navy,'LineWidth',1.7); hold on;
yline(0,'--','Color',gray,'LineWidth',1.0);
xline(.5,':','Color',red,'LineWidth',1.2);
plot(yOuterClone(2),(2*pi)^2-(pi/yOuterClone(2))^2,'o', ...
    'Color',teal,'MarkerFaceColor',teal,'MarkerSize',6);
plot(xThree(5),(2*pi)^2-(pi/xThree(5))^2,'s', ...
    'Color',yellow,'MarkerFaceColor',yellow,'MarkerSize',6);
grid on; xlabel('groove width w/a'); ylabel('\beta_1^2 a^2');
title('(b) First transverse cutoff');
legend('q=1 dispersion','\beta_1^2=0','two-groove limit', ...
    'outer small groove','three-groove center','Location','southeast');
style_axes(gca);

nexttile;
semilogy(Klist,cloneOdd.sigma_ratio,'-o','Color',navy, ...
    'MarkerFaceColor',navy,'LineWidth',1.5); hold on;
semilogy(Klist,cloneEvenSigma,'-s','Color',red, ...
    'MarkerFaceColor',red,'LineWidth',1.4);
semilogy(Klist,cross{3}.sigma_ratio,'--^','Color',teal,'LineWidth',1.3);
semilogy(Klist,cross{4}.sigma_ratio,'--v','Color',green,'LineWidth',1.3);
yline(1e-8,':','Color',gray);
grid on; xlabel('transverse truncation K');
ylabel('\sigma_{min}/\sigma_{max}'); title('(c) Strict compatibility');
legend('outer clone: odd','outer clone: even','odd search point', ...
    'even search point','10^{-8} guide','Location','southwest');
style_axes(gca);

nexttile;
semilogy(Klist,cloneOdd.full_raw_residual,'-o','Color',navy, ...
    'MarkerFaceColor',navy,'LineWidth',1.5); hold on;
semilogy(Klist,cloneEvenRaw,'-s','Color',red, ...
    'MarkerFaceColor',red,'LineWidth',1.4);
semilogy(Klist,cross{3}.full_raw_residual,'--^','Color',teal,'LineWidth',1.3);
semilogy(Klist,cross{4}.full_raw_residual,'--v','Color',green,'LineWidth',1.3);
yline(1e-8,':','Color',gray);
grid on; xlabel('transverse truncation K');
ylabel('||Fz||_2/||z||_2'); title('(d) Full-equation residual');
legend('outer clone: odd','outer clone: even','odd search point', ...
    'even search point','10^{-8} guide','Location','best');
style_axes(gca);

nexttile;
qPlot=(0:min(5,numel(outerPressure)-1)).';
dataP=[outerPressure(qPlot+1),p3(qPlot+1,1),p3(qPlot+1,2)];
b=bar(qPlot,dataP,.86);
b(1).FaceColor=navy; b(2).FaceColor=teal; b(3).FaceColor=yellow;
grid on; xlabel('transverse order q'); ylabel('surface-pressure fraction');
title('(e) Aperture-pressure composition');
legend('two-groove approximate mode','three-groove outer', ...
    'three-groove center','Location','northoutside');
style_axes(gca);

nexttile;
dataV=[outerVelocity(qPlot+1),v3(qPlot+1,1),v3(qPlot+1,2)];
b=bar(qPlot,dataV,.86);
b(1).FaceColor=navy; b(2).FaceColor=teal; b(3).FaceColor=yellow;
grid on; xlabel('transverse order q'); ylabel('surface-velocity fraction');
title('(f) Aperture-velocity composition');
legend('two-groove approximate mode','three-groove outer', ...
    'three-groove center','Location','northoutside');
style_axes(gca);

allLegend=findobj(fig,'Type','Legend');
set(allLegend,'Color','w','TextColor',[.10 .10 .12], ...
    'EdgeColor',[.55 .55 .58]);

figureFile=fullfile(outputDir, ...
    'TwoSymmetricGrooveHighMode_feasibility.png');
exportgraphics(fig,figureFile,'Resolution',220);

csvFile=fullfile(outputDir, ...
    'TwoSymmetricGrooveHighMode_feasibility.csv');
writetable(comparisonTable,csvFile);
matFile=fullfile(outputDir, ...
    'TwoSymmetricGrooveHighMode_feasibility.mat');
save(matFile,'xThree','yOuterClone','yOddSearch','yEvenSearch', ...
    'Klist','Nlist','comparisonTable','cross','cloneEvenSigma', ...
    'cloneEvenRaw','outerPressure','outerVelocity','p3','v3', ...
    'outerGamma','outerDecay','centralBeta');

fprintf('\nTwo identical side-groove high-mode feasibility\n');
fprintf('  outer-clone y=[%.12g %.12g %.12g], fill=%.6f\n', ...
    yOuterClone,2*yOuterClone(2)+yOuterClone(3));
fprintf('  q1 outer: gamma*a=%.6f, exp(-gamma*d)=%.3e\n', ...
    outerGamma,outerDecay);
fprintf('  q1 three-groove center: beta*a=%.6f (propagating)\n',centralBeta);
fprintf('  outer clone K=31 odd:  sigma %.3e, raw %.3e\n', ...
    cloneOdd.sigma_ratio(end),cloneOdd.full_raw_residual(end));
fprintf('  outer clone K=31 even: sigma %.3e, raw %.3e\n', ...
    cloneEvenSigma(end),cloneEvenRaw(end));
fprintf('  strict removed orders: -1 0 1\n');
fprintf('  verdict: no converged two-identical-rectangular-groove root found.\n');
fprintf('  figure: %s\n  data:   %s\n  table:  %s\n', ...
    figureFile,matFile,csvFile);

assignin('base','TwoSymmetricGrooveHighModeFeasibility',struct( ...
    'table',comparisonTable,'cross',{cross},'figure',figureFile, ...
    'mat',matFile,'csv',csvFile));

function cfg=make_cfg(y,N,K)
cfg=struct('a',1,'lambda',1,'theta_i_deg',0, ...
    'depths',[y(1) y(1)],'widths',[y(2) y(2)],'gaps',y(3), ...
    'N',N,'K',K,'solve_scattering',false);
end

function R=evaluate_parity(cfg,parity)
if strcmp(parity,'odd')
    R=ni2019_two_groove_odd_strict_operator(cfg);
elseif strcmp(parity,'even')
    R=ni2019_two_groove_even_strict_operator(cfg);
else
    error('Unknown parity %s.',parity);
end
end

function style_axes(ax)
ax.FontName='Helvetica'; ax.FontSize=10; ax.LineWidth=.8;
ax.TickDir='out'; ax.Color='w'; ax.XColor=[.10 .10 .12];
ax.YColor=[.10 .10 .12]; ax.GridColor=[.68 .68 .70];
ax.MinorGridColor=[.82 .82 .84]; ax.Title.Color=[.10 .10 .12];
box(ax,'on');
end
