%% High-truncation convergence audit for the strict two-groove root
% This audit separates three questions:
%   1. Does a strict root exist after reoptimization at each truncation?
%   2. Do the optimized geometric parameters converge?
%   3. Does one fixed high-truncation geometry remain close to singular
%      when evaluated at neighboring truncations?
clear; close all; clc;

thisFile = mfilename('fullpath');
analysisDir = fileparts(thisFile);
repoDir = fileparts(analysisDir);
solverDir = fullfile(repoDir, 'Ni2019_MATLAB');
outDir = fullfile(analysisDir, 'extended_convergence');
if ~exist(outDir, 'dir'), mkdir(outDir); end
addpath(solverDir);

designFile = fullfile(solverDir, 'results', ...
    'StrictRayleighBIC_200kHz_min1mm.mat');
D = load(designFile);
x = D.xFinal(:).';

fTarget = D.fTarget;
cWater = D.cWater;
minimumWidth = 1e-3;
kappaRange = [.04, .13];
OmegaMinimum = 1 - kappaRange(2);
widthMinimumOverA = fTarget * minimumWidth / (cWater * OmegaMinimum);
lb = [kappaRange(1), .10, .10, .20, widthMinimumOverA, .08];
ub = [kappaRange(2), .80, .80, .68, .40, .28];
fillMax = .96;

% Continue beyond the published (313,39) root while preserving N=8K+1.
rootTruncations = [D.rootTruncations; 345,43; 377,47];
nOld = size(D.rootTruncations, 1);
nRoot = size(rootTruncations, 1);
xSequence = nan(nRoot, 6);
strictSigma = nan(nRoot, 1);
strictRaw = nan(nRoot, 1);
fullSigma = nan(nRoot, 1);
fullRaw = nan(nRoot, 1);
thresholdError = nan(nRoot, 1);
finiteOpenError = nan(nRoot, 1);

xSequence(1:nOld,:) = D.xSequence;
strictSigma(1:nOld) = D.rootSigma;
strictRaw(1:nOld) = D.rootRaw;

for j = 1:nRoot
    NK = rootTruncations(j,:);
    if j > nOld
        cfgSeed = make_cfg(x, NK(1), NK(2));
        refined = ni2019_refine_strict_rayleigh_bic(cfgSeed, x, ...
            'Truncation', NK, 'LowerBounds', lb, 'UpperBounds', ub, ...
            'FillMax', fillMax, 'MaxFunctionEvaluations', 1200, ...
            'MaxIterations', 160, 'Display', 'off');
        x = refined.x;
        xSequence(j,:) = x;
        strictSigma(j) = refined.sigma_ratio;
        strictRaw(j) = refined.strict_operator.strict_residual;
    else
        x = xSequence(j,:);
    end
    [fullSigma(j), fullRaw(j), thresholdError(j), finiteOpenError(j)] = ...
        full_diagnostic(xSequence(j,:), NK(1), NK(2));
    fprintf(['root N=%3d K=%2d: strict %.3e, full %.3e, ' ...
        '|A_RA|/||C|| %.3e, |A_open|/||C|| %.3e\n'], ...
        NK(1), NK(2), strictSigma(j), fullSigma(j), ...
        thresholdError(j), finiteOpenError(j));
end

xHighest = xSequence(end,:);
fixedTruncations = [281,35; 313,39; 345,43; 377,47; 409,51];
nFixed = size(fixedTruncations,1);
fixedStrictSigma = nan(nFixed,1);
fixedFullSigma = nan(nFixed,1);
fixedThresholdError = nan(nFixed,1);
fixedFiniteOpenError = nan(nFixed,1);
for j = 1:nFixed
    NK = fixedTruncations(j,:);
    R = ni2019_strict_rayleigh_operator(make_cfg(xHighest,NK(1),NK(2)), ...
        'TargetOrder', -1);
    fixedStrictSigma(j) = R.sigma_ratio;
    [fixedFullSigma(j),~,fixedThresholdError(j),fixedFiniteOpenError(j)] = ...
        full_diagnostic(xHighest,NK(1),NK(2));
    fprintf(['fixed-high root N=%3d K=%2d: strict %.3e, full %.3e, ' ...
        '|A_RA|/||C|| %.3e\n'], NK(1), NK(2), fixedStrictSigma(j), ...
        fixedFullSigma(j), fixedThresholdError(j));
end

parameterStep = [nan; max(abs(diff(xSequence,1,1)),[],2)];
rootTable = table(rootTruncations(:,1), rootTruncations(:,2), ...
    xSequence(:,1), xSequence(:,2), xSequence(:,3), xSequence(:,4), ...
    xSequence(:,5), xSequence(:,6), strictSigma, strictRaw, fullSigma, ...
    fullRaw, thresholdError, finiteOpenError, parameterStep, ...
    'VariableNames', {'N','K','kappa','d1_over_a','d2_over_a', ...
    'w1_over_a','w2_over_a','g_over_a','strict_sigma_ratio', ...
    'strict_raw_residual','full_sigma_ratio','full_raw_residual', ...
    'threshold_amplitude_over_Cnorm','finite_open_amplitude_over_Cnorm', ...
    'max_parameter_step'});
fixedTable = table(fixedTruncations(:,1), fixedTruncations(:,2), ...
    fixedStrictSigma, fixedFullSigma, fixedThresholdError, ...
    fixedFiniteOpenError, 'VariableNames', {'N','K', ...
    'strict_sigma_ratio','full_sigma_ratio', ...
    'threshold_amplitude_over_Cnorm','finite_open_amplitude_over_Cnorm'});
writetable(rootTable, fullfile(outDir,'extended_root_sequence.csv'));
writetable(fixedTable, fullfile(outDir,'fixed_high_root_crosscheck.csv'));
save(fullfile(outDir,'extended_convergence.mat'), 'D','lb','ub','fillMax', ...
    'rootTruncations','xSequence','strictSigma','strictRaw','fullSigma', ...
    'fullRaw','thresholdError','finiteOpenError','parameterStep', ...
    'fixedTruncations','fixedStrictSigma','fixedFullSigma', ...
    'fixedThresholdError','fixedFiniteOpenError','rootTable','fixedTable');

%% Restrained PRL-style audit figure
blue = [0.10,0.31,0.55]; orange = [0.86,0.38,0.08];
red = [0.76,0.15,0.18]; gray = [0.43,0.47,0.51];
fig = figure('Color','w','Units','inches','Position',[.5,.5,7.0,5.3]);
set(fig,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica', ...
    'DefaultAxesFontSize',8,'DefaultTextFontSize',8, ...
    'DefaultAxesLineWidth',.65,'DefaultLineLineWidth',1.15);
tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

nexttile;
plot(rootTruncations(:,2),xSequence(:,1),'-o','Color',blue); hold on;
plot(rootTruncations(:,2),xSequence(:,2),'-s','Color',orange);
plot(rootTruncations(:,2),xSequence(:,3),'-^','Color',red);
xlabel('groove truncation K'); ylabel('normalized parameter');
legend('\kappa','d_1/a','d_2/a','Location','best','Box','off');
title('(a) root-parameter sequence');

nexttile;
semilogy(rootTruncations(:,2),max(strictSigma,1e-18),'-o','Color',blue); hold on;
semilogy(rootTruncations(:,2),max(fullSigma,1e-18),'-s','Color',red);
yline(1e-8,'--','Color',gray);
xlabel('groove truncation K'); ylabel('relative singular value');
legend('radiation-reduced','full operator','10^{-8} guide', ...
    'Location','best','Box','off');
title('(b) reoptimized roots');

nexttile;
semilogy(fixedTruncations(:,2),max(fixedStrictSigma,1e-18),'-o','Color',blue); hold on;
semilogy(fixedTruncations(:,2),max(fixedFullSigma,1e-18),'-s','Color',red);
yline(1e-8,'--','Color',gray);
xlabel('groove truncation K'); ylabel('relative singular value');
legend('radiation-reduced','full operator','10^{-8} guide', ...
    'Location','best','Box','off');
title('(c) fixed highest-K geometry');

nexttile;
semilogy(rootTruncations(:,2),max(thresholdError,1e-18),'-o','Color',orange); hold on;
semilogy(rootTruncations(:,2),max(finiteOpenError,1e-18),'-s','Color',blue);
yline(1e-8,'--','Color',gray);
xlabel('groove truncation K'); ylabel('amplitude / ||C||_2');
legend('|A_{-1}|/||C||','|A_0|/||C||','10^{-8} guide', ...
    'Location','best','Box','off');
title('(d) unconstrained full-mode radiation');

axs = findall(fig,'Type','axes');
for j = 1:numel(axs)
    set(axs(j),'Box','on','Layer','top','TickDir','out', ...
        'XColor',[.12,.14,.18],'YColor',[.12,.14,.18], ...
        'GridColor',[.78,.80,.83],'GridAlpha',.34, ...
        'TitleFontWeight','normal');
    grid(axs(j),'on');
end
exportgraphics(fig,fullfile(outDir,'figS_extended_convergence.pdf'), ...
    'ContentType','vector','BackgroundColor','white');
close(fig);

function cfg = make_cfg(x,N,K)
Omega = 1-x(1);
cfg = struct('a',1,'lambda',1/Omega, ...
    'theta_i_deg',asind(x(1)/Omega),'depths',x(2:3), ...
    'widths',x(4:5),'gaps',x(6),'N',N,'K',K, ...
    'solve_scattering',false);
end

function [sigmaRatio,rawResidual,targetError,finiteError] = ...
    full_diagnostic(x,N,K)
cfg = make_cfg(x,N,K);
op = ni2019_full_eigen_operator(cfg);
[~,S,V] = svd(op.Fscaled,'econ');
sv = diag(S);
z = V(:,end)./transpose(op.column_scale);
z = z/max(norm(z),eps);
A = z(1:op.N);
C = z(op.N+1:end);
cn = max(norm(C),eps);
target = op.orders==-1;
tol = 1e-8*max(abs(op.k0),1);
finite = abs(imag(op.ky))<=tol & real(op.ky)>tol & ~target;
sigmaRatio = sv(end)/max(sv(1),eps);
rawResidual = norm(op.F*z)/max(norm(z),eps);
targetError = norm(A(target))/cn;
finiteError = norm(A(finite))/cn;
end
