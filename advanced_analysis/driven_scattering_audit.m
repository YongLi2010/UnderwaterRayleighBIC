%DRIVEN_SCATTERING_AUDIT Re-run the driven-scattering checks at fixed geometry.
%
% This audit is deliberately independent of the cached PRL-gallery values.
% Every numerical row is obtained from ni2019_modal_solver at the final
% manufacturable two-groove geometry.  The finite-truncation route and the
% separate off-Rayleigh point are kept as distinct calculations.
%
% Outputs (all under advanced_analysis/driven_scattering/):
%   driven_scattering_data.csv          route and off-Rayleigh audit rows
%   driven_scattering_audit.mat          machine-readable MATLAB data
%   fig_driven_scattering.pdf            four-panel vector figure
%   README.txt                           data/figure provenance
%
% Run from the repository root with:
%   matlab -batch "run('advanced_analysis/driven_scattering_audit.m')"

clear; close all; clc;

scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);
solverDir = fullfile(repoDir,'Ni2019_MATLAB');
outDir = fullfile(scriptDir,'driven_scattering');
if ~exist(outDir,'dir'), mkdir(outDir); end
addpath(solverDir);

% Final fixed geometry (SI units) from the independently audited design file.
designPath = fullfile(solverDir,'results','StrictRayleighBIC_200kHz_min1mm.mat');
design = load(designPath);
a = design.aPhysical;
c = design.cWater;
fRoute = 200216.874212;       % Hz; finite-truncation route operating point
thetaRoute = 6.0;             % degree; signed partner is -6 degree
fOff = 202430.0;              % Hz; separate ordinary off-Rayleigh point
thetaOff = 32.6328;           % degree
widths = design.widthsMm*1e-3;
depths = design.depthsMm*1e-3;
gaps = design.gapMm*1e-3;

routeNK = [81 11; 121 15; 201 25; 313 39; 401 50];
offNK = [81 11; 121 15; 161 20; 201 25; 313 39; 401 50];

% The target order is the -1 order for positive incidence and +1 for
% negative incidence.  The ordinary off-Rayleigh point is evaluated at +theta.
routePlus = repmat(empty_record(),size(routeNK,1),1);
routeMinus = repmat(empty_record(),size(routeNK,1),1);
offRows = repmat(empty_record(),size(offNK,1),1);

for j = 1:size(routeNK,1)
    N = routeNK(j,1); K = routeNK(j,2);
    routePlus(j) = driven_row('route_plus',N,K,thetaRoute,fRoute,-1, ...
        a,c,widths,depths,gaps);
    routeMinus(j) = driven_row('route_minus',N,K,-thetaRoute,fRoute,+1, ...
        a,c,widths,depths,gaps);
end

for j = 1:size(offNK,1)
    N = offNK(j,1); K = offNK(j,2);
    offRows(j) = driven_row('off_rayleigh',N,K,thetaOff,fOff,-1, ...
        a,c,widths,depths,gaps);
end

% The requested reciprocity comparison is made at (N,K)=(313,39), while
% (401,50) remains in the fixed-geometry sweep as a sensitivity check.
routeRef = find(routeNK(:,1)==313 & routeNK(:,2)==39,1);

% Explicit acceptance checks.  These should fail loudly if a future solver
% change makes energy closure or the reciprocal signed-order partner unreliable.
allRows = [routePlus; routeMinus; offRows];
maxEnergyError = max([allRows.energy_error]);
if maxEnergyError > 1e-8
    error('Driven-scattering audit failed: max energy error %.3g > 1e-8.',maxEnergyError);
end
recipTargetDiff = abs(routePlus(routeRef).eta_target-routeMinus(routeRef).eta_target);
recipZeroDiff = abs(routePlus(routeRef).eta_zero-routeMinus(routeRef).eta_zero);
% The reciprocal output directions have opposite signed angles, so compare
% theta_out(+theta) with -theta_out(-theta), not with the raw signed values.
recipAngleDiff = abs(routePlus(routeRef).theta_out+routeMinus(routeRef).theta_out);
if recipTargetDiff > 1e-8 || recipZeroDiff > 1e-8 || recipAngleDiff > 1e-8
    error(['Driven-scattering audit failed: reciprocal signed-order partner ', ...
        'differs beyond tolerance (eta %.3g, eta0 %.3g, angle %.3g deg).'], ...
        recipTargetDiff,recipZeroDiff,recipAngleDiff);
end

% Convert to a compact, solver-traceable table.  The signs and target orders
% remain explicit so that the signed-order reciprocity statement is auditable.
rows = [routePlus; routeMinus; offRows];
T = struct2table(rows);
writetable(T,fullfile(outDir,'driven_scattering_data.csv'));

% Save the complete records plus the exact design/provenance metadata.
audit = struct();
audit.design_path = designPath;
audit.a_m = a;
audit.c_m_per_s = c;
audit.widths_m = widths;
audit.depths_m = depths;
audit.gaps_m = gaps;
audit.route_frequency_Hz = fRoute;
audit.route_theta_deg = thetaRoute;
audit.off_frequency_Hz = fOff;
audit.off_theta_deg = thetaOff;
audit.route_truncations = routeNK;
audit.off_truncations = offNK;
audit.route_plus = routePlus;
audit.route_minus = routeMinus;
audit.off_rayleigh = offRows;
audit.max_energy_error = maxEnergyError;
audit.reciprocity_target_difference = recipTargetDiff;
audit.reciprocity_zero_difference = recipZeroDiff;
audit.reciprocity_angle_difference_deg = recipAngleDiff;
save(fullfile(outDir,'driven_scattering_audit.mat'),'-struct','audit','-v7');

% Four-panel PRL-style audit figure.  Panels answer complementary questions:
% (a) Is the practical route converged at fixed geometry? (No.)
% (b) Are the two signed incidences reciprocal? (Yes, within solver error.)
% (c) Does the separate off-Rayleigh point converge? (It is a distinct state.)
% (d) Are closure and conditioning reported alongside efficiencies? (Yes.)
make_figure(routePlus,routeMinus,offRows,routeNK,offNK,routeRef, ...
    recipTargetDiff,recipZeroDiff,recipAngleDiff,maxEnergyError, ...
    fullfile(outDir,'fig_driven_scattering.pdf'));

fid = fopen(fullfile(outDir,'README.txt'),'w');
if fid < 0, error('Could not create driven-scattering README.'); end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'Driven-scattering audit generated by advanced_analysis/driven_scattering_audit.m\n');
fprintf(fid,'Geometry: final fixed two-groove design loaded from StrictRayleighBIC_200kHz_min1mm.mat\n');
fprintf(fid,'Route: theta=+/-%.12g deg, f=%.12g Hz, fixed geometry, N/K sweep.\n',thetaRoute,fRoute);
fprintf(fid,'Off-Rayleigh: theta=+%.12g deg, f=%.12g Hz, fixed geometry, N/K sweep.\n',thetaOff,fOff);
fprintf(fid,'All efficiencies, closure, condition numbers, and angles are direct solver outputs.\n');
fprintf(fid,'The route is labeled truncation-sensitive; the off-Rayleigh point is not BIC-enhanced.\n');
fprintf(fid,'Reciprocity checks at N/K=%d/%d: target diff %.16g, zero-order diff %.16g, angle diff %.16g deg.\n', ...
    routeNK(routeRef,1),routeNK(routeRef,2),recipTargetDiff,recipZeroDiff,recipAngleDiff);
fprintf(fid,'Maximum energy closure error over all rows: %.16g.\n',maxEnergyError);

fprintf('Driven-scattering audit complete.\n');
fprintf('  route (+theta) eta_target at N/K=313/39: %.16g\n',routePlus(routeRef).eta_target);
fprintf('  route (-theta) eta_target at N/K=313/39: %.16g\n',routeMinus(routeRef).eta_target);
fprintf('  off-Rayleigh eta_target at N/K=401/50: %.16g\n',offRows(end).eta_target);
fprintf('  max energy error: %.3g; reciprocity target difference: %.3g\n', ...
    maxEnergyError,recipTargetDiff);

function row = empty_record()
row = struct('case_name','','N',0,'K',0,'theta_i_deg',0,'frequency_Hz',0, ...
    'lambda_m',0,'target_order',0,'eta_target',0,'eta_zero',0, ...
    'A_target_abs',0,'A_zero_abs',0,'theta_out',NaN,'energy_error',0, ...
    'condition_number',0,'total_efficiency',0,'target_is_propagating',false);
end

function row = driven_row(caseName,N,K,theta,freq,targetOrder,a,c,widths,depths,gaps)
cfg = struct();
cfg.lambda = c/freq;
cfg.a = a;
cfg.theta_i_deg = theta;
cfg.widths = widths;
cfg.depths = depths;
cfg.gaps = gaps;
cfg.N = N;
cfg.K = K;
cfg.solve_scattering = true;
R = ni2019_modal_solver(cfg);
it = find(R.orders==targetOrder,1);
i0 = find(R.orders==0,1);
if isempty(it) || isempty(i0)
    error('Expected target/zero order is absent at N/K=%d/%d.',N,K);
end
row = empty_record();
row.case_name = caseName;
row.N = N; row.K = K;
row.theta_i_deg = theta;
row.frequency_Hz = freq;
row.lambda_m = cfg.lambda;
row.target_order = targetOrder;
row.eta_target = R.eta(it);
row.eta_zero = R.eta(i0);
row.A_target_abs = abs(R.A(it));
row.A_zero_abs = abs(R.A(i0));
row.theta_out = R.theta_deg(it);
row.energy_error = R.energy_error;
row.condition_number = R.condition_number;
row.total_efficiency = R.total_efficiency;
row.target_is_propagating = R.is_propagating(it);
if ~row.target_is_propagating
    error('Target order is not propagating at %s N/K=%d/%d.',caseName,N,K);
end
end

function make_figure(routePlus,routeMinus,offRows,routeNK,offNK,routeRef, ...
    recipTargetDiff,recipZeroDiff,recipAngleDiff,maxEnergyError,pdfPath)
% MATLAB is used for both plotting and vector export in this solver-native
% audit so the source data and the rendered panels share one deterministic run.
set(groot,'defaultAxesFontName','Helvetica');
set(groot,'defaultTextFontName','Helvetica');
set(groot,'defaultAxesFontSize',8);
set(groot,'defaultTextFontSize',8);
set(groot,'defaultLegendFontName','Helvetica');
set(groot,'defaultLegendFontSize',7);
fig = figure('Color','w','Units','inches','Position',[1 1 6.95 5.20], ...
    'PaperUnits','inches','PaperPosition',[0 0 6.95 5.20]);
tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

blue = [0.08 0.28 0.57]; orange = [0.87 0.35 0.08];
teal = [0.05 0.52 0.49]; gray = [0.35 0.35 0.35];
routeLabel = arrayfun(@(r)sprintf('%d/%d',r.N,r.K),routePlus,'UniformOutput',false);
offLabel = arrayfun(@(r)sprintf('%d/%d',r.N,r.K),offRows,'UniformOutput',false);
xRoute = 1:numel(routePlus); xOff = 1:numel(offRows);

% (a) Fixed-geometry route convergence.
nexttile;
plot(xRoute,[routePlus.eta_target],'-o','Color',blue,'LineWidth',1.25, ...
    'MarkerSize',4.5,'MarkerFaceColor',blue); hold on;
plot(xRoute,[routePlus.eta_zero],'--s','Color',gray,'LineWidth',1.0, ...
    'MarkerSize',3.8,'MarkerFaceColor','w');
ylabel('\eta'); xlabel('Floquet/groove truncation (N/K)');
set(gca,'XTick',xRoute,'XTickLabel',routeLabel,'YLim',[0 1.05]);
legend({'target order','0 order'},'Location','best','Box','off');
title('(a) Fixed-geometry Rayleigh route','FontWeight','normal');
text(0.04,0.10,'not converged across N/K', 'Units','normalized', ...
    'Color',[0.65 0.12 0.08],'FontSize',7);
text(0.04,0.19,'f=200.216874 kHz; fixed geometry','Units','normalized', ...
    'Color',gray,'FontSize',7);
box off; grid on; set(gca,'GridAlpha',.14);

% (b) Signed-angle reciprocity at the final requested truncation.
nexttile;
vals = [routePlus(routeRef).eta_target,routeMinus(routeRef).eta_target; ...
    routePlus(routeRef).eta_zero,routeMinus(routeRef).eta_zero];
b = bar(vals,'grouped','BarWidth',.72); b(1).FaceColor=blue; b(2).FaceColor=orange;
b(1).EdgeColor='none'; b(2).EdgeColor='none';
set(gca,'XTick',1:2,'XTickLabel',{'target order','0 order'},'YLim',[0 .8]);
ylabel('\eta at N/K=313/39');
legend({'+6^\circ incidence','-6^\circ incidence'}, ...
    'Location','northoutside','Orientation','horizontal','Box','off');
title('(b) Reciprocal signed-order pair','FontWeight','normal'); box off; grid on; set(gca,'GridAlpha',.14);
text(.04,.08,sprintf('%s%s_{target}=%.1e',char(92),'eta',recipTargetDiff), ...
    'Units','normalized','Color',gray,'FontSize',7);

% (c) Separate off-Rayleigh convergence.
nexttile;
plot(xOff,[offRows.eta_target],'-o','Color',teal,'LineWidth',1.25, ...
    'MarkerSize',4.5,'MarkerFaceColor',teal); hold on;
plot(xOff,[offRows.eta_zero],'--s','Color',gray,'LineWidth',1.0, ...
    'MarkerSize',3.8,'MarkerFaceColor','w');
ylabel('\eta'); xlabel('Floquet/groove truncation (N/K)');
set(gca,'XTick',xOff,'XTickLabel',offLabel,'YLim',[0 1.05]);
legend({'n=-1 target','0 order'},'Location','best','Box','off');
title('(c) Separate off-Rayleigh point','FontWeight','normal');
text(0.04,0.10,'ordinary; not BIC-enhanced','Units','normalized', ...
    'Color',[0.05 0.35 0.32],'FontSize',7);
text(0.04,0.19,'f=202.430 kHz; fixed geometry','Units','normalized', ...
    'Color',gray,'FontSize',7);
box off; grid on; set(gca,'GridAlpha',.14);

% (d) Closure and conditioning reported on the same route sweep.
nexttile;
yyaxis left;
plot(xRoute,[routePlus.energy_error],'-o','Color',orange,'LineWidth',1.1, ...
    'MarkerSize',4,'MarkerFaceColor',orange); hold on;
set(gca,'YScale','log','YLim',[1e-16 1e-8]);
ylabel('|1-\Sigma\eta|');
yyaxis right;
plot(xRoute,[routePlus.condition_number],'-s','Color',blue,'LineWidth',1.1, ...
    'MarkerSize',4,'MarkerFaceColor',blue);
set(gca,'YScale','log'); ylabel('condition number');
xlabel('Floquet/groove truncation (N/K)');
set(gca,'XTick',xRoute,'XTickLabel',routeLabel);
title('(d) Numerical diagnostics','FontWeight','normal'); box off; grid on; set(gca,'GridAlpha',.14);
text(0.04,0.10,sprintf('max closure error %.1e',maxEnergyError), ...
    'Units','normalized','Color',gray,'FontSize',7);

% Keep panel labels compact and export editable vector text/paths.
axesHandles = findall(fig,'Type','axes');
for ax = axesHandles(:).'
    ax.LineWidth=0.65; ax.TickDir='out'; ax.FontName='Helvetica'; ax.FontSize=8;
    ax.XColor=[.15 .15 .15]; ax.YColor=[.15 .15 .15];
end
% Use an explicit physical paper size so the vector PDF has no letter-page
% whitespace around the compact 2x2 panel layout.
set(fig,'PaperUnits','inches','PaperSize',[6.95 5.20], ...
    'PaperPosition',[0 0 6.95 5.20],'PaperPositionMode','manual');
print(fig,pdfPath,'-dpdf','-painters');
close(fig);
end
