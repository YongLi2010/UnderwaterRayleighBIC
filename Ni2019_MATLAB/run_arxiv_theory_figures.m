function run_arxiv_theory_figures()
%RUN_ARXIV_THEORY_FIGURES  Generate five theory figures for the arXiv draft.
%
% All numerical panels are read from the verified v5 cache.  The script does
% not call the solver, so the geometry and truncation labels cannot drift from
% the cached calculations.  exportgraphics(...,'ContentType','vector') keeps
% axes, text, lines, and markers as vector objects; heat-map samples remain
% embedded image data where MATLAB requires them.

rootDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(rootDir);
resultDir = fullfile(rootDir, 'results');
outDir = fullfile(projectDir, 'arxiv_theory_paper', 'figures');
if ~exist(outDir, 'dir'), mkdir(outDir); end

cacheFile = fullfile(resultDir, 'PRL_article_figure_data_v5.mat');
designFile = fullfile(resultDir, 'StrictRayleighBIC_200kHz_min1mm.mat');
poleFile = fullfile(resultDir, 'StrictRayleighBIC_200kHz_min1mm_poles.mat');
singleFile = fullfile(resultDir, 'SingleGroove_strict_Rayleigh_search.mat');
assert(exist(cacheFile, 'file') == 2, 'Required v5 cache is missing: %s', cacheFile);
assert(exist(designFile, 'file') == 2, 'Required final geometry cache is missing: %s', designFile);
assert(exist(poleFile, 'file') == 2, 'Required pole cache is missing: %s', poleFile);
assert(exist(singleFile, 'file') == 2, 'Required ablation cache is missing: %s', singleFile);

C = load(cacheFile);
D = load(designFile);
P0 = load(poleFile);
single = load(singleFile);
assert_cache_compatibility(C, D, P0, single);

p = figure_palette();
kappa0 = D.xFinal(1);
Omega0 = D.OmegaFinal;
theta0 = D.thetaFinal;
f0 = D.fTarget;
aPhysical = D.aPhysical;
c0 = D.cWater;
P = P0.P;

% No supertitles: the manuscript captions carry interpretation, while each
% panel carries only its local label and variables.
fig = make_figure([7.0, 5.45]);
t = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(t); draw_unit_cell(D, p); title('(a) unit cell and final geometry');
nexttile(t); draw_topology(kappa0, Omega0, p); title('(b) Rayleigh-channel topology');
nexttile(t); draw_channel_count(kappa0, Omega0, p); title('(c) channels at the strict point');
nexttile(t); draw_bic_location(theta0, f0, aPhysical, c0, p); title('(d) strict BIC location');
format_figure(fig); export_vector(fig, fullfile(outDir, 'fig1_structure_topology.pdf')); close(fig);

fig = make_figure([7.0, 5.45]);
t = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(t); draw_pole_branch(P, kappa0, Omega0, p); title('(a) pole meets the Rayleigh branch point');
nexttile(t); draw_quadratic_q(P, p); title('(b) local finite-model Q scaling');
nexttile(t); draw_strict_zeros(C.strictMode, C, p); title('(c) strict open-channel zeros');
nexttile(t); draw_convergence(D, p); title('(d) two-groove root convergence');
format_figure(fig); export_vector(fig, fullfile(outDir, 'fig2_pole_linewidth_q.pdf')); close(fig);

fig = make_figure([7.0, 5.75]);
t = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(t); draw_scattering_map(C, 1, D, p); title('(a) anomalous-order efficiency');
nexttile(t); draw_scattering_map(C, 2, D, p); title('(b) specular-order efficiency');
nexttile(t); draw_spectral_cuts(C, p); title('(c) fixed-angle cuts through the anomaly');
nexttile(t); draw_representative_field(C, p); title('(d) driven field at the practical route point');
format_figure(fig); export_vector(fig, fullfile(outDir, 'fig3_scattering_fields.pdf')); close(fig);

fig = make_figure([7.0, 5.75]);
t = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(t); draw_two_cavity(D, p); title('(a) two-cavity radiation-control unit cell');
nexttile(t); draw_phasor(C.phasors(1, :), 'A_0', p); title('(b) cancellation of A_0');
nexttile(t); draw_phasor(C.phasors(2, :), 'A_{-1}', p); title('(c) cancellation at the Rayleigh order');
nexttile(t); draw_ablation(single, D, p); title('(d) single-cavity ablation');
nexttile(t); draw_dimensions(D, p); title('(e) final dimensions');
nexttile(t); draw_tolerance(C, p); title('(f) fixed-root fabrication tolerance');
format_figure(fig); export_vector(fig, fullfile(outDir, 'fig4_cancellation_tolerance.pdf')); close(fig);

fig = make_figure([7.0, 5.75]);
t = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(t); draw_signed_routing(C, p); title('(a) reciprocal signed-order routing');
nexttile(t); draw_scattering_matrices(C, p); title('(b) power-normalized S matrices');
nexttile(t); draw_no_go_audit(C, p); title('(c) reciprocity no-go audit');
nexttile(t); draw_regimes(C, D, p); title('(d) distinct operating regimes');
nexttile(t); draw_offrayleigh_convergence(C, p); title('(e) off-Rayleigh convergence');
nexttile(t); draw_offrayleigh_spectrum(C, p); title('(f) separate off-Rayleigh device point');
format_figure(fig); export_vector(fig, fullfile(outDir, 'fig5_routing_device.pdf')); close(fig);

% Prospective supplements are deliberately separate from the five main
% theory figures.  Every panel is labeled as proposed/theory-only so that a
% reader cannot mistake the protocol drawings for measurements.
fig = make_figure([7.0, 5.35]);
t = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(t); draw_proposed_tank(p); title('(a) proposed water-tank geometry (theory plan)');
nexttile(t); draw_proposed_scan(p); title('(b) proposed coarse/fine scan (no measurements)');
nexttile(t); draw_proposed_extraction(p); title('(c) proposed complex-field extraction (theory plan)');
nexttile(t); draw_proposed_ringdown(p); title('(d) proposed ring-down check (schematic only)');
format_figure(fig); export_vector(fig, fullfile(outDir, 'figS1_experiment_protocol.pdf')); close(fig);

fig = make_figure([7.0, 5.35]);
t = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(t); draw_proposed_controls(p); title('(a) proposed control set (theory only)');
nexttile(t); draw_proposed_finite_size(p); title('(b) proposed M=20/40/60 finite-size plan');
nexttile(t); draw_theory_tolerance(C, p); title('(c) theory-only tolerance prediction');
nexttile(t); draw_sound_speed_prediction(D, p); title('(d) theory-only sound-speed prediction');
format_figure(fig); export_vector(fig, fullfile(outDir, 'figS2_controls_tolerance.pdf')); close(fig);

write_manifest(outDir, cacheFile, D, P, C);
fprintf('Generated five vector PDFs and figure_manifest.txt in %s\n', outDir);
end

function assert_cache_compatibility(C, D, P0, single)
requiredCache = {'eta', 'cutEtaAnom', 'cutA0Phase', 'drivenFields', ...
    'fieldResults', 'phasors', 'Splus', 'Sminus', 'matrixReciprocityError', ...
    'matrixUnitarityError', 'toleranceSigma', 'offRayleighNK', ...
    'offRayleighEta', 'offRayleighEtaFrequency', 'offRayleighPoint', ...
    'directionTheta', 'etaPlus', 'etaMinus'};
for j = 1:numel(requiredCache)
    assert(isfield(C, requiredCache{j}), 'v5 cache lacks required field %s', requiredCache{j});
end
requiredDesign = {'xFinal', 'OmegaFinal', 'thetaFinal', 'fTarget', 'aPhysical', ...
    'cWater', 'depthsMm', 'widthsMm', 'gapMm', 'crossTruncations', ...
    'crossSigma', 'rootTruncations', 'rootSigma'};
for j = 1:numel(requiredDesign)
    assert(isfield(D, requiredDesign{j}), 'geometry cache lacks required field %s', requiredDesign{j});
end
assert(isfield(P0, 'P') && isfield(P0.P, 'kappa') && isfield(P0.P, 'Q'), ...
    'pole cache is incompatible');
assert(isfield(single, 'rootTruncations') && isfield(single, 'rootSigma'), ...
    'single-cavity cache is incompatible');
assert(numel(C.drivenFields) >= 3 && numel(C.fieldResults) >= 3, ...
    'v5 cache has no representative practical field');
assert(numel(D.depthsMm) == 2 && numel(D.widthsMm) == 2, ...
    'final geometry is not the expected two-groove design');
end

function fig = make_figure(sizeInches)
fig = figure('Color', 'w', 'Units', 'inches', ...
    'Position', [0.5, 0.5, sizeInches(1), sizeInches(2)], ...
    'PaperUnits', 'inches', 'PaperPosition', [0, 0, sizeInches(1), sizeInches(2)]);
set(fig, 'DefaultAxesFontName', 'Helvetica', 'DefaultTextFontName', 'Helvetica', ...
    'DefaultAxesFontSize', 8, 'DefaultTextFontSize', 8, ...
    'DefaultAxesLineWidth', 0.65, 'DefaultLineLineWidth', 1.15);
end

function format_figure(fig)
axesList = findall(fig, 'Type', 'axes');
for j = 1:numel(axesList)
    ax = axesList(j);
    set(ax, 'FontName', 'Helvetica', 'FontSize', 8, 'LineWidth', 0.65, ...
        'Box', 'on', 'Layer', 'top', 'TickDir', 'out', 'TickLength', [0.018, 0.018], ...
        'XColor', [0.12, 0.14, 0.18], 'YColor', [0.12, 0.14, 0.18], ...
        'TitleFontWeight', 'normal');
    grid(ax, 'on'); ax.GridColor = [0.78, 0.80, 0.83]; ax.GridAlpha = 0.34;
    ax.MinorGridAlpha = 0.16;
    set(get(ax, 'Title'), 'FontName', 'Helvetica', 'FontSize', 8.3, 'FontWeight', 'normal');
    labels = [get(ax, 'XLabel'), get(ax, 'YLabel')];
    set(labels, 'FontName', 'Helvetica', 'FontSize', 8, 'FontWeight', 'normal');
end
legends = findall(fig, 'Type', 'legend');
for j = 1:numel(legends)
    set(legends(j), 'FontName', 'Helvetica', 'FontSize', 7.1, 'Box', 'off', ...
        'TextColor', [0.12, 0.14, 0.18]);
end
colorbars = findall(fig, 'Type', 'colorbar');
for j = 1:numel(colorbars)
    set(colorbars(j), 'FontName', 'Helvetica', 'FontSize', 7.2, 'LineWidth', 0.5);
end
end

function export_vector(fig, fileName)
drawnow;
exportgraphics(fig, fileName, 'ContentType', 'vector', 'BackgroundColor', 'white');
assert(exist(fileName, 'file') == 2, 'Vector export failed: %s', fileName);
end

function p = figure_palette()
p.ink = [0.12, 0.14, 0.18]; p.blue = [0.10, 0.31, 0.55];
p.teal = [0.08, 0.48, 0.46]; p.orange = [0.86, 0.38, 0.08];
p.red = [0.76, 0.15, 0.18]; p.purple = [0.42, 0.29, 0.61];
p.gray = [0.43, 0.47, 0.51]; p.lightGray = [0.86, 0.88, 0.90];
p.lightBlue = [0.91, 0.96, 0.98]; p.white = [1, 1, 1];
end

function draw_unit_cell(D, p)
a = 1e3 * D.aPhysical;
d = D.depthsMm(:).';
w = D.widthsMm(:).';
g = D.gapMm;
occupied = sum(w) + g;
x1 = (a - occupied) / 2;
x2 = x1 + w(1) + g;
rectangle('Position', [0, -max(d) - 0.55, a, max(d) + 0.55], ...
    'FaceColor', p.lightBlue, 'EdgeColor', p.ink, 'LineWidth', 0.8); hold on;
rectangle('Position', [x1, -d(1), w(1), d(1)], 'FaceColor', p.white, ...
    'EdgeColor', p.blue, 'LineWidth', 1.4);
rectangle('Position', [x2, -d(2), w(2), d(2)], 'FaceColor', p.white, ...
    'EdgeColor', p.orange, 'LineWidth', 1.4);
plot([0, a], [0, 0], '-', 'Color', p.ink, 'LineWidth', 1.0);
text(x1 + w(1) / 2, -d(1) / 2, sprintf('d_1=%.2f', d(1)), ...
    'HorizontalAlignment', 'center', 'Color', p.blue, 'FontSize', 7.4);
text(x2 + w(2) / 2, -d(2) / 2, sprintf('d_2=%.2f', d(2)), ...
    'HorizontalAlignment', 'center', 'Color', p.orange, 'FontSize', 7.4);
text(x1 + w(1) / 2, 0.24, sprintf('w_1=%.2f', w(1)), ...
    'HorizontalAlignment', 'center', 'Color', p.blue, 'FontSize', 7.2);
text(x2 + w(2) / 2, 0.24, sprintf('w_2=%.2f', w(2)), ...
    'HorizontalAlignment', 'center', 'Color', p.orange, 'FontSize', 7.2);
text((x1 + w(1) + x2) / 2, 0.56, sprintf('g=%.2f', g), ...
    'HorizontalAlignment', 'center', 'Color', p.ink, 'FontSize', 7.2);
text(a / 2, 0.72, sprintf('a=%.3f mm', a), ...
    'HorizontalAlignment', 'center', 'Color', p.ink, 'FontSize', 7.4);
xlim([-0.15 * a, 1.15 * a]); ylim([-max(d) - 0.65, 1.00]);
xlabel('x (mm)'); ylabel('y (mm)');
end

function draw_topology(kappa0, Omega0, p)
k = linspace(-0.32, 0.32, 500);
plot(k, abs(k - 1), '-', 'Color', p.orange, 'LineWidth', 1.35); hold on;
plot(k, abs(k), '-', 'Color', p.gray, 'LineWidth', 1.0);
plot(k, abs(k + 1), '-', 'Color', p.blue, 'LineWidth', 1.35);
plot(kappa0, Omega0, 'o', 'MarkerSize', 7.5, 'MarkerFaceColor', p.red, ...
    'MarkerEdgeColor', p.white, 'LineWidth', 0.8);
text(kappa0 + 0.014, Omega0 + 0.026, 'strict BIC', 'Color', p.red, 'FontSize', 7.5);
xlabel('\kappa'); ylabel('\Omega=fa/c');
xlim([-0.32, 0.32]); ylim([0.70, 1.15]);
legend('n=-1 Rayleigh', 'n=0 light line', 'n=+1 Rayleigh', 'BIC', ...
    'Location', 'north', 'NumColumns', 2);
end

function draw_channel_count(kappa0, Omega0, p)
orders = -1:1;
kyPlus = real(sqrt(complex(Omega0^2 - (kappa0 + orders).^2)));
kyMinus = real(sqrt(complex(Omega0^2 - (-kappa0 + orders).^2)));
bar(orders - 0.17, kyPlus, 0.32, 'FaceColor', p.blue, 'EdgeColor', 'none'); hold on;
bar(orders + 0.17, kyMinus, 0.32, 'FaceColor', p.orange, 'EdgeColor', 'none');
yline(0, '-', 'Color', p.ink, 'LineWidth', 0.8);
plot([-1, 1], [0, 0], 'o', 'Color', p.red, 'MarkerFaceColor', p.red, 'MarkerSize', 4.5);
text(-0.94, max([kyPlus, kyMinus]) * 0.88, '+\theta', 'Color', p.blue, 'FontSize', 7.3);
text(0.46, max([kyPlus, kyMinus]) * 0.88, '-\theta', 'Color', p.orange, 'FontSize', 7.3);
text(-1, 0.055, 'grazing', 'HorizontalAlignment', 'center', 'Color', p.red, 'FontSize', 7.0);
text(1, 0.055, 'grazing', 'HorizontalAlignment', 'center', 'Color', p.red, 'FontSize', 7.0);
set(gca, 'XTick', orders, 'XTickLabel', {'n=-1', 'n=0', 'n=+1'});
xlabel('Floquet order'); ylabel('Re(k_y a/2\pi)');
legend('+\theta', '-\theta', 'threshold', 'Location', 'northwest');
ylim([-0.02, max([kyPlus, kyMinus]) * 1.25]);
end

function draw_bic_location(theta0, f0, aPhysical, c0, p)
theta = linspace(0, 12, 300);
fRA = c0 ./ (aPhysical * (1 + sind(theta))) / 1e3;
plot(theta, fRA, '-', 'Color', p.blue, 'LineWidth', 1.35); hold on;
plot(theta0, f0 / 1e3, 'o', 'MarkerSize', 7.5, 'MarkerFaceColor', p.red, ...
    'MarkerEdgeColor', p.white, 'LineWidth', 0.8);
plot([theta0, theta0], [min(fRA), f0 / 1e3], ':', 'Color', p.red, 'LineWidth', 0.8);
text(theta0 + 0.25, f0 / 1e3 + 0.75, sprintf('(%.3f deg, %.3f kHz)', theta0, f0 / 1e3), ...
    'Color', p.red, 'FontSize', 7.3);
xlabel('\theta_i (deg)'); ylabel('f (kHz)');
legend('n=-1 Rayleigh threshold', 'strict BIC', 'Location', 'southwest');
xlim([0, 12]); ylim([min(fRA) - 1.2, max(fRA) + 1.2]);
end

function draw_pole_branch(P, kappa0, Omega0, p)
plot(P.kappa, real(P.Omega), 'o-', 'Color', p.blue, 'MarkerFaceColor', p.blue, ...
    'MarkerSize', 3.2, 'LineWidth', 1.15); hold on;
plot(P.kappa, 1 - P.kappa, '--', 'Color', p.orange, 'LineWidth', 1.2);
plot(kappa0, Omega0, 'o', 'MarkerSize', 7.5, 'MarkerFaceColor', p.red, ...
    'MarkerEdgeColor', p.white, 'LineWidth', 0.8);
plot([kappa0 - 0.00036, kappa0 - 0.00003], [Omega0 - 0.00012, Omega0 - 0.00002], ':', ...
    'Color', p.red, 'LineWidth', 0.8);
text(kappa0 - 0.00078, Omega0 - 0.00015, 'BIC', 'Color', p.red, 'FontSize', 7.5);
xlabel('\kappa'); ylabel('Re \Omega');
xlim([min(P.kappa) - 0.00025, max(P.kappa) + 0.00015]);
ylim([min([real(P.Omega(:)); Omega0]) - 0.00003, max([real(P.Omega(:)); Omega0]) + 0.00004]);
legend('leaky pole', 'n=-1 Rayleigh line', 'strict point', 'Location', 'northwest');
end

function draw_quadratic_q(P, p)
valid = isfinite(P.Q) & P.Q > 0 & abs(P.delta_kappa) > 0;
dk = abs(P.delta_kappa(valid)); linewidth = max(abs(imag(P.Omega(valid))), realmin); q = P.Q(valid);
[dk, order] = sort(dk); linewidth = linewidth(order); q = q(order);
yyaxis left;
loglog(dk, linewidth, 'o-', 'Color', p.blue, 'MarkerFaceColor', p.blue, ...
    'MarkerSize', 3.1, 'LineWidth', 1.05); ylabel('-Im \Omega');
yyaxis right;
loglog(dk, q, 's-', 'Color', p.orange, 'MarkerFaceColor', p.orange, ...
    'MarkerSize', 3.1, 'LineWidth', 1.05); ylabel('Q');
xlabel('|\delta\kappa|'); grid on; box on;
fit = polyfit(log(dk), log(q), 1);
text(0.000013, 2.8e9, sprintf('slope = %.2f', fit(1)), 'Color', p.orange, 'FontSize', 7.4);
text(0.000013, 7.5e7, 'local finite model: Q \propto |\delta\kappa|^{-2}', ...
    'Color', p.ink, 'FontSize', 7.2);
text(0.05, 0.92, 'local cached pole track: N/K=121/15', 'Units', 'normalized', ...
    'Color', p.ink, 'FontSize', 7.0);
legend('-Im \Omega', 'Q', 'Location', 'southwest');
end

function draw_strict_zeros(S, C, p)
orders = [-2, -1, 0, 1, 2]; idx = zeros(size(orders));
for j = 1:numel(orders), idx(j) = find(S.full_operator.orders == orders(j), 1); end
values = max(abs(S.mode.A(idx)), 1e-16);
semilogy(orders, values, 'o', 'Color', p.blue, 'MarkerFaceColor', p.blue, ...
    'MarkerSize', 5.2, 'LineWidth', 1.0); hold on;
semilogy([-1, 0], values(2:3), 'o', 'Color', p.red, 'MarkerFaceColor', p.red, 'MarkerSize', 6.8);
yline(1e-8, '--', 'Color', p.gray, 'LineWidth', 0.9);
text(-0.98, 2e-7, 'A_{-1}=0', 'Color', p.red, 'FontSize', 7.2, 'HorizontalAlignment', 'center');
text(0.02, 2e-7, 'A_0=0', 'Color', p.red, 'FontSize', 7.2, 'HorizontalAlignment', 'center');
set(gca, 'XTick', orders, 'XTickLabel', {'-2', '-1', '0', '+1', '+2'});
xlabel('order n'); ylabel('|A_n^{BIC}|'); ylim([1e-16, max(1e-1, 2 * max(values))]);
legend('homogeneous mode', 'strict open channels', '10^{-8} guide', 'Location', 'southwest');
text(0.55, 2e-2, sprintf('E_{ext}(20a)=%.3g; finite', C.exteriorEnergy(end)), ...
    'Color', p.teal, 'FontSize', 7.0);
text(0.55, 7e-3, 'square-integrable; grazing diverges', ...
    'Color', p.ink, 'FontSize', 7.0);
end

function draw_convergence(D, p)
semilogy(D.crossTruncations(:, 2), max(D.crossSigma, 1e-18), 'o-', ...
    'Color', p.gray, 'MarkerFaceColor', p.gray, 'MarkerSize', 3.8, 'LineWidth', 1.0); hold on;
semilogy(D.rootTruncations(:, 2), max(D.rootSigma, 1e-18), 's-', ...
    'Color', p.red, 'MarkerFaceColor', p.red, 'MarkerSize', 4.2, 'LineWidth', 1.1);
yline(1e-8, '--', 'Color', p.orange, 'LineWidth', 0.9);
xlabel('groove truncation K'); ylabel('\sigma_{min}/\sigma_{max}');
set(gca, 'XTick', unique([D.crossTruncations(:, 2); D.rootTruncations(:, 2)])); xtickangle(35);
ylim([1e-17, 0.4]); legend('fixed geometry', 'reoptimized root', 'strict criterion', 'Location', 'southwest');
text(0.07, 0.88, 'two independent open-channel zeros', 'Units', 'normalized', ...
    'Color', p.ink, 'FontSize', 7.2);
end

function draw_scattering_map(C, orderIndex, D, p)
data = C.eta(:, :, orderIndex).';
imagesc(C.thetaValues, C.frequencyValues / 1e3, data); axis xy; hold on;
colormap(gca, parula(256)); clim([0, 1]); cb = colorbar; cb.Label.String = '\eta';
fRA = D.cWater ./ (D.aPhysical * (1 + sind(C.thetaValues))) / 1e3;
plot(C.thetaValues, fRA, '--', 'Color', p.ink, 'LineWidth', 1.1);
plot(C.hotPoint(1), C.hotPoint(2) / 1e3, 'o', 'Color', p.orange, ...
    'MarkerFaceColor', p.orange, 'MarkerEdgeColor', p.white, 'MarkerSize', 5.8);
plot(C.routePoint(1), C.routePoint(2) / 1e3, 's', 'Color', p.blue, ...
    'MarkerFaceColor', p.blue, 'MarkerEdgeColor', p.white, 'MarkerSize', 5.4);
plot(D.thetaFinal, D.fTarget / 1e3, 'o', 'Color', p.red, ...
    'MarkerFaceColor', p.red, 'MarkerEdgeColor', p.white, 'MarkerSize', 5.8);
text(D.thetaFinal + 0.10, D.fTarget / 1e3 + 0.30, 'BIC', 'Color', p.red, 'FontSize', 7.1);
xlabel('\theta_i (deg)'); ylabel('f (kHz)');
if orderIndex == 1
    text(0.04, 0.91, 'dashed: n=-1 Rayleigh line', 'Units', 'normalized', ...
        'Color', p.ink, 'FontSize', 7.1);
end
end

function draw_spectral_cuts(C, p)
cutColors = [p.blue; p.teal; p.orange; p.purple; p.red];
for j = 1:numel(C.cutTheta)
    plot(C.cutFrequency / 1e3, C.cutEtaAnom(j, :), 'Color', cutColors(j, :), ...
        'LineWidth', 1.05); hold on;
    xline(C.cutRayleigh(j) / 1e3, ':', 'Color', cutColors(j, :), 'LineWidth', 0.75, ...
        'HandleVisibility', 'off');
end
xlabel('f (kHz)'); ylabel('\eta_{-1}');
xlim([min(C.cutFrequency), max(C.cutFrequency)] / 1e3); ylim([0, 0.36]);
legend(compose('%.2f deg', C.cutTheta), 'Location', 'northeast', 'NumColumns', 1);
text(0.05, 0.90, 'dotted: Rayleigh thresholds', 'Units', 'normalized', ...
    'Color', p.ink, 'FontSize', 7.1);
end

function draw_representative_field(C, p)
F = C.drivenFields{3}; R = C.fieldResults{3};
level = log10(max(abs(F.p), 1e-6));
imagesc(F.x, F.y, level); axis xy; hold on;
colormap(gca, parula(256)); clim([max(-6, max(level(:)) - 5), max(level(:))]);
cb = colorbar; cb.Label.String = 'log_{10}|p|'; draw_grooves(R, p.white);
strideY = 30; strideX = 27;
X = F.X(1:strideY:end, 1:strideX:end); Y = F.Y(1:strideY:end, 1:strideX:end);
Ix = F.Ix(1:strideY:end, 1:strideX:end); Iy = F.Iy(1:strideY:end, 1:strideX:end);
M = max(sqrt(Ix.^2 + Iy.^2), 1e-12);
quiver(X, Y, Ix ./ M, Iy ./ M, 0.33, 'Color', p.white, 'LineWidth', 0.45, 'MaxHeadSize', 0.55);
idx = order_ids(R); etaTarget = R.eta(idx(1));
fField = C.fieldPoints(3, 2) / 1e3;
text(0.04, 0.92, sprintf('f=%.3f kHz, eta_-1=%.3f', fField, etaTarget), ...
    'Units', 'normalized', 'Color', p.white, 'FontSize', 7.1);
xlabel('x/a'); ylabel('y/a');
end

function draw_two_cavity(D, p)
a = 1; d = D.xFinal(2:3); w = D.xFinal(4:5); g = D.xFinal(6);
x1 = (a - sum(w) - g) / 2; x2 = x1 + w(1) + g;
rectangle('Position', [0, -0.86, 1, 0.86], 'FaceColor', p.lightBlue, ...
    'EdgeColor', p.ink, 'LineWidth', 0.8); hold on;
rectangle('Position', [x1, -d(1), w(1), d(1)], 'FaceColor', p.white, ...
    'EdgeColor', p.blue, 'LineWidth', 1.3);
rectangle('Position', [x2, -d(2), w(2), d(2)], 'FaceColor', p.white, ...
    'EdgeColor', p.orange, 'LineWidth', 1.3);
text(x1 + w(1) / 2, -0.38, 'storage', 'HorizontalAlignment', 'center', ...
    'Color', p.blue, 'FontSize', 7.3);
text(x2 + w(2) / 2, -0.10, 'phase trim', 'HorizontalAlignment', 'center', ...
    'Color', p.orange, 'FontSize', 7.2);
quiver(0.50, 0.18, -0.27, -0.20, 0, 'Color', p.ink, 'LineWidth', 1.1, 'MaxHeadSize', 0.65);
text(0.56, 0.21, 'incident n=0', 'Color', p.ink, 'FontSize', 7.2);
quiver(0.46, 0.02, -0.24, 0.00, 0, 'Color', p.red, 'LineWidth', 1.1, 'MaxHeadSize', 0.65);
text(0.55, 0.03, 'n=-1 grazing', 'Color', p.red, 'FontSize', 7.2);
xlim([-0.08, 1.08]); ylim([-0.93, 0.34]); xlabel('x/a'); ylabel('y/a');
end

function draw_phasor(z, labelText, p)
scale = max(abs(z)); if scale == 0, scale = 1; end; z = z / scale;
quiver(0, 0, real(z(1)), imag(z(1)), 0, 'Color', p.blue, 'LineWidth', 1.65, ...
    'MaxHeadSize', 0.42); hold on;
quiver(real(z(1)), imag(z(1)), real(z(2)), imag(z(2)), 0, 'Color', p.orange, ...
    'LineWidth', 1.65, 'MaxHeadSize', 0.42);
plot([0, real(sum(z))], [0, imag(sum(z))], '--', 'Color', p.red, 'LineWidth', 1.0);
plot(0, 0, 'o', 'Color', p.ink, 'MarkerFaceColor', p.ink, 'MarkerSize', 3.5);
plot(real(sum(z)), imag(sum(z)), 'o', 'Color', p.red, 'MarkerFaceColor', p.red, 'MarkerSize', 4.0);
axis equal; xlim([-1.25, 1.25]); ylim([-1.25, 1.25]);
xlabel('Re contribution'); ylabel('Im contribution');
legend('groove 1', 'groove 2', 'sum', 'Location', 'southwest');
text(0.05, 0.90, sprintf('%s: |sum|/|z_1|=%.1e', labelText, abs(sum(z))), ...
    'Units', 'normalized', 'Color', p.red, 'FontSize', 7.0);
end

function draw_ablation(single, D, p)
semilogy(single.rootTruncations(:, 2), max(single.rootSigma, 1e-8), 'o-', ...
    'Color', p.gray, 'MarkerFaceColor', p.gray, 'MarkerSize', 3.7, 'LineWidth', 1.0); hold on;
semilogy(D.rootTruncations(:, 2), max(D.rootSigma, 1e-18), 's-', ...
    'Color', p.red, 'MarkerFaceColor', p.red, 'MarkerSize', 3.9, 'LineWidth', 1.1);
yline(1e-8, '--', 'Color', p.orange, 'LineWidth', 0.9);
xlabel('groove truncation K'); ylabel('\sigma_{min}/\sigma_{max}'); ylim([1e-17, 0.3]);
legend('single cavity', 'two cavities', 'strict criterion', 'Location', 'southwest');
text(0.08, 0.88, 'one cavity retains a finite residual', 'Units', 'normalized', ...
    'Color', p.gray, 'FontSize', 7.1);
end

function draw_dimensions(D, p)
vals = [D.depthsMm(:), D.widthsMm(:)];
bar(vals, 'grouped'); hold on; yline(1, '--', 'Color', p.red, 'LineWidth', 0.9);
set(gca, 'XTick', 1:2, 'XTickLabel', {'groove 1', 'groove 2'});
xlabel('cavity'); ylabel('dimension (mm)'); ylim([0, 5.2]);
legend('depth d', 'width w', '1 mm lower bound', 'Location', 'northwest');
text(0.05, 0.16, sprintf('a=%.3f mm,  g=%.3f mm', D.aPhysical * 1e3, D.gapMm), ...
    'Units', 'normalized', 'Color', p.ink, 'FontSize', 7.2);
end

function draw_tolerance(C, p)
score = -log10(max(C.toleranceSigma, 1e-16));
imagesc(C.errorValuesMm, C.errorValuesMm, score); axis xy; hold on;
colormap(gca, parula(256)); clim([1.5, 16]); cb = colorbar; cb.Label.String = '-log_{10}(\sigma ratio)';
contour(C.errorValuesMm, C.errorValuesMm, score, [4, 6, 10], ...
    'Color', p.ink, 'LineWidth', 0.45);
plot(0, 0, 'o', 'Color', p.red, 'MarkerFaceColor', 'none', 'MarkerSize', 5.2);
xlabel('\Delta w_2 (mm)'); ylabel('\Delta d_2 (mm)');
text(0.06, 0.91, 'N/K=121/15; fixed operating point', 'Units', 'normalized', ...
    'Color', p.ink, 'FontSize', 7.0);
text(0.06, 0.83, sprintf('center: %.1e', C.toleranceSigma(16, 16)), 'Units', 'normalized', ...
    'Color', p.red, 'FontSize', 7.0);
end

function draw_signed_routing(C, p)
plot(C.directionTheta, C.etaPlus, '-', 'Color', p.blue, 'LineWidth', 1.15); hold on;
plot(C.directionTheta, C.etaMinus, '--', 'Color', p.orange, 'LineWidth', 1.15);
plot(C.directionTheta, C.eta0Plus, ':', 'Color', p.gray, 'LineWidth', 1.0);
plot(C.directionTheta, C.eta0Minus, ':', 'Color', p.gray, 'LineWidth', 1.0, 'HandleVisibility', 'off');
xlabel('|\theta_i| (deg)'); ylabel('normal-power efficiency'); ylim([0, 0.34]);
legend('+\theta\rightarrow n=-1', '-\theta\rightarrow n=+1', 'n=0', 'Location', 'south');
text(0.05, 0.90, sprintf('f=%.3f kHz', C.directionFrequency / 1e3), ...
    'Units', 'normalized', 'Color', p.ink, 'FontSize', 7.2);
text(0.05, 0.80, 'reciprocal: direction swaps signed order', 'Units', 'normalized', ...
    'Color', p.red, 'FontSize', 7.1);
end

function draw_scattering_matrices(C, p)
M = [abs(C.Splus), nan(2, 1), abs(C.Sminus.')];
imagesc(M); axis image; clim([0, 1]); hold on; colormap(gca, parula(256)); colorbar;
set(gca, 'XTick', [1, 2, 4, 5], 'XTickLabel', {'S+ (0)', 'S+ (-1)', 'S-^T (0)', 'S-^T (+1)'}, ...
    'YTick', [1, 2], 'YTickLabel', {'out 1', 'out 2'});
for r = 1:2
    for c = [1, 2, 4, 5]
        value = M(r, c);
        text(c, r, sprintf('%.3f', value), 'HorizontalAlignment', 'center', ...
            'Color', ternary_color(value), 'FontSize', 7.0, 'FontWeight', 'bold');
    end
end
xlabel('input port representation'); ylabel('output port');
text(0.05, 0.05, sprintf('||S_+-S_-^T||_F=%.1e', C.matrixReciprocityError), ...
    'Units', 'normalized', 'Color', p.ink, 'FontSize', 7.0);
end

function draw_no_go_audit(C, p)
residuals = [C.matrixReciprocityError, C.matrixUnitarityError(:).'];
bars = -log10(max(residuals, 1e-18));
bar(1:3, bars, 0.56, 'FaceColor', p.blue, 'EdgeColor', 'none'); hold on;
for j = 1:3
    text(j, bars(j) + 0.35, sprintf('10^{-%d}', round(bars(j))), ...
        'HorizontalAlignment', 'center', 'Color', p.ink, 'FontSize', 7.0);
end
set(gca, 'XTick', 1:3, 'XTickLabel', {'reciprocity', 'unitarity +', 'unitarity -'});
xlabel('audit'); ylabel('-log_{10}(residual)'); ylim([0, max(bars) + 2.2]);
text(0.05, 0.17, 'reciprocal port basis', 'Units', 'normalized', ...
    'Color', p.red, 'FontSize', 7.3);
text(0.05, 0.09, 'no one-way isolation', 'Units', 'normalized', ...
    'Color', p.ink, 'FontSize', 7.2);
end

function draw_regimes(C, D, p)
theta = linspace(0, 40, 500);
fRA = D.cWater ./ (D.aPhysical * (1 + sind(theta))) / 1e3;
plot(theta, fRA, '--', 'Color', p.gray, 'LineWidth', 1.0); hold on;
plot(D.thetaFinal, D.fTarget / 1e3, 'o', 'Color', p.red, 'MarkerFaceColor', p.red, ...
    'MarkerEdgeColor', p.white, 'MarkerSize', 6.0);
plot(C.routePoint(1), C.routePoint(2) / 1e3, 's', 'Color', p.blue, 'MarkerFaceColor', p.blue, ...
    'MarkerEdgeColor', p.white, 'MarkerSize', 5.8);
plot(C.offRayleighPoint(1), C.offRayleighPoint(2) / 1e3, 'd', 'Color', p.teal, ...
    'MarkerFaceColor', p.teal, 'MarkerEdgeColor', p.white, 'MarkerSize', 6.0);
offEta = nearest_value(C.offRayleighFrequency, C.offRayleighEtaFrequency, C.offRayleighPoint(2));
text(D.thetaFinal + 0.8, D.fTarget / 1e3 + 3.0, 'strict BIC', 'Color', p.red, 'FontSize', 7.0);
text(C.routePoint(1) + 0.8, C.routePoint(2) / 1e3 - 4.2, 'near-RA route', 'Color', p.blue, 'FontSize', 7.0);
text(C.offRayleighPoint(1) - 10.2, C.offRayleighPoint(2) / 1e3 + 3.2, ...
    sprintf('off-RA, eta=%.4f', offEta), 'Color', p.teal, 'FontSize', 7.0);
xlabel('|\theta_i| (deg)'); ylabel('f (kHz)'); xlim([0, 40]); ylim([130, 220]);
legend('n=-1 Rayleigh threshold', 'strict BIC', 'near-RA route', 'off-RA point', ...
    'Location', 'southwest');
end

function draw_offrayleigh_convergence(C, p)
x = 1:size(C.offRayleighNK, 1);
plot(x, C.offRayleighEta, 'o-', 'Color', p.teal, 'MarkerFaceColor', p.teal, ...
    'MarkerSize', 4.0, 'LineWidth', 1.05); hold on;
plot(x, C.offRayleighEta0, 's-', 'Color', p.gray, 'MarkerFaceColor', p.gray, ...
    'MarkerSize', 3.7, 'LineWidth', 1.0);
yline(0.99, '--', 'Color', p.orange, 'LineWidth', 0.9);
set(gca, 'XTick', x, 'XTickLabel', compose('%d/%d', C.offRayleighNK(:, 1), C.offRayleighNK(:, 2)));
xtickangle(35); xlabel('modal truncation N/K'); ylabel('normal-power efficiency'); ylim([0.6, 1.04]);
legend('target anomalous order', 'specular n=0', '99% guide', 'Location', 'southwest');
text(0.05, 0.90, sprintf('point value eta=%.4f', C.offRayleighEta(end)), ...
    'Units', 'normalized', 'Color', p.teal, 'FontSize', 7.2);
end

function draw_offrayleigh_spectrum(C, p)
plot(C.offRayleighFrequency / 1e3, C.offRayleighEtaFrequency, '-', ...
    'Color', p.teal, 'LineWidth', 1.15); hold on;
plot(C.offRayleighFrequency / 1e3, C.offRayleighEta0Frequency, '-', ...
    'Color', p.gray, 'LineWidth', 1.0);
yline(0.99, '--', 'Color', p.orange, 'LineWidth', 0.8);
idx = nearest_index(C.offRayleighFrequency, C.offRayleighPoint(2));
plot(C.offRayleighPoint(2) / 1e3, C.offRayleighEtaFrequency(idx), 'd', ...
    'Color', p.red, 'MarkerFaceColor', p.red, 'MarkerSize', 5.8);
xlabel('f (kHz)'); ylabel('normal-power efficiency'); ylim([0.6, 1.04]);
legend('target n=-1', 'specular n=0', '99% guide', 'device point', 'Location', 'southwest');
text(0.05, 0.90, sprintf('eta=%.4f at %.3f kHz', C.offRayleighEtaFrequency(idx), C.offRayleighPoint(2) / 1e3), ...
    'Units', 'normalized', 'Color', p.red, 'FontSize', 7.0);
text(0.05, 0.80, 'off-Rayleigh; not BIC-enhanced', 'Units', 'normalized', ...
    'Color', p.ink, 'FontSize', 7.1);
end

function draw_grooves(R, color)
yline(0, '-', 'Color', color, 'LineWidth', 0.9);
for ell = 1:numel(R.widths)
    rectangle('Position', [R.xleft(ell), -R.depths(ell), R.widths(ell), R.depths(ell)], ...
        'EdgeColor', color, 'LineWidth', 0.9);
end
end

function ids = order_ids(R)
ids = [find(R.orders == -1, 1), find(R.orders == 0, 1), find(R.orders == 1, 1)];
end

function value = nearest_value(x, y, target)
value = y(nearest_index(x, target));
end

function idx = nearest_index(x, target)
[~, idx] = min(abs(x - target));
end

function color = ternary_color(value)
if value > 0.55, color = 'k'; else, color = 'w'; end
end

function write_manifest(outDir, cacheFile, D, P, C)
manifest = fullfile(outDir, 'figure_manifest.txt');
fid = fopen(manifest, 'w');
assert(fid >= 0, 'Could not create figure manifest: %s', manifest);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Theory figure manifest (vector PDF export)\n');
fprintf(fid, 'Data source: %s\n', cacheFile);
fprintf(fid, 'Geometry: a=%.6g mm, f0=%.6g kHz, theta0=%.6g deg, N/K field=313/39\n', ...
    D.aPhysical * 1e3, D.fTarget / 1e3, D.thetaFinal);
fprintf(fid, 'Numerical truncation: maps N/K=61/9; route fields and reciprocal audit N/K=313/39;\n');
fprintf(fid, '  strict homogeneous mode N/K=313/39; pole track uses cached continuation (N/K=121/15);\n');
fprintf(fid, '  tolerance uses the cached (121,15) strict root; off-Rayleigh spectra use N/K=201/25;\n');
fprintf(fid, '  off-Rayleigh convergence reports N/K = %s.\n', ...
    strjoin(compose('%d/%d', C.offRayleighNK(:, 1), C.offRayleighNK(:, 2)), ', '));
fprintf(fid, '\nfig1_structure_topology.pdf\n');
fprintf(fid, '  (a) final two-groove unit cell and dimensions; (b) n=-1,0,+1 Rayleigh topology;\n');
fprintf(fid, '  (c) open/threshold channel count for opposite incidence signs; (d) strict point on the Rayleigh curve.\n');
fprintf(fid, '\nfig2_pole_linewidth_q.pdf\n');
fprintf(fid, '  (a) cached leaky pole approaching the n=-1 branch point; (b) local finite-model Q scaling;\n');
fprintf(fid, '  (c) homogeneous strict open-channel zeros with finite exterior-energy/square-integrability audit;\n');
fprintf(fid, '  (d) fixed-geometry versus reoptimized convergence.\n');
fprintf(fid, '\nfig3_scattering_fields.pdf\n');
fprintf(fid, '  (a,b) cached angle-frequency eta maps with Rayleigh threshold; (c) five fixed-angle cuts;\n');
fprintf(fid, '  (d) cached driven pressure field at the practical Rayleigh route point.\n');
fprintf(fid, '\nfig4_cancellation_tolerance.pdf\n');
fprintf(fid, '  (a) two-cavity geometry; (b,c) complex radiation phasor cancellation for A0 and A-1;\n');
fprintf(fid, '  (d) single-cavity ablation; (e) dimensions and 1 mm bound; (f) cached fixed-root tolerance map.\n');
fprintf(fid, '\nfig5_routing_device.pdf\n');
fprintf(fid, '  (a) reciprocal signed-order routing; (b) power-normalized S matrices; (c) reciprocity/energy no-go audit;\n');
fprintf(fid, '  (d) strict BIC, near-RA route, and separate off-Rayleigh point; (e,f) off-Rayleigh convergence and spectrum.\n');
fprintf(fid, '  The off-Rayleigh point is labeled as not BIC-enhanced. No experimental panels are included.\n');
fprintf(fid, '\nfigS1_experiment_protocol.pdf\n');
fprintf(fid, '  (a) proposed tank geometry; (b) proposed scan sequence; (c) proposed complex-field to Floquet extraction;\n');
fprintf(fid, '  (d) proposed ring-down check. All panels are schematics/theory plans; no experimental data are shown.\n');
fprintf(fid, '\nfigS2_controls_tolerance.pdf\n');
fprintf(fid, '  (a) proposed full/single/blocked-cavity controls; (b) proposed M=20/40/60 finite-size plan;\n');
fprintf(fid, '  (c) cached theory-only tolerance prediction; (d) theory-only sound-speed prediction.\n');
fprintf(fid, '\nPole diagnostics: cached finite-point local finite-model log-log fit gives Q slope approximately %.4f.\n', pole_q_slope(P));
end

function slope = pole_q_slope(P)
valid = isfinite(P.Q) & P.Q > 0 & abs(P.delta_kappa) > 0;
fit = polyfit(log(abs(P.delta_kappa(valid))), log(P.Q(valid)), 1);
slope = fit(1);
end

function draw_proposed_tank(p)
axis off; hold on;
rectangle('Position', [0.04, 0.08, 0.92, 0.82], 'FaceColor', p.lightBlue, ...
    'EdgeColor', p.blue, 'LineWidth', 1.0);
rectangle('Position', [0.35, 0.29, 0.38, 0.07], 'FaceColor', p.lightGray, ...
    'EdgeColor', p.ink, 'LineWidth', 0.8);
for j = 0:7
    rectangle('Position', [0.37 + 0.043 * j, 0.30, 0.020, 0.045], ...
        'FaceColor', p.white, 'EdgeColor', 'none');
end
rectangle('Position', [0.10, 0.52, 0.10, 0.18], 'Curvature', [0.2, 0.2], ...
    'FaceColor', p.orange, 'EdgeColor', p.ink, 'LineWidth', 0.8);
quiver(0.20, 0.61, 0.19, -0.16, 0, 'Color', p.blue, 'LineWidth', 1.4, 'MaxHeadSize', 0.55);
plot(0.62 + 0.23 * cosd(linspace(10, 170, 80)), 0.43 + 0.23 * sind(linspace(10, 170, 80)), ...
    '--', 'Color', p.red, 'LineWidth', 1.0);
plot(0.62, 0.66, 'o', 'Color', p.red, 'MarkerFaceColor', p.red, 'MarkerSize', 4.5);
text(0.07, 0.77, 'proposed source', 'Color', p.blue, 'FontSize', 7.2);
text(0.49, 0.22, 'periodic sample', 'Color', p.ink, 'FontSize', 7.2);
text(0.66, 0.78, 'proposed hydrophone arc', 'Color', p.red, 'FontSize', 7.0);
text(0.06, 0.13, 'schematic; no tank data are shown', 'Color', p.ink, 'FontSize', 7.2);
xlim([0, 1]); ylim([0, 1]);
end

function draw_proposed_scan(p)
axis off; hold on;
steps = {'1  proposed ToF / sound-speed calibration', ...
    '2  proposed coarse f--theta map', '3  proposed fine scan near Rayleigh', ...
    '4  proposed complex p(x,y) raster'};
for j = 1:4
    y = 0.86 - (j - 1) * 0.20;
    rectangle('Position', [0.07, y - 0.09, 0.86, 0.12], 'Curvature', 0.03, ...
        'EdgeColor', p.gray, 'FaceColor', [0.97, 0.97, 0.97], 'LineWidth', 0.7);
    text(0.11, y - 0.034, steps{j}, 'Color', p.ink, 'FontSize', 7.1);
    if j < 4
        quiver(0.50, y - 0.10, 0, -0.065, 0, 'Color', p.red, 'LineWidth', 1.1, 'MaxHeadSize', 0.8);
    end
end
text(0.50, 0.05, 'proposed sequence only; no measured scan is implied', ...
    'HorizontalAlignment', 'center', 'Color', p.red, 'FontSize', 7.1);
xlim([0, 1]); ylim([0, 1]);
end

function draw_proposed_extraction(p)
axis off; hold on;
boxX = [0.05, 0.36, 0.67];
labels = {'complex p(x,y)', 'Floquet fit', 'A_n, eta_n'};
edge = [p.blue; p.orange; p.teal];
for j = 1:3
    rectangle('Position', [boxX(j), 0.42, 0.25, 0.22], 'Curvature', 0.05, ...
        'FaceColor', [0.97, 0.98, 0.99], 'EdgeColor', edge(j, :), 'LineWidth', 1.0);
    text(boxX(j) + 0.125, 0.53, labels{j}, 'HorizontalAlignment', 'center', ...
        'Color', p.ink, 'FontSize', 7.1);
    if j < 3
        quiver(boxX(j) + 0.25, 0.53, 0.055, 0, 0, 'Color', p.red, 'LineWidth', 1.1, 'MaxHeadSize', 0.8);
    end
end
text(0.50, 0.27, 'proposed least-squares reconstruction of complex fields', ...
    'HorizontalAlignment', 'center', 'Color', p.blue, 'FontSize', 7.2);
text(0.50, 0.16, 'eta_n = Re(k_{y,n}) |A_n|^2 / k_{y,inc}', ...
    'HorizontalAlignment', 'center', 'Color', p.ink, 'FontSize', 7.2);
text(0.50, 0.08, 'theory protocol; no experimental A_n are plotted', ...
    'HorizontalAlignment', 'center', 'Color', p.red, 'FontSize', 7.0);
xlim([0, 1]); ylim([0, 1]);
end

function draw_proposed_ringdown(p)
t = linspace(0, 5, 160);
trace = exp(-0.75 * t) .* cos(2 * pi * 1.2 * t);
plot(t, trace, 'Color', p.blue, 'LineWidth', 1.15); hold on;
plot(t, exp(-0.75 * t), '--', 'Color', p.orange, 'LineWidth', 0.9);
plot(t, -exp(-0.75 * t), '--', 'Color', p.orange, 'LineWidth', 0.9);
xlabel('proposed time window (arb. units)'); ylabel('normalized p(t) (schematic)');
legend('proposed ring-down trace', 'envelope', 'Location', 'northeast');
text(0.05, 0.12, 'theory/schematic only; no ring-down data', 'Units', 'normalized', ...
    'Color', p.red, 'FontSize', 7.1);
end

function draw_proposed_controls(p)
axis off; hold on;
names = {'full two-cavity', 'single cavity', 'blocked small cavity'};
edge = [p.red; p.gray; p.orange];
for j = 1:3
    y = 0.71 - (j - 1) * 0.27;
    rectangle('Position', [0.07, y, 0.86, 0.15], 'FaceColor', [0.97, 0.97, 0.97], ...
        'EdgeColor', edge(j, :), 'LineWidth', 1.1);
    rectangle('Position', [0.18, y + 0.03, 0.39, 0.08], 'FaceColor', p.white, ...
        'EdgeColor', p.blue, 'LineWidth', 0.8);
    if j ~= 2
        rectangle('Position', [0.65, y + 0.05, 0.12, 0.06], 'FaceColor', p.white, ...
            'EdgeColor', p.orange, 'LineWidth', 0.8);
    end
    text(0.82, y + 0.075, names{j}, 'HorizontalAlignment', 'right', ...
        'Color', edge(j, :), 'FontSize', 7.0);
end
text(0.50, 0.10, 'proposed controls isolate storage, phase trim, and cavity removal', ...
    'HorizontalAlignment', 'center', 'Color', p.ink, 'FontSize', 7.1);
xlim([0, 1]); ylim([0, 1]);
end

function draw_proposed_finite_size(p)
theta = linspace(-90, -55, 700);
M = [20, 40, 60];
cols = [p.blue; p.orange; p.teal];
for j = 1:3
    psi = 0.12 * (sind(theta) - sind(-72));
    den = sin(psi / 2); num = sin(M(j) * psi / 2);
    af = ones(size(psi)); idx = abs(den) > 1e-12; af(idx) = abs(num(idx) ./ (M(j) * den(idx))).^2;
    plot(theta, af, 'Color', cols(j, :), 'LineWidth', 1.05); hold on;
end
xlabel('proposed output angle (deg)'); ylabel('normalized array factor (theory)'); ylim([0, 1.05]);
legend('M=20', 'M=40', 'M=60', 'Location', 'northwest');
text(0.05, 0.10, 'proposed finite-size check; curves are theoretical guides', ...
    'Units', 'normalized', 'Color', p.red, 'FontSize', 7.0);
end

function draw_theory_tolerance(C, p)
score = -log10(max(C.toleranceSigma, 1e-16));
imagesc(C.errorValuesMm, C.errorValuesMm, score); axis xy; hold on;
colormap(gca, parula(256)); clim([1.5, 16]); cb = colorbar; cb.Label.String = '-log_{10}(\sigma ratio)';
contour(C.errorValuesMm, C.errorValuesMm, score, [4, 6, 10], ...
    'Color', p.ink, 'LineWidth', 0.45);
plot(0, 0, 'o', 'Color', p.red, 'MarkerFaceColor', 'none', 'MarkerSize', 5.0);
xlabel('\Delta w_2 (mm)'); ylabel('\Delta d_2 (mm)');
text(0.05, 0.91, 'cached theoretical prediction, N/K=121/15', 'Units', 'normalized', ...
    'Color', p.ink, 'FontSize', 7.0);
text(0.05, 0.82, 'not a fabrication dataset', 'Units', 'normalized', ...
    'Color', p.red, 'FontSize', 7.1);
end

function draw_sound_speed_prediction(D, p)
cValues = linspace(1450, 1530, 161);
fShift = D.fTarget * cValues / D.cWater / 1e3;
plot(cValues, fShift, '-', 'Color', p.blue, 'LineWidth', 1.15); hold on;
plot(D.cWater, D.fTarget / 1e3, 'o', 'Color', p.red, 'MarkerFaceColor', p.red, 'MarkerSize', 5.0);
xlabel('water sound speed (m/s)'); ylabel('predicted f_{BIC} (kHz)');
legend('theory scaling f\propto c', 'reference medium', 'Location', 'northwest');
text(0.05, 0.11, 'theory-only environmental tuning; no measurements', 'Units', 'normalized', ...
    'Color', p.red, 'FontSize', 7.0);
end
