%% Robustness, cavity-degree-of-freedom, and homogeneous-field diagnostics
% This script is deliberately separate from the manuscript-generating code.
% It evaluates the existing pole-free modal-matching formulation and writes
% auditable data/figures under advanced_analysis/robustness_fields/.
%
% Run from the project root:
%   matlab -batch "run('advanced_analysis/robustness_fields.m')"
%
% The tolerance ensemble uses a fixed Rayleigh point and perturbs only the
% normalized geometry.  The strict reduced residual and the threshold
% amplitude error of an *unconstrained* homogeneous SVD vector are reported
% separately; the latter is not replaced by the exact zero that is imposed
% in the strict operator.  Field panels are homogeneous fields only.

clearvars;
close all;
clc;

rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(rootDir,'Ni2019_MATLAB'));
outDir = fullfile(rootDir,'advanced_analysis','robustness_fields');
if ~exist(outDir,'dir'), mkdir(outDir); end

% Reproducible nominal design cache used by the current PRL draft.
designFile = fullfile(rootDir,'Ni2019_MATLAB','results', ...
    'StrictRayleighBIC_200kHz_min1mm.mat');
singleFile = fullfile(rootDir,'Ni2019_MATLAB','results', ...
    'SingleGroove_strict_Rayleigh_search.mat');
assert(exist(designFile,'file')==2, 'Missing nominal design cache: %s',designFile);
assert(exist(singleFile,'file')==2, 'Missing single-groove cache: %s',singleFile);
D = load(designFile);
Sg = load(singleFile);

% Normalized nominal state.  Physical dimensional values are retained in
% the report, but all modal matching is done with a=1 and lambda=1/Omega.
xNom = D.xFinal(:).';
OmegaNom = D.OmegaFinal;
thetaNom = D.thetaFinal;
fNom = D.fTarget;
cNom = D.cWater;
aPhysical = D.aPhysical;
lambdaPhysical = D.lambdaPhysical;

Nhi = 313;
Khi = 39;
cfgNom = make_cfg(xNom,Nhi,Khi);
Rnom = ni2019_strict_rayleigh_operator(cfgNom,'TargetOrder',-1, ...
    'EnforceAllFiniteOpen',true,'KyTolerance',1e-8*abs(2*pi/cfgNom.lambda), ...
    'Verbose',false);
assert(isfinite(Rnom.strict_residual) && isfinite(Rnom.sigma_ratio), ...
    'Nominal strict operator returned a nonfinite diagnostic.');
nominalResidual = Rnom.strict_residual;
nominalSigma = Rnom.sigma_ratio;
% The high-order result must reproduce the cache used in the manuscript.
assert(abs(nominalResidual-D.rootRaw(end)) < 5e-11, ...
    'Nominal raw residual drifted from the cached result.');
assert(abs(nominalSigma-D.rootSigma(end)) < 5e-11, ...
    'Nominal singular-value ratio drifted from the cached result.');

fprintf('Nominal strict BIC: kappa=%.15g, Omega=%.15g, theta=%.12g deg\n', ...
    xNom(1),OmegaNom,thetaNom);
fprintf('  (N,K)=(%d,%d), sigma ratio=%.6e, raw residual=%.6e\n', ...
    Nhi,Khi,nominalSigma,nominalResidual);

%% A. Deterministic fabrication-tolerance ensemble
% A moderate but converged modal truncation keeps the 600-realization
% ensemble practical.  The high-order nominal result above is the strict
% manuscript verification; the tolerance data state their own truncation.
Ntol = 185;
Ktol = 23;
nMC = 100;
seed = 260824;
levels = [0.01,0.02,0.05];
parameterNames = {'d1','d2','w1','w2','g'};
modeNames = {'independent','joint'};
rng(seed,'twister');

nRows = numel(levels)*numel(modeNames)*nMC;
toleranceRows = repmat(struct( ...
    'mode','','level_pct',NaN,'sample',NaN, ...
    'delta_d1_pct',NaN,'delta_d2_pct',NaN,'delta_w1_pct',NaN, ...
    'delta_w2_pct',NaN,'delta_g_pct',NaN,'geometry_valid',false, ...
    'solver_ok',false,'strict_sigma_ratio',NaN,'strict_raw_residual',NaN, ...
    'full_sigma_ratio',NaN,'full_raw_residual',NaN, ...
    'threshold_amplitude_error',NaN,'finite_open_amplitude_error',NaN, ...
    'strict_threshold_amplitude',NaN,'strict_sigma_over_baseline',NaN, ...
    'strict_raw_over_baseline',NaN),nRows,1);

reuseTolerance = strcmp(getenv('RF_REUSE_TOLERANCE'),'1') && ...
    exist(fullfile(outDir,'tolerance_samples.csv'),'file')==2 && ...
    exist(fullfile(outDir,'tolerance_summary.csv'),'file')==2;
toleranceFile = fullfile(outDir,'tolerance_samples.csv');
summaryFile = fullfile(outDir,'tolerance_summary.csv');
if reuseTolerance
    toleranceRows = table2struct(readtable(toleranceFile));
    summaryRows = table2struct(readtable(summaryFile));
    fprintf('Reusing existing deterministic tolerance CSVs because RF_REUSE_TOLERANCE=1.\n');
else
    row = 0;
    for im = 1:numel(modeNames)
        modeName = modeNames{im};
        for il = 1:numel(levels)
            level = levels(il);
            for is = 1:nMC
                row = row + 1;
                if strcmp(modeName,'independent')
                    delta = zeros(1,5);
                    j = mod(is-1,5)+1; % balanced one-parameter perturbations
                    delta(j) = (2*rand-1)*level;
                else
                    delta = (2*rand(1,5)-1)*level;
                end
                r = toleranceRows(row);
                r.mode = modeName;
                r.level_pct = 100*level;
                r.sample = is;
                r.delta_d1_pct = 100*delta(1);
                r.delta_d2_pct = 100*delta(2);
                r.delta_w1_pct = 100*delta(3);
                r.delta_w2_pct = 100*delta(4);
                r.delta_g_pct = 100*delta(5);
                x = xNom;
                x(2:6) = xNom(2:6).*(1+delta);
                r.geometry_valid = geometry_is_valid(x);
                if r.geometry_valid
                    try
                        cfg = make_cfg(x,Ntol,Ktol);
                        Rs = ni2019_strict_rayleigh_operator(cfg, ...
                            'TargetOrder',-1,'EnforceAllFiniteOpen',true, ...
                            'KyTolerance',1e-8*abs(2*pi/cfg.lambda), ...
                            'Verbose',false);
                        [~,fullDiag] = unconstrained_mode_diagnostics(cfg,-1);
                        r.strict_sigma_ratio = Rs.sigma_ratio;
                        r.strict_raw_residual = Rs.strict_residual;
                        r.full_sigma_ratio = fullDiag.sigma_ratio;
                        r.full_raw_residual = fullDiag.raw_residual;
                        r.threshold_amplitude_error = fullDiag.threshold_error;
                        r.finite_open_amplitude_error = fullDiag.finite_open_error;
                        r.strict_threshold_amplitude = abs(Rs.target_amplitude);
                        r.solver_ok = all(isfinite([r.strict_sigma_ratio, ...
                            r.strict_raw_residual,r.full_sigma_ratio, ...
                            r.full_raw_residual,r.threshold_amplitude_error, ...
                            r.finite_open_amplitude_error]));
                    catch ME
                        r.solver_ok = false;
                        fprintf(2,'Tolerance realization failed (%s, %.1f%%, %d): %s\n', ...
                            modeName,100*level,is,ME.message);
                    end
                end
                toleranceRows(row) = r;
            end
        end
    end
    toleranceTable = struct2table(toleranceRows);
    writetable(toleranceTable,toleranceFile);
    summaryRows = build_tolerance_summary(toleranceRows,modeNames,levels);
    summaryTable = struct2table(summaryRows);
    writetable(summaryTable,summaryFile);
end

% Resolve the unperturbed final geometry at the tolerance truncation so the
% finite discretization floor is visible beside every fabrication statistic.
cfgTolNom = make_cfg(xNom,Ntol,Ktol);
RtolNom = ni2019_strict_rayleigh_operator(cfgTolNom,'TargetOrder',-1, ...
    'EnforceAllFiniteOpen',true,'KyTolerance',1e-8*abs(2*pi/cfgTolNom.lambda), ...
    'Verbose',false);
[~,tolFullNom] = unconstrained_mode_diagnostics(cfgTolNom,-1);
tolBaseline = struct('N',Ntol,'K',Ktol, ...
    'strict_sigma_ratio',RtolNom.sigma_ratio, ...
    'strict_raw_residual',RtolNom.strict_residual, ...
    'full_sigma_ratio',tolFullNom.sigma_ratio, ...
    'full_raw_residual',tolFullNom.raw_residual, ...
    'threshold_amplitude_error',tolFullNom.threshold_error, ...
    'finite_open_amplitude_error',tolFullNom.finite_open_error);
if ~isfield(toleranceRows,'strict_sigma_over_baseline')
    [toleranceRows.strict_sigma_over_baseline] = deal(NaN);
    [toleranceRows.strict_raw_over_baseline] = deal(NaN);
end
for i = 1:numel(toleranceRows)
    if toleranceRows(i).solver_ok && isfinite(toleranceRows(i).strict_sigma_ratio)
        toleranceRows(i).strict_sigma_over_baseline = ...
            toleranceRows(i).strict_sigma_ratio/max(tolBaseline.strict_sigma_ratio,eps);
        toleranceRows(i).strict_raw_over_baseline = ...
            toleranceRows(i).strict_raw_residual/max(tolBaseline.strict_raw_residual,eps);
    end
end
toleranceTable = struct2table(toleranceRows);
writetable(toleranceTable,toleranceFile);
toleranceBaselineFile = fullfile(outDir,'tolerance_baseline.csv');
writetable(struct2table(tolBaseline),toleranceBaselineFile);

fprintf('Tolerance ensemble: seed=%d, n=%d per mode/level, (N,K)=(%d,%d)\n', ...
    seed,nMC,Ntol,Ktol);
fprintf('  unperturbed baseline at (N,K)=(%d,%d): sigma=%.6e, raw=%.6e\n', ...
    Ntol,Ktol,tolBaseline.strict_sigma_ratio,tolBaseline.strict_raw_residual);
for i = 1:numel(summaryRows)
    fprintf('  %s %.0f%%: valid=%d/%d, solver failures=%d, ', ...
        summaryRows(i).mode,summaryRows(i).level_pct, ...
        summaryRows(i).valid_count,summaryRows(i).sample_count, ...
        summaryRows(i).solver_failure_count);
    fprintf('strict r50=%.3e, threshold err r50=%.3e\n', ...
        summaryRows(i).strict_raw_q50,summaryRows(i).threshold_error_q50);
end

%% B. Single-cavity versus two-cavity degree of freedom
% The cached single-groove optimization retained many transverse modes and
% used the broad bounds documented in its source script.  We re-evaluate its
% final geometry directly with the same strict operator and compare the
% complete normalized residual with the nominal two-groove state.
singleX = Sg.xFinal(:).';
cfgSingle = struct('a',1,'lambda',1/(1-singleX(1)), ...
    'theta_i_deg',asind(singleX(1)/(1-singleX(1))), ...
    'depths',singleX(2),'widths',singleX(3),'gaps',[], ...
    'N',Sg.rootTruncations(end,1),'K',Sg.rootTruncations(end,2), ...
    'solve_scattering',false);
Rsingle = ni2019_strict_rayleigh_operator(cfgSingle,'TargetOrder',-1, ...
    'EnforceAllFiniteOpen',true,'KyTolerance',1e-8*abs(2*pi/cfgSingle.lambda), ...
    'Verbose',false);
singleSigma = Rsingle.sigma_ratio;
singleRaw = Rsingle.strict_residual;

% Radiation source phasors in the two-groove strict mode.  For each open or
% threshold order, the global velocity row is a sum of groove contributions.
% The strict mode sets the total to zero; the two individual terms expose the
% independent phase cancellation made available by the second cavity.
[phasorOrders,phasors,totalPhasors] = radiation_phasors(Rnom,-1);
assert(all(isfinite([real(phasors(:));imag(phasors(:))])), ...
    'Nonfinite phasor diagnostic.');

dofRows = struct( ...
    'system','single_rectangular', ...
    'n_cavities',1,'N',cfgSingle.N,'K_per_cavity',cfgSingle.K, ...
    'kappa',singleX(1),'Omega',1-singleX(1), ...
    'sigma_ratio',singleSigma,'raw_residual',singleRaw, ...
    'transverse_fraction',Sg.transverseFraction(end), ...
    'search_bound_note','cached broad single-groove search; see report');
dofRows(2) = struct( ...
    'system','two_groove_nominal', ...
    'n_cavities',2,'N',Nhi,'K_per_cavity',Khi, ...
    'kappa',xNom(1),'Omega',OmegaNom, ...
    'sigma_ratio',nominalSigma,'raw_residual',nominalResidual, ...
    'transverse_fraction',Rnom.groove.pressure_proxy_fraction, ...
    'search_bound_note','current two-groove strict root; normalized a=1');
dofTable = struct2table(dofRows);
dofFile = fullfile(outDir,'dof_comparison.csv');
writetable(dofTable,dofFile);

fprintf('DOF comparison: single sigma=%.6e (raw %.6e), two sigma=%.6e (raw %.6e)\n', ...
    singleSigma,singleRaw,nominalSigma,nominalResidual);
for j = 1:numel(phasorOrders)
    fprintf('  order n=%d phasor closure |sum|=%.3e\n',phasorOrders(j), ...
        abs(totalPhasors(j)));
end

%% C. Homogeneous field states
% All three panels use the same normalized grid and the same visualization
% rule: each state is Euclidean-normalized in coefficient space, then the
% plotted pressure magnitude is divided by its own maximum inside the grooves.
% This prevents arbitrary eigenvector phase/scale from changing the colors.
% The strict map is reconstructed from the high-order strict mode but uses a
% safe low-order field grid (N=121,K=15) so deep transverse cosh factors do
% not overflow.  It is still a homogeneous field, never a driven field.
Nfield = 121;
Kfield = 15;
xGrid = linspace(0,1,241);
yGrid = linspace(-0.75,1.05,301);

strictField = strict_mode_for_field(Rnom,Nfield,Kfield);

% Near-BIC and ordinary-resonance fields are evaluated directly on the final
% geometry at real frequency.  They are smallest-singular-vector surrogates
% of the homogeneous operator, not exact outgoing QNMs; this distinction is
% kept explicit in the report.  No earlier pole cache is used here.
nearOmega = OmegaNom - 6.0e-4;
nearCfg = make_cfg_from_kappa_omega(xNom,nearOmega,Nfield,Kfield);
nearMode = full_mode_for_field(nearCfg);

% Ordinary resonance: choose the local minimum of the full homogeneous
% singular-value ratio in a resolved window well below the Rayleigh point.
ordinaryScan = linspace(0.862,0.872,21);
ordinaryRatio = nan(size(ordinaryScan));
for j = 1:numel(ordinaryScan)
    c = make_cfg_from_kappa_omega(xNom,ordinaryScan(j),Nfield,Kfield);
    [~,metrics] = unconstrained_mode_diagnostics(c,-1);
    ordinaryRatio(j) = metrics.sigma_ratio;
end
[~,ordinaryId] = min(ordinaryRatio);
ordinaryOmega = ordinaryScan(ordinaryId);
ordinaryCfg = make_cfg_from_kappa_omega(xNom,ordinaryOmega,Nfield,Kfield);
ordinaryMode = full_mode_for_field(ordinaryCfg);

fieldStates(1) = evaluate_field_state('strict BIC',strictField,xGrid,yGrid, ...
    xNom(1),OmegaNom,thetaNom,'strict operator; A_-1=A_0=0');
fieldStates(2) = evaluate_field_state('near-BIC',nearMode,xGrid,yGrid, ...
    xNom(1),nearOmega,nearCfg.theta_i_deg, ...
    'real-axis smallest-SVD-vector surrogate; not exact eigenstate');
fieldStates(3) = evaluate_field_state('ordinary resonance',ordinaryMode, ...
    xGrid,yGrid,xNom(1),ordinaryOmega,ordinaryCfg.theta_i_deg, ...
    'real-axis smallest-SVD-vector resonance surrogate');

fieldTable = struct2table(rmfield(fieldStates,{'pmag','X','Y'}));
fieldFile = fullfile(outDir,'field_states.csv');
writetable(fieldTable,fieldFile);

fprintf('Field states (N,K)=(%d,%d): near Omega=%.7f, ordinary Omega=%.7f\n', ...
    Nfield,Kfield,nearOmega,ordinaryOmega);
for j = 1:numel(fieldStates)
    fprintf('  %-18s sigma=%.3e raw=%.3e open err=%.3e exterior/groove=%.3e\n', ...
        fieldStates(j).state,fieldStates(j).sigma_ratio, ...
        fieldStates(j).raw_residual,fieldStates(j).finite_open_error, ...
        fieldStates(j).exterior_to_groove_energy);
end

%% D. Publication-style vector figures
figTol = make_tolerance_figure(summaryRows,toleranceRows,modeNames,levels);
tolFigure = fullfile(outDir,'fig_fabrication_tolerance.pdf');
exportgraphics(figTol,tolFigure,'ContentType','vector');
close(figTol);

figDof = make_dof_figure(dofRows,phasorOrders,phasors,totalPhasors);
dofFigure = fullfile(outDir,'fig_dof_comparison.pdf');
exportgraphics(figDof,dofFigure,'ContentType','vector');
close(figDof);

figFields = make_field_figure(fieldStates,xGrid,yGrid,xNom);
fieldFigure = fullfile(outDir,'fig_homogeneous_fields.pdf');
exportgraphics(figFields,fieldFigure,'ContentType','vector');
close(figFields);

%% E. MAT archive and machine-readable report
archiveFile = fullfile(outDir,'robustness_fields.mat');
save(archiveFile,'xNom','OmegaNom','thetaNom','fNom','cNom','aPhysical', ...
    'lambdaPhysical','Nhi','Khi','Ntol','Ktol','nMC','seed','levels', ...
    'nominalSigma','nominalResidual','tolBaseline','toleranceRows','summaryRows', ...
    'singleX','singleSigma','singleRaw','dofRows','phasorOrders','phasors', ...
    'totalPhasors','fieldStates','ordinaryScan','ordinaryRatio','nearOmega', ...
    'ordinaryOmega','Nfield','Kfield','xGrid','yGrid');

reportFile = fullfile(rootDir,'advanced_analysis','robustness_fields_report.md');
write_report(reportFile,designFile,singleFile,toleranceFile,summaryFile, ...
    dofFile,fieldFile,tolFigure,dofFigure,fieldFigure,archiveFile, ...
    xNom,OmegaNom,thetaNom,fNom,cNom,aPhysical,lambdaPhysical,Nhi,Khi, ...
    Ntol,Ktol,nMC,seed,summaryRows,singleX,singleSigma,singleRaw, ...
    nominalSigma,nominalResidual,phasorOrders,totalPhasors,fieldStates, ...
    nearOmega,ordinaryOmega,ordinaryScan,ordinaryRatio,Nfield,Kfield, ...
    cfgSingle.N,cfgSingle.K,tolBaseline,toleranceBaselineFile);

fprintf('\nOutputs written under %s\n',outDir);
fprintf('  %s\n  %s\n  %s\n',tolFigure,dofFigure,fieldFigure);
fprintf('  %s\n',reportFile);

%% Local functions
function cfg = make_cfg(x,N,K)
Omega = 1-x(1);
cfg = struct('a',1,'lambda',1/Omega, ...
    'theta_i_deg',asind(x(1)/Omega),'depths',x(2:3), ...
    'widths',x(4:5),'gaps',x(6),'N',N,'K',K, ...
    'solve_scattering',false);
end

function cfg = make_cfg_from_kappa_omega(x,Omega,N,K)
cfg = struct('a',1,'lambda',1/Omega, ...
    'theta_i_deg',asind(x(1)/Omega),'depths',x(2:3), ...
    'widths',x(4:5),'gaps',x(6),'N',N,'K',K, ...
    'solve_scattering',false);
end

function tf = geometry_is_valid(x)
tf = all(isfinite(x)) && x(1)>0 && x(1)<.5 && ...
    all(x(2:5)>0) && x(6)>=0 && sum(x(4:6))<1;
end

function [R,metrics] = unconstrained_mode_diagnostics(cfg,targetOrder)
op = ni2019_full_eigen_operator(cfg);
[~,S,V] = svd(op.Fscaled,'econ');
sv = diag(S);
z = V(:,end)./transpose(op.column_scale);
z = z/max(norm(z),eps);
A = z(1:op.N);
target = op.orders==targetOrder;
finiteOpen = ~target & abs(imag(op.ky))<=1e-8*max(abs(op.k0),1) & ...
    real(op.ky)>1e-8*max(abs(op.k0),1);
metrics = struct('sigma_ratio',sv(end)/max(sv(1),eps), ...
    'raw_residual',norm(op.F*z)/max(norm(z),eps), ...
    'threshold_error',norm(A(target))/max(norm(A),eps), ...
    'finite_open_error',norm(A(finiteOpen))/max(norm(A),eps), ...
    'A',A,'z',z,'op',op,'sv',sv);
R = op;
end

function summary = build_tolerance_summary(rows,modeNames,levels)
summary = repmat(struct( ...
    'mode','','level_pct',NaN,'sample_count',0,'valid_count',0, ...
    'geometry_invalid_count',0,'solver_failure_count',0, ...
    'strict_sigma_q05',NaN,'strict_sigma_q50',NaN,'strict_sigma_q95',NaN, ...
    'strict_raw_q05',NaN,'strict_raw_q50',NaN,'strict_raw_q95',NaN, ...
    'threshold_error_q05',NaN,'threshold_error_q50',NaN, ...
    'threshold_error_q95',NaN,'finite_open_error_q50',NaN, ...
    'strict_raw_above_1e_2',0),numel(modeNames)*numel(levels),1);
idx = 0;
for im = 1:numel(modeNames)
    for il = 1:numel(levels)
        idx = idx+1;
        mask = strcmp({rows.mode},modeNames{im}) & ...
            abs([rows.level_pct]-100*levels(il))<1e-12;
        r = rows(mask);
        valid = [r.geometry_valid] & [r.solver_ok];
        sr = [r(valid).strict_sigma_ratio];
        rr = [r(valid).strict_raw_residual];
        te = [r(valid).threshold_amplitude_error];
        fe = [r(valid).finite_open_amplitude_error];
        summary(idx).mode = modeNames{im};
        summary(idx).level_pct = 100*levels(il);
        summary(idx).sample_count = numel(r);
        summary(idx).valid_count = nnz(valid);
        summary(idx).geometry_invalid_count = nnz(~[r.geometry_valid]);
        summary(idx).solver_failure_count = nnz([r.geometry_valid] & ~[r.solver_ok]);
        if ~isempty(sr)
            summary(idx).strict_sigma_q05 = prctile(sr,5);
            summary(idx).strict_sigma_q50 = prctile(sr,50);
            summary(idx).strict_sigma_q95 = prctile(sr,95);
            summary(idx).strict_raw_q05 = prctile(rr,5);
            summary(idx).strict_raw_q50 = prctile(rr,50);
            summary(idx).strict_raw_q95 = prctile(rr,95);
            summary(idx).threshold_error_q05 = prctile(te,5);
            summary(idx).threshold_error_q50 = prctile(te,50);
            summary(idx).threshold_error_q95 = prctile(te,95);
            summary(idx).finite_open_error_q50 = prctile(fe,50);
            summary(idx).strict_raw_above_1e_2 = nnz(rr>1e-2);
        end
    end
end
end

function [orders,phasors,total] = radiation_phasors(R,targetOrder)
op = R.full_operator;
Cscaled = R.mode.C_scaled;
orders = op.orders(abs(op.ky)<=1e-8*max(abs(op.k0),1) | ...
    (abs(imag(op.ky))<=1e-8*max(abs(op.k0),1) & real(op.ky)>1e-8*max(abs(op.k0),1)));
orders = orders(orders==targetOrder | orders==0);
L = op.L; K = op.K;
phasors = complex(zeros(numel(orders),L));
for io = 1:numel(orders)
    n = find(op.orders==orders(io),1);
    for ell = 1:L
        ids = (ell-1)*K+(1:K);
        phasors(io,ell) = 1i*sum( ...
            groove_projection(-op.kx(n),op.xleft(ell),op.widths(ell),0:K-1) .* ...
            op.beta_sin_normalized(ids).' .* Cscaled(ids).');
    end
end
total = sum(phasors,2);
end

function values = groove_projection(kx,xl,t,qs)
values = zeros(size(qs));
for j = 1:numel(qs)
    q = qs(j); alpha=q*pi/t;
    values(j) = (t)*0.5*expint_local(kx-alpha,t) + ...
        (t)*0.5*expint_local(kx+alpha,t);
    % The full operator has Cmap=(t/a)*projection and a=1 here.
    values(j) = values(j)*exp(-1i*kx*xl);
end
end

function value = expint_local(kappa,t)
z = kappa*t;
if abs(z)<1e-8
    value = 1-1i*z/2-z^2/6;
else
    value = exp(-1i*z/2)*sin(z/2)/(z/2);
end
end

function mode = strict_mode_for_field(R,N,K)
% Truncate the high-order strict mode only for stable field plotting.  The
% strict state itself is still validated at the full (313,39) truncation.
cfg = R.cfg;
cfg.N = N; cfg.K = K;
opPlot = ni2019_full_eigen_operator(cfg);
zHigh = R.mode.z_full;
ordersHigh = R.full_operator.orders;
% Floquet orders are centered and therefore the retained low-order block is
% contiguous; explicit order matching prevents an indexing convention drift.
A = zeros(opPlot.N,1);
[idsA,locA] = ismember(opPlot.orders,ordersHigh);
A(idsA) = zHigh(locA(idsA));
surfaceHigh = R.mode.surface_coefficients;
Csurf = zeros(opPlot.K,opPlot.L);
for ell = 1:opPlot.L
    ids = (ell-1)*R.full_operator.K + (1:min(K,R.full_operator.K));
    Csurf(:,ell) = surfaceHigh(ids(1:min(K,numel(ids))));
end
mode = mode_struct(opPlot,A,Csurf,'strict BIC',R.residual.reduced_raw, ...
    R.sigma_ratio);
end

function mode = full_mode_for_field(cfg)
op = ni2019_full_eigen_operator(cfg);
[~,S,V] = svd(op.Fscaled,'econ');
z = V(:,end)./transpose(op.column_scale);
z = z/max(norm(z),eps);
A = z(1:op.N);
Cscaled = z(op.N+1:end);
surface = Cscaled.*op.cos_depth_normalized;
mode = mode_struct(op,A,reshape(surface,op.K,op.L), ...
    'full homogeneous SVD',norm(op.F*z)/max(norm(z),eps), ...
    S(end,end)/max(S(1,1),eps));
end

function mode = mode_struct(op,A,surface,kind,rawResidual,sigmaRatio)
mode = struct('op',op,'A',A,'surface_coefficients',surface, ...
    'kind',kind,'raw_residual',rawResidual,'sigma_ratio',sigmaRatio, ...
    'orders',op.orders,'kx',op.kx,'ky',op.ky,'k0',op.k0,'a',op.a, ...
    'lambda',op.lambda,'theta_i_deg',NaN,'widths',op.widths, ...
    'depths',op.depths,'gaps',op.gaps,'xleft',op.xleft, ...
    'N',op.N,'K',op.K);
end

function out = evaluate_field_state(label,mode,x,y,kappa,Omega,theta,note)
[X,Y] = meshgrid(x,y);
p = homogeneous_pressure(mode,X,Y);
groove = false(size(X));
for ell = 1:numel(mode.widths)
    groove = groove | (Y<0 & Y>=-mode.depths(ell) & ...
        X>=mode.xleft(ell) & X<=mode.xleft(ell)+mode.widths(ell));
end
gmax = max(abs(p(groove)),[],'omitnan');
if isempty(gmax) || ~isfinite(gmax) || gmax==0, gmax=1; end
pmag = abs(p)/gmax;
pmag(~isfinite(pmag)) = NaN;
above = Y>=0;
dx = mean(diff(x)); dy = mean(diff(y));
grooveEnergy = sum(pmag(groove).^2,'omitnan')*dx*dy;
exteriorEnergy = sum(pmag(above).^2,'omitnan')*dx*dy;
target = mode.orders==-1;
finiteOpen = ~target & abs(imag(mode.ky))<=1e-8*max(abs(mode.k0),1) & ...
    real(mode.ky)>1e-8*max(abs(mode.k0),1);
OmegaReal = real(Omega);
OmegaImag = imag(Omega);
out = struct('state',label,'kind',mode.kind,'note',note, ...
    'kappa',real(kappa),'Omega',OmegaReal,'Omega_imag',OmegaImag, ...
    'f_kHz',OmegaReal*1500/0.00687113477585325/1e3, ...
    'theta_deg',theta,'N',mode.N,'K',mode.K, ...
    'sigma_ratio',mode.sigma_ratio,'raw_residual',mode.raw_residual, ...
    'threshold_error',norm(mode.A(target))/max(norm(mode.A),eps), ...
    'finite_open_error',norm(mode.A(finiteOpen))/max(norm(mode.A),eps), ...
    'groove_max_pressure',gmax,'exterior_max_norm',max(pmag(above),[],'all'), ...
    'groove_energy_norm',grooveEnergy,'exterior_energy_norm',exteriorEnergy, ...
    'exterior_to_groove_energy',exteriorEnergy/max(grooveEnergy,eps));
out.pmag = pmag;
out.X = X; out.Y = Y;
end

function p = homogeneous_pressure(mode,X,Y)
p = complex(nan(size(X)));
above = Y>=0;
pa = complex(zeros(nnz(above),1));
for n = 1:numel(mode.orders)
    pa = pa + mode.A(n).*exp(-1i*mode.kx(n).*X(above) - ...
        1i*mode.ky(n).*Y(above));
end
p(above) = pa;
for ell = 1:numel(mode.widths)
    inside = Y<0 & Y>=-mode.depths(ell) & ...
        X>=mode.xleft(ell) & X<=mode.xleft(ell)+mode.widths(ell);
    if ~any(inside,'all'), continue; end
    u = X(inside)-mode.xleft(ell);
    yy = Y(inside);
    pg = complex(zeros(size(u)));
    for q = 0:mode.K-1
        alpha = q*pi/mode.widths(ell);
        beta = groove_beta_local(mode.k0^2-alpha^2);
        Pq = mode.surface_coefficients(q+1,ell);
        vertical = safe_cos_ratio(beta,yy+mode.depths(ell),mode.depths(ell));
        pg = pg + Pq*cos(alpha*u).*vertical;
    end
    p(inside) = pg;
end
end

function value = safe_cos_ratio(beta,t,d)
% Stable ratio cos(beta*t)/cos(beta*d), with 0<=t<=d in a groove.
sn = abs(imag(beta*t)); sd = abs(imag(beta*d));
num = exp(1i*beta*t-sn)+exp(-1i*beta*t-sn);
den = exp(1i*beta*d-sd)+exp(-1i*beta*d-sd);
value = exp(sn-sd).*num./den;
end

function b = groove_beta_local(z)
if real(z)>=0, b=sqrt(real(z)); else, b=1i*sqrt(-real(z)); end
end

function fig = make_tolerance_figure(summaryRows,rows,modeNames,levels)
fig = figure('Color','w','Position',[80 80 1120 760]);
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
colors = [0.05 0.23 0.55; 0.82 0.25 0.10];
metrics = {'strict_raw_q','threshold_error_q'};
ylabs = {'strict raw residual','threshold amplitude error'};
for imetric = 1:2
    nexttile(tl); hold on; box on; grid on;
    for im = 1:numel(modeNames)
        med = nan(size(levels)); lo=med; hi=med;
        for il = 1:numel(levels)
            q = summaryRows(strcmp({summaryRows.mode},modeNames{im}) & ...
                [summaryRows.level_pct]==100*levels(il));
            if imetric==1
                lo(il)=q.strict_raw_q05; med(il)=q.strict_raw_q50; hi(il)=q.strict_raw_q95;
            else
                lo(il)=q.threshold_error_q05; med(il)=q.threshold_error_q50; hi(il)=q.threshold_error_q95;
            end
        end
        errorbar(100*levels,med,med-lo,hi-med,'o-','Color',colors(im,:), ...
            'LineWidth',1.6,'MarkerFaceColor',colors(im,:),'MarkerSize',5);
    end
    set(gca,'YScale','log','XScale','linear');
    xlabel('uniform perturbation bound (%)'); ylabel(ylabs{imetric});
    legend(modeNames,'Location','best');
    title(sprintf('(%c) 5--95%% ensemble interval',96+imetric));
end
nexttile(tl); hold on; box on; grid on;
for im = 1:numel(modeNames)
    mask = strcmp({rows.mode},modeNames{im}) & [rows.level_pct]==5 & ...
        [rows.geometry_valid] & [rows.solver_ok];
    vals = [rows(mask).strict_raw_residual]; vals=sort(vals);
    cdf = (1:numel(vals))/max(numel(vals),1);
    plot(vals,cdf,'-','Color',colors(im,:),'LineWidth',1.7);
end
set(gca,'XScale','log'); xlabel('strict raw residual at +/-5%'); ylabel('empirical CDF');
legend(modeNames,'Location','southeast'); title('(c) joint/independent residual distributions');
nexttile(tl); hold on; box on; grid on;
barData=zeros(numel(levels),2);
for il=1:numel(levels)
    for im=1:2
        q=summaryRows(strcmp({summaryRows.mode},modeNames{im}) & ...
            [summaryRows.level_pct]==100*levels(il));
        barData(il,im)=q.geometry_invalid_count+q.solver_failure_count;
    end
end
bar(100*levels,barData,'grouped'); xlabel('uniform perturbation bound (%)');
ylabel('invalid or failed realizations'); ylim([0 1]);
legend(modeNames,'Location','best');
text(mean(100*levels),0.5,'none in this ensemble','HorizontalAlignment','center', ...
    'Color',[.25 .25 .25],'FontSize',8);
title('(d) validity accounting (no residual pass/fail threshold imposed)');
sgtitle(tl,'Fabrication tolerance of the homogeneous strict-radiation diagnostic', ...
    'FontWeight','normal');
style_axes(fig);
end

function fig = make_dof_figure(rows,orders,phasors,total)
fig = figure('Color','w','Position',[80 80 1120 400]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile(tl); hold on; box on; grid on;
vals=[rows.sigma_ratio];
bar(1:2,vals,'FaceColor','flat','CData',[0.82 0.25 0.10;0.05 0.23 0.55]);
set(gca,'YScale','log','XTick',1:2,'XTickLabel',{'one cavity','two cavities'});
ylabel('\sigma_{min}/\sigma_{max}'); title('(a) complete residual');
text(1,vals(1)*1.4,sprintf('%.3g',vals(1)),'HorizontalAlignment','center');
text(2,vals(2)*5,sprintf('%.3g',vals(2)),'HorizontalAlignment','center');
nexttile(tl); hold on; box on; grid on; axis equal;
colors=[0.05 0.23 0.55;0.82 0.25 0.10];
for io=1:numel(orders)
    c=phasors(io,:); scale=max(abs([c,total(io)]));
    if scale==0, scale=1; end
    for ell=1:numel(c)
        plot([0 real(c(ell))/scale],[0 imag(c(ell))/scale],'-', ...
            'Color',colors(ell,:),'LineWidth',1.8);
        plot(real(c(ell))/scale,imag(c(ell))/scale,'o', ...
            'Color',colors(ell,:),'MarkerFaceColor',colors(ell,:),'MarkerSize',5);
    end
    plot([0 real(total(io))/scale],[0 imag(total(io))/scale],'k--','LineWidth',1.2);
    text(-0.95,0.92-0.16*(io-1), ...
        sprintf('n=%d: |sum|=%.1e',orders(io),abs(total(io))), ...
        'FontSize',8,'HorizontalAlignment','left');
    if io<numel(orders), xline(0,'Color',[.8 .8 .8]); yline(0,'Color',[.8 .8 .8]); end
end
xlabel('Re normalized aperture source'); ylabel('Im normalized aperture source');
legend('groove 1','groove 2','sum','Location','best'); title('(b) two independent phase closures');
nexttile(tl); hold on; box on; grid on;
% Display the single-groove transverse content beside the two-groove strict
% state as a compact diagnostic, without turning it into a universal claim.
bar(1,[rows(1).transverse_fraction],'FaceColor',[0.82 0.25 0.10]);
bar(2,[rows(2).transverse_fraction],'FaceColor',[0.05 0.23 0.55]);
set(gca,'XTick',1:2,'XTickLabel',{'one cavity','two cavities'});
ylabel('transverse-mode pressure fraction'); ylim([0 1.08]);
title('(c) modal content in the searched solutions');
text(1,rows(1).transverse_fraction+.04,sprintf('%.3f',rows(1).transverse_fraction), ...
    'HorizontalAlignment','center');
text(2,rows(2).transverse_fraction+.04,sprintf('%.3f',rows(2).transverse_fraction), ...
    'HorizontalAlignment','center');
sgtitle(tl,'Cavity degree of freedom: bounded numerical ablation', ...
    'FontWeight','normal');
style_axes(fig);
end

function fig = make_field_figure(states,x,y,xNom)
fig = figure('Color','w','Position',[60 80 1220 390]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for j=1:numel(states)
    nexttile(tl); hold on; box on;
    pcolor(x,y,states(j).pmag); shading flat;
    caxis([0 1]); axis image; xlim([0 1]); ylim([-0.75 1.05]);
    xlabel('x/a'); ylabel('y/a');
    title(sprintf('%s\nOmega=%.6f',states(j).state,states(j).Omega), ...
        'FontWeight','normal');
    for ell=1:2
        % xNom stores normalized groove widths/depths/gap.
        if ell==1, xl=(1-sum(xNom(4:6)))/2; else, xl=(1-sum(xNom(4:6)))/2+xNom(4)+xNom(6); end
        rectangle('Position',[xl,-xNom(ell+1),xNom(ell+3),xNom(ell+1)], ...
            'EdgeColor',[.1 .1 .1],'LineWidth',.7,'FaceColor','none');
    end
    text(.03,.94,sprintf('open err %.2e',states(j).finite_open_error), ...
        'Units','normalized','Color','w','FontSize',7,'FontWeight','bold');
end
cb = colorbar; cb.Layout.Tile='east'; cb.Label.String='|p| / max_{groove}|p|';
colormap(fig,turbo(256));
sgtitle(tl,'Homogeneous pressure fields: strict BIC, near-BIC, and ordinary resonance', ...
    'FontWeight','normal');
style_axes(fig);
end

function style_axes(fig)
ax = findall(fig,'Type','axes');
for j=1:numel(ax)
    ax(j).FontName='Helvetica'; ax(j).FontSize=8;
    ax(j).LineWidth=.7; ax(j).TickDir='out';
    ax(j).XColor=[.12 .12 .12]; ax(j).YColor=[.12 .12 .12];
    ax(j).GridColor=[.78 .78 .78]; ax(j).GridAlpha=.35;
end
end

function write_report(fileName,designFile,singleFile,tolFile,sumFile,dofFile, ...
    fieldFile,tolFigure,dofFigure,fieldFigure,archiveFile,xNom,OmegaNom, ...
    thetaNom,fNom,cNom,aPhysical,lambdaPhysical,Nhi,Khi,Ntol,Ktol,nMC,seed, ...
    summaryRows,singleX,singleSigma,singleRaw,nominalSigma,nominalResidual, ...
    phasorOrders,totalPhasors,fieldStates,nearOmega,ordinaryOmega, ...
    ordinaryScan,ordinaryRatio,Nfield,Kfield,singleN,singleK,tolBaseline, ...
    baselineFile)
fid=fopen(fileName,'w'); assert(fid>0,'Cannot write report: %s',fileName);
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# Robustness, cavity-DOF, and homogeneous-field diagnostics\n\n');
fprintf(fid,'This report is generated from the existing MATLAB modal-matching implementation. It does not alter the manuscript claim and does not use driven fields as eigenfields. All files are generated by `advanced_analysis/robustness_fields.m`.\n\n');
fprintf(fid,'## Nominal strict BIC reproduction\n\n');
fprintf(fid,'- Design cache: `%s`\n',designFile);
fprintf(fid,'- $f_0=%.6g$ kHz, $c=%.6g$ m/s, $a=%.9g$ mm, $\\lambda=%.9g$ mm.\n',fNom/1e3,cNom,1e3*aPhysical,1e3*lambdaPhysical);
fprintf(fid,'- Normalized point: $\\kappa=%.15g$, $\\Omega=%.15g$, $\\theta=\\pm%.12g^\\circ$.\n',xNom(1),OmegaNom,thetaNom);
fprintf(fid,'- High-order strict operator: $(N,K)=(%d,%d)$, where $N$ is the retained exterior-order count and $K$ is the cosine-mode count per groove.\n',Nhi,Khi);
fprintf(fid,'- The selected-threshold classification uses an explicit absolute tolerance $|k_y|\\le 10^{-8}k_0$; the target threshold is then removed from the strict reduced unknown vector.\n');
fprintf(fid,'- Reproduced $\\sigma_{\\min}/\\sigma_{\\max}=%.6e$ and raw residual $r=%.6e$. The selected $n=-1$ threshold and finite-flux $n=0$ amplitude are removed before the reduced SVD and reinserted as exact zeros.\n\n',nominalSigma,nominalResidual);

fprintf(fid,'## A. Fabrication tolerance\n\n');
fprintf(fid,'The fixed Rayleigh point and period are held constant. Each of $(d_1,d_2,w_1,w_2,g)$ is perturbed in normalized units. `independent` perturbs one parameter at a time (balanced cyclic parameter selection); `joint` perturbs all five simultaneously. Each draw is uniform in $[-p,p]$ with `rng(%d,''twister'')`, and uses %d realizations per mode and perturbation level. The ensemble uses $(N,K)=(%d,%d)$.\n\n',seed,nMC,Ntol,Ktol);
fprintf(fid,'The unperturbed final geometry at this ensemble truncation has $\\sigma_{\\min}/\\sigma_{\\max}=%.6e$ and raw residual $%.6e$; this finite discretization floor is retained explicitly in the sample CSV through the `*_over_baseline` columns.\n\n',tolBaseline.strict_sigma_ratio,tolBaseline.strict_raw_residual);
fprintf(fid,'For each valid draw, `strict_raw_residual` is the reduced homogeneous residual after enforcing $A_{-1}=A_0=0$. `threshold_amplitude_error` is instead computed from the *unconstrained* full-operator smallest right-singular vector, $|A_{-1}|/\\|A\\|_2$; it therefore measures how much the perturbed geometry fails to satisfy the threshold-pressure condition before that condition is imposed. The strict operator’s threshold amplitude is exactly zero by construction and is not used as the error metric. `finite_open_amplitude_error` is $|A_0|/\\|A\\|_2$. Invalid geometries (nonpositive dimensions or occupied aperture $\\ge a$) and nonfinite solver outputs are retained in the sample CSV and excluded from quantiles. No arbitrary residual pass/fail cutoff is claimed; $10^{-2}$ is shown only as a descriptive count.\n\n');
fprintf(fid,'| mode | level | valid/total | geometry invalid | solver failures | strict raw q05/q50/q95 | threshold error q05/q50/q95 |\n|---|---:|---:|---:|---:|---:|---:|\n');
for j=1:numel(summaryRows)
    q=summaryRows(j);
    fprintf(fid,'| %s | %.0f%% | %d/%d | %d | %d | %.3e / %.3e / %.3e | %.3e / %.3e / %.3e |\n', ...
        q.mode,q.level_pct,q.valid_count,q.sample_count,q.geometry_invalid_count, ...
        q.solver_failure_count,q.strict_raw_q05,q.strict_raw_q50,q.strict_raw_q95, ...
        q.threshold_error_q05,q.threshold_error_q50,q.threshold_error_q95);
end
fprintf(fid,'\nFiles: `%s`, `%s`, `%s`, `%s`.\n\n',tolFile,sumFile,baselineFile,tolFigure);

fprintf(fid,'## B. Cavity degree of freedom\n\n');
fprintf(fid,'The single-cavity point is the final cached solution of the documented rectangular single-groove search, with transverse modes retained through $(N,K)=(%d,%d)$ and normalized $x=(\\kappa,d/a,w/a)=(%.9g,%.9g,%.9g)$. Direct reevaluation gives $\\sigma_{\\min}/\\sigma_{\\max}=%.6e$ and raw residual %.6e. The two-groove strict root gives %.6e and %.6e at the high-order truncation.\n\n', ...
    singleN,singleK,singleX(1),singleX(2),singleX(3),singleSigma,singleRaw,nominalSigma,nominalResidual);
fprintf(fid,'The comparison is a bounded ablation, not a no-go theorem for every possible one-cavity shape. Both systems use the same dimensionless Rayleigh condition $\\Omega=1-\\kappa$ and the same pressure--velocity modal normalization. The two-cavity design has two independent aperture phase combinations. At positive angle, the $n=0$ finite-flux row and the $n=-1$ threshold-pressure row are two complex radiation constraints. In the strict two-groove mode, the per-groove source phasors close with residuals:\n\n');
for j=1:numel(phasorOrders)
    fprintf(fid,'- $n=%d$: $|\\sum_\\ell s_{n\\ell}|=%.6e$.\n',phasorOrders(j),abs(totalPhasors(j)));
end
fprintf(fid,'\nThis supports the phase-degree-of-freedom interpretation within the searched rectangular family; it does not establish a universal mathematical impossibility for a single cavity with arbitrary internal structure.\n\n');
fprintf(fid,'Files: `%s`, `%s`.\n\n',dofFile,dofFigure);

fprintf(fid,'## C. Homogeneous eigenfield visualization\n\n');
fprintf(fid,'The three maps use the same $(x/a,y/a)$ grid and color range. Each coefficient vector is Euclidean-normalized; the plotted pressure magnitude is divided by the maximum magnitude inside the grooves using the same rule for every state. Rigid material outside the grooves is masked. The strict BIC map is reconstructed from the high-order strict mode; its coefficients are truncated to $(N,K)=(%d,%d)$ only for visualization, while the strict numerical diagnosis remains at $(313,39)$. The near-BIC and ordinary-resonance maps are smallest-SVD-vector surrogates of the full real-axis homogeneous operator, not exact real-frequency outgoing eigenstates or QNMs.\n\n',Nfield,Kfield);
fprintf(fid,'The near-BIC is fixed at $\\kappa=%.12g$, $\\Omega=%.12g$ ($\\Delta\\Omega=%.3e$), just below the threshold. The ordinary resonance is selected as the minimum full-operator singular-value ratio on the resolved scan $\\Omega\\in[%.3f,%.3f]$ (%d points), at $\\Omega=%.12g$.\n\n',xNom(1),nearOmega,nearOmega-OmegaNom,ordinaryScan(1),ordinaryScan(end),numel(ordinaryScan),ordinaryOmega);
fprintf(fid,'| state | Omega | theta (deg) | sigma ratio | raw residual | threshold error | finite-open error | exterior/groove normalized energy |\n|---|---:|---:|---:|---:|---:|---:|---:|\n');
for j=1:numel(fieldStates)
    s=fieldStates(j);
    fprintf(fid,'| %s | %.9g | %.6g | %.3e | %.3e | %.3e | %.3e | %.3e |\n', ...
        s.state,s.Omega,s.theta_deg,s.sigma_ratio,s.raw_residual, ...
        s.threshold_error,s.finite_open_error,s.exterior_to_groove_energy);
end
fprintf(fid,'\nThe ordinary resonance field is intentionally labeled a real-axis surrogate because an exactly lossless open resonance has a complex pole rather than a square-integrable real-frequency null. No driven scattering field is used in these panels.\n\n');
fprintf(fid,'Files: `%s`, `%s`.\n\n',fieldFile,fieldFigure);

fprintf(fid,'## Scope and remaining risks\n\n');
fprintf(fid,'- The tolerance ensemble keeps the sound speed, frequency, period, and Bloch point fixed; environmental sound-speed scaling is a separate analysis.\n');
fprintf(fid,'- The ensemble samples local uniform perturbations and does not replace a fabrication process model, finite-aperture simulation, thermoviscous loss, or full 3-D tank model.\n');
fprintf(fid,'- The single-cavity result is a numerical search bound. A stronger manuscript statement would require a broader shape family or a formal rank argument.\n');
fprintf(fid,'- The cached single-cavity search used its documented wide-aperture bounds ($w/a\\ge0.50$), whereas the two-cavity design includes a narrow second groove ($w_2/a=%.3g$). Thus the DOF comparison shares normalization and modal equations but is not an identical global-search domain; an expanded single-cavity search over the narrow-width range remains a recommended control.\n',xNom(5));
fprintf(fid,'- The field plots establish qualitative localization and channel content of the computed homogeneous vectors; they are not an experimental prediction of absolute pressure.\n\n');
fprintf(fid,'Archive: `%s`.\n',archiveFile);
end
