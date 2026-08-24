%% Scaling and environmental robustness of the strict Rayleigh BIC
% This script is deliberately self-contained at the analysis level.  It
% calls the existing pole-free modal-matching operator and does not replace
% that solver by a fitted or analytical surrogate.
%
% Run from the repository root with
%   matlab -batch "run('advanced_analysis/scaling_environment.m')"
%
% The two scaling experiments distinguish exact similarity retuning from a
% fixed physical device.  The former keeps all normalized quantities and
% therefore the Rayleigh-BIC mechanism invariant.  The latter is a control:
% changing c at fixed geometry, frequency, and angle changes a/lambda and
% moves the device away from the calibrated normalized BIC point.

clear; close all; clc;

repoRoot = pwd;
if ~isfile(fullfile(repoRoot,'Ni2019_MATLAB','ni2019_strict_rayleigh_operator.m'))
    thisFile = mfilename('fullpath');
    if isempty(thisFile)
        error('Run this script from the repository root.');
    end
    repoRoot = fileparts(fileparts(thisFile));
end
solverDir = fullfile(repoRoot,'Ni2019_MATLAB');
addpath(solverDir);
outputDir = fullfile(repoRoot,'advanced_analysis','scaling_environment');
if ~exist(outputDir,'dir'), mkdir(outputDir); end

% Nominal physical design produced by the existing 200 kHz / 1 mm search.
parameterFile = fullfile(solverDir,'results', ...
    'StrictRayleighBIC_200kHz_min1mm.mat');
if ~isfile(parameterFile)
    error('Nominal parameter file is missing: %s',parameterFile);
end
D0 = load(parameterFile);

c0 = D0.cWater;                    % nominal water sound speed (m/s)
f0 = D0.fTarget;                   % nominal operating frequency (Hz)
N = 313;                           % reported strict-operator truncation
K = 39;
a0 = D0.aPhysical;                 % m
lambda0 = D0.lambdaPhysical;       % m
Omega0 = D0.OmegaFinal;
kappa0 = D0.xFinal(1);
theta0 = D0.thetaFinal;
depths0 = a0*D0.xFinal(2:3);
widths0 = a0*D0.xFinal(4:5);
gap0 = a0*D0.xFinal(6);

% Internal consistency checks make accidental use of an older optimized
% geometry fail loudly instead of silently producing a second design.
assert(abs(a0/lambda0 - Omega0) < 2e-12, ...
    'The nominal a/lambda is inconsistent with Omega.');
assert(abs(c0/lambda0 - f0) < 2e-6, ...
    'The nominal wavelength is inconsistent with f0 and c0.');
assert(abs(kappa0 - Omega0*sind(theta0)) < 2e-12, ...
    'The nominal kappa is inconsistent with the incidence angle.');
assert(all(depths0 > 0) && all(widths0 > 0) && gap0 > 0, ...
    'Nominal geometry must be positive.');

% Normalized geometry is the invariant object in the similarity analysis.
normalizedGeometry = struct('d_over_a',depths0/a0, ...
    'w_over_a',widths0/a0,'g_over_a',gap0/a0);

fprintf('Scaling/environment analysis of the strict Rayleigh BIC\n');
fprintf('  nominal: a = %.9f mm, lambda = %.9f mm, f = %.6f kHz\n', ...
    1e3*a0,1e3*lambda0,1e-3*f0);
fprintf('  kappa = %.12f, Omega = %.12f, theta = %.9f deg\n', ...
    kappa0,Omega0,theta0);
fprintf('  normalized geometry d/a=[%.12f %.12f], w/a=[%.12f %.12f], g/a=%.12f\n', ...
    normalizedGeometry.d_over_a,normalizedGeometry.w_over_a,normalizedGeometry.g_over_a);

% Reproduce the nominal strict homogeneous cancellation before any sweep.
cfg0 = make_cfg(a0,lambda0,theta0,depths0,widths0,gap0,N,K);
nominal = ni2019_strict_rayleigh_operator(cfg0,'TargetOrder',-1, ...
    'EnforceAllFiniteOpen',true,'NormalizeMode',true,'Verbose',false, ...
    'KyTolerance',1e-7*max(abs(2*pi/lambda0),1));
nominalSummary = struct('sigma_ratio',nominal.sigma_ratio, ...
    'strict_residual',nominal.strict_residual, ...
    'full_residual',nominal.residual.full_raw, ...
    'finite_open_power',nominal.radiation.finite_open_power, ...
    'threshold_amplitude_norm',nominal.radiation.grazing_amplitude_norm, ...
    'target_is_rayleigh',nominal.target_is_rayleigh, ...
    'removed_orders',nominal.removed_orders(:).');
if ~nominal.target_is_rayleigh
    error('Nominal target order is not classified as Rayleigh.');
end
fprintf('  nominal strict operator: sigma ratio %.6e, raw residual %.6e\n', ...
    nominalSummary.sigma_ratio,nominalSummary.strict_residual);
fprintf('  nominal removed orders: [%s]\n',strtrim(sprintf('%d ',nominalSummary.removed_orders)));

%% A. Exact dimensional scaling
% For fixed normalized geometry and fixed (kappa,Omega),
%       Omega0 = a f_BIC / c,
%       f_BIC(c,a) = Omega0 c/a.
% We evaluate the full modal operator at each dimensional point, while all
% normalized parameters are generated from the same nominal design.

cSweep = (1400:10:1600).';          % requested environmental range (m/s)
scaleFactors = linspace(0.50,2.00,16).';

fBIC_c = Omega0*cSweep/a0;
retunedResidual = nan(size(cSweep));
retunedSigma = nan(size(cSweep));
retunedFinitePower = nan(size(cSweep));
retunedOmega = nan(size(cSweep));
retunedKappa = nan(size(cSweep));
for j = 1:numel(cSweep)
    % c changes, but a and all physical dimensions stay at the nominal
    % scale; f is retuned so that lambda=a/Omega0 exactly.
    lambda = cSweep(j)/fBIC_c(j);
    R = ni2019_strict_rayleigh_operator( ...
        make_cfg(a0,lambda,theta0,depths0,widths0,gap0,N,K), ...
        'TargetOrder',-1,'EnforceAllFiniteOpen',true, ...
        'NormalizeMode',true,'Verbose',false, ...
        'KyTolerance',1e-7*max(abs(2*pi/lambda),1));
    retunedResidual(j) = R.strict_residual;
    retunedSigma(j) = R.sigma_ratio;
    retunedFinitePower(j) = R.radiation.finite_open_power;
    retunedOmega(j) = a0/lambda;
    retunedKappa(j) = retunedOmega(j)*sind(theta0);
end

fBIC_scale = Omega0*c0./(scaleFactors*a0);
scaleResidual = nan(size(scaleFactors));
scaleSigma = nan(size(scaleFactors));
scaleFinitePower = nan(size(scaleFactors));
scaleOmega = nan(size(scaleFactors));
scaleKappa = nan(size(scaleFactors));
for j = 1:numel(scaleFactors)
    s = scaleFactors(j);
    a = s*a0;
    lambda = c0/fBIC_scale(j);     % = s*a0/Omega0
    R = ni2019_strict_rayleigh_operator( ...
        make_cfg(a,lambda,theta0,s*depths0,s*widths0,s*gap0,N,K), ...
        'TargetOrder',-1,'EnforceAllFiniteOpen',true, ...
        'NormalizeMode',true,'Verbose',false, ...
        'KyTolerance',1e-7*max(abs(2*pi/lambda),1));
    scaleResidual(j) = R.strict_residual;
    scaleSigma(j) = R.sigma_ratio;
    scaleFinitePower(j) = R.radiation.finite_open_power;
    scaleOmega(j) = a/lambda;
    scaleKappa(j) = scaleOmega(j)*sind(theta0);
end

%% B. Environmental robustness and fixed-device control
% Similarity-retuned results above isolate the exact scale invariance.  For
% the control, keep the actual device, f=200 kHz, and theta=theta_BIC fixed.
% Then lambda=c/f changes and the n=-1 order is no longer at the nominal
% threshold except at c0.  The full homogeneous operator is used here:
% its smallest singular vector supplies an ordinary real-frequency mode,
% and its minimum singular residual plus open-channel amplitude quantify
% the loss of the strict cancellation.  We also report the forced-zero
% strict operator residual, but label it as an over-constrained diagnostic,
% not as a BIC test away from a Rayleigh point.

fixedFrequency = f0*ones(size(cSweep));
fixedOmega = a0.*fixedFrequency./cSweep;
fixedKappa = fixedOmega*sind(theta0);
fixedRayleighDetuning = abs(abs(fixedKappa-1)-fixedOmega);
fixedTargetKyOverK0 = complex(nan(size(cSweep)));
fixedFullSigma = nan(size(cSweep));
fixedFullRawResidual = nan(size(cSweep));
fixedOpenAmplitudeRatio = nan(size(cSweep));
fixedOpenPower = nan(size(cSweep));
fixedOpenPowerFraction = nan(size(cSweep));
fixedForcedResidual = nan(size(cSweep));
fixedTargetThreshold = false(size(cSweep));
fixedOpenOrders = cell(size(cSweep));

for j = 1:numel(cSweep)
    lambda = cSweep(j)/f0;
    cfg = make_cfg(a0,lambda,theta0,depths0,widths0,gap0,N,K);
    op = ni2019_full_eigen_operator(cfg);
    D = full_mode_diagnostic(op);
    fixedTargetKyOverK0(j) = D.target_ky/abs(op.k0);
    fixedFullSigma(j) = D.sigma_ratio;
    fixedFullRawResidual(j) = D.raw_residual;
    fixedOpenAmplitudeRatio(j) = D.open_amplitude_ratio;
    fixedOpenPower(j) = D.open_power;
    fixedOpenPowerFraction(j) = D.open_power_fraction;
    fixedTargetThreshold(j) = D.target_is_threshold;
    fixedOpenOrders{j} = D.open_orders(:).';
    % This call is intentionally only a control.  At c ~= c0 the target
    % order is not grazing, so its residual cannot establish a BIC.
    Rforced = ni2019_strict_rayleigh_operator(cfg,'TargetOrder',-1, ...
        'EnforceAllFiniteOpen',true,'NormalizeMode',true,'Verbose',false, ...
        'KyTolerance',1e-7*max(abs(2*pi/lambda),1));
    fixedForcedResidual(j) = Rforced.strict_residual;
end

%% Tables and machine-readable output
scalingBySpeed = table(cSweep,fBIC_c, ...
    (fBIC_c-f0),retunedOmega,retunedKappa,retunedSigma, ...
    retunedResidual,retunedFinitePower, ...
    'VariableNames',{'sound_speed_m_per_s','f_BIC_Hz','delta_f_Hz', ...
    'Omega','kappa','sigma_ratio','strict_residual','finite_open_power'});
scalingByScale = table(scaleFactors,fBIC_scale, ...
    (fBIC_scale-f0),scaleOmega,scaleKappa,scaleSigma, ...
    scaleResidual,scaleFinitePower, ...
    'VariableNames',{'scale_factor','f_BIC_Hz','delta_f_Hz', ...
    'Omega','kappa','sigma_ratio','strict_residual','finite_open_power'});
environmentControl = table(cSweep,fixedFrequency,fixedOmega,fixedKappa, ...
    fixedRayleighDetuning,fixedTargetKyOverK0,fixedTargetThreshold, ...
    fixedFullSigma,fixedFullRawResidual,fixedOpenAmplitudeRatio, ...
    fixedOpenPower,fixedOpenPowerFraction,fixedForcedResidual, ...
    'VariableNames',{'sound_speed_m_per_s','fixed_frequency_Hz','Omega', ...
    'kappa','rayleigh_detuning','target_ky_over_k0','target_is_threshold', ...
    'full_operator_sigma_ratio','full_operator_raw_residual', ...
    'open_amplitude_ratio','open_power','open_power_fraction', ...
    'forced_zero_control_residual'});

writetable(scalingBySpeed,fullfile(outputDir,'fBIC_vs_sound_speed.csv'));
writetable(scalingByScale,fullfile(outputDir,'fBIC_vs_scale_factor.csv'));
writetable(environmentControl,fullfile(outputDir,'fixed_device_sound_speed_control.csv'));

summary = struct();
summary.nominal = nominalSummary;
summary.constants = struct('c0_m_per_s',c0,'f0_Hz',f0,'a0_m',a0, ...
    'lambda0_m',lambda0,'Omega0',Omega0,'kappa0',kappa0, ...
    'theta0_deg',theta0,'N',N,'K',K, ...
    'd_over_a',normalizedGeometry.d_over_a, ...
    'w_over_a',normalizedGeometry.w_over_a,'g_over_a',normalizedGeometry.g_over_a);
summary.scalingBySpeed = scalingBySpeed;
summary.scalingByScale = scalingByScale;
summary.environmentControl = environmentControl;
summary.paths = struct('speedCsv',fullfile(outputDir,'fBIC_vs_sound_speed.csv'), ...
    'scaleCsv',fullfile(outputDir,'fBIC_vs_scale_factor.csv'), ...
    'controlCsv',fullfile(outputDir,'fixed_device_sound_speed_control.csv'));
save(fullfile(outputDir,'scaling_environment_data.mat'), ...
    'summary','nominal','scalingBySpeed','scalingByScale','environmentControl', ...
    'cSweep','scaleFactors','fixedOpenOrders','fixedTargetKyOverK0');

%% C. Publication-quality vector figures
navy = [0.035 0.105 0.300];
orange = [0.900 0.280 0.070];
teal = [0.000 0.460 0.430];
gray = [0.420 0.420 0.420];
lightGray = [0.78 0.78 0.78];

% Scaling: the two independent dimensional axes collapse onto the same
% similarity law without hiding the exact numerical modal evaluation.
figScaling = figure('Visible','off','Color','w','Units','inches', ...
    'Position',[1 1 7.05 3.30]);
tl = tiledlayout(figScaling,1,2,'TileSpacing','compact','Padding','compact');
nexttile(tl,1);
plot(cSweep,fBIC_c/1e3,'-','Color',navy,'LineWidth',1.8); hold on;
plot(c0,f0/1e3,'o','Color',orange,'MarkerFaceColor',orange,'MarkerSize',6);
yline(f0/1e3,'--','Color',lightGray,'LineWidth',1.0);
xlabel('$c$ (m s$^{-1}$)','Interpreter','latex');
ylabel('$f_{\rm BIC}$ (kHz)','Interpreter','latex');
title('(a) Sound-speed scaling','FontWeight','normal');
legend({'$f_{\rm BIC}=\Omega_0c/a_0$','nominal'}, ...
    'Interpreter','latex','Location','northwest','Box','off');
format_axes(gca);
nexttile(tl,2);
plot(scaleFactors,fBIC_scale/1e3,'-','Color',teal,'LineWidth',1.8); hold on;
plot(1,f0/1e3,'o','Color',orange,'MarkerFaceColor',orange,'MarkerSize',6);
yline(f0/1e3,'--','Color',lightGray,'LineWidth',1.0);
xlabel('$s=a/a_0$','Interpreter','latex');
ylabel('$f_{\rm BIC}$ (kHz)','Interpreter','latex');
title('(b) Geometric scaling','FontWeight','normal');
legend({'$f_{\rm BIC}=f_0/s$','nominal'}, ...
    'Interpreter','latex','Location','northeast','Box','off');
format_axes(gca);
exportgraphics(figScaling,fullfile(outputDir,'fig_scaling.pdf'), ...
    'ContentType','vector');
close(figScaling);

% Environment: exact retuning remains dark; a fixed device leaves the
% Rayleigh point and develops ordinary open-channel radiation.
figEnv = figure('Visible','off','Color','w','Units','inches', ...
    'Position',[1 1 7.05 6.45]);
tl = tiledlayout(figEnv,2,2,'TileSpacing','compact','Padding','compact');
nexttile(tl,1);
plot(cSweep,fBIC_c/1e3,'-','Color',navy,'LineWidth',1.8); hold on;
plot(c0,f0/1e3,'o','Color',orange,'MarkerFaceColor',orange,'MarkerSize',6);
xlabel('$c$ (m s$^{-1}$)','Interpreter','latex');
ylabel('$f_{\rm BIC}$ (kHz)','Interpreter','latex');
title('(a) Retuned BIC frequency','FontWeight','normal');
legend({'similarity retuning','nominal'},'Location','northwest','Box','off');
format_axes(gca);
nexttile(tl,2);
semilogy(cSweep,max(retunedResidual,realmin),'-','Color',navy,'LineWidth',1.8); hold on;
semilogy(cSweep,max(retunedSigma,realmin),'--','Color',teal,'LineWidth',1.4);
xline(c0,':','Color',gray,'LineWidth',1.0);
xlabel('$c$ (m s$^{-1}$)','Interpreter','latex');
ylabel('dimensionless residual','Interpreter','latex');
title('(b) Retuned homogeneous cancellation','FontWeight','normal');
legend({'raw residual','singular-value ratio'},'Location','southwest','Box','off');
format_axes(gca);
nexttile(tl,3);
plot(cSweep,fixedFullRawResidual,'-','Color',orange,'LineWidth',1.8); hold on;
plot(cSweep,fixedOpenAmplitudeRatio,'--','Color',teal,'LineWidth',1.5);
xline(c0,':','Color',gray,'LineWidth',1.0);
xlabel('$c$ (m s$^{-1}$)','Interpreter','latex');
ylabel('dimensionless control metric','Interpreter','latex');
title('(c) Fixed device and frequency','FontWeight','normal');
legend({'full-operator residual','open-amplitude ratio'}, ...
    'Location','northwest','Box','off');
format_axes(gca);
nexttile(tl,4);
plot(cSweep,fixedOmega,'-','Color',navy,'LineWidth',1.8); hold on;
plot(cSweep,abs(fixedKappa-1),'-','Color',orange,'LineWidth',1.6);
yline(1,'--','Color',lightGray,'LineWidth',1.0);
xline(c0,':','Color',gray,'LineWidth',1.0);
xlabel('$c$ (m s$^{-1}$)','Interpreter','latex');
ylabel('normalized frequency / threshold','Interpreter','latex');
title('(d) Leaving the Rayleigh condition','FontWeight','normal');
legend({'$\Omega=a f/c$','$|\kappa-1|$'}, ...
    'Interpreter','latex','Location','best','Box','off');
format_axes(gca);
exportgraphics(figEnv,fullfile(outputDir,'fig_environment.pdf'), ...
    'ContentType','vector');
close(figEnv);

fprintf('\nExact similarity retuning\n');
fprintf('  c sweep: f_BIC = %.6f--%.6f kHz; residual range %.3e--%.3e\n', ...
    min(fBIC_c)/1e3,max(fBIC_c)/1e3,min(retunedResidual),max(retunedResidual));
fprintf('  scale sweep: f_BIC = %.6f--%.6f kHz; residual range %.3e--%.3e\n', ...
    min(fBIC_scale)/1e3,max(fBIC_scale)/1e3,min(scaleResidual),max(scaleResidual));
fprintf('  max |Omega-Omega0| over both sweeps = %.3e\n', ...
    max([max(abs(retunedOmega-Omega0)),max(abs(scaleOmega-Omega0))]));
fprintf('\nFixed device/frequency control\n');
fprintf('  fixed-c residual range %.3e--%.3e; open-amplitude ratio range %.3e--%.3e\n', ...
    min(fixedFullRawResidual),max(fixedFullRawResidual), ...
    min(fixedOpenAmplitudeRatio),max(fixedOpenAmplitudeRatio));
fprintf('  nominal target threshold flag = %d; off-nominal flags = %d of %d\n', ...
    fixedTargetThreshold(cSweep==c0),nnz(~fixedTargetThreshold),numel(cSweep));
fprintf('\nOutputs\n  %s\n  %s\n  %s\n  %s\n  %s\n', ...
    fullfile(outputDir,'fBIC_vs_sound_speed.csv'), ...
    fullfile(outputDir,'fBIC_vs_scale_factor.csv'), ...
    fullfile(outputDir,'fixed_device_sound_speed_control.csv'), ...
    fullfile(outputDir,'scaling_environment_data.mat'), ...
    fullfile(outputDir,'fig_scaling.pdf'));
fprintf('  %s\n',fullfile(outputDir,'fig_environment.pdf'));

function cfg = make_cfg(a,lambda,theta,depths,widths,gap,N,K)
cfg = struct('a',a,'lambda',lambda,'theta_i_deg',theta, ...
    'depths',depths,'widths',widths,'gaps',gap, ...
    'N',N,'K',K,'solve_scattering',false);
end

function D = full_mode_diagnostic(op)
% Smallest full-operator singular vector, with no radiation amplitudes
% deleted.  This is an ordinary real-frequency control diagnostic.
[~,S,V] = svd(op.Fscaled,'econ');
y = V(:,end);
u = y./op.column_scale(:);
u = u/max(norm(u),eps);
A = u(1:op.N);
tol = 1e-7*max(abs(op.k0),1);
finiteOpen = abs(imag(op.ky))<=tol & real(op.ky)>tol;
threshold = abs(op.ky)<=tol;
powerByOrder = max(real(op.ky),0)./max(abs(op.k0),eps).*abs(A).^2;
openPower = sum(powerByOrder(finiteOpen));
totalPower = sum(powerByOrder);
D = struct('sigma_ratio',S(end,end)/max(S(1,1),eps), ...
    'raw_residual',norm(op.F*u)/max(norm(u),eps), ...
    'target_ky',op.ky(op.orders==-1), ...
    'target_is_threshold',any(threshold & op.orders==-1), ...
    'open_orders',op.orders(finiteOpen), ...
    'open_amplitude_ratio',norm(A(finiteOpen))/max(norm(A),eps), ...
    'open_power',openPower, ...
    'open_power_fraction',openPower/max(totalPower,eps));
end

function format_axes(ax)
set(ax,'Color','w','FontName','Helvetica','FontSize',9, ...
    'LineWidth',0.8,'TickDir','out','Box','on','XColor',[.1 .1 .1], ...
    'YColor',[.1 .1 .1],'GridColor',[.82 .82 .82]);
grid(ax,'on');
ax.Title.Color=[.1 .1 .1];
ax.XLabel.Color=[.1 .1 .1];
ax.YLabel.Color=[.1 .1 .1];
end
