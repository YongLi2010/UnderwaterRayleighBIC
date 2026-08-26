%% Figure 4 source data: experimentally accessible driven response
% Data only.  All figure drawing and export are performed by the selected
% Python backend.  The calculation follows the final 180-kHz eigenpole with
% an adaptive frequency grid measured in local linewidths, then evaluates
% the full N=345, K=43 driven modal-matching system at every point.
clear; clc;
rootDir = fileparts(mfilename('fullpath'));
outDir = fullfile(rootDir,'results','fig4_experimental_observables_180k');
if ~exist(outDir,'dir'), mkdir(outDir); end
addpath(rootDir);

S = load(fullfile(rootDir,'results', ...
    'StrictRayleighBIC_180kHz_7p10deg_final.mat'));
Tall = readtable(fullfile(rootDir,'results', ...
    'fig3_radiation_contributions_180k','radiation_contributions.csv'));

cWater = 1500;
aPhysical = S.a;
kBIC = S.x(1);
OBIC = 1-kBIC;
thetaBIC = asind(kBIC/OBIC);

% Keep the driven far-field analysis on the physical open side of the
% Rayleigh point, where both n=0 and n=-1 carry normal power.
thetaAll = asind(Tall.kappa./Tall.Omega_real);
openMask = Tall.kappa >= kBIC-1e-13 & thetaAll <= thetaBIC+1.45;
T = Tall(openMask,:);

% A wide normalized-detuning stencil resolves the complex line shape and
% also provides a local off-resonant background on both sides of the pole.
u = unique([linspace(-32,-8,9),linspace(-8,-2,9), ...
    linspace(-2,2,49),linspace(2,8,9),linspace(8,32,9)]).';
nK = height(T); nU = numel(u); nP = nK*nU;

kappa = T.kappa;
OmegaReal = T.Omega_real;
OmegaImag = abs(T.Omega_imag);
thetaDeg = asind(kappa./OmegaReal);
poleFrequencyHz = OmegaReal*cWater/aPhysical;
linewidthHz = 2*OmegaImag*cWater/aPhysical;
rayleighFrequencyHz = cWater./(aPhysical*(1+sind(thetaDeg)));

A0flat = complex(nan(nP,1)); Am1flat = complex(nan(nP,1));
eta0flat = nan(nP,1); etam1flat = nan(nP,1);
conditionFlat = nan(nP,1); energyErrorFlat = nan(nP,1);
minDrivenLinewidthHz = 1e-5;

cfgBase = S.cfg;
cfgBase.N = 345;
cfgBase.K = 43;
cfgBase.solve_scattering = true;

fprintf('Fig. 4 adaptive driven response: %d angles x %d detunings.\n', ...
    nK,nU);
tic;
parfor idx = 1:nP
    [ik,iu] = ind2sub([nK,nU],idx);
    if linewidthHz(ik)<minDrivenLinewidthHz
        continue;
    end
    frequencyHz = poleFrequencyHz(ik) + 0.5*u(iu)*linewidthHz(ik);
    Omega = frequencyHz*aPhysical/cWater;
    cfg = cfgBase;
    cfg.lambda = 1/Omega;
    cfg.theta_i_deg = thetaDeg(ik);
    R = ni2019_modal_solver(cfg);
    id0 = find(R.orders==0,1);
    idm = find(R.orders==-1,1);
    A0flat(idx) = R.A(id0);
    Am1flat(idx) = R.A(idm);
    eta0flat(idx) = R.eta(id0);
    etam1flat(idx) = R.eta(idm);
    conditionFlat(idx) = R.condition_number;
    energyErrorFlat(idx) = R.energy_error;
end
elapsedSeconds = toc;

A0 = reshape(A0flat,[nK,nU]);
Am1 = reshape(Am1flat,[nK,nU]);
eta0 = reshape(eta0flat,[nK,nU]);
etam1 = reshape(etam1flat,[nK,nU]);
conditionNumber = reshape(conditionFlat,[nK,nU]);
energyError = reshape(energyErrorFlat,[nK,nU]);

% Background-subtracted resonant amplitudes and their integrated spectral
% weights.  Both quantities are directly obtainable from frequency sweeps.
A0Background = 0.5*(A0(:,1)+A0(:,end));
Am1Background = 0.5*(Am1(:,1)+Am1(:,end));
deltaA0 = A0-A0Background;
deltaAm1 = Am1-Am1Background;
frequencyGridHz = poleFrequencyHz + 0.5*linewidthHz.*u.';
spectralWeight0 = nan(nK,1); spectralWeightM1 = nan(nK,1);
peakContrast0 = nan(nK,1); peakContrastM1 = nan(nK,1);
for ik = 1:nK
    valid = isfinite(real(deltaA0(ik,:))) & isfinite(real(deltaAm1(ik,:)));
    if nnz(valid)<3, continue; end
    spectralWeight0(ik) = trapz(frequencyGridHz(ik,valid), ...
        abs(deltaA0(ik,valid)).^2);
    spectralWeightM1(ik) = trapz(frequencyGridHz(ik,valid), ...
        abs(deltaAm1(ik,valid)).^2);
    peakContrast0(ik) = max(abs(deltaA0(ik,valid)));
    peakContrastM1(ik) = max(abs(deltaAm1(ik,valid)));
end

% Long-form source table.  The exact BIC row is intentionally retained with
% NaN driven amplitudes: an ideal BIC is not a forced-scattering solution.
[IK,IU] = ndgrid((1:nK).',(1:nU).');
frequencyHz = poleFrequencyHz(IK) + 0.5*u(IU).*linewidthHz(IK);
response = table(kappa(IK(:)),thetaDeg(IK(:)), ...
    poleFrequencyHz(IK(:)),rayleighFrequencyHz(IK(:)), ...
    linewidthHz(IK(:)),u(IU(:)),frequencyHz(:), ...
    real(A0(:)),imag(A0(:)),real(Am1(:)),imag(Am1(:)), ...
    eta0(:),etam1(:),conditionNumber(:),energyError(:), ...
    'VariableNames',{'kappa','theta_deg','pole_frequency_hz', ...
    'rayleigh_frequency_hz','linewidth_hz','normalized_detuning', ...
    'frequency_hz','A0_real','A0_imag','Am1_real','Am1_imag', ...
    'eta0','etam1','condition_number','energy_error'});
writetable(response,fullfile(outDir,'adaptive_scattering_response.csv'));

summary = table(kappa,thetaDeg,poleFrequencyHz,rayleighFrequencyHz, ...
    linewidthHz,T.Q,T.b0_sum_abs,T.bm1_sum_abs, ...
    spectralWeight0,spectralWeightM1,peakContrast0,peakContrastM1, ...
    'VariableNames',{'kappa','theta_deg','pole_frequency_hz', ...
    'rayleigh_frequency_hz','linewidth_hz','Q', ...
    'radiation_A0','radiation_Am1','spectral_weight_A0', ...
    'spectral_weight_Am1','peak_contrast_A0','peak_contrast_Am1'});
writetable(summary,fullfile(outDir,'pole_and_radiation_summary.csv'));

% A representative open-side driven state for the measurement-emulation
% field panel.  It is deliberately not the unexcitable strict BIC.
[~,fieldIndex] = min(abs(thetaDeg-(thetaBIC+0.55)));
fieldFrequencyHz = poleFrequencyHz(fieldIndex);
fieldOmega = fieldFrequencyHz*aPhysical/cWater;
fieldThetaDeg = thetaDeg(fieldIndex);
cfgField = cfgBase;
cfgField.lambda = 1/fieldOmega;
cfgField.theta_i_deg = fieldThetaDeg;
RF = ni2019_modal_solver(cfgField);

xOverA = linspace(-5,5,801);
yOverA = linspace(0.12,3.2,321);
[X,Y] = meshgrid(xOverA,yOverA);
pIncident = exp(-1i*(RF.k0*sind(fieldThetaDeg))*X + ...
    1i*RF.ky_incident*Y);
pScattered = complex(zeros(size(X)));
for n = 1:numel(RF.orders)
    pScattered = pScattered + RF.A(n).* ...
        exp(-1i*RF.kx(n).*X - 1i*RF.ky(n).*Y);
end
pTotal = pIncident+pScattered;
fieldScale = max(abs(pTotal(:)));
pTotal = pTotal/max(fieldScale,eps);
pScattered = pScattered/max(fieldScale,eps);

writematrix(xOverA(:),fullfile(outDir,'field_x_over_a.csv'));
writematrix(yOverA(:),fullfile(outDir,'field_y_over_a.csv'));
writematrix(real(pTotal),fullfile(outDir,'field_total_real.csv'));
writematrix(imag(pTotal),fullfile(outDir,'field_total_imag.csv'));
writematrix(real(pScattered),fullfile(outDir,'field_scattered_real.csv'));
writematrix(imag(pScattered),fullfile(outDir,'field_scattered_imag.csv'));

id0 = find(RF.orders==0,1); idm = find(RF.orders==-1,1);
fieldMetadata = table(fieldThetaDeg,fieldFrequencyHz,kappa(fieldIndex), ...
    linewidthHz(fieldIndex),real(RF.A(id0)),imag(RF.A(id0)), ...
    real(RF.A(idm)),imag(RF.A(idm)),RF.eta(id0),RF.eta(idm), ...
    'VariableNames',{'theta_deg','frequency_hz','kappa', ...
    'linewidth_hz','A0_real','A0_imag','Am1_real','Am1_imag', ...
    'eta0','etam1'});
writetable(fieldMetadata,fullfile(outDir,'field_metadata.csv'));

save(fullfile(outDir,'fig4_experimental_observables.mat'), ...
    'S','T','u','response','summary','A0','Am1','eta0','etam1', ...
    'deltaA0','deltaAm1','frequencyGridHz','spectralWeight0', ...
    'spectralWeightM1','peakContrast0','peakContrastM1', ...
    'conditionNumber','energyError','xOverA','yOverA','pTotal', ...
    'pScattered','fieldMetadata','elapsedSeconds','-v7.3');

fid = fopen(fullfile(outDir,'README.txt'),'w');
fprintf(fid,'Figure 4 theory-only source data; no experimental observations.\n');
fprintf(fid,'Final design: f_BIC=180 kHz, theta_BIC=%.9f deg.\n',thetaBIC);
fprintf(fid,'Full driven truncation: N=345, K=43.\n');
fprintf(fid,'Adaptive grid: %d angles x %d normalized detunings.\n',nK,nU);
fprintf(fid,'Open-side angular window: %.9f to %.9f deg.\n', ...
    min(thetaDeg),max(thetaDeg));
fprintf(fid,'Exact BIC driven row is NaN because the ideal forced problem is singular/dark.\n');
fprintf(fid,'Discrete symbols in the figure are sampled theory placeholders, not experiment.\n');
fprintf(fid,'Runtime: %.3f s.\n',elapsedSeconds);
fclose(fid);

fprintf('Completed in %.2f s.  Max energy error %.3e.\n', ...
    elapsedSeconds,max(energyError(:),[],'omitnan'));
fprintf('Field state: theta %.6f deg, f %.6f kHz, eta0 %.6f, eta-1 %.6f.\n', ...
    fieldThetaDeg,fieldFrequencyHz/1e3,RF.eta(id0),RF.eta(idm));
