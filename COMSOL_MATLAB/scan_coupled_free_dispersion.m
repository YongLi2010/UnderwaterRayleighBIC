%% Dense physical-sheet dispersion scan for the free stainless-steel plate
% Each kx point is solved as a coupled pressure-acoustics/solid-mechanics
% eigenproblem with a fluid PML.  The target acoustic branch is connected by
% pressure-field overlap, never by raw eigenvalue index.

clear; clc;
import com.comsol.model.*
import com.comsol.model.util.*

rootDir = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(rootDir,'COMSOL_MATLAB','rounded_180k_results');
modelFile = fullfile(outDir,'RayleighBIC_rounded_180k_coupled_free.mph');
model = mphload(modelFile,'CoupledDispersion');
model.study('std1').feature('eig').set('neigs',16);
model.study('std1').feature('eig').set('shift','f0');

a = model.param.evaluate('a');
c0 = model.param.evaluate('c0');
kappaSeed = model.param.evaluate('kxB')*a/(2*pi);
[coord,pStored] = pressureSamples(model,10);
fStoredAll = mphglobal(model,'freq');
fStored = fStoredAll(10);

dk = 0.002;
kappa = kappaSeed + (-35:35)*dk;
nK = numel(kappa);
seedIndex = 36;
branchF = complex(nan(nK,1));
branchQ = nan(nK,1);
branchOverlap = nan(nK,1);
branchModeIndex = nan(nK,1);
branchP = cell(nK,1);
allRows = cell(nK,1);

fprintf('Coupled free-plate scan: kappa %.6f to %.6f, dk %.4f\n', ...
    kappa(1),kappa(end),dk);

% Seed point, referenced to the already solved coupled mode.
[branchF(seedIndex),branchQ(seedIndex),branchOverlap(seedIndex), ...
    branchModeIndex(seedIndex),branchP{seedIndex},allRows{seedIndex}] = ...
    solvePoint(model,kappa(seedIndex),coord,pStored,fStored,a,c0);

% Continue toward Gamma.
for i = seedIndex-1:-1:1
    [branchF(i),branchQ(i),branchOverlap(i),branchModeIndex(i), ...
        branchP{i},allRows{i}] = solvePoint(model,kappa(i),coord, ...
        branchP{i+1},branchF(i+1),a,c0);
    fprintf('down %3d/%3d: kappa=%.6f, f=%.6f kHz, Q=%.4g, O=%.4f\n', ...
        seedIndex-i,nK-1,kappa(i),real(branchF(i))/1e3, ...
        branchQ(i),branchOverlap(i));
end

% Continue away from Gamma.
for i = seedIndex+1:nK
    [branchF(i),branchQ(i),branchOverlap(i),branchModeIndex(i), ...
        branchP{i},allRows{i}] = solvePoint(model,kappa(i),coord, ...
        branchP{i-1},branchF(i-1),a,c0);
    fprintf('up   %3d/%3d: kappa=%.6f, f=%.6f kHz, Q=%.4g, O=%.4f\n', ...
        i-seedIndex,nK-1,kappa(i),real(branchF(i))/1e3, ...
        branchQ(i),branchOverlap(i));
end

branchOmega = real(branchF)*a/c0;
rayleighOmega = 1-kappa(:);
deltaOmega = branchOmega-rayleighOmega;
branchTable = table(kappa(:),real(branchF),imag(branchF),branchQ, ...
    branchOverlap,branchModeIndex,branchOmega,rayleighOmega,deltaOmega, ...
    'VariableNames',{'kappa','frequency_real_hz','frequency_imag_hz','Q', ...
    'pressure_overlap','local_mode_index','Omega','Omega_Rayleigh', ...
    'delta_Omega'});
writetable(branchTable,fullfile(outDir, ...
    'coupled_free_target_branch_dense.csv'));

allTable = vertcat(allRows{:});
writetable(allTable,fullfile(outDir,'coupled_free_all_eigenvalues_dense.csv'));

% Locate the real-frequency crossing with the n=-1 Rayleigh line.
crossIndex = find(deltaOmega(1:end-1).*deltaOmega(2:end)<=0,1,'first');
if ~isempty(crossIndex)
    t = -deltaOmega(crossIndex)/(deltaOmega(crossIndex+1)-deltaOmega(crossIndex));
    kappaCross = kappa(crossIndex)+t*(kappa(crossIndex+1)-kappa(crossIndex));
    omegaCross = branchOmega(crossIndex)+t*(branchOmega(crossIndex+1)-branchOmega(crossIndex));
    fCross = omegaCross*c0/a;
    fprintf('Rayleigh crossing (linear interpolation): kappa=%.9f, Omega=%.9f, f=%.6f kHz\n', ...
        kappaCross,omegaCross,fCross/1e3);
else
    kappaCross=nan; omegaCross=nan; fCross=nan;
    fprintf('No Rayleigh crossing found inside scan window.\n');
end

% Direct negative-k checks of reciprocity; the full plot may then mirror the
% positive-k data without pretending the symmetry was untested.
checkIndices = [6,seedIndex,nK-5];
reciprocityRows = cell(numel(checkIndices),1);
for q = 1:numel(checkIndices)
    i = checkIndices(q);
    [fNeg,qNeg,oNeg,idxNeg] = solvePoint(model,-kappa(i),coord, ...
        conj(branchP{i}),branchF(i),a,c0);
    reciprocityRows{q} = table(kappa(i),real(branchF(i)),real(fNeg), ...
        real(fNeg-branchF(i)),qNeg,oNeg,idxNeg, ...
        'VariableNames',{'abs_kappa','positive_frequency_hz', ...
        'negative_frequency_hz','difference_hz','negative_Q', ...
        'pressure_overlap','negative_mode_index'});
end
reciprocityTable = vertcat(reciprocityRows{:});
writetable(reciprocityTable,fullfile(outDir, ...
    'coupled_free_reciprocity_checks.csv'));

save(fullfile(outDir,'coupled_free_dispersion_dense.mat'), ...
    'kappa','branchF','branchQ','branchOverlap','branchModeIndex', ...
    'branchOmega','rayleighOmega','deltaOmega','kappaCross','omegaCross', ...
    'fCross','allTable','reciprocityTable');
fprintf('Dispersion scan completed.\n');


function [fSelected,qSelected,oSelected,idxSelected,pSelected,allTable] = ...
    solvePoint(model,kappa,coord,pPrevious,fPrevious,a,c0)
model.param.set('kxB',sprintf('2*pi*(%.16g)/a',kappa));
model.study('std1').run;
freq = mphglobal(model,'freq');
freq = freq(:);
n = numel(freq);
omega = real(freq)*a/c0;
qAll = abs(real(freq)./(2*imag(freq)));
allTable = table(repmat(kappa,n,1),(1:n).',real(freq),imag(freq), ...
    qAll,omega,'VariableNames',{'kappa','local_mode_index', ...
    'frequency_real_hz','frequency_imag_hz','Q','Omega'});

candidate = find(abs(real(freq)-real(fPrevious))<10e3 & ...
    real(freq)>100e3 & real(freq)<260e3);
if isempty(candidate)
    candidate = find(real(freq)>100e3 & real(freq)<260e3);
end
overlap = -inf(n,1);
pCache = cell(n,1);
for j = reshape(candidate,1,[])
    pj = mphinterp(model,'acpr.p_t','coord',coord,'solnum',j);
    valid = isfinite(pj(:)) & isfinite(pPrevious(:));
    v = pj(valid); r = pPrevious(valid);
    denom = norm(v)*norm(r);
    if denom>0
        overlap(j)=abs(r(:)'*v(:))/denom;
    end
    pCache{j}=pj;
end
[oSelected,idxSelected] = max(overlap);
if ~isfinite(oSelected) || oSelected<0.35
    error('Branch tracking failed at kappa %.9f (best overlap %.5g).', ...
        kappa,oSelected);
end
fSelected = freq(idxSelected);
qSelected = qAll(idxSelected);
pSelected = pCache{idxSelected};
end


function [coord,p] = pressureSamples(model,solnum)
x1 = model.param.evaluate('x1');
x2 = model.param.evaluate('x2');
w1 = model.param.evaluate('w1');
w2 = model.param.evaluate('w2');
d1 = model.param.evaluate('d1');
d2 = model.param.evaluate('d2');
a = model.param.evaluate('a');
[X1,Y1] = meshgrid(linspace(x1+0.03*w1,x1+0.97*w1,28), ...
    linspace(-0.97*d1,-0.03*d1,22));
[X2,Y2] = meshgrid(linspace(x2+0.05*w2,x2+0.95*w2,10), ...
    linspace(-0.94*d2,-0.06*d2,8));
[XE,YE] = meshgrid(linspace(0.02*a,0.98*a,28), ...
    linspace(0.08e-3,5e-3,12));
coord = [X1(:).',X2(:).',XE(:).';Y1(:).',Y2(:).',YE(:).'];
p = mphinterp(model,'acpr.p_t','coord',coord,'solnum',solnum);
end
