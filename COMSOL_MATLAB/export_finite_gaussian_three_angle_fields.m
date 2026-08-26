%% Three representative finite-sample Gaussian-beam scattering states
% Each state follows the physical eigenpole branch: Gamma, the rounded-cell
% Rayleigh BIC, and a point above the BIC.  The same Gaussian source is solved
% with and without the sample so that p_ref = p_total - p_inc is obtained on
% an identical grid.  Complex scan-line data are exported for angular spectra.
clear; clc;
import com.comsol.model.*
import com.comsol.model.util.*

rootDir=fileparts(fileparts(mfilename('fullpath')));
outDir=fullfile(rootDir,'COMSOL_MATLAB','rounded_180k_results', ...
    'finite_20period_experiment');
sampleFile=fullfile(outDir, ...
    'RayleighBIC_finite_20period_80mm_source_hard_boundary.mph');
referenceFile=fullfile(outDir, ...
    'RayleighBIC_finite_80mm_gaussian_incident_reference.mph');
assert(isfile(sampleFile),'Missing finite-sample base model.');
assert(isfile(referenceFile), ...
    'Run build_finite_incident_reference.m before this script.');

sample=mphload(sampleFile,'FiniteGaussianThreeStateSample');
reference=mphload(referenceFile,'FiniteGaussianThreeStateReference');

sourceExpression=[ ...
    'aSrc*exp(-((x-W/2)/wBeam)^2)*' ...
    'exp(-i*2*pi*freq/c0*sin(thetaSrc)*(x-W/2))'];
for model={sample,reference}
    model{1}.param.set('wBeam','25[mm]');
    model{1}.component('comp1').physics('acpr').feature('src1').set( ...
        'nacc',sourceExpression);
end

% Frequencies are the real part of the physical outgoing eigenpole branch.
% State 2 uses the independently refined rounded-cell COMSOL Rayleigh point.
stateLabel={'Gamma','Rayleigh_BIC','above_BIC'};
thetaDeg=[0,7.05,10.05404];
frequencyHz=[178780.2005,180046.681327,180498.322];

x=readmatrix(fullfile(outDir,'finite_bic_field_x_m.csv')).';
y=readmatrix(fullfile(outDir,'finite_bic_field_y_m.csv')).';
[X,Y]=meshgrid(x,y);
waterMask=Y>=0;
xScan=linspace(0,sample.param.evaluate('W'),2001);
yScan=20e-3*ones(size(xScan));

summary=table('Size',[numel(thetaDeg),8], ...
    'VariableTypes',repmat({'double'},1,8), ...
    'VariableNames',{'state','theta_deg','frequency_hz','waist_mm', ...
    'incident_rms_pa','reflected_rms_pa','incident_max_pa','reflected_max_pa'});

for state=1:numel(thetaDeg)
    thetaString=sprintf('%.10g[deg]',thetaDeg(state));
    frequencyString=sprintf('%.12g[Hz]',frequencyHz(state));
    for model={sample,reference}
        model{1}.param.set('thetaSrc',thetaString);
        model{1}.study('std1').feature('freq').set('plist',frequencyString);
        model{1}.study('std1').run;
    end

    pTotal=mphinterp(sample,'acpr.p_t', ...
        'coord',[X(waterMask).';Y(waterMask).']);
    pIncident=mphinterp(reference,'acpr.p_t', ...
        'coord',[X(waterMask).';Y(waterMask).']);
    pReflected=pTotal-pIncident;

    incidentField=complex(nan(size(X)));
    reflectedField=complex(nan(size(X)));
    incidentField(waterMask)=pIncident;
    reflectedField(waterMask)=pReflected;

    prefix=sprintf('gaussian_state_%d_%s',state,lower(stateLabel{state}));
    writematrix(real(incidentField),fullfile(outDir,[prefix '_incident_real_pa.csv']));
    writematrix(imag(incidentField),fullfile(outDir,[prefix '_incident_imag_pa.csv']));
    writematrix(real(reflectedField),fullfile(outDir,[prefix '_reflected_real_pa.csv']));
    writematrix(imag(reflectedField),fullfile(outDir,[prefix '_reflected_imag_pa.csv']));

    pTotalScan=mphinterp(sample,'acpr.p_t','coord',[xScan;yScan]);
    pIncidentScan=mphinterp(reference,'acpr.p_t','coord',[xScan;yScan]);
    pReflectedScan=pTotalScan-pIncidentScan;
    scanTable=table(xScan(:),yScan(:),real(pIncidentScan(:)), ...
        imag(pIncidentScan(:)),real(pReflectedScan(:)), ...
        imag(pReflectedScan(:)),abs(pIncidentScan(:)), ...
        abs(pReflectedScan(:)), ...
        'VariableNames',{'x_m','y_m','incident_real_pa','incident_imag_pa', ...
        'reflected_real_pa','reflected_imag_pa','incident_abs_pa', ...
        'reflected_abs_pa'});
    writetable(scanTable,fullfile(outDir,[prefix '_scanline.csv']));

    summary{state,:}=[state,thetaDeg(state),frequencyHz(state),25, ...
        sqrt(mean(abs(pIncidentScan).^2)),sqrt(mean(abs(pReflectedScan).^2)), ...
        max(abs(pIncident)),max(abs(pReflected))];
    fprintf('%s: theta %.5f deg, f %.6f kHz, reflected RMS %.6g Pa\n', ...
        stateLabel{state},thetaDeg(state),frequencyHz(state)/1e3, ...
        summary.reflected_rms_pa(state));
end

writetable(summary,fullfile(outDir,'gaussian_three_state_summary.csv'));
mphsave(sample,fullfile(outDir, ...
    'RayleighBIC_finite_20period_gaussian_three_states.mph'));
mphsave(reference,fullfile(outDir, ...
    'RayleighBIC_finite_gaussian_reference_three_states.mph'));
disp(summary);
