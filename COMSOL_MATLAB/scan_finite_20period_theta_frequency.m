%% Refined finite-sample angle-frequency scan around the rounded quasi-BIC
% Loads the solved 20-period, 80-mm-aperture model and resolves the narrow
% response that a 20-Hz exploratory sweep can miss.

clear; clc;
import com.comsol.model.*
import com.comsol.model.util.*

rootDir=fileparts(fileparts(mfilename('fullpath')));
outDir=fullfile(rootDir,'COMSOL_MATLAB','rounded_180k_results', ...
    'finite_20period_experiment');
modelFile=fullfile(outDir, ...
    'RayleighBIC_finite_20period_80mm_source_hard_boundary.mph');
model=mphload(modelFile,'Finite20PeriodScan');

thetaList=[6.80 6.95 7.05 7.15 7.30];
freqList=(179.960e3:2:180.120e3).';
model.study('std1').feature('freq').set('plist', ...
    sprintf('range(%.12g[Hz],2[Hz],%.12g[Hz])',freqList(1),freqList(end)));

a=model.param.evaluate('a'); x1=model.param.evaluate('x1');
w1=model.param.evaluate('w1'); d1=model.param.evaluate('d1');
probeCells=7:12;
probeCoord=zeros(2,numel(probeCells));
for j=1:numel(probeCells)
    probeCoord(:,j)=[probeCells(j)*a+x1+w1/2;-0.65*d1];
end

allRows=[];
for it=1:numel(thetaList)
    th=thetaList(it);
    model.param.set('thetaSrc',sprintf('%.8g[deg]',th));
    fprintf('Finite scan %d/%d: theta = %.3f deg, %d frequencies\n', ...
        it,numel(thetaList),th,numel(freqList));
    model.study('std1').run;
    f=real(mphglobal(model,'freq')); f=f(:);
    p=mphinterp(model,'acpr.p_t','coord',probeCoord);
    if size(p,1)==numel(probeCells), pByState=p.'; else, pByState=p; end
    probeRms=sqrt(mean(abs(pByState).^2,2));
    probeCoherent=mean(pByState,2);
    rows=table(repmat(th,numel(f),1),f,probeRms, ...
        real(probeCoherent),imag(probeCoherent),abs(probeCoherent), ...
        angle(probeCoherent),'VariableNames',{'theta_deg','frequency_hz', ...
        'central_six_probe_rms_pa','coherent_probe_real_pa', ...
        'coherent_probe_imag_pa','coherent_probe_abs_pa','coherent_probe_phase_rad'});
    allRows=[allRows;rows]; %#ok<AGROW>
end

writetable(allRows,fullfile(outDir,'finite_theta_frequency_refined.csv'));

% Save the model with the final angular cut and its complete frequency set.
mphsave(model,fullfile(outDir, ...
    'RayleighBIC_finite_20period_80mm_source_refined_scan.mph'));

% Report one maximum per angular cut without interpreting it as a strict BIC.
summary=table('Size',[numel(thetaList),4],'VariableTypes', ...
    {'double','double','double','double'},'VariableNames', ...
    {'theta_deg','peak_frequency_hz','peak_probe_rms_pa','sampled_fwhm_hz'});
for it=1:numel(thetaList)
    use=allRows.theta_deg==thetaList(it);
    f=allRows.frequency_hz(use); y=allRows.central_six_probe_rms_pa(use).^2;
    [ym,im]=max(y); baseline=prctile(y,10); half=baseline+(ym-baseline)/2;
    above=find(y>=half);
    if numel(above)>1, width=f(above(end))-f(above(1)); else, width=NaN; end
    summary{it,:}=[thetaList(it),f(im),sqrt(ym),width];
end
writetable(summary,fullfile(outDir,'finite_theta_frequency_summary.csv'));
disp(summary);
