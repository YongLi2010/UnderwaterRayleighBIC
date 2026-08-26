%% Export finite-sample fields around the rounded periodic quasi-BIC
clear; clc;
import com.comsol.model.*
import com.comsol.model.util.*

rootDir=fileparts(fileparts(mfilename('fullpath')));
outDir=fullfile(rootDir,'COMSOL_MATLAB','rounded_180k_results', ...
    'finite_20period_experiment');
model=mphload(fullfile(outDir, ...
    'RayleighBIC_finite_20period_80mm_source_hard_boundary.mph'), ...
    'Finite20PeriodFields');

theta=7.05;
f0=180046.681327;
freqList=f0+[-5 0 5];
model.param.set('thetaSrc',sprintf('%.8g[deg]',theta));
model.study('std1').feature('freq').set('plist',sprintf('%.12g[Hz] %.12g[Hz] %.12g[Hz]',freqList));
model.study('std1').run;

W=model.param.evaluate('W'); H=model.param.evaluate('H');
a=model.param.evaluate('a'); x1=model.param.evaluate('x1'); x2=model.param.evaluate('x2');
w1=model.param.evaluate('w1'); w2=model.param.evaluate('w2');
d1=model.param.evaluate('d1'); d2=model.param.evaluate('d2');
x=linspace(0,W,1001); y=linspace(-6e-3,H,601); [X,Y]=meshgrid(x,y);
mask=Y>=0;
for cellId=0:19
    gx1=cellId*a+x1; gx2=cellId*a+x2;
    mask=mask | (X>=gx1 & X<=gx1+w1 & Y>=-d1) | ...
        (X>=gx2 & X<=gx2+w2 & Y>=-d2);
end
writematrix(x(:),fullfile(outDir,'finite_bic_field_x_m.csv'));
writematrix(y(:),fullfile(outDir,'finite_bic_field_y_m.csv'));

xScan=linspace(0,W,1001); yScan=20e-3*ones(size(xScan));
fieldSummary=table('Size',[3,4],'VariableTypes',repmat({'double'},1,4), ...
    'VariableNames',{'state','frequency_hz','max_pressure_pa','scanline_rms_pa'});
for state=1:3
    p=mphinterp(model,'acpr.p_t','coord',[X(mask).';Y(mask).'],'solnum',state);
    field=complex(nan(size(X))); field(mask)=p;
    writematrix(real(field),fullfile(outDir,sprintf('finite_bic_state_%d_real_pa.csv',state)));
    writematrix(imag(field),fullfile(outDir,sprintf('finite_bic_state_%d_imag_pa.csv',state)));
    ps=mphinterp(model,'acpr.p_t','coord',[xScan;yScan],'solnum',state);
    scanTable=table(xScan(:),yScan(:),real(ps(:)),imag(ps(:)),abs(ps(:)), ...
        'VariableNames',{'x_m','y_m','pressure_real_pa','pressure_imag_pa','pressure_abs_pa'});
    writetable(scanTable,fullfile(outDir,sprintf('finite_bic_state_%d_scanline.csv',state)));
    fieldSummary{state,:}=[state,freqList(state),max(abs(p)),sqrt(mean(abs(ps).^2))];
end
writetable(fieldSummary,fullfile(outDir,'finite_bic_field_summary.csv'));
mphsave(model,fullfile(outDir,'RayleighBIC_finite_20period_80mm_source_bic_fields.mph'));
disp(fieldSummary);
