%% Finite 20-period rigid sample with an 80-mm phased source aperture
% Two-dimensional numerical experiment.  The explicit stainless-steel
% domain is solved in the sound-hard limit to match the manuscript theory.

clear; clc;
import com.comsol.model.*
import com.comsol.model.util.*

rootDir=fileparts(fileparts(mfilename('fullpath')));
outDir=fullfile(rootDir,'COMSOL_MATLAB','rounded_180k_results', ...
    'finite_20period_experiment');
if ~exist(outDir,'dir'), mkdir(outDir); end
ModelUtil.remove('Finite20Period');
model=ModelUtil.create('Finite20Period');
model.label('20-period rigid underwater sample - 80 mm source aperture');

% Geometry and experiment parameters.
model.param.set('a','7.42[mm]');
model.param.set('w1','4.41[mm]'); model.param.set('d1','5.15[mm]');
model.param.set('w2','1.00[mm]'); model.param.set('d2','1.32[mm]');
model.param.set('g12','1.00[mm]'); model.param.set('g21','1.01[mm]');
model.param.set('x1','g21/2'); model.param.set('x2','x1+w1+g12');
model.param.set('Ncell','20'); model.param.set('W','Ncell*a');
model.param.set('tPlate','10[mm]'); model.param.set('H','80[mm]');
model.param.set('Wsrc','80[mm]'); model.param.set('xsrcL','(W-Wsrc)/2');
model.param.set('xsrcR','(W+Wsrc)/2');
model.param.set('c0','1500[m/s]'); model.param.set('rho0','1000[kg/m^3]');
model.param.set('thetaSrc','7.10[deg]'); model.param.set('aSrc','1[m/s^2]');

model.component.create('comp1',true);
model.component('comp1').geom.create('geom1',2);
g=model.component('comp1').geom('geom1'); g.lengthUnit('m');

% Stainless plate and forty grooves.
g.feature.create('plate','Rectangle');
g.feature('plate').set('size',{'W','tPlate'});
g.feature('plate').set('pos',{'0','-tPlate'});
cutGrooveTags=cell(1,40);
fluidGrooveTags=cell(1,40);
for cellId=0:19
    c1=sprintf('cut1_%02d',cellId+1); c2=sprintf('cut2_%02d',cellId+1);
    f1=sprintf('fill1_%02d',cellId+1); f2=sprintf('fill2_%02d',cellId+1);
    cutGrooveTags{2*cellId+1}=c1; cutGrooveTags{2*cellId+2}=c2;
    fluidGrooveTags{2*cellId+1}=f1; fluidGrooveTags{2*cellId+2}=f2;
    g.feature.create(c1,'Rectangle');
    g.feature(c1).set('size',{'w1','d1'});
    g.feature(c1).set('pos',{sprintf('%d*a+x1',cellId),'-d1'});
    g.feature.create(c2,'Rectangle');
    g.feature(c2).set('size',{'w2','d2'});
    g.feature(c2).set('pos',{sprintf('%d*a+x2',cellId),'-d2'});
    g.feature.create(f1,'Rectangle');
    g.feature(f1).set('size',{'w1','d1'});
    g.feature(f1).set('pos',{sprintf('%d*a+x1',cellId),'-d1'});
    g.feature.create(f2,'Rectangle');
    g.feature(f2).set('size',{'w2','d2'});
    g.feature(f2).set('pos',{sprintf('%d*a+x2',cellId),'-d2'});
end
g.feature.create('steel','Difference');
g.feature('steel').selection('input').set({'plate'});
g.feature('steel').selection('input2').set(cutGrooveTags);
g.feature('steel').set('selresult',true);

% Three adjacent water rectangles retain vertices at the 80-mm aperture
% endpoints after the union, allowing an exact source-boundary selection.
g.feature.create('waterL','Rectangle');
g.feature('waterL').set('size',{'xsrcL','H'}); g.feature('waterL').set('pos',{'0','0'});
g.feature.create('waterM','Rectangle');
g.feature('waterM').set('size',{'Wsrc','H'}); g.feature('waterM').set('pos',{'xsrcL','0'});
g.feature.create('waterR','Rectangle');
g.feature('waterR').set('size',{'W-xsrcR','H'}); g.feature('waterR').set('pos',{'xsrcR','0'});
g.feature.create('fluid','Union');
g.feature('fluid').selection('input').set([{'waterL','waterM','waterR'},fluidGrooveTags]);
g.feature('fluid').set('intbnd',false); g.feature('fluid').set('selresult',true);
g.feature.create('all','Union');
g.feature('all').selection('input').set({'steel','fluid'});
g.feature('all').set('intbnd',true);
g.run;

fluid=mphgetselection(model,'geom1_fluid_dom'); fluidDom=fluid.entities;
steel=mphgetselection(model,'geom1_steel_dom'); steelDom=steel.entities;
W=model.param.evaluate('W'); H=model.param.evaluate('H');
Wsrc=model.param.evaluate('Wsrc'); xL=0.5*(W-Wsrc); xR=0.5*(W+Wsrc); tol=1e-9;
sourceB=mphselectbox(model,'geom1',[xL-tol,xR+tol;H-tol,H+tol],'boundary');
topB=mphselectbox(model,'geom1',[-tol,W+tol;H-tol,H+tol],'boundary');
leftB=mphselectbox(model,'geom1',[-tol,tol;-tol,H+tol],'boundary');
rightB=mphselectbox(model,'geom1',[W-tol,W+tol;-tol,H+tol],'boundary');
radiationB=setdiff(unique([topB,leftB,rightB]),sourceB);
fprintf('Finite selections: fluid=%s steel=%s source=%s radiation=%s\n', ...
    mat2str(fluidDom),mat2str(steelDom),mat2str(sourceB),mat2str(radiationB));
if isempty(sourceB), error('The 80-mm top source boundary was not resolved.'); end

model.component('comp1').material.create('water','Common');
model.component('comp1').material('water').selection.set(fluidDom);
model.component('comp1').material('water').propertyGroup('def').set('density','rho0');
model.component('comp1').material('water').propertyGroup('def').set('soundspeed','c0');
model.component('comp1').material.create('steelmat','Common');
model.component('comp1').material('steelmat').label('Stainless steel');
model.component('comp1').material('steelmat').selection.set(steelDom);
model.component('comp1').material('steelmat').propertyGroup('def').set('density','7850[kg/m^3]');
model.component('comp1').material('steelmat').propertyGroup('def').set('youngsmodulus','200[GPa]');
model.component('comp1').material('steelmat').propertyGroup('def').set('poissonsratio','0.30');

model.component('comp1').physics.create('acpr','PressureAcoustics','geom1');
model.component('comp1').physics('acpr').selection.set(fluidDom);
model.component('comp1').physics('acpr').create('src1','NormalAcceleration',1);
model.component('comp1').physics('acpr').feature('src1').selection.set(sourceB);
model.component('comp1').physics('acpr').feature('src1').set('nacc', ...
    'aSrc*exp(-i*2*pi*freq/c0*sin(thetaSrc)*(x-W/2))');
model.component('comp1').physics('acpr').create('pwr1','PlaneWaveRadiation',1);
model.component('comp1').physics('acpr').feature('pwr1').selection.set(radiationB);

model.component('comp1').mesh.create('mesh1');
model.component('comp1').mesh('mesh1').feature('size').set('custom','on');
model.component('comp1').mesh('mesh1').feature('size').set('hmax','0.55[mm]');
model.component('comp1').mesh('mesh1').feature('size').set('hmin','0.035[mm]');
model.component('comp1').mesh('mesh1').feature('size').set('hgrad',1.25);
model.component('comp1').mesh('mesh1').feature('size').set('hcurve',0.25);
model.component('comp1').mesh('mesh1').feature('size').set('hnarrow',0.75);
model.component('comp1').mesh('mesh1').create('ftri1','FreeTri');
model.component('comp1').mesh('mesh1').run;

model.study.create('std1'); model.study('std1').create('freq','Frequency');
model.study('std1').feature('freq').set('plist', ...
    'range(179.60[kHz],0.020[kHz],180.40[kHz])');
model.study('std1').setGenPlots(false);
model.study('std1').run;

freq=real(mphglobal(model,'freq')); freq=freq(:);
centralCell=10; probeX=centralCell*model.param.evaluate('a')+ ...
    model.param.evaluate('x1+w1/2');
probe=mphinterp(model,'acpr.p_t','coord',[probeX;-0.65*model.param.evaluate('d1')]);
probe=probe(:); [~,peakIndex]=max(abs(probe));

% Export the peak field and a hydrophone scan plane 20 mm above the sample.
x=linspace(0,W,1001); y=linspace(-6e-3,H,601);
[X,Y]=meshgrid(x,y);
mask=Y>=0;
for cellId=0:19
    gx1=cellId*model.param.evaluate('a')+model.param.evaluate('x1');
    gx2=cellId*model.param.evaluate('a')+model.param.evaluate('x2');
    mask=mask | (X>=gx1 & X<=gx1+model.param.evaluate('w1') & ...
        Y>=-model.param.evaluate('d1')) | ...
        (X>=gx2 & X<=gx2+model.param.evaluate('w2') & ...
        Y>=-model.param.evaluate('d2'));
end
p=mphinterp(model,'acpr.p_t','coord',[X(mask).';Y(mask).'],'solnum',peakIndex);
field=complex(nan(size(X))); field(mask)=p;
writematrix(x(:),fullfile(outDir,'finite_field_x_m.csv'));
writematrix(y(:),fullfile(outDir,'finite_field_y_m.csv'));
writematrix(real(field),fullfile(outDir,'finite_field_real_pa.csv'));
writematrix(imag(field),fullfile(outDir,'finite_field_imag_pa.csv'));

xScan=linspace(0,W,1001); yScan=20e-3*ones(size(xScan));
pScan=mphinterp(model,'acpr.p_t','coord',[xScan;yScan],'solnum',peakIndex);
scanTable=table(xScan(:),yScan(:),real(pScan(:)),imag(pScan(:)),abs(pScan(:)), ...
    'VariableNames',{'x_m','y_m','pressure_real_pa','pressure_imag_pa','pressure_abs_pa'});
writetable(scanTable,fullfile(outDir,'hydrophone_scanline_peak.csv'));

spectrum=table(freq,real(probe),imag(probe),abs(probe),angle(probe), ...
    'VariableNames',{'frequency_hz','probe_real_pa','probe_imag_pa', ...
    'probe_abs_pa','probe_phase_rad'});
writetable(spectrum,fullfile(outDir,'finite_sample_frequency_spectrum.csv'));

model.result.create('pg1','PlotGroup2D'); model.result('pg1').set('looplevel',peakIndex);
model.result('pg1').feature.create('surf1','Surface');
model.result('pg1').feature('surf1').set('expr','real(acpr.p_t)');
model.result('pg1').feature('surf1').set('colortable','Wave');
mphsave(model,fullfile(outDir, ...
    'RayleighBIC_finite_20period_80mm_source_hard_boundary.mph'));

summary=table(20,W*1e3,Wsrc*1e3,model.param.evaluate('thetaSrc')*180/pi, ...
    freq(peakIndex),abs(probe(peakIndex)),peakIndex, ...
    'VariableNames',{'number_of_periods','sample_width_mm','source_width_mm', ...
    'source_angle_deg','peak_frequency_hz','peak_probe_abs_pa','peak_solution'});
writetable(summary,fullfile(outDir,'finite_model_summary.csv'));
disp(summary);
