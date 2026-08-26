%% Incident-field reference for finite-aperture field separation
% The sample is removed and the lower boundary is made radiating.  With the
% same 80-mm phased source, this gives p_inc on the exact grid used by the
% finite-sample model.  The reflected/scattered field is p_total-p_inc.
clear; clc;
import com.comsol.model.*
import com.comsol.model.util.*

rootDir=fileparts(fileparts(mfilename('fullpath')));
outDir=fullfile(rootDir,'COMSOL_MATLAB','rounded_180k_results', ...
    'finite_20period_experiment');
ModelUtil.remove('FiniteIncidentReference');
model=ModelUtil.create('FiniteIncidentReference');
model.label('80 mm finite-aperture incident-field reference');

model.param.set('W','148.4[mm]'); model.param.set('H','80[mm]');
model.param.set('Wsrc','80[mm]'); model.param.set('xsrcL','(W-Wsrc)/2');
model.param.set('xsrcR','(W+Wsrc)/2');
model.param.set('c0','1500[m/s]'); model.param.set('rho0','1000[kg/m^3]');
model.param.set('thetaSrc','7.05[deg]'); model.param.set('aSrc','1[m/s^2]');
model.param.set('wBeam','25[mm]');

model.component.create('comp1',true);
model.component('comp1').geom.create('geom1',2);
g=model.component('comp1').geom('geom1'); g.lengthUnit('m');
g.feature.create('waterL','Rectangle');
g.feature('waterL').set('size',{'xsrcL','H'}); g.feature('waterL').set('pos',{'0','0'});
g.feature.create('waterM','Rectangle');
g.feature('waterM').set('size',{'Wsrc','H'}); g.feature('waterM').set('pos',{'xsrcL','0'});
g.feature.create('waterR','Rectangle');
g.feature('waterR').set('size',{'W-xsrcR','H'}); g.feature('waterR').set('pos',{'xsrcR','0'});
g.feature.create('fluid','Union');
g.feature('fluid').selection('input').set({'waterL','waterM','waterR'});
g.feature('fluid').set('intbnd',true); g.feature('fluid').set('selresult',true);
g.run;

fluid=mphgetselection(model,'geom1_fluid_dom'); fluidDom=fluid.entities;
W=model.param.evaluate('W'); H=model.param.evaluate('H');
xL=model.param.evaluate('xsrcL'); xR=model.param.evaluate('xsrcR'); tol=1e-9;
sourceB=mphselectbox(model,'geom1',[xL-tol,xR+tol;H-tol,H+tol],'boundary');
topB=mphselectbox(model,'geom1',[-tol,W+tol;H-tol,H+tol],'boundary');
bottomB=mphselectbox(model,'geom1',[-tol,W+tol;-tol,tol],'boundary');
leftB=mphselectbox(model,'geom1',[-tol,tol;-tol,H+tol],'boundary');
rightB=mphselectbox(model,'geom1',[W-tol,W+tol;-tol,H+tol],'boundary');
radiationB=setdiff(unique([topB,bottomB,leftB,rightB]),sourceB);
fprintf('Reference selections: domain=%s source=%s radiation=%s\n', ...
    mat2str(fluidDom),mat2str(sourceB),mat2str(radiationB));

model.component('comp1').material.create('water','Common');
model.component('comp1').material('water').selection.set(fluidDom);
model.component('comp1').material('water').propertyGroup('def').set('density','rho0');
model.component('comp1').material('water').propertyGroup('def').set('soundspeed','c0');
model.component('comp1').physics.create('acpr','PressureAcoustics','geom1');
model.component('comp1').physics('acpr').selection.set(fluidDom);
model.component('comp1').physics('acpr').create('src1','NormalAcceleration',1);
model.component('comp1').physics('acpr').feature('src1').selection.set(sourceB);
model.component('comp1').physics('acpr').feature('src1').set('nacc', ...
    'aSrc*exp(-((x-W/2)/wBeam)^2)*exp(-i*2*pi*freq/c0*sin(thetaSrc)*(x-W/2))');
model.component('comp1').physics('acpr').create('pwr1','PlaneWaveRadiation',1);
model.component('comp1').physics('acpr').feature('pwr1').selection.set(radiationB);

model.component('comp1').mesh.create('mesh1');
model.component('comp1').mesh('mesh1').feature('size').set('custom','on');
model.component('comp1').mesh('mesh1').feature('size').set('hmax','0.55[mm]');
model.component('comp1').mesh('mesh1').feature('size').set('hmin','0.05[mm]');
model.component('comp1').mesh('mesh1').feature('size').set('hgrad',1.25);
model.component('comp1').mesh('mesh1').create('ftri1','FreeTri');
model.component('comp1').mesh('mesh1').run;

f0=180046.681327;
model.study.create('std1'); model.study('std1').create('freq','Frequency');
model.study('std1').feature('freq').set('plist',sprintf('%.12g[Hz]',f0));
model.study('std1').setGenPlots(false); model.study('std1').run;

x=readmatrix(fullfile(outDir,'finite_bic_field_x_m.csv')).';
y=readmatrix(fullfile(outDir,'finite_bic_field_y_m.csv')).';
[X,Y]=meshgrid(x,y); mask=Y>=0;
p=mphinterp(model,'acpr.p_t','coord',[X(mask).';Y(mask).']);
field=complex(nan(size(X))); field(mask)=p;
writematrix(real(field),fullfile(outDir,'finite_gaussian_incident_real_pa.csv'));
writematrix(imag(field),fullfile(outDir,'finite_gaussian_incident_imag_pa.csv'));
mphsave(model,fullfile(outDir,'RayleighBIC_finite_80mm_gaussian_incident_reference.mph'));

fprintf('Incident reference complete: max |p_inc| = %.6g Pa.\n',max(abs(p)));
