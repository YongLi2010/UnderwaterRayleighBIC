%% Build and solve two COMSOL 6.4 models for the rounded 180-kHz design
% Requires a running COMSOL Multiphysics Server and LiveLink for MATLAB.
% The stainless-steel plate is represented explicitly with a 10-mm total
% thickness.  The acoustic calculation uses the sound-hard limit on the
% water--steel boundary, matching the modal-matching theory.  The steel
% material is retained in the MPH files for a later acoustic--structure
% coupling extension.

clear; clc;
import com.comsol.model.*
import com.comsol.model.util.*

rootDir = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(rootDir,'COMSOL_MATLAB','rounded_180k_results');
if ~exist(outDir,'dir'), mkdir(outDir); end

P = struct;
P.a = 7.42e-3;
P.w1 = 4.41e-3;
P.d1 = 5.15e-3;
P.w2 = 1.00e-3;
P.d2 = 1.32e-3;
P.g12 = 1.00e-3;
P.g21 = 1.01e-3;
P.x1 = 0.5*P.g21;
P.x2 = P.x1+P.w1+P.g12;
P.tPlate = 10.00e-3;
P.H = 25.0e-3;
P.tPML = 8.0e-3;
P.c = 1500;
P.rho = 1000;
P.f0 = 180e3;
P.theta = 7.10;

assert(abs(P.x2+P.w2+0.5*P.g21-P.a)<1e-12, ...
    'Rounded unit-cell dimensions do not close.');

fprintf('Building eigenfrequency model ...\n');
[eigenModel,eigenMeta] = buildEigenModel(P,outDir);
mphsave(eigenModel,fullfile(outDir,'RayleighBIC_rounded_180k_eigenfrequency.mph'));
writetable(struct2table(eigenMeta), ...
    fullfile(outDir,'eigenfrequency_summary.csv'));

fprintf('Building driven scattering model ...\n');
[scatModel,scatMeta] = buildScatteringModel(P,outDir);
mphsave(scatModel,fullfile(outDir,'RayleighBIC_rounded_180k_scattering.mph'));
writetable(struct2table(scatMeta), ...
    fullfile(outDir,'scattering_summary.csv'));

save(fullfile(outDir,'rounded_180k_comsol_metadata.mat'),'P','eigenMeta','scatMeta');
fprintf('Completed. Models written to %s\n',outDir);


function [model,meta] = buildEigenModel(P,outDir)
import com.comsol.model.*
import com.comsol.model.util.*
ModelUtil.remove('RayleighBIC_eigen');
model = ModelUtil.create('RayleighBIC_eigen');
model.label('Rounded 180-kHz Rayleigh BIC - eigenfrequency');
setParameters(model,P);
model.component.create('comp1',true);
model.component('comp1').geom.create('geom1',2);
buildGeometry(model,true);
model.component('comp1').geom('geom1').run;

% Use geometry-generated domain selections.  mphselectcoords searches for
% vertices (not arbitrary interior points), so it must not be used to select
% a domain from its centre coordinate.
a=P.a; H=P.H; tp=P.tPML; tol=1e-9;
fluidDom = namedEntities(model,'geom1_fluid_dom');
pmlDom = namedEntities(model,'geom1_pmlRect_dom');
steelDom = namedEntities(model,'geom1_steel_dom');
leftB = mphselectbox(model,'geom1',[-tol,tol;-tol,H+tp+tol], ...
    'boundary');
rightB = mphselectbox(model,'geom1',[a-tol,a+tol;-tol,H+tp+tol], ...
    'boundary');
fprintf('Eigen selections: fluid=%s, pml=%s, steel=%s, left=%s, right=%s\n', ...
    mat2str(fluidDom),mat2str(pmlDom),mat2str(steelDom), ...
    mat2str(leftB),mat2str(rightB));

addMaterials(model,unique([fluidDom,pmlDom]),steelDom);
model.component('comp1').physics.create('acpr','PressureAcoustics','geom1');
model.component('comp1').physics('acpr').selection.set(unique([fluidDom,pmlDom]));
model.component('comp1').physics('acpr').create('pc1','PeriodicCondition',1);
model.component('comp1').physics('acpr').feature('pc1').selection.set([leftB,rightB]);
model.component('comp1').physics('acpr').feature('pc1').set('PeriodicType','Floquet');
model.component('comp1').physics('acpr').feature('pc1') ...
    .set('kFloquet',{'kxB','0','0'});

model.component('comp1').coordSystem.create('pml1','PML');
model.component('comp1').coordSystem('pml1').selection.set(pmlDom);
model.component('comp1').coordSystem('pml1').set('ScalingType','Cartesian');
model.component('comp1').coordSystem('pml1').set('wavelengthSourceType','userDefined');
model.component('comp1').coordSystem('pml1').set('typicalWavelength','lambda0');

addMesh(model);
model.study.create('std1');
model.study('std1').create('eig','Eigenfrequency');
model.study('std1').feature('eig').set('neigs',8);
model.study('std1').feature('eig').set('shift','f0');
model.study('std1').setGenPlots(false);
model.study('std1').run;

freqs = mphglobal(model,'freq');
pL = mphinterp(model,'acpr.p_t','coord',[P.x1+P.w1/2;-0.65*P.d1]);
pS = mphinterp(model,'acpr.p_t','coord',[P.x2+P.w2/2;-0.55*P.d2]);
pE = mphinterp(model,'acpr.p_t','coord',[P.a/2;8e-3]);
score = (abs(pL)+abs(pS))./max(abs(pE),1e-15);
window = abs(real(freqs)-P.f0)<35e3;
score(~window) = -inf;
[~,modeIndex] = max(score);

model.result.create('pg1','PlotGroup2D');
model.result('pg1').label('Localized eigenpressure');
model.result('pg1').set('looplevel',modeIndex);
model.result('pg1').feature.create('surf1','Surface');
model.result('pg1').feature('surf1').set('expr','real(acpr.p_t)');
model.result('pg1').feature('surf1').set('colortable','Wave');

exportEigenField(model,P,modeIndex,outDir);
meta = struct('selected_mode',modeIndex, ...
    'eigenfrequency_real_hz',real(freqs(modeIndex)), ...
    'eigenfrequency_imag_hz',imag(freqs(modeIndex)), ...
    'Q_estimate',abs(real(freqs(modeIndex))/(2*imag(freqs(modeIndex)))), ...
    'localization_score',score(modeIndex), ...
    'kx_rad_per_m',2*pi*P.f0/P.c*sind(P.theta));
end


function [model,meta] = buildScatteringModel(P,outDir)
import com.comsol.model.*
import com.comsol.model.util.*
ModelUtil.remove('RayleighBIC_scattering');
model = ModelUtil.create('RayleighBIC_scattering');
model.label('Rounded 180-kHz Rayleigh BIC - periodic-port scattering');
setParameters(model,P);
model.component.create('comp1',true);
model.component('comp1').geom.create('geom1',2);
buildGeometry(model,false);
model.component('comp1').geom('geom1').run;

a=P.a; H=P.H; tol=1e-9;
fluidDom = namedEntities(model,'geom1_fluid_dom');
steelDom = namedEntities(model,'geom1_steel_dom');
leftB = mphselectbox(model,'geom1',[-tol,tol;-tol,H+tol], ...
    'boundary');
rightB = mphselectbox(model,'geom1',[a-tol,a+tol;-tol,H+tol], ...
    'boundary');
topB = mphselectbox(model,'geom1',[-tol,a+tol;H-tol,H+tol], ...
    'boundary');
fprintf('Scattering selections: fluid=%s, steel=%s, left=%s, right=%s, top=%s\n', ...
    mat2str(fluidDom),mat2str(steelDom),mat2str(leftB), ...
    mat2str(rightB),mat2str(topB));

addMaterials(model,fluidDom,steelDom);
model.component('comp1').physics.create('acpr','PressureAcoustics','geom1');
model.component('comp1').physics('acpr').selection.set(fluidDom);
model.component('comp1').physics('acpr').create('pport1','PeriodicPort',1);
model.component('comp1').physics('acpr').feature('pport1').selection.set(topB);
model.component('comp1').physics('acpr').feature('pport1').set('IncidentWave','PowerPerLength');
model.component('comp1').physics('acpr').feature('pport1').set('P0','1[W/m]');
model.component('comp1').physics('acpr').feature('pport1').set('thetai','theta');
for n = [-3,-2,-1,1,2,3]
    tag = sprintf('dport%s',strrep(sprintf('%+d',n),'-','m'));
    tag = strrep(tag,'+','p');
    model.component('comp1').physics('acpr').feature('pport1') ...
        .create(tag,'DiffractionOrderPort',1);
    model.component('comp1').physics('acpr').feature('pport1') ...
        .feature(tag).setIndex('m',n,0);
end
model.component('comp1').physics('acpr').create('pc1','PeriodicCondition',1);
model.component('comp1').physics('acpr').feature('pc1').selection.set([leftB,rightB]);
model.component('comp1').physics('acpr').feature('pc1').set('PeriodicType','Floquet');
model.component('comp1').physics('acpr').feature('pc1') ...
    .set('kFloquet_src','root.comp1.acpr.pport1.kitx');

addMesh(model);
model.study.create('std1');
model.study('std1').create('freq','Frequency');
model.study('std1').feature('freq').set('plist', ...
    'range(179.60[kHz],0.02[kHz],180.40[kHz])');
model.study('std1').setGenPlots(false);
model.study('std1').run;

freqs = mphglobal(model,'freq');
pProbe = mphinterp(model,'acpr.p_t','coord', ...
    [P.x1+P.w1/2;-0.65*P.d1]);
[~,fieldIndex] = max(abs(pProbe));
model.result.create('pg1','PlotGroup2D');
model.result('pg1').label('Driven total pressure');
model.result('pg1').set('looplevel',fieldIndex);
model.result('pg1').feature.create('surf1','Surface');
model.result('pg1').feature('surf1').set('expr','real(acpr.p_t)');
model.result('pg1').feature('surf1').set('colortable','Wave');

exportScatteringField(model,P,fieldIndex,outDir);
meta = struct('selected_solution',fieldIndex, ...
    'selected_frequency_hz',real(freqs(fieldIndex)), ...
    'probe_pressure_abs_pa',abs(pProbe(fieldIndex)), ...
    'rayleigh_frequency_hz',P.c/(P.a*(1+sind(P.theta))), ...
    'theta_deg',P.theta);
end


function setParameters(model,P)
model.param.set('a',sprintf('%.12g[m]',P.a));
model.param.set('w1',sprintf('%.12g[m]',P.w1));
model.param.set('d1',sprintf('%.12g[m]',P.d1));
model.param.set('w2',sprintf('%.12g[m]',P.w2));
model.param.set('d2',sprintf('%.12g[m]',P.d2));
model.param.set('g12',sprintf('%.12g[m]',P.g12));
model.param.set('g21',sprintf('%.12g[m]',P.g21));
model.param.set('x1','g21/2');
model.param.set('x2','x1+w1+g12');
model.param.set('tPlate',sprintf('%.12g[m]',P.tPlate));
model.param.set('H',sprintf('%.12g[m]',P.H));
model.param.set('tPML',sprintf('%.12g[m]',P.tPML));
model.param.set('c0',sprintf('%.12g[m/s]',P.c));
model.param.set('rho0',sprintf('%.12g[kg/m^3]',P.rho));
model.param.set('f0',sprintf('%.12g[Hz]',P.f0));
model.param.set('theta',sprintf('%.12g[deg]',P.theta));
model.param.set('lambda0','c0/f0');
model.param.set('kxB','2*pi*f0/c0*sin(theta)');
end


function buildGeometry(model,withPML)
g = model.component('comp1').geom('geom1');
g.lengthUnit('m');
g.feature.create('plate','Rectangle');
g.feature('plate').set('size',{'a','tPlate'});
g.feature('plate').set('pos',{'0','-tPlate'});
g.feature.create('groove1','Rectangle');
g.feature('groove1').set('size',{'w1','d1'});
g.feature('groove1').set('pos',{'x1','-d1'});
g.feature.create('groove2','Rectangle');
g.feature('groove2').set('size',{'w2','d2'});
g.feature('groove2').set('pos',{'x2','-d2'});
g.feature.create('steel','Difference');
g.feature('steel').selection('input').set({'plate'});
g.feature('steel').selection('input2').set({'groove1','groove2'});
g.feature('steel').set('selresult',true);

g.feature.create('waterTop','Rectangle');
g.feature('waterTop').set('size',{'a','H'});
g.feature('waterTop').set('pos',{'0','0'});
g.feature.create('waterG1','Rectangle');
g.feature('waterG1').set('size',{'w1','d1'});
g.feature('waterG1').set('pos',{'x1','-d1'});
g.feature.create('waterG2','Rectangle');
g.feature('waterG2').set('size',{'w2','d2'});
g.feature('waterG2').set('pos',{'x2','-d2'});
g.feature.create('fluid','Union');
g.feature('fluid').selection('input').set({'waterTop','waterG1','waterG2'});
g.feature('fluid').set('intbnd',false);
g.feature('fluid').set('selresult',true);

inputs = {'steel','fluid'};
if withPML
    g.feature.create('pmlRect','Rectangle');
    g.feature('pmlRect').set('size',{'a','tPML'});
    g.feature('pmlRect').set('pos',{'0','H'});
    g.feature('pmlRect').set('selresult',true);
    inputs = {'steel','fluid','pmlRect'};
end
g.feature.create('all','Union');
g.feature('all').selection('input').set(inputs);
g.feature('all').set('intbnd',true);
end


function addMaterials(model,waterDomains,steelDomains)
model.component('comp1').material.create('water','Common');
model.component('comp1').material('water').label('Water');
model.component('comp1').material('water').selection.set(waterDomains);
model.component('comp1').material('water').propertyGroup('def').set('density','rho0');
model.component('comp1').material('water').propertyGroup('def').set('soundspeed','c0');
model.component('comp1').material.create('steelmat','Common');
model.component('comp1').material('steelmat').label('Stainless steel');
model.component('comp1').material('steelmat').selection.set(steelDomains);
model.component('comp1').material('steelmat').propertyGroup('def').set('density','7850[kg/m^3]');
model.component('comp1').material('steelmat').propertyGroup('def').set('youngsmodulus','200[GPa]');
model.component('comp1').material('steelmat').propertyGroup('def').set('poissonsratio','0.30');
end


function addMesh(model)
model.component('comp1').mesh.create('mesh1');
model.component('comp1').mesh('mesh1').feature('size').set('custom','on');
model.component('comp1').mesh('mesh1').feature('size').set('hmax','0.32[mm]');
model.component('comp1').mesh('mesh1').feature('size').set('hmin','0.025[mm]');
model.component('comp1').mesh('mesh1').feature('size').set('hgrad',1.25);
model.component('comp1').mesh('mesh1').feature('size').set('hcurve',0.25);
model.component('comp1').mesh('mesh1').feature('size').set('hnarrow',0.8);
model.component('comp1').mesh('mesh1').create('ftri1','FreeTri');
model.component('comp1').mesh('mesh1').run;
end


function exportEigenField(model,P,solnum,outDir)
x = linspace(0,P.a,321);
y = linspace(-P.tPlate,P.H,501);
[X,Y] = meshgrid(x,y);
mask = Y>=0 | ...
    (X>=P.x1 & X<=P.x1+P.w1 & Y>=-P.d1) | ...
    (X>=P.x2 & X<=P.x2+P.w2 & Y>=-P.d2);
coord = [X(mask).';Y(mask).'];
val = mphinterp(model,'acpr.p_t','coord',coord,'solnum',solnum);
field = complex(nan(size(X))); field(mask)=val;
scale = max(abs(field(mask))); field = field/max(scale,eps);
writematrix(x(:),fullfile(outDir,'eigen_field_x_m.csv'));
writematrix(y(:),fullfile(outDir,'eigen_field_y_m.csv'));
writematrix(real(field),fullfile(outDir,'eigen_field_real.csv'));
writematrix(imag(field),fullfile(outDir,'eigen_field_imag.csv'));
end


function exportScatteringField(model,P,solnum,outDir)
x = linspace(0,P.a,321);
y = linspace(-P.tPlate,P.H,501);
[X,Y] = meshgrid(x,y);
mask = Y>=0 | ...
    (X>=P.x1 & X<=P.x1+P.w1 & Y>=-P.d1) | ...
    (X>=P.x2 & X<=P.x2+P.w2 & Y>=-P.d2);
coord = [X(mask).';Y(mask).'];
val = mphinterp(model,'acpr.p_t','coord',coord,'solnum',solnum);
field = complex(nan(size(X))); field(mask)=val;
writematrix(x(:),fullfile(outDir,'scattering_field_x_m.csv'));
writematrix(y(:),fullfile(outDir,'scattering_field_y_m.csv'));
writematrix(real(field),fullfile(outDir,'scattering_field_real.csv'));
writematrix(imag(field),fullfile(outDir,'scattering_field_imag.csv'));
end


function ids = namedEntities(model,tag)
info = mphgetselection(model,tag);
ids = reshape(info.entities,1,[]);
if isempty(ids)
    error('Generated geometry selection %s is empty.',tag);
end
end
