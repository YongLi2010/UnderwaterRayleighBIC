%% Compare rigid-wall and finite-elastic stainless-steel eigenmodes
% The two coupled limits differ only in the lower-face mechanical condition:
% (1) traction-free lower face and (2) fully fixed lower face.

clear; clc;
import com.comsol.model.*
import com.comsol.model.util.*

rootDir = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(rootDir,'COMSOL_MATLAB','rounded_180k_results');
baseFile = fullfile(outDir,'RayleighBIC_rounded_180k_eigenfrequency.mph');

rigid = mphload(baseFile,'RigidReference');
fRigidAll = mphglobal(rigid,'freq');
fRigid = fRigidAll(4);
[sampleCoord,pRigid] = referenceSamples(rigid,4);

fprintf('Rigid reference: %.9f %+.9fi kHz\n', ...
    real(fRigid)/1e3,imag(fRigid)/1e3);

[freeModel,freeResult] = solveCoupled(baseFile,'CoupledFree',false, ...
    sampleCoord,pRigid,real(fRigid));
mphsave(freeModel,fullfile(outDir, ...
    'RayleighBIC_rounded_180k_coupled_free.mph'));

[fixedModel,fixedResult] = solveCoupled(baseFile,'CoupledFixed',true, ...
    sampleCoord,pRigid,real(fRigid));
mphsave(fixedModel,fullfile(outDir, ...
    'RayleighBIC_rounded_180k_coupled_fixed.mph'));

rigidResult = struct( ...
    'case_name',"rigid_water_only", ...
    'frequency_real_hz',real(fRigid), ...
    'frequency_imag_hz',imag(fRigid), ...
    'Q_estimate',abs(real(fRigid)/(2*imag(fRigid))), ...
    'shift_from_rigid_hz',0, ...
    'relative_shift_percent',0, ...
    'pressure_overlap_with_rigid',1);

T = struct2table([rigidResult,freeResult,fixedResult]);
writetable(T,fullfile(outDir,'rigid_vs_coupled_eigen_comparison.csv'));
disp(T);


function [model,result] = solveCoupled(baseFile,tag,fixBottom,coord,pRef,fRef)
model = mphload(baseFile,tag);
g = 'geom1'; c = 'comp1'; tol = 1e-9;
a = model.param.evaluate('a');
tPlate = model.param.evaluate('tPlate');

steelInfo = mphgetselection(model,'geom1_steel_dom');
fluidInfo = mphgetselection(model,'geom1_fluid_dom');
steelDom = reshape(steelInfo.entities,1,[]);
fluidDom = reshape(fluidInfo.entities,1,[]);
steelB = mphgetadj(model,g,'boundary','domain',steelDom);
fluidB = mphgetadj(model,g,'boundary','domain',fluidDom);
interfaceB = intersect(steelB,fluidB);
leftSteelB = mphselectbox(model,g,[-tol,tol;-tPlate-tol,tol], ...
    'boundary','adjnumber',steelDom);
rightSteelB = mphselectbox(model,g,[a-tol,a+tol;-tPlate-tol,tol], ...
    'boundary','adjnumber',steelDom);
bottomB = mphselectbox(model,g,[-tol,a+tol;-tPlate-tol,-tPlate+tol], ...
    'boundary','adjnumber',steelDom);

fprintf('%s: steel=%s, interface=%s, left=%s, right=%s, bottom=%s\n', ...
    tag,mat2str(steelDom),mat2str(interfaceB),mat2str(leftSteelB), ...
    mat2str(rightSteelB),mat2str(bottomB));

model.component(c).physics.create('solid','SolidMechanics',g);
model.component(c).physics('solid').selection.set(steelDom);
model.component(c).physics('solid').create('pc1','PeriodicCondition',1);
model.component(c).physics('solid').feature('pc1').selection.set( ...
    [leftSteelB,rightSteelB]);
model.component(c).physics('solid').feature('pc1').set( ...
    'PeriodicType','Floquet');
model.component(c).physics('solid').feature('pc1').set( ...
    'kFloquet',{'kxB','0','0'});

if fixBottom
    model.component(c).physics('solid').create('fix1','Fixed',1);
    model.component(c).physics('solid').feature('fix1').selection.set(bottomB);
end

model.component(c).multiphysics.create('asb1', ...
    'AcousticStructureBoundary',1);
model.component(c).multiphysics('asb1').set('Acoustics_physics','acpr');
model.component(c).multiphysics('asb1').set('Structure_physics','solid');
model.component(c).multiphysics('asb1').selection.set(interfaceB);

model.study('std1').feature('eig').set('neigs',20);
model.study('std1').feature('eig').set('shift','f0');
model.study('std1').run;

freq = mphglobal(model,'freq');
nmode = numel(freq);
overlap = zeros(nmode,1);
for j = 1:nmode
    pj = mphinterp(model,'acpr.p_t','coord',coord,'solnum',j);
    valid = isfinite(pj(:)) & isfinite(pRef(:));
    v = pj(valid); r = pRef(valid);
    denom = norm(v)*norm(r);
    if denom > 0
        overlap(j) = abs(r(:)'*v(:))/denom;
    end
end

window = abs(real(freq)-fRef)<35e3;
overlap(~window) = -inf;
[bestOverlap,idx] = max(overlap);
f = freq(idx);
if fixBottom
    caseName = "coupled_fixed_bottom";
else
    caseName = "coupled_free_bottom";
end
result = struct( ...
    'case_name',caseName, ...
    'frequency_real_hz',real(f), ...
    'frequency_imag_hz',imag(f), ...
    'Q_estimate',abs(real(f)/(2*imag(f))), ...
    'shift_from_rigid_hz',real(f)-fRef, ...
    'relative_shift_percent',100*(real(f)-fRef)/fRef, ...
    'pressure_overlap_with_rigid',bestOverlap);
fprintf('%s selected mode %d: %.9f %+.9fi kHz, Q=%.6g, overlap=%.6f\n', ...
    tag,idx,real(f)/1e3,imag(f)/1e3,result.Q_estimate,bestOverlap);
end


function [coord,p] = referenceSamples(model,solnum)
% Sample the large groove, small groove, and near-surface exterior.  Modal
% overlap on this shared physical grid tracks the acoustic branch through
% the dense set of elastic plate modes without relying on eigenvalue order.
x1 = model.param.evaluate('x1');
x2 = model.param.evaluate('x2');
w1 = model.param.evaluate('w1');
w2 = model.param.evaluate('w2');
d1 = model.param.evaluate('d1');
d2 = model.param.evaluate('d2');
a = model.param.evaluate('a');

[X1,Y1] = meshgrid(linspace(x1+0.02*w1,x1+0.98*w1,36), ...
    linspace(-0.98*d1,-0.02*d1,30));
[X2,Y2] = meshgrid(linspace(x2+0.04*w2,x2+0.96*w2,14), ...
    linspace(-0.96*d2,-0.04*d2,10));
[XE,YE] = meshgrid(linspace(0.01*a,0.99*a,40),linspace(0.05e-3,6e-3,18));
coord = [X1(:).',X2(:).',XE(:).';Y1(:).',Y2(:).',YE(:).'];
p = mphinterp(model,'acpr.p_t','coord',coord,'solnum',solnum);
end
