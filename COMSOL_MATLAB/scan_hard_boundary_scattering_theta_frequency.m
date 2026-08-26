%% COMSOL numerical scattering experiment for the rigid rounded sample
% Direct periodic-port calculations versus incidence angle and frequency.
% All diffraction efficiencies are port powers divided by incident power.

clear; clc;
import com.comsol.model.*
import com.comsol.model.util.*

rootDir=fileparts(fileparts(mfilename('fullpath')));
inFile=fullfile(rootDir,'COMSOL_MATLAB','rounded_180k_results', ...
    'RayleighBIC_rounded_180k_scattering.mph');
outDir=fullfile(rootDir,'COMSOL_MATLAB','rounded_180k_results', ...
    'hard_boundary_scattering_scan');
if ~exist(outDir,'dir'), mkdir(outDir); end
model=mphload(inFile,'RigidScatteringScan');

a=model.param.evaluate('a'); c0=model.param.evaluate('c0');
x1=model.param.evaluate('x1'); w1=model.param.evaluate('w1');
d1=model.param.evaluate('d1');

%% Coarse two-dimensional experimental map
thetaMap=(6.20:0.05:8.20).';
frequencyMap=(179.40e3:20:180.60e3).';
nT=numel(thetaMap); nF=numel(frequencyMap);
eta0=nan(nT,nF); etaM1=nan(nT,nF); etaP1=nan(nT,nF);
etaTotal=nan(nT,nF); knM1=complex(nan(nT,nF));
probeAbs=nan(nT,nF);
plistCoarse='range(179.40[kHz],0.020[kHz],180.60[kHz])';

fprintf('COMSOL coarse map: %d angles x %d frequencies.\n',nT,nF);
tic;
for it=1:nT
    model.param.set('theta',sprintf('%.12g[deg]',thetaMap(it)));
    model.study('std1').feature('freq').set('plist',plistCoarse);
    model.study('std1').run;
    f=real(mphglobal(model,'freq'));
    assert(numel(f)==nF && max(abs(f(:)-frequencyMap))<1e-4);
    Pin=real(mphglobal(model,'acpr.pport1.P_in'));
    eta0(it,:)=reshape(real(mphglobal(model, ...
        'acpr.pport1.P_out'))./Pin,1,[]);
    etaM1(it,:)=reshape(real(mphglobal(model, ...
        'acpr.pport1.dportm1.P_out'))./Pin,1,[]);
    etaP1(it,:)=reshape(real(mphglobal(model, ...
        'acpr.pport1.dportp1.P_out'))./Pin,1,[]);
    etaTotal(it,:)=reshape(real(mphglobal(model, ...
        'acpr.pport1.P_out_tot'))./Pin,1,[]);
    knM1(it,:)=reshape(mphglobal(model, ...
        'acpr.pport1.dportm1.kn'),1,[]);
    p=mphinterp(model,'acpr.p_t','coord',[x1+w1/2;-0.65*d1]);
    probeAbs(it,:)=reshape(abs(p),1,[]);
    fprintf('  angle %2d/%2d: %.3f deg\n',it,nT,thetaMap(it));
end
elapsedCoarse=toc;
rayleighFrequency=c0./(a*(1+sind(thetaMap)));
[IT,JF]=ndgrid(1:nT,1:nF);
coarseTable=table(thetaMap(IT(:)),frequencyMap(JF(:)), ...
    rayleighFrequency(IT(:)),eta0(:),etaM1(:),etaP1(:), ...
    etaTotal(:),real(knM1(:)),imag(knM1(:)),probeAbs(:), ...
    'VariableNames',{'theta_deg','frequency_hz','rayleigh_frequency_hz', ...
    'eta_0','eta_m1','eta_p1','eta_total','kn_m1_real','kn_m1_imag', ...
    'large_groove_probe_abs_pa'});
writetable(coarseTable,fullfile(outDir,'comsol_coarse_angle_frequency_map.csv'));

%% Sub-hertz frequency sweeps at experimentally selected angles
cutAngles=[6.70,6.90,7.10,7.30,7.50].';
cutRows=cell(numel(cutAngles),1);
fprintf('COMSOL fine frequency cuts: %d angles.\n',numel(cutAngles));
tic;
for ia=1:numel(cutAngles)
    theta=cutAngles(ia);
    fRA=c0/(a*(1+sind(theta)));
    % Exclude exactly zero normal wave number; a finite sample would
    % regularize that single mathematical threshold point.
    plist=sprintf(['range(%.12g[Hz],0.5[Hz],%.12g[Hz]) ', ...
        'range(%.12g[Hz],0.5[Hz],%.12g[Hz])'], ...
        fRA-250,fRA-0.5,fRA+0.5,fRA+250);
    model.param.set('theta',sprintf('%.12g[deg]',theta));
    model.study('std1').feature('freq').set('plist',plist);
    model.study('std1').run;
    f=real(mphglobal(model,'freq')); f=f(:);
    Pin=real(mphglobal(model,'acpr.pport1.P_in')); Pin=Pin(:);
    p0=real(mphglobal(model,'acpr.pport1.P_out'))./Pin;
    pm=real(mphglobal(model,'acpr.pport1.dportm1.P_out'))./Pin;
    pp=real(mphglobal(model,'acpr.pport1.dportp1.P_out'))./Pin;
    pt=real(mphglobal(model,'acpr.pport1.P_out_tot'))./Pin;
    kn=mphglobal(model,'acpr.pport1.dportm1.kn'); kn=kn(:);
    probe=mphinterp(model,'acpr.p_t','coord',[x1+w1/2;-0.65*d1]);
    probe=probe(:);
    n=numel(f);
    cutRows{ia}=table(repmat(theta,n,1),repmat(fRA,n,1),f,f-fRA, ...
        p0(:),pm(:),pp(:),pt(:),real(kn),imag(kn),real(probe), ...
        imag(probe),abs(probe),angle(probe), ...
        'VariableNames',{'theta_deg','rayleigh_frequency_hz', ...
        'frequency_hz','detuning_from_rayleigh_hz','eta_0','eta_m1', ...
        'eta_p1','eta_total','kn_m1_real','kn_m1_imag', ...
        'probe_real_pa','probe_imag_pa','probe_abs_pa','probe_phase_rad'});
    fprintf('  fine cut %d/%d: %.3f deg, %d frequencies\n', ...
        ia,numel(cutAngles),theta,n);
end
elapsedFine=toc;
fineTable=vertcat(cutRows{:});
writetable(fineTable,fullfile(outDir,'comsol_fine_frequency_cuts.csv'));

%% Three COMSOL field states around the strongest theta=7.10 deg response
target=fineTable(abs(fineTable.theta_deg-7.10)<1e-10,:);
[~,peakIndex]=max(target.probe_abs_pa);
peakFrequency=target.frequency_hz(peakIndex);
fieldFrequency=peakFrequency+[-5,0,5];
model.param.set('theta','7.10[deg]');
model.study('std1').feature('freq').set('plist',sprintf( ...
    '%.12g[Hz] %.12g[Hz] %.12g[Hz]',fieldFrequency));
model.study('std1').run;

x=linspace(0,a,321); y=linspace(-10e-3,25e-3,501);
[X,Y]=meshgrid(x,y);
mask=Y>=0 | ...
    (X>=model.param.evaluate('x1') & X<=model.param.evaluate('x1+w1') & ...
    Y>=-model.param.evaluate('d1')) | ...
    (X>=model.param.evaluate('x2') & X<=model.param.evaluate('x2+w2') & ...
    Y>=-model.param.evaluate('d2'));
coord=[X(mask).';Y(mask).'];
fieldRows=cell(3,1);
for state=1:3
    p=mphinterp(model,'acpr.p_t','coord',coord,'solnum',state);
    field=complex(nan(size(X))); field(mask)=p;
    writematrix(real(field),fullfile(outDir,sprintf( ...
        'field_state_%d_real_pa.csv',state)));
    writematrix(imag(field),fullfile(outDir,sprintf( ...
        'field_state_%d_imag_pa.csv',state)));
    Pin=real(mphglobal(model,'acpr.pport1.P_in','solnum',state));
    p0=real(mphglobal(model,'acpr.pport1.P_out','solnum',state))/Pin;
    pm=real(mphglobal(model,'acpr.pport1.dportm1.P_out','solnum',state))/Pin;
    pt=real(mphglobal(model,'acpr.pport1.P_out_tot','solnum',state))/Pin;
    fieldRows{state}=table(state,fieldFrequency(state),p0,pm,pt, ...
        max(abs(field(:)),[],'omitnan'),'VariableNames', ...
        {'state','frequency_hz','eta_0','eta_m1','eta_total', ...
        'max_pressure_pa'});
end
writematrix(x(:),fullfile(outDir,'field_x_m.csv'));
writematrix(y(:),fullfile(outDir,'field_y_m.csv'));
fieldTable=vertcat(fieldRows{:});
writetable(fieldTable,fullfile(outDir,'comsol_field_states.csv'));

model.label('Rigid rounded 180-kHz sample - theta/f scattering scan');
mphsave(model,fullfile(outDir, ...
    'RayleighBIC_rounded_180k_hard_boundary_scattering_scan.mph'));
save(fullfile(outDir,'comsol_hard_boundary_scattering_scan.mat'), ...
    'thetaMap','frequencyMap','rayleighFrequency','eta0','etaM1', ...
    'etaP1','etaTotal','knM1','probeAbs','fineTable','fieldTable', ...
    'elapsedCoarse','elapsedFine','-v7.3');

fprintf('COMSOL scattering scan completed: coarse %.2f s, fine %.2f s.\n', ...
    elapsedCoarse,elapsedFine);
fprintf('theta=7.10 deg strongest sampled field: %.9f kHz.\n', ...
    peakFrequency/1e3);
disp(fieldTable);
