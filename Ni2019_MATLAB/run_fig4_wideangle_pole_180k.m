%% Wide-angle physical eigenpole continuation for Fig. 4
% Resolves the resonance frequency and linewidth from approximately 0 to
% 20 degrees.  The dense logarithmic stencil around the strict BIC prevents
% the vanishing linewidth from being hidden by a uniform angular grid.
clear; clc;
rootDir=fileparts(mfilename('fullpath'));
outDir=fullfile(rootDir,'results','fig4_wideangle_180k');
if ~exist(outDir,'dir'), mkdir(outDir); end
addpath(rootDir);

S=load(fullfile(rootDir,'results','StrictRayleighBIC_180kHz_7p10deg_final.mat'));
cfg=S.cfg; cfg.solve_scattering=false;
kB=S.x(1); OB=1-kB; c0=1500; a=S.a;

near=logspace(-7,-3,31).';
leftOuter=linspace(0.002,kB-1e-4,57).';
rightOuter=linspace(0.002,0.225,113).';
leftK=kB-unique([near;leftOuter]); leftK=sort(leftK,'descend');
rightK=kB+unique([near;rightOuter]); rightK=sort(rightK,'ascend');

fprintf('Wide-angle pole continuation: %d left + endpoint + %d right.\n', ...
    numel(leftK),numel(rightK));
left=continue_side(leftK,complex(OB,1e-11),cfg,c0,a);
right=continue_side(rightK,complex(OB,1e-11),cfg,c0,a);

endpoint=table(kB,OB,0,asind(kB/OB),OB*c0/a,0,Inf, ...
    'VariableNames',{'kappa','Omega_real','Omega_imag','theta_deg', ...
    'frequency_hz','linewidth_hz','Q'});
data=[flipud(left);endpoint;right];
data=sortrows(data,'theta_deg');
data=data(data.theta_deg>=0 & data.theta_deg<=20.05,:);
writetable(data,fullfile(outDir,'wideangle_physical_pole.csv'));
save(fullfile(outDir,'wideangle_physical_pole.mat'),'data','S','cfg','-v7.3');

fprintf('Saved %d points: theta %.6f to %.6f deg.\n',height(data), ...
    min(data.theta_deg),max(data.theta_deg));
for target=[0 4 6 7 7.0996629 8 10 15 20]
    [~,j]=min(abs(data.theta_deg-target));
    fprintf('theta=%8.5f deg f=%10.6f kHz linewidth=%12.5g Hz Q=%12.5g\n', ...
        data.theta_deg(j),data.frequency_hz(j)/1e3,data.linewidth_hz(j),data.Q(j));
end

function T=continue_side(kList,seed,cfg,c0,a)
rows=cell(numel(kList),1);
for j=1:numel(kList)
    p=ni2019_refine_outgoing_pole_kappa(cfg,kList(j),seed, ...
        'OuterIterations',12,'Display','off', ...
        'FunctionTolerance',5e-12,'StepTolerance',5e-12);
    seed=p.Omega;
    theta=asind(kList(j)/real(p.Omega));
    linewidth=2*abs(imag(p.Omega))*c0/a;
    q=real(p.Omega)/max(2*abs(imag(p.Omega)),eps);
    rows{j}=table(kList(j),real(p.Omega),imag(p.Omega),theta, ...
        real(p.Omega)*c0/a,linewidth,q,'VariableNames', ...
        {'kappa','Omega_real','Omega_imag','theta_deg', ...
        'frequency_hz','linewidth_hz','Q'});
    if mod(j,10)==0 || j==numel(kList)
        fprintf('  %d/%d: theta %.4f deg, linewidth %.4g Hz\n', ...
            j,numel(kList),theta,linewidth);
    end
end
T=vertcat(rows{:});
end
