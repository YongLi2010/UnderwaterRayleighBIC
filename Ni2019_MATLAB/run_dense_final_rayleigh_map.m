%% Dense driven-scattering map around the audited 200 kHz Rayleigh BIC
% This script only recomputes the local angle-frequency atlas.  It does not
% run the legacy gallery code, whose router-convergence constants are not
% part of the audited final-root data set.
clear; clc;

rootDir=fileparts(mfilename('fullpath'));
resultDir=fullfile(rootDir,'results');
paperDir=fullfile(fileparts(rootDir),'arxiv_theory_paper');
dataDir=fullfile(paperDir,'figure_data');
if ~exist(dataDir,'dir'), mkdir(dataDir); end

D=load(fullfile(resultDir,'StrictRayleighBIC_200kHz_min1mm.mat'));
x=D.xFinal;
aPhysical=D.aPhysical;
cWater=D.cWater;

% Focused PRL map: 0.0075 deg and 9.375 Hz increments.
thetaValues=linspace(4.4,6.2,241);
frequencyValues=linspace(198.5e3,201.5e3,321);
nTheta=numel(thetaValues);
nFrequency=numel(frequencyValues);
nPoint=nTheta*nFrequency;
mapN=61;
mapK=9;

absAFlat=nan(nPoint,3);
phaseAFlat=nan(nPoint,3);
etaFlat=nan(nPoint,3);
grooveNormFlat=nan(nPoint,1);
conditionFlat=nan(nPoint,1);

fprintf('Dense final-geometry scan: %d x %d = %d points (N=%d,K=%d)\n', ...
    nTheta,nFrequency,nPoint,mapN,mapK);
tic;
parfor index=1:nPoint
    [it,jf]=ind2sub([nTheta,nFrequency],index);
    R=solve_physical_dense(thetaValues(it),frequencyValues(jf),x, ...
        aPhysical,cWater,mapN,mapK);
    ids=order_ids_dense(R);
    absAFlat(index,:)=abs(R.A(ids));
    phaseAFlat(index,:)=angle(R.A(ids));
    etaFlat(index,:)=R.eta(ids);
    grooveNormFlat(index)=norm(R.groove_surface_coefficients(:));
    conditionFlat(index)=R.condition_number;
end
elapsedSeconds=toc;

absA=reshape(absAFlat,[nTheta,nFrequency,3]);
phaseA=reshape(phaseAFlat,[nTheta,nFrequency,3]);
eta=reshape(etaFlat,[nTheta,nFrequency,3]);
grooveNorm=reshape(grooveNormFlat,[nTheta,nFrequency]);
conditionNumber=reshape(conditionFlat,[nTheta,nFrequency]);

kappa=D.xFinal(1);
Omega=1-kappa;
thetaBIC=asind(kappa/Omega);
frequencyBIC=D.fTarget;
rayleighFrequency=cWater./(aPhysical*(1+sind(thetaValues)));

outFile=fullfile(resultDir,'DenseFinalRayleighMap_241x321.mat');
save(outFile,'thetaValues','frequencyValues','absA','phaseA','eta', ...
    'grooveNorm','conditionNumber','rayleighFrequency','thetaBIC', ...
    'frequencyBIC','mapN','mapK','elapsedSeconds','x','aPhysical','cWater','-v7.3');

% Preserve the orientation expected by the Python PRL figure script:
% rows are angle samples and columns are frequency samples.
writematrix(thetaValues(:),fullfile(dataDir,'fig2_map_theta_deg.csv'));
writematrix(frequencyValues(:),fullfile(dataDir,'fig2_map_frequency_hz.csv'));
writematrix(absA(:,:,1),fullfile(dataDir,'fig2_map_absAm1.csv'));
writematrix(absA(:,:,2),fullfile(dataDir,'fig2_map_absA0.csv'));
writematrix(absA(:,:,3),fullfile(dataDir,'fig2_map_absAp1.csv'));
writematrix(eta(:,:,1),fullfile(dataDir,'fig2_map_eta_m1.csv'));
writematrix(eta(:,:,2),fullfile(dataDir,'fig2_map_eta_0.csv'));
writematrix(eta(:,:,3),fullfile(dataDir,'fig2_map_eta_p1.csv'));
writematrix(grooveNorm,fullfile(dataDir,'fig2_map_groove_norm.csv'));
writematrix(conditionNumber,fullfile(dataDir,'fig2_map_condition_number.csv'));

fprintf('Completed in %.2f s. Saved %s\n',elapsedSeconds,outFile);
fprintf('Grid increments: dtheta=%.6f deg, df=%.3f Hz\n', ...
    thetaValues(2)-thetaValues(1),frequencyValues(2)-frequencyValues(1));

function R=solve_physical_dense(theta,f,x,aPhysical,c,N,K)
Omega=f*aPhysical/c;
cfg=struct('a',1,'lambda',1/Omega,'theta_i_deg',theta, ...
    'depths',x(2:3),'widths',x(4:5),'gaps',x(6), ...
    'N',N,'K',K,'solve_scattering',true);
R=ni2019_modal_solver(cfg);
end

function ids=order_ids_dense(R)
ids=[find(R.orders==-1,1),find(R.orders==0,1),find(R.orders==1,1)];
end
