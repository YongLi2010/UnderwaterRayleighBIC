%% Independent multistart pole census on the positive half Brillouin zone
% This avoids silently carrying a pole across a Rayleigh square-root branch
% point.  Reciprocity supplies the negative half only after independent
% spot checks at negative kappa.
clear; clc;
rootDir=fileparts(mfilename('fullpath'));
outDir=fullfile(rootDir,'results','full_complex_band_180k');
S=load(fullfile(rootDir,'results','StrictRayleighBIC_180kHz_7p10deg_final.mat'));
G=load(fullfile(outDir,'gamma_discovery.mat'));
cfg=S.cfg; cfg.N=89; cfg.K=11; cfg.solve_scattering=false;
kappa=(0:0.01:0.5).';
rootsByK=cell(size(kappa)); rootsByK{1}=G.roots;
fixedReal=0.055:0.10:1.155;
fixedImag=[0.002,0.035,0.16];

fprintf('Independent positive-k pole census: %d k-points.\n',numel(kappa));
for m=2:numel(kappa)
    previous=rootsByK{m-1};
    previousSeeds=complex([]);
    if ~isempty(previous), previousSeeds=[previous.Omega]; end
    rootsByK{m}=ni2019_discover_outgoing_poles_kappa(cfg,kappa(m), ...
        'OmegaRange',[0.015,1.205],'RealSeeds',fixedReal, ...
        'ImagSeeds',fixedImag,'ComplexSeeds',previousSeeds, ...
        'MergeTolerance',5e-4,'Verbose',false);
    fprintf('  %3d/%3d  kappa=%.2f  roots=%d\n',m,numel(kappa), ...
        kappa(m),numel(rootsByK{m}));
    save(fullfile(outDir,'positive_census_checkpoint.mat'), ...
        'rootsByK','kappa','cfg','fixedReal','fixedImag');
end

% Independent negative-k spot checks; these are not mirrored inputs.
kappaCheck=[-0.11,-0.25,-0.40,-0.50];
negativeCheck=cell(size(kappaCheck));
for m=1:numel(kappaCheck)
    positiveSeeds=complex([]);
    [~,id]=min(abs(kappa-abs(kappaCheck(m))));
    if ~isempty(rootsByK{id}), positiveSeeds=[rootsByK{id}.Omega]; end
    negativeCheck{m}=ni2019_discover_outgoing_poles_kappa(cfg,kappaCheck(m), ...
        'OmegaRange',[0.015,1.205],'RealSeeds',fixedReal, ...
        'ImagSeeds',fixedImag,'ComplexSeeds',positiveSeeds, ...
        'MergeTolerance',5e-4,'Verbose',false);
end

tracks=connect_roots(kappa,rootsByK);
kappaBIC=S.x(1); OmegaBIC=1-kappaBIC;
save(fullfile(outDir,'positive_pole_census.mat'),'rootsByK','kappa','tracks', ...
    'negativeCheck','kappaCheck','cfg','fixedReal','fixedImag', ...
    'kappaBIC','OmegaBIC','-v7.3');

fprintf('\nCensus complete: %d connected track segments.\n',numel(tracks));
for j=1:numel(tracks)
    fprintf('  track %2d: %2d points, kappa %.2f--%.2f, ReOmega %.3f--%.3f\n', ...
        j,numel(tracks(j).kappa),min(tracks(j).kappa),max(tracks(j).kappa), ...
        min(real(tracks(j).Omega)),max(real(tracks(j).Omega)));
end

% Assignment-independent reciprocity spot-check errors.
for m=1:numel(kappaCheck)
    [~,id]=min(abs(kappa-abs(kappaCheck(m))));
    po=sort([rootsByK{id}.Omega]); pn=sort([negativeCheck{m}.Omega]);
    n=min(numel(po),numel(pn));
    if n==0, err=NaN; else, err=max(abs(po(1:n)-pn(1:n))); end
    fprintf('  reciprocity kappa=%+.2f: n+/%d n-/%d error %.3e\n', ...
        kappaCheck(m),numel(po),numel(pn),err);
end

function tracks=connect_roots(kappa,rootsByK)
tracks=struct('kappa',{},'Omega',{},'Q',{},'sigma_ratio',{},'raw_residual',{});
active=[];
for m=1:numel(kappa)
    roots=rootsByK{m}; nR=numel(roots);
    if nR==0, active=[]; continue; end
    om=[roots.Omega];
    if isempty(active)
        for r=1:nR, [tracks,active]=start_track(tracks,active,kappa(m),roots(r)); end
        continue;
    end
    pairs=[];
    for ia=1:numel(active)
        last=tracks(active(ia)).Omega(end);
        for r=1:nR
            pairs(end+1,:)=[abs(last-om(r)),ia,r]; %#ok<AGROW>
        end
    end
    pairs=sortrows(pairs,1); usedA=false(1,numel(active)); usedR=false(1,nR);
    for p=1:size(pairs,1)
        cost=pairs(p,1); ia=pairs(p,2); r=pairs(p,3);
        if cost>0.16 || usedA(ia) || usedR(r), continue; end
        id=active(ia); tracks(id)=append_track(tracks(id),kappa(m),roots(r));
        usedA(ia)=true; usedR(r)=true;
    end
    newActive=active(usedA);
    for r=find(~usedR)
        [tracks,newId]=start_track(tracks,[],kappa(m),roots(r));
        newActive(end+1)=newId; %#ok<AGROW>
    end
    active=newActive;
end
end

function [tracks,active]=start_track(tracks,active,kappa,root)
id=numel(tracks)+1;
tracks(id)=struct('kappa',kappa,'Omega',root.Omega,'Q',root.Q, ...
    'sigma_ratio',root.sigma_ratio,'raw_residual',root.raw_residual);
active(end+1)=id;
end

function track=append_track(track,kappa,root)
track.kappa(end+1)=kappa; track.Omega(end+1)=root.Omega;
track.Q(end+1)=root.Q; track.sigma_ratio(end+1)=root.sigma_ratio;
track.raw_residual(end+1)=root.raw_residual;
end
