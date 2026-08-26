%% Discover all groove-dominated outgoing poles at Gamma in 0<Omega<1.18
clear; clc;
rootDir=fileparts(mfilename('fullpath'));
outDir=fullfile(rootDir,'results','full_complex_band_180k');
if ~exist(outDir,'dir'), mkdir(outDir); end
S=load(fullfile(rootDir,'results','StrictRayleighBIC_180kHz_7p10deg_final.mat'));

% Moderate discovery truncation; accepted roots are subsequently refined at
% the production truncation in the full-band driver.
cfg=S.cfg; cfg.N=89; cfg.K=11; cfg.solve_scattering=false;
realSeeds=0.06:0.035:1.18;
imagSeeds=[0.001,0.012,0.05,0.16];
roots=struct([]); attempts=0;
fprintf('Gamma-pole discovery: %d seeds, (N,K)=(%d,%d)\n', ...
    numel(realSeeds)*numel(imagSeeds),cfg.N,cfg.K);

for ir=1:numel(realSeeds)
    for ii=1:numel(imagSeeds)
        attempts=attempts+1;
        seed=complex(realSeeds(ir),imagSeeds(ii));
        try
            p=ni2019_refine_outgoing_pole_kappa(cfg,0,seed, ...
                'OuterIterations',7,'Display','off');
        catch ME
            fprintf('  failed seed %.4f%+.4fi: %s\n',real(seed),imag(seed),ME.message);
            continue;
        end
        om=p.Omega;
        valid=isfinite(real(om)) && isfinite(imag(om)) && ...
            real(om)>0.015 && real(om)<1.205 && imag(om)>-1e-7 && ...
            imag(om)<0.8 && p.sigma_ratio<2e-7 && p.raw_residual<2e-4;
        if ~valid, continue; end
        if isempty(roots)
            roots=p;
        else
            distance=arrayfun(@(r)abs(r.Omega-om),roots);
            [dmin,id]=min(distance);
            if dmin<2e-4
                if p.sigma_ratio<roots(id).sigma_ratio, roots(id)=p; end
            else
                roots(end+1)=p; %#ok<SAGROW>
            end
        end
    end
    if mod(ir,5)==0
        fprintf('  %3d/%3d real seeds; %d unique candidates\n', ...
            ir,numel(realSeeds),numel(roots));
    end
end

if isempty(roots), error('No Gamma poles were discovered.'); end
[~,order]=sort(real([roots.Omega])); roots=roots(order);

% Merge roots that converge more tightly after sorting.
keep=true(size(roots));
for j=2:numel(roots)
    previous=find(keep(1:j-1),1,'last');
    if ~isempty(previous) && abs(roots(j).Omega-roots(previous).Omega)<5e-4
        if roots(j).sigma_ratio<roots(previous).sigma_ratio
            keep(previous)=false;
        else
            keep(j)=false;
        end
    end
end
roots=roots(keep);

fprintf('\nAccepted Gamma poles:\n');
for j=1:numel(roots)
    fprintf('  %2d  Omega=% .10f%+.10fi  Q=%10.3g  sigma=%.3e raw=%.3e\n', ...
        j,real(roots(j).Omega),imag(roots(j).Omega),roots(j).Q, ...
        roots(j).sigma_ratio,roots(j).raw_residual);
end
save(fullfile(outDir,'gamma_discovery.mat'),'roots','cfg','realSeeds', ...
    'imagSeeds','attempts');
