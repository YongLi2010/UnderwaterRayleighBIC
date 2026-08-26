%% Selective truncation validation of the full complex-band census
clear; clc;
rootDir=fileparts(mfilename('fullpath'));
dataDir=fullfile(rootDir,'results','full_complex_band_180k');
S=load(fullfile(rootDir,'results','StrictRayleighBIC_180kHz_7p10deg_final.mat'));
D=load(fullfile(dataDir,'positive_pole_census.mat'));
cfg=S.cfg; cfg.N=121; cfg.K=15; cfg.solve_scattering=false;
kCheck=[0,.11,.25,.50];
rows=[];
for m=1:numel(kCheck)
    [~,id]=min(abs(D.kappa-kCheck(m))); roots=D.rootsByK{id};
    for j=1:numel(roots)
        try
            p=ni2019_refine_outgoing_pole_kappa(cfg,D.kappa(id),roots(j).Omega, ...
                'OuterIterations',9,'Display','off');
            rows=[rows;D.kappa(id),j,real(roots(j).Omega),imag(roots(j).Omega), ...
                real(p.Omega),imag(p.Omega),abs(p.Omega-roots(j).Omega), ...
                p.sigma_ratio,p.raw_residual]; %#ok<AGROW>
        catch
            rows=[rows;D.kappa(id),j,real(roots(j).Omega),imag(roots(j).Omega), ...
                NaN,NaN,NaN,NaN,NaN]; %#ok<AGROW>
        end
    end
end
validation=array2table(rows,'VariableNames',{'kappa','root_id', ...
    'Omega_low_real','Omega_low_imag','Omega_high_real','Omega_high_imag', ...
    'absolute_shift','sigma_high','raw_high'});
writetable(validation,fullfile(dataDir,'full_band_truncation_validation.csv'));
save(fullfile(dataDir,'full_band_truncation_validation.mat'),'validation','cfg');
fprintf('Full-band selective truncation validation (89/11 -> 121/15):\n');
fprintf('  roots checked: %d\n',height(validation));
fprintf('  median |Delta Omega|: %.3e\n',median(validation.absolute_shift,'omitnan'));
fprintf('  max |Delta Omega|: %.3e\n',max(validation.absolute_shift,[],'omitnan'));
fprintf('  max high sigma: %.3e\n',max(validation.sigma_high,[],'omitnan'));
fprintf('  max high raw residual: %.3e\n',max(validation.raw_high,[],'omitnan'));
