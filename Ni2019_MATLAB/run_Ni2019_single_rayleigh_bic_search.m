%% Off-Gamma single-Rayleigh BIC search (stage 1 only)
clear; close all; clc;

a=1;
cfg=struct('lambda',1,'a',a,'theta_i_deg',0, ...
    'widths',[.085 .527]*a,'depths',[.487 .041]*a,'gaps',.101*a, ...
    'N',41,'K',7,'solve_scattering',false);

coarse=ni2019_optimize_single_rayleigh_bic(cfg,'Starts',16,'Display','off');
fprintf('Single-Rayleigh coarse search\n');
print_result(coarse);

% Re-optimize after increasing the modal truncation.
cfgFine=coarse.cfg; cfgFine.N=61; cfgFine.K=10;
fine=ni2019_optimize_single_rayleigh_bic(cfgFine,'Starts',4, ...
    'KappaRange',sort(coarse.kappa+[-.08 .08]), ...
    'DepthRange',[max(.005,min(coarse.depths_over_a)-.12), ...
                  min(.90,max(coarse.depths_over_a)+.12)], ...
    'WidthRange',[max(.02,min(coarse.widths_over_a)-.12), ...
                  min(.75,max(coarse.widths_over_a)+.12)], ...
    'GapRange',[max(.01,coarse.gap_over_a-.10), ...
                min(.40,coarse.gap_over_a+.10)],'Display','off');
fprintf('\nSingle-Rayleigh fine search\n');
print_result(fine);

outputDir=fullfile(pwd,'results');
if ~exist(outputDir,'dir'), mkdir(outputDir); end
save(fullfile(outputDir,'SingleRayleighBIC_search.mat'),'coarse','fine');
assignin('base','SingleRayleighBIC',fine);

function print_result(r)
fprintf('  kappa = %.12f, Omega = %.12f, theta = %.6f deg\n', ...
    r.kappa,r.Omega,r.theta_deg);
fprintf('  [d1 d2]/a = [%.12f %.12f]\n',r.depths_over_a);
fprintf('  [w1 w2 g]/a = [%.12f %.12f %.12f]\n', ...
    r.widths_over_a,r.gap_over_a);
fprintf('  sigma_min/sigma_max = %.3e\n',r.diagnostics.sigma_ratio);
fprintf('  finite-channel radiation fraction = %.3e\n', ...
    r.diagnostics.radiation_fraction);
fprintf('  grazing-order fraction = %.3e\n',r.diagnostics.grazing_fraction);
end
