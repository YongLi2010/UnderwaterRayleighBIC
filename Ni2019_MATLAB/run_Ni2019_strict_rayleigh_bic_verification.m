%% Three-part strict verification of the proposed single-Rayleigh candidate
clear; close all; clc;
navy=[.035 .10 .30]; orange=[1 .36 .08]; gray=[.45 .45 .45]; red=[.90 .05 .06];

kappa0=0.240035321859; Omega0=1-kappa0;
cfg=struct('a',1,'lambda',1/Omega0, ...
    'theta_i_deg',asind(kappa0/Omega0), ...
    'depths',[0.682459207769 0.122072875277], ...
    'widths',[0.658186442030 0.052212626051], ...
    'gaps',0.144800546474,'N',101,'K',14,'solve_scattering',false);

% Part 1: independently locate a leaky pole away from q=0, then track it.
poleTrack=ni2019_track_leaky_pole_to_rayleigh(cfg,kappa0, ...
    'DeltaStart',-1e-2,'DeltaEnd',-1e-6,'NumSteps',31,'Verbose',true);
part1Pass=poleTrack.same_branch && ~any(poleTrack.jump_detected) && ...
    min(poleTrack.mode_overlap(2:end-1))>.99 && ...
    poleTrack.independent_start.groove_content>.5;

% Part 2: separate the homogeneous denominator from aperture radiation.
radiationZero=ni2019_radiation_zero_diagnostic(cfg, ...
    'Kappa',kappa0,'Omega',Omega0,'TargetOrder',-1, ...
    'NList',[41 61 81 101],'KList',[6 8 10 12 14 16], ...
    'GridN',101,'GridK',14, ...
    'KappaOffsets',[-.02 -.01 -.005 0 .005 .01 .02], ...
    'OmegaOffsets',[-.02 -.01 0 .01 .02],'Verbose',true);
part2Pass=radiationZero.converged_intersection;

% Part 3: exact modal energy integrals and a reconstructed eigenfield.
energyDiagnostic=ni2019_eigenmode_energy_diagnostic(cfg, ...
    'NList',[41 61 81 101],'KList',[6 10 14], ...
    'HeightList',[1 2 5 10 20],'ReturnField',true, ...
    'X',linspace(0,1,321),'Y',linspace(-max(cfg.depths),2,361), ...
    'Verbose',true);
part3Pass=energyDiagnostic.base.strict_square_integrable;
overallPass=part1Pass && part2Pass && part3Pass;

fprintf('\nSTRICT THREE-PART VERDICT\n');
fprintf('  Part 1, nontrivial pole branch:       %s\n',passfail(part1Pass));
fprintf('  Part 2, converged radiation zero:     %s\n',passfail(part2Pass));
fprintf('  Part 3, square-integrable eigenfield: %s\n',passfail(part3Pass));
fprintf('  Overall Rayleigh-BIC verification:    %s\n',passfail(overallPass));
if ~overallPass
    fprintf('  Classification: Rayleigh-threshold pole candidate, not a verified BIC.\n');
end

fig=figure('Color','w','Position',[60 50 1160 870]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% (a) Independently located leaky pole reaches the Rayleigh endpoint.
nexttile;
plot(poleTrack.kappa,real(poleTrack.Omega),'-','Color',navy,'LineWidth',1.9); hold on;
plot(poleTrack.kappa,abs(poleTrack.kappa-1),'-','Color',orange,'LineWidth',1.5);
plot(poleTrack.kappa(1),real(poleTrack.Omega(1)),'s', ...
    'MarkerFaceColor',navy,'MarkerEdgeColor','w','MarkerSize',8);
plot(kappa0,Omega0,'o','MarkerFaceColor',red,'MarkerEdgeColor','w','MarkerSize',9);
xlabel('\kappa'); ylabel('Re \Omega'); grid on; box on;
lg=legend('independently tracked pole','n=-1 Rayleigh anomaly', ...
    'nonzero-q start','endpoint','Location','best'); paper_legend(lg);
title(sprintf('(a) Pole branch: %s',passfail(part1Pass))); paper_axes(gca);

% (b) The apparent numerator/denominator zero is not K converged.
nexttile;
Kvals=cellfun(@(s)s.K,radiationZero.center_by_K);
poleByK=cellfun(@(s)s.pole_denominator,radiationZero.center_by_K);
numByK=cellfun(@(s)s.finite_channel_numerator_residual, ...
    radiationZero.center_by_K);
semilogy(Kvals,poleByK,'o-','Color',navy,'LineWidth',1.7,'MarkerFaceColor',navy); hold on;
semilogy(Kvals,numByK,'s-','Color',orange,'LineWidth',1.7,'MarkerFaceColor',orange);
xlabel('groove truncation K'); ylabel('normalized residual'); grid on; box on;
lg=legend('pole denominator','radiation numerator','Location','best'); paper_legend(lg);
title(sprintf('(b) Pole-zero convergence: %s',passfail(part2Pass))); paper_axes(gca);

% (c) The grazing harmonic gives an exterior energy linear in H.
nexttile;
H=energyDiagnostic.base.energy.height_list;
Eg=energyDiagnostic.base.energy.exterior_target_grazing;
plot(H,Eg,'o-','Color',orange,'LineWidth',1.9,'MarkerFaceColor',orange); hold on;
yline(energyDiagnostic.base.energy.groove_total,'--','Color',navy,'LineWidth',1.5);
xlabel('exterior integration height H/a'); ylabel('dimensionless modal energy');
grid on; box on;
lg=legend('n=-1 grazing exterior energy','total groove energy','Location','best');
paper_legend(lg);
title(sprintf('(c) Square integrability: %s',passfail(part3Pass))); paper_axes(gca);

% (d) Homogeneous field: non-decaying grazing exterior component is visible.
nexttile;
fld=energyDiagnostic.base.mode.field;
level=log10(abs(fld.p)/max(abs(fld.p),[],'all')+1e-8);
imagesc(fld.x,fld.y,level); axis xy; clim([-5 0]); hold on;
for ell=1:numel(cfg.widths)
    xl=energyDiagnostic.base.operator.xleft(ell);
    rectangle('Position',[xl,-cfg.depths(ell),cfg.widths(ell),cfg.depths(ell)], ...
        'EdgeColor','k','LineWidth',1.1);
end
yline(0,'k-','LineWidth',1.1); xlabel('x/a'); ylabel('y/a');
colormap(ni2019_viridis(256)); cb=colorbar; cb.Label.String='log_{10}|p/p_{max}|';
title('(d) Threshold eigenfield'); paper_axes(gca);

outputDir=fullfile(pwd,'results');
if ~exist(outputDir,'dir'), mkdir(outputDir); end
exportgraphics(fig,fullfile(outputDir,'StrictRayleighBIC_three_part_verification.png'), ...
    'Resolution',220);
save(fullfile(outputDir,'StrictRayleighBIC_three_part_verification.mat'), ...
    'cfg','kappa0','Omega0','poleTrack','radiationZero','energyDiagnostic', ...
    'part1Pass','part2Pass','part3Pass','overallPass');

function s=passfail(tf)
if tf, s='PASS'; else, s='FAIL'; end
end

function paper_axes(ax)
set(ax,'Color','w','XColor','k','YColor','k','GridColor',[.75 .75 .75], ...
    'MinorGridColor',[.86 .86 .86]);
ax.Title.Color='k'; ax.XLabel.Color='k'; ax.YLabel.Color='k';
end

function paper_legend(lg)
set(lg,'Color','w','TextColor','k','EdgeColor',[.25 .25 .25]);
end
