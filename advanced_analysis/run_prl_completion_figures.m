%% PRL completion figures for the strict underwater Rayleigh BIC
% This script redraws the four main figures and four theory supplements from
% audited solver outputs.  It never reads the superseded intermediate-root
% pole cache and never uses the invalid hard-coded router-convergence table.
%
% Run from the project root:
%   matlab -batch "run('advanced_analysis/run_prl_completion_figures.m')"

clearvars;
close all;
clc;

analysisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(analysisDir);
solverDir = fullfile(rootDir,'Ni2019_MATLAB');
outDir = fullfile(rootDir,'arxiv_theory_paper','figures');
if ~exist(outDir,'dir'), mkdir(outDir); end
addpath(solverDir);

designFile = fullfile(solverDir,'results', ...
    'StrictRayleighBIC_200kHz_min1mm.mat');
cacheFile = fullfile(solverDir,'results','PRL_article_figure_data_v5.mat');
poleFile = fullfile(analysisDir,'final_root_pole','final_root_pole_track.mat');
convFile = fullfile(analysisDir,'extended_convergence','extended_convergence.mat');
robustFile = fullfile(analysisDir,'robustness_fields','robustness_fields.mat');
scalingFile = fullfile(analysisDir,'scaling_environment', ...
    'scaling_environment_data.mat');
required = {designFile,cacheFile,poleFile,convFile,robustFile,scalingFile};
for j = 1:numel(required)
    assert(exist(required{j},'file')==2,'Missing audited input: %s',required{j});
end

D = load(designFile);
C = load(cacheFile,'strictMode','heightValues','exteriorEnergy', ...
    'grazingControl','phasors','phasorOrders');
P = load(poleFile);
E = load(convFile);
R = load(robustFile);
S = load(scalingFile);
p = palette();

%% Figure 1: geometry and channel topology
fig = new_figure([7.0 4.55]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl); draw_geometry(ax,D,p); panel_label(ax,'a');
ax = nexttile(tl); draw_rayleigh_map(ax,D,p); panel_label(ax,'b');
ax = nexttile(tl); draw_channel_status(ax,D,p); panel_label(ax,'c');
ax = nexttile(tl); draw_physical_rayleigh(ax,D,p); panel_label(ax,'d');
finish_figure(fig);
export_pdf(fig,fullfile(outDir,'fig1_structure_topology.pdf'));

%% Figure 2: strict null, final-root pole, and modal convergence
fig = new_figure([7.0 4.75]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl); draw_final_pole(ax,P,D,p); panel_label(ax,'a');
ax = nexttile(tl); draw_final_q(ax,P,p); panel_label(ax,'b');
ax = nexttile(tl); draw_strict_amplitudes(ax,C,D,P,p); panel_label(ax,'c');
ax = nexttile(tl); draw_root_sequence(ax,E,p); panel_label(ax,'d');
finish_figure(fig);
export_pdf(fig,fullfile(outDir,'fig2_strict_bic_evidence.pdf'));

%% Figure 3: homogeneous fields and quantitative localization
fig = new_figure([7.0 5.15]);
tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
fieldAxes = gobjects(1,3);
for j = 1:3
    ax = nexttile(tl,j);
    fieldAxes(j) = ax;
    draw_homogeneous_field(ax,R.fieldStates(j),R.xNom,p);
    panel_label(ax,char('a'+j-1));
end
ax = nexttile(tl,4,[1 2]);
draw_exterior_norm(ax,C,p); panel_label(ax,'d');
ax = nexttile(tl,6);
draw_field_metrics(ax,R.fieldStates,p); panel_label(ax,'e');
colormap(fig,ni2019_viridis(256));
cb = colorbar(fieldAxes(3),'eastoutside');
cb.Layout.Tile = 'east';
cb.Label.String = '|p| / max_{groove}|p|';
cb.FontName = 'Helvetica'; cb.FontSize = 7.2;
finish_figure(fig);
export_pdf(fig,fullfile(outDir,'fig3_homogeneous_fields.pdf'));

%% Figure 4: two-cavity phase degrees of freedom and fabrication sensitivity
fig = new_figure([7.0 5.25]);
tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl); draw_dof_ablation(ax,R,p); panel_label(ax,'a');
ax = nexttile(tl); draw_one_phasor(ax,R,-1,p); panel_label(ax,'b');
ax = nexttile(tl); draw_one_phasor(ax,R,0,p); panel_label(ax,'c');
ax = nexttile(tl); draw_tolerance_quantile(ax,R,'strict',p); panel_label(ax,'d');
ax = nexttile(tl); draw_tolerance_quantile(ax,R,'threshold',p); panel_label(ax,'e');
ax = nexttile(tl); draw_parameter_sensitivity(ax,R,p); panel_label(ax,'f');
finish_figure(fig);
export_pdf(fig,fullfile(outDir,'fig4_dof_tolerance.pdf'));

%% Supplement S1: convergence details
fig = new_figure([7.0 5.0]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl); draw_parameter_drift(ax,E,p); panel_label(ax,'a');
ax = nexttile(tl); draw_root_sequence(ax,E,p); panel_label(ax,'b');
ax = nexttile(tl); draw_fixed_crosscheck(ax,E,p); panel_label(ax,'c');
ax = nexttile(tl); draw_fixed_channel_errors(ax,E,p); panel_label(ax,'d');
finish_figure(fig);
export_pdf(fig,fullfile(outDir,'figS1_extended_convergence.pdf'));

%% Supplement S2: exact similarity and fixed-device sound-speed control
fig = new_figure([7.0 5.0]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl); draw_scaling_speed(ax,S,p); panel_label(ax,'a');
ax = nexttile(tl); draw_scaling_size(ax,S,p); panel_label(ax,'b');
ax = nexttile(tl); draw_retuned_residual(ax,S,p); panel_label(ax,'c');
ax = nexttile(tl); draw_fixed_environment(ax,S,p); panel_label(ax,'d');
finish_figure(fig);
export_pdf(fig,fullfile(outDir,'figS2_scaling_environment.pdf'));

%% Supplement S3: complete tolerance distributions
fig = new_figure([7.0 5.0]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl); draw_tolerance_quantile(ax,R,'strict',p); panel_label(ax,'a');
ax = nexttile(tl); draw_tolerance_quantile(ax,R,'threshold',p); panel_label(ax,'b');
ax = nexttile(tl); draw_tolerance_cdf(ax,R,p); panel_label(ax,'c');
ax = nexttile(tl); draw_tolerance_scatter(ax,R,p); panel_label(ax,'d');
finish_figure(fig);
export_pdf(fig,fullfile(outDir,'figS3_tolerance_statistics.pdf'));

%% Supplement S4: proposed water-tank protocol (schematic only)
fig = new_figure([7.0 4.65]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl); draw_tank_schematic(ax,p); panel_label(ax,'a');
ax = nexttile(tl); draw_scan_plan(ax,p); panel_label(ax,'b');
ax = nexttile(tl); draw_floquet_extraction(ax,p); panel_label(ax,'c');
ax = nexttile(tl); draw_ringdown_plan(ax,p); panel_label(ax,'d');
finish_figure(fig);
export_pdf(fig,fullfile(outDir,'figS5_proposed_experiment.pdf'));

write_manifest(outDir,required,P,E,R,S);
close all;
fprintf('Generated four main and four supplemental PRL figures in %s\n',outDir);

%% Figure helpers
function fig = new_figure(sizeInches)
fig = figure('Visible','off','Color','w','Units','inches', ...
    'Position',[0.5 0.5 sizeInches(1) sizeInches(2)], ...
    'PaperUnits','inches','PaperPosition',[0 0 sizeInches(1) sizeInches(2)]);
set(fig,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica', ...
    'DefaultAxesFontSize',8,'DefaultTextFontSize',8, ...
    'DefaultAxesLineWidth',0.65,'DefaultLineLineWidth',1.15);
end

function finish_figure(fig)
axs = findall(fig,'Type','axes');
for j = 1:numel(axs)
    ax = axs(j);
    set(ax,'FontName','Helvetica','FontSize',8,'LineWidth',0.65, ...
        'Box','on','Layer','top','TickDir','out','TickLength',[0.018 0.018], ...
        'XColor',[0.12 0.13 0.15],'YColor',[0.12 0.13 0.15], ...
        'GridColor',[0.82 0.83 0.85],'GridAlpha',0.25, ...
        'MinorGridAlpha',0.12);
    ax.TitleFontWeight = 'normal';
    ax.TitleFontSizeMultiplier = 1;
end
legs = findall(fig,'Type','legend');
for j = 1:numel(legs)
    set(legs(j),'FontName','Helvetica','FontSize',7.0,'Box','off');
end
end

function export_pdf(fig,fileName)
drawnow;
exportgraphics(fig,fileName,'ContentType','vector','BackgroundColor','white');
assert(exist(fileName,'file')==2,'Failed to export %s',fileName);
close(fig);
end

function panel_label(ax,label)
text(ax,-0.12,1.07,sprintf('(%s)',label),'Units','normalized', ...
    'FontName','Helvetica','FontSize',9,'FontWeight','bold', ...
    'Color',[0.08 0.09 0.10],'Clipping','off');
end

function p = palette()
p.blue=[0.00 0.35 0.62]; p.orange=[0.84 0.32 0.00];
p.green=[0.00 0.50 0.38]; p.red=[0.75 0.12 0.16];
p.purple=[0.55 0.28 0.58]; p.gray=[0.40 0.43 0.47];
p.light=[0.92 0.95 0.97]; p.dark=[0.10 0.11 0.13];
p.gold=[0.90 0.62 0.00];
end

function draw_geometry(ax,D,p)
hold(ax,'on');
a = 1e3*D.aPhysical; d = D.depthsMm(:).';
w = D.widthsMm(:).'; g = D.gapMm;
x1=(a-sum(w)-g)/2; x2=x1+w(1)+g;
patch(ax,[0 a a 0],[-max(d)-0.45 -max(d)-0.45 0 0],p.light, ...
    'EdgeColor',p.dark,'LineWidth',0.8);
rectangle(ax,'Position',[x1 -d(1) w(1) d(1)],'FaceColor','w', ...
    'EdgeColor',p.blue,'LineWidth',1.4);
rectangle(ax,'Position',[x2 -d(2) w(2) d(2)],'FaceColor','w', ...
    'EdgeColor',p.orange,'LineWidth',1.4);
plot(ax,[0 a],[0 0],'-','Color',p.dark,'LineWidth',0.9);
quiver(ax,0.39*a,0.94,0.72,-0.70,0,'Color',p.dark,'LineWidth',1.0, ...
    'MaxHeadSize',0.7);
quiver(ax,0.40*a,0.12,-1.25,0,0,'Color',p.red,'LineWidth',1.0, ...
    'MaxHeadSize',0.7);
text(ax,0.43*a,0.85,'n=0','Color',p.dark,'FontSize',7.3);
text(ax,0.43*a,0.10,'n=-1, k_y=0','Color',p.red,'FontSize',7.3);
text(ax,x1+w(1)/2,-0.54*d(1),sprintf('d_1=%.3f',d(1)), ...
    'HorizontalAlignment','center','Color',p.blue,'FontSize',7.2);
text(ax,x2+w(2)/2,-0.53*d(2),sprintf('d_2=%.3f',d(2)), ...
    'HorizontalAlignment','center','Color',p.orange,'FontSize',7.2);
text(ax,x1+w(1)/2,0.28,sprintf('w_1=%.3f',w(1)), ...
    'HorizontalAlignment','center','Color',p.blue,'FontSize',7.1);
text(ax,x2+w(2)/2,0.28,sprintf('w_2=%.3f',w(2)), ...
    'HorizontalAlignment','center','Color',p.orange,'FontSize',7.1);
text(ax,a/2,1.18,sprintf('a=%.6f mm;  g=%.3f mm',a,g), ...
    'HorizontalAlignment','center','Color',p.dark,'FontSize',7.2);
xlabel(ax,'x (mm)'); ylabel(ax,'y (mm)');
xlim(ax,[-0.08*a 1.08*a]); ylim(ax,[-max(d)-0.55 1.40]);
title(ax,'two-groove water-loaded metagrating');
end

function draw_rayleigh_map(ax,D,p)
hold(ax,'on');
k=linspace(-0.26,0.26,401);
plot(ax,k,abs(k-1),'-','Color',p.orange,'LineWidth',1.35);
plot(ax,k,abs(k+1),'-','Color',p.blue,'LineWidth',1.35);
k0=D.xFinal(1); O0=D.OmegaFinal;
plot(ax,k0,O0,'o','MarkerSize',6.5,'MarkerFaceColor',p.red, ...
    'MarkerEdgeColor','w','LineWidth',0.7);
plot(ax,-k0,O0,'o','MarkerSize',6.5,'MarkerFaceColor',p.red, ...
    'MarkerEdgeColor','w','LineWidth',0.7);
text(ax,k0+0.014,O0-0.035,'-1 threshold','Color',p.orange,'FontSize',7.0);
text(ax,-k0-0.014,O0-0.035,'+1 threshold','Color',p.blue, ...
    'FontSize',7.0,'HorizontalAlignment','right');
text(ax,0,1.075,'n=0 propagating','HorizontalAlignment','center', ...
    'Color',p.dark,'FontSize',7.2);
xlabel(ax,'Bloch number  \kappa'); ylabel(ax,'normalized frequency  \Omega');
xlim(ax,[-0.26 0.26]); ylim(ax,[0.74 1.12]); grid(ax,'on');
legend(ax,'n=-1 Rayleigh line','n=+1 Rayleigh line','BIC pair', ...
    'Location','south','NumColumns',2);
title(ax,'single-Rayleigh topology');
end

function draw_channel_status(ax,D,p)
hold(ax,'on');
n=-1:1; k=D.xFinal(1); O=D.OmegaFinal;
qplus=sqrt(max(O^2-(k+n).^2,0));
qminus=sqrt(max(O^2-(-k+n).^2,0));
b1=bar(ax,n-0.17,qplus,0.32,'FaceColor',p.blue,'EdgeColor','none');
b2=bar(ax,n+0.17,qminus,0.32,'FaceColor',p.orange,'EdgeColor','none');
plot(ax,-1,0,'o','Color',p.red,'MarkerFaceColor','w','MarkerSize',5,'LineWidth',1.2);
plot(ax,+1,0,'o','Color',p.red,'MarkerFaceColor','w','MarkerSize',5,'LineWidth',1.2);
text(ax,-1,0.055,'threshold','HorizontalAlignment','center','Color',p.red,'FontSize',7);
text(ax,+1,0.055,'threshold','HorizontalAlignment','center','Color',p.red,'FontSize',7);
text(ax,-1,0.22,'closed at +\theta','HorizontalAlignment','center','Color',p.gray,'FontSize',6.8,'Rotation',90);
text(ax,+1,0.22,'closed at -\theta','HorizontalAlignment','center','Color',p.gray,'FontSize',6.8,'Rotation',90);
set(ax,'XTick',n,'XTickLabel',{'-1','0','+1'});
xlabel(ax,'Floquet order n'); ylabel(ax,'Re(k_{y,n}a/2\pi)');
ylim(ax,[-0.02 1.12*max([qplus qminus])]);
legend(ax,[b1 b2],{'+\theta','-\theta'},'Location','northwest');
title(ax,'open, threshold, and evanescent orders');
end

function draw_physical_rayleigh(ax,D,p)
hold(ax,'on');
theta=linspace(0,12,400);
fRA=D.cWater./(D.aPhysical*(1+sind(theta)))/1e3;
plot(ax,theta,fRA,'-','Color',p.blue,'LineWidth',1.35);
plot(ax,D.thetaFinal,D.fTarget/1e3,'o','MarkerSize',6.5, ...
    'MarkerFaceColor',p.red,'MarkerEdgeColor','w');
xline(ax,D.thetaFinal,':','Color',p.gray,'LineWidth',0.8);
yline(ax,D.fTarget/1e3,':','Color',p.gray,'LineWidth',0.8);
text(ax,D.thetaFinal+0.28,D.fTarget/1e3+1.1, ...
    sprintf('%.6f^{\\circ}, 200 kHz',D.thetaFinal),'Color',p.red,'FontSize',7.1);
xlabel(ax,'|\theta| (deg)'); ylabel(ax,'Rayleigh frequency (kHz)');
xlim(ax,[0 12]); grid(ax,'on');
title(ax,'dimensional operating point');
end

function draw_final_pole(ax,P,D,p)
tr=P.track; finite=tr.delta_kappa~=0 & isfinite(tr.Q) & tr.Q>0;
hold(ax,'on');
plot(ax,tr.kappa(finite),real(tr.Omega(finite)),'o-', ...
    'Color',p.blue,'MarkerFaceColor',p.blue,'MarkerSize',3.2);
plot(ax,tr.kappa(finite),1-tr.kappa(finite),'--','Color',p.orange,'LineWidth',1.1);
plot(ax,D.xFinal(1),D.OmegaFinal,'o','MarkerSize',6.5, ...
    'MarkerFaceColor',p.red,'MarkerEdgeColor','w');
text(ax,D.xFinal(1)-1.7e-4,D.OmegaFinal-1.9e-5,'strict endpoint', ...
    'HorizontalAlignment','right','Color',p.red,'FontSize',7.1);
xlabel(ax,'\kappa'); ylabel(ax,'Re \Omega_p'); grid(ax,'on');
legend(ax,'continued pole','n=-1 Rayleigh line','full-operator null', ...
    'Location','northwest');
title(ax,'same final geometry, N/K=313/39');
end

function draw_final_q(ax,P,p)
tr=P.track; finite=tr.delta_kappa~=0 & isfinite(tr.Q) & tr.Q>0 & tr.sigma_ratio<1e-7;
dk=abs(tr.delta_kappa(finite)); q=tr.Q(finite);
[dk,id]=sort(dk); q=q(id);
loglog(ax,dk,q,'o-','Color',p.blue,'MarkerFaceColor',p.blue,'MarkerSize',3.2); hold(ax,'on');
fitBand=P.fitBand; fitDk=sort(abs(tr.delta_kappa(fitBand)));
loglog(ax,fitDk,exp(P.fitQ(2))*fitDk.^P.fitQ(1),'--','Color',p.orange,'LineWidth',1.1);
text(ax,0.08,0.14,sprintf('local slope %.3f',P.fitQ(1)), ...
    'Units','normalized','Color',p.orange,'FontSize',7.2);
text(ax,0.08,0.06,sprintf('Q_{max}=%.2e',max(q)), ...
    'Units','normalized','Color',p.dark,'FontSize',7.2);
xlabel(ax,'|\Delta\kappa|'); ylabel(ax,'Q'); grid(ax,'on');
legend(ax,'complex pole','local fit','Location','southwest');
title(ax,'finite-model pole linewidth');
end

function draw_strict_amplitudes(ax,C,D,P,p)
orders=-3:3; opOrders=C.strictMode.full_operator.orders;
A=C.strictMode.mode.A; vals=nan(size(orders));
for j=1:numel(orders), vals(j)=abs(A(opOrders==orders(j))); end
vals=max(vals,1e-16);
semilogy(ax,orders,vals,'o','Color',p.blue,'MarkerFaceColor',p.blue, ...
    'MarkerSize',4.8,'LineWidth',0.9); hold(ax,'on');
semilogy(ax,[-1 0],vals(orders==-1 | orders==0),'o','Color',p.red, ...
    'MarkerFaceColor',p.red,'MarkerSize',6.0);
text(ax,-1,3e-15,'threshold','HorizontalAlignment','center','Color',p.red,'FontSize',6.9);
text(ax,0,3e-15,'open','HorizontalAlignment','center','Color',p.red,'FontSize',6.9);
text(ax,0.04,0.93,sprintf('reduced  %.2e',D.rootSigma(end)), ...
    'Units','normalized','Color',p.dark,'FontSize',7.0);
text(ax,0.04,0.84,sprintf('full       %.2e',P.fullEndpointSigma), ...
    'Units','normalized','Color',p.dark,'FontSize',7.0);
set(ax,'XTick',orders,'XTickLabel',{'-3','-2','-1','0','+1','+2','+3'});
xlabel(ax,'order n'); ylabel(ax,'|A_n^{hom}|');
ylim(ax,[5e-17 1]); grid(ax,'on');
title(ax,'strict homogeneous channel amplitudes');
end

function draw_root_sequence(ax,E,p)
K=E.rootTable.K;
semilogy(ax,K,E.rootTable.strict_sigma_ratio,'o-','Color',p.blue, ...
    'MarkerFaceColor',p.blue,'MarkerSize',3.8); hold(ax,'on');
semilogy(ax,K,E.rootTable.full_sigma_ratio,'s--','Color',p.orange, ...
    'MarkerFaceColor','w','MarkerSize',3.8);
semilogy(ax,K,max(E.rootTable.threshold_amplitude_over_Cnorm,1e-17),'^-', ...
    'Color',p.green,'MarkerFaceColor',p.green,'MarkerSize',3.8);
semilogy(ax,K,max(E.rootTable.finite_open_amplitude_over_Cnorm,1e-17),'v-', ...
    'Color',p.purple,'MarkerFaceColor',p.purple,'MarkerSize',3.8);
text(ax,0.05,0.09,sprintf('max |\\Delta x| at K=47: %.1e', ...
    E.rootTable.max_parameter_step(end)),'Units','normalized', ...
    'Color',p.dark,'FontSize',6.9);
xlabel(ax,'groove truncation K'); ylabel(ax,'dimensionless diagnostic');
ylim(ax,[5e-17 1e-10]); grid(ax,'on');
legend(ax,'strict operator','full operator','threshold amplitude', ...
    'finite-open amplitude','Location','northwest','NumColumns',2);
title(ax,'reoptimized root sequence');
end

function draw_homogeneous_field(ax,state,xNom,p)
levels=linspace(0,1,33);
contourf(ax,state.X,state.Y,min(state.pmag,1),levels,'LineStyle','none'); hold(ax,'on');
mask_solid_regions(ax,xNom,-0.72);
draw_grooves(ax,xNom,p.dark);
axis(ax,'equal'); xlim(ax,[0 1]); ylim(ax,[-0.72 1.02]);
xlabel(ax,'x/a'); ylabel(ax,'y/a'); caxis(ax,[0 1]);
name=state.state;
if contains(name,'near'), name='near-BIC surrogate'; end
if contains(name,'ordinary'), name='ordinary-resonance surrogate'; end
title(ax,sprintf('%s,  \\Omega=%.6f',name,real(state.Omega)));
end

function mask_solid_regions(ax,xNom,ymin)
w=xNom(4:5); d=xNom(2:3); g=xNom(6);
x1=(1-sum(w)-g)/2; x2=x1+w(1)+g;
rects=[0,ymin,x1,-ymin; ...
    x1,ymin,w(1),max(d(1)+ymin,0); ...
    x1+w(1),ymin,x2-(x1+w(1)),-ymin; ...
    x2,ymin,w(2),max(d(2)+ymin,0); ...
    x2+w(2),ymin,1-(x2+w(2)),-ymin];
for ir=1:size(rects,1)
    if rects(ir,3)>0 && rects(ir,4)>0
        rectangle(ax,'Position',rects(ir,:),'FaceColor','w', ...
            'EdgeColor','none');
    end
end
end

function draw_grooves(ax,xNom,color)
w=xNom(4:5); d=xNom(2:3); g=xNom(6);
x1=(1-sum(w)-g)/2; x2=x1+w(1)+g;
plot(ax,[0 x1 x1 x1+w(1) x1+w(1) x2 x2 x2+w(2) x2+w(2) 1], ...
    [0 0 -d(1) -d(1) 0 0 -d(2) -d(2) 0 0],'-', ...
    'Color',color,'LineWidth',0.8);
end

function draw_exterior_norm(ax,C,p)
semilogy(ax,C.heightValues,C.exteriorEnergy,'-','Color',p.blue,'LineWidth',1.4); hold(ax,'on');
semilogy(ax,C.heightValues,C.grazingControl,'--','Color',p.red,'LineWidth',1.25);
text(ax,10.5,C.exteriorEnergy(end)*1.13,'saturates','Color',p.blue,'FontSize',7.0);
text(ax,13,C.grazingControl(end)*0.80,'linear threshold tail','Color',p.red,'FontSize',7.0);
xlabel(ax,'integration height H/a'); ylabel(ax,'exterior norm proxy');
xlim(ax,[0 max(C.heightValues)]); grid(ax,'on');
legend(ax,'strict BIC','same evanescent field + nonzero grazing amplitude', ...
    'Location','southeast');
title(ax,'square-integrability test');
end

function draw_field_metrics(ax,states,p)
v=[max([states.threshold_error],1e-16); ...
   max([states.finite_open_error],1e-16); ...
   max([states.exterior_to_groove_energy],1e-16)].';
cols=[p.red;p.orange;p.blue]; hold(ax,'on');
for j=1:3
    plot(ax,1:3,v(:,j),'o-','Color',cols(j,:), ...
        'MarkerFaceColor',cols(j,:),'MarkerSize',4.1,'LineWidth',1.05);
end
set(ax,'YScale','log','XTick',1:3, ...
    'XTickLabel',{'strict','near','ordinary'});
xtickangle(ax,18); ylabel(ax,'normalized metric'); grid(ax,'on');
legend(ax,'|A_{-1}|/||A||','|A_0|/||A||','exterior/groove norm', ...
    'Location','southwest');
title(ax,'channel leakage and localization');
end

function draw_dof_ablation(ax,R,p)
vals=[R.singleSigma R.nominalSigma];
semilogy(ax,1:2,vals,'o-','Color',p.dark,'LineWidth',1.0, ...
    'MarkerSize',7,'MarkerFaceColor',p.blue); hold(ax,'on');
plot(ax,1,vals(1),'o','Color',p.gray,'MarkerFaceColor',p.gray,'MarkerSize',7);
set(ax,'YScale','log','XTick',1:2,'XTickLabel',{'one cavity','two cavities'});
ylim(ax,[5e-17 1]); ylabel(ax,'\sigma_{min}/\sigma_{max}'); grid(ax,'on');
text(ax,1,vals(1)*2.0,sprintf('%.3f',vals(1)),'HorizontalAlignment','center','FontSize',7.1);
text(ax,2,vals(2)*9,sprintf('%.2e',vals(2)),'HorizontalAlignment','center','FontSize',7.1);
text(ax,0.05,0.05,'searched rectangular family','Units','normalized','Color',p.dark,'FontSize',6.9);
title(ax,'bounded cavity ablation');
end

function draw_one_phasor(ax,R,targetOrder,p)
id=find(R.phasorOrders==targetOrder,1); z=R.phasors(id,:);
scale=max(abs(z)); if scale==0, scale=1; end; z=z/scale;
hold(ax,'on');
quiver(ax,0,0,real(z(1)),imag(z(1)),0,'Color',p.blue,'LineWidth',1.55,'MaxHeadSize',0.45);
quiver(ax,real(z(1)),imag(z(1)),real(z(2)),imag(z(2)),0, ...
    'Color',p.orange,'LineWidth',1.55,'MaxHeadSize',0.45);
plot(ax,real(sum(z)),imag(sum(z)),'o','Color',p.red,'MarkerFaceColor',p.red,'MarkerSize',4.5);
plot(ax,0,0,'o','Color',p.dark,'MarkerFaceColor',p.dark,'MarkerSize',3.5);
axis(ax,'equal'); lim=1.18; xlim(ax,[-lim lim]); ylim(ax,[-lim lim]);
xline(ax,0,':','Color',[.75 .76 .78]); yline(ax,0,':','Color',[.75 .76 .78]);
xlabel(ax,'Re aperture source'); ylabel(ax,'Im aperture source');
text(ax,0.04,0.92,sprintf('n=%d; closure %.1e',targetOrder,abs(sum(z))), ...
    'Units','normalized','Color',p.red,'FontSize',7.0);
legend(ax,'groove 1','groove 2','sum','Location','southwest');
title(ax,sprintf('A_{%d}^{hom}=0',targetOrder));
end

function draw_tolerance_quantile(ax,R,metric,p)
hold(ax,'on'); modes={'independent','joint'}; cols=[p.blue;p.orange];
for im=1:2
    med=nan(size(R.levels)); lo=med; hi=med;
    for il=1:numel(R.levels)
        q=R.summaryRows(strcmp({R.summaryRows.mode},modes{im}) & ...
            [R.summaryRows.level_pct]==100*R.levels(il));
        if strcmp(metric,'strict')
            lo(il)=q.strict_raw_q05; med(il)=q.strict_raw_q50; hi(il)=q.strict_raw_q95;
        else
            lo(il)=q.threshold_error_q05; med(il)=q.threshold_error_q50; hi(il)=q.threshold_error_q95;
        end
    end
    errorbar(ax,100*R.levels,med,med-lo,hi-med,'o-', ...
        'Color',cols(im,:),'MarkerFaceColor',cols(im,:), ...
        'MarkerSize',3.8,'LineWidth',1.05,'CapSize',4);
end
set(ax,'YScale','log','XTick',100*R.levels);
xlabel(ax,'perturbation bound (%)'); grid(ax,'on');
if strcmp(metric,'strict')
    ylabel(ax,'strict raw residual'); title(ax,'fixed-point radiation residual');
    yline(ax,R.tolBaseline.strict_raw_residual,':','Color',p.gray,'LineWidth',0.8, ...
        'HandleVisibility','off');
else
    ylabel(ax,'|A_{-1}|/||A||'); title(ax,'unconstrained threshold error');
    yline(ax,R.tolBaseline.threshold_amplitude_error,':','Color',p.gray,'LineWidth',0.8, ...
        'HandleVisibility','off');
end
legend(ax,'one parameter','all five parameters','Location','northwest');
end

function draw_parameter_sensitivity(ax,R,p)
names={'d_1','d_2','w_1','w_2','g'};
fields={'delta_d1_pct','delta_d2_pct','delta_w1_pct','delta_w2_pct','delta_g_pct'};
rows=R.toleranceRows;
mask=strcmp({rows.mode},'independent') & [rows.level_pct]==1 & [rows.solver_ok];
rows=rows(mask); med=nan(1,5); q95=med;
for j=1:5
    active=arrayfun(@(x) abs(x.(fields{j}))>0,rows);
    vals=[rows(active).threshold_amplitude_error];
    med(j)=median(vals); q95(j)=prctile(vals,95);
end
bar(ax,1:5,med,0.62,'FaceColor',p.blue,'EdgeColor','none'); hold(ax,'on');
for j=1:5
    plot(ax,[j j],[med(j) q95(j)],'-','Color',p.dark,'LineWidth',0.9);
    plot(ax,j,q95(j),'_','Color',p.dark,'MarkerSize',6);
end
set(ax,'YScale','log','XTick',1:5,'XTickLabel',names);
xlabel(ax,'perturbed dimension'); ylabel(ax,'median / 95th percentile'); grid(ax,'on');
title(ax,'one-at-a-time threshold error, \pm1%');
end

function draw_parameter_drift(ax,E,p)
K=E.rootTable.K; x=E.xSequence; ref=x(end,:); dx=(x-ref)*1e5;
cols=[p.red;p.blue;p.orange;p.green;p.purple;p.gray];
hold(ax,'on');
for j=1:6
    plot(ax,K,dx(:,j),'o-','Color',cols(j,:),'MarkerSize',3,'LineWidth',0.9);
end
yline(ax,0,':','Color',p.dark);
xlabel(ax,'K'); ylabel(ax,'10^5(x_K-x_{K=47})'); grid(ax,'on');
legend(ax,'\kappa','d_1/a','d_2/a','w_1/a','w_2/a','g/a', ...
    'Location','eastoutside');
title(ax,'reoptimized parameter drift');
end

function draw_fixed_crosscheck(ax,E,p)
semilogy(ax,E.fixedTable.K,E.fixedTable.strict_sigma_ratio,'o-', ...
    'Color',p.gray,'MarkerFaceColor',p.gray,'MarkerSize',4); hold(ax,'on');
semilogy(ax,E.fixedTable.K,E.fixedTable.full_sigma_ratio,'s--', ...
    'Color',p.orange,'MarkerFaceColor','w','MarkerSize',4);
xlabel(ax,'K'); ylabel(ax,'\sigma_{min}/\sigma_{max}'); grid(ax,'on');
legend(ax,'strict','full','Location','southwest');
title(ax,'fixed K=47 geometry across truncations');
end

function draw_fixed_channel_errors(ax,E,p)
semilogy(ax,E.fixedTable.K,max(E.fixedTable.threshold_amplitude_over_Cnorm,1e-17), ...
    'o-','Color',p.red,'MarkerFaceColor',p.red,'MarkerSize',4); hold(ax,'on');
semilogy(ax,E.fixedTable.K,max(E.fixedTable.finite_open_amplitude_over_Cnorm,1e-17), ...
    's-','Color',p.blue,'MarkerFaceColor',p.blue,'MarkerSize',4);
xlabel(ax,'K'); ylabel(ax,'amplitude / ||C||'); grid(ax,'on');
legend(ax,'threshold n=-1','finite-open n=0','Location','southwest');
title(ax,'fixed-geometry channel mismatch');
end

function draw_scaling_speed(ax,S,p)
t=S.scalingBySpeed;
plot(ax,t.sound_speed_m_per_s,t.f_BIC_Hz/1e3,'o-', ...
    'Color',p.blue,'MarkerFaceColor',p.blue,'MarkerSize',3.5); hold(ax,'on');
plot(ax,[1400 1600],[1400 1600]*(200/1500),'--','Color',p.orange);
xlabel(ax,'sound speed c (m/s)'); ylabel(ax,'retuned f_{BIC} (kHz)'); grid(ax,'on');
legend(ax,'modal operator','f_0c/c_0','Location','northwest');
title(ax,'fixed geometry, exact frequency retuning');
end

function draw_scaling_size(ax,S,p)
t=S.scalingByScale;
plot(ax,t.scale_factor,t.f_BIC_Hz/1e3,'o-', ...
    'Color',p.blue,'MarkerFaceColor',p.blue,'MarkerSize',3.5); hold(ax,'on');
plot(ax,t.scale_factor,200./t.scale_factor,'--','Color',p.orange);
xlabel(ax,'uniform scale factor s'); ylabel(ax,'retuned f_{BIC} (kHz)'); grid(ax,'on');
legend(ax,'modal operator','f_0/s','Location','northeast');
title(ax,'geometric similarity');
end

function draw_retuned_residual(ax,S,p)
t=S.scalingBySpeed;
semilogy(ax,t.sound_speed_m_per_s,t.strict_residual,'o-', ...
    'Color',p.green,'MarkerFaceColor',p.green,'MarkerSize',3.4); hold(ax,'on');
semilogy(ax,t.sound_speed_m_per_s,t.sigma_ratio,'s--', ...
    'Color',p.purple,'MarkerFaceColor','w','MarkerSize',3.4);
xlabel(ax,'sound speed c (m/s)'); ylabel(ax,'homogeneous diagnostic'); grid(ax,'on');
legend(ax,'raw residual','singular-value ratio','Location','northwest');
title(ax,'retuned strict null');
end

function draw_fixed_environment(ax,S,p)
t=S.environmentControl;
yyaxis(ax,'left');
semilogy(ax,t.sound_speed_m_per_s,t.full_operator_raw_residual,'o-', ...
    'Color',p.red,'MarkerFaceColor',p.red,'MarkerSize',3.4);
ylabel(ax,'full-operator raw residual');
yyaxis(ax,'right');
plot(ax,t.sound_speed_m_per_s,t.open_amplitude_ratio,'s-', ...
    'Color',p.blue,'MarkerFaceColor',p.blue,'MarkerSize',3.4);
ylabel(ax,'open-amplitude ratio');
xline(ax,1500,':','Color',p.dark,'LineWidth',0.9);
xlabel(ax,'sound speed c (m/s)'); grid(ax,'on');
legend(ax,'full residual','open amplitude','nominal c','Location','northwest');
title(ax,'fixed 200-kHz device: not a BIC off nominal c');
end

function draw_tolerance_cdf(ax,R,p)
hold(ax,'on'); modes={'independent','joint'}; cols=[p.blue;p.orange];
for im=1:2
    mask=strcmp({R.toleranceRows.mode},modes{im}) & ...
        [R.toleranceRows.level_pct]==5 & [R.toleranceRows.solver_ok];
    vals=sort([R.toleranceRows(mask).strict_raw_residual]);
    plot(ax,vals,(1:numel(vals))/numel(vals),'-','Color',cols(im,:),'LineWidth',1.25);
end
set(ax,'XScale','log'); xlabel(ax,'strict raw residual at \pm5%');
ylabel(ax,'empirical CDF'); grid(ax,'on');
legend(ax,'one parameter','all five parameters','Location','northwest');
title(ax,'complete residual distributions');
end

function draw_tolerance_scatter(ax,R,p)
hold(ax,'on'); marks={'o','s','^'};
for il=1:numel(R.levels)
    mask=strcmp({R.toleranceRows.mode},'joint') & ...
        [R.toleranceRows.level_pct]==100*R.levels(il) & [R.toleranceRows.solver_ok];
    x=[R.toleranceRows(mask).threshold_amplitude_error];
    y=[R.toleranceRows(mask).finite_open_amplitude_error];
    loglog(ax,x,y,marks{il},'Color',p.blue,'MarkerSize',3.2, ...
        'MarkerFaceColor','none','LineWidth',0.65);
end
xlabel(ax,'|A_{-1}|/||A||'); ylabel(ax,'|A_0|/||A||'); grid(ax,'on');
legend(ax,'\pm1%','\pm2%','\pm5%','Location','southeast');
title(ax,'joint perturbations reopen both constraints');
end

function draw_tank_schematic(ax,p)
axis(ax,'off'); hold(ax,'on');
rectangle(ax,'Position',[0.03 0.08 0.94 0.82],'FaceColor',[0.91 0.96 0.98], ...
    'EdgeColor',p.blue,'LineWidth',0.9);
rectangle(ax,'Position',[0.36 0.25 0.40 0.07],'FaceColor',[0.72 0.74 0.77], ...
    'EdgeColor',p.dark,'LineWidth',0.7);
for j=0:8
    rectangle(ax,'Position',[0.38+0.040*j 0.26 0.017 0.045], ...
        'FaceColor','w','EdgeColor','none');
end
rectangle(ax,'Position',[0.09 0.56 0.10 0.17],'Curvature',0.15, ...
    'FaceColor',p.orange,'EdgeColor',p.dark);
quiver(ax,0.19,0.64,0.22,-0.24,0,'Color',p.dark,'LineWidth',1.2,'MaxHeadSize',0.6);
plot(ax,0.65+0.23*cosd(linspace(15,165,80)), ...
    0.40+0.23*sind(linspace(15,165,80)),'--','Color',p.red,'LineWidth',1.0);
plot(ax,0.65,0.63,'o','Color',p.red,'MarkerFaceColor',p.red,'MarkerSize',4);
text(ax,0.08,0.78,'source','Color',p.orange);
text(ax,0.48,0.19,'finite sample + rigid reference','Color',p.dark,'FontSize',7);
text(ax,0.68,0.77,'hydrophone arc / raster','Color',p.red,'FontSize',7);
text(ax,0.04,0.11,'PROPOSED; no measurements','Color',p.red,'FontSize',7.2);
xlim(ax,[0 1]); ylim(ax,[0 1]); title(ax,'water-tank geometry');
end

function draw_scan_plan(ax,p)
hold(ax,'on');
theta=linspace(3.7,7.0,120); f=linspace(196,205,120);
[T,F]=meshgrid(theta,f);
feature=exp(-((F-(DUMMY_RA(T)))./0.13).^2).*exp(-((T-5.25)./0.9).^2);
contourf(ax,theta,f,feature,20,'LineStyle','none');
plot(ax,theta,DUMMY_RA(theta),'--','Color','w','LineWidth',1.0);
plot(ax,5.251216,200,'o','Color',p.red,'MarkerFaceColor',p.red,'MarkerSize',5);
xlabel(ax,'|\theta| (deg)'); ylabel(ax,'f (kHz)');
text(ax,0.04,0.91,'coarse map -> fine threshold scan','Units','normalized', ...
    'Color','w','FontSize',7);
title(ax,'proposed angle-frequency scan');
end

function f=DUMMY_RA(theta)
% Schematic line passing through the final BIC; used only in the proposed
% experiment panel, not as a numerical data source.
f=200*(1+sind(5.251216))./(1+sind(theta));
end

function draw_floquet_extraction(ax,p)
axis(ax,'off'); hold(ax,'on');
x=[0.06 0.38 0.70]; labels={'complex p(x,y_0)','Floquet least squares','A_n, \eta_n'};
cols=[p.blue;p.orange;p.green];
for j=1:3
    rectangle(ax,'Position',[x(j) 0.40 0.24 0.22],'Curvature',0.04, ...
        'FaceColor',[0.97 0.98 0.99],'EdgeColor',cols(j,:),'LineWidth',1.0);
    text(ax,x(j)+0.12,0.51,labels{j},'HorizontalAlignment','center','FontSize',7.2);
    if j<3
        quiver(ax,x(j)+0.24,0.51,0.065,0,0,'Color',p.dark,'LineWidth',1.0,'MaxHeadSize',0.7);
    end
end
text(ax,0.50,0.25,'flat-plate calibration + measured incident field', ...
    'HorizontalAlignment','center','Color',p.dark,'FontSize',7.1);
text(ax,0.50,0.14,'separate grazing pressure from normal power', ...
    'HorizontalAlignment','center','Color',p.red,'FontSize',7.1);
xlim(ax,[0 1]); ylim(ax,[0 1]); title(ax,'proposed diffraction-amplitude extraction');
end

function draw_ringdown_plan(ax,p)
t=linspace(0,6,320); trace=exp(-0.58*t).*cos(2*pi*1.35*t);
plot(ax,t,trace,'-','Color',p.blue,'LineWidth',1.05); hold(ax,'on');
plot(ax,t,exp(-0.58*t),'--','Color',p.orange,'LineWidth',0.9);
plot(ax,t,-exp(-0.58*t),'--','Color',p.orange,'LineWidth',0.9);
xlabel(ax,'time after switch-off (arb. units)'); ylabel(ax,'normalized pressure');
text(ax,0.05,0.10,'Q_{ring}=\pi f\tau_A; schematic only', ...
    'Units','normalized','Color',p.red,'FontSize',7.1);
grid(ax,'on'); title(ax,'proposed finite-sample ring-down');
end

function write_manifest(outDir,inputs,P,E,R,S)
fileName=fullfile(outDir,'figure_manifest.txt');
fid=fopen(fileName,'w'); assert(fid>0,'Cannot write %s',fileName);
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'Audited PRL figure manifest\nGenerated by advanced_analysis/run_prl_completion_figures.m\n\n');
fprintf(fid,'Inputs:\n');
for j=1:numel(inputs), fprintf(fid,'  %s\n',inputs{j}); end
fprintf(fid,'\nMain figures:\n');
fprintf(fid,'fig1_structure_topology.pdf -- final geometry and exact channel topology.\n');
fprintf(fid,'fig2_strict_bic_evidence.pdf -- final-root pole, Q, strict amplitudes, reoptimized sequence.\n');
fprintf(fid,'fig3_homogeneous_fields.pdf -- strict/near/ordinary homogeneous fields and localization metrics.\n');
fprintf(fid,'fig4_dof_tolerance.pdf -- bounded cavity ablation, two phasor closures, tolerance statistics.\n');
fprintf(fid,'\nSupplementary figures:\n');
fprintf(fid,'figS1_extended_convergence.pdf -- root drift and fixed-geometry crosschecks.\n');
fprintf(fid,'figS2_scaling_environment.pdf -- exact similarity and fixed-device sound-speed control.\n');
fprintf(fid,'figS3_tolerance_statistics.pdf -- full one-parameter/joint tolerance statistics.\n');
fprintf(fid,'figS4_driven_scattering.pdf -- fixed-geometry driven audit; generated separately.\n');
fprintf(fid,'figS5_proposed_experiment.pdf -- proposed protocol; no experimental data.\n');
fprintf(fid,'\nPole slope %.9f; endpoint full sigma %.6e.\n',P.fitQ(1),P.fullEndpointSigma);
fprintf(fid,'Highest reoptimized root N/K=%d/%d; max parameter step %.6e.\n', ...
    E.rootTable.N(end),E.rootTable.K(end),E.rootTable.max_parameter_step(end));
fprintf(fid,'Tolerance seed %d; %d samples per mode/level; ensemble N/K=%d/%d.\n', ...
    R.seed,R.nMC,R.Ntol,R.Ktol);
fprintf(fid,'Similarity speed samples %d; scale samples %d.\n', ...
    height(S.scalingBySpeed),height(S.scalingByScale));
fprintf(fid,'The superseded intermediate-root pole cache and invalid router table are not used.\n');
end
