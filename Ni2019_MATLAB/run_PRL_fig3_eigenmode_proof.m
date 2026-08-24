%% PRL Figure 3 -- strict eigenmode proof of the acoustic Rayleigh BIC
clear; close all; clc;

thisFile = mfilename('fullpath');
solverDir = fileparts(thisFile);
repoDir = fileparts(solverDir);
figureDir = fullfile(repoDir,'arxiv_theory_paper','figures');
dataDir = fullfile(repoDir,'arxiv_theory_paper','figure_data');
if ~exist(figureDir,'dir'), mkdir(figureDir); end
addpath(solverDir);

S = load(fullfile(solverDir,'results','StrictRayleighBIC_200kHz_min1mm.mat'));
Rdata = load(fullfile(repoDir,'advanced_analysis','robustness_fields', ...
    'robustness_fields.mat'));
Tnull = readtable(fullfile(dataDir,'fig2_strict_rayleigh_line.csv'));
x = S.xFinal(:).';
kappaBIC = x(1); OmegaBIC = 1-kappaBIC;

% Recompute the audited high-order strict eigenvector used for the channel
% and square-integrability diagnostics.
cfg = struct('a',1,'lambda',1/OmegaBIC, ...
    'theta_i_deg',asind(kappaBIC/OmegaBIC), ...
    'depths',x(2:3),'widths',x(4:5),'gaps',x(6), ...
    'N',313,'K',39,'solve_scattering',false);
Rstrict = ni2019_strict_rayleigh_operator(cfg,'TargetOrder',-1, ...
    'KyTolerance',1e-6);
assert(Rstrict.sigma_ratio<2e-15);

% Analytic exterior-energy integral over one period.  Orthogonality of the
% Floquet factors eliminates cross terms after integration over x/a in [0,1].
A = Rstrict.mode.z_full(1:Rstrict.full_operator.N);
ky = Rstrict.full_operator.ky;
orders = Rstrict.full_operator.orders;
threshold = orders==-1;
finiteOpen = orders==0;
evanescent = ~(threshold|finiteOpen);
alpha = abs(imag(ky(evanescent)));
amp2 = abs(A(evanescent)).^2;
valid = alpha>1e-12;
alpha = alpha(valid); amp2 = amp2(valid);
H = linspace(0,20,301);
Estrict = zeros(size(H));
for j=1:numel(H)
    Estrict(j) = sum(amp2.*(1-exp(-2*alpha*H(j)))./(2*alpha));
end
controlAmplitude = 0.05;
Econtrol = Estrict + controlAmplitude^2*H;
writetable(table(H(:),Estrict(:),Econtrol(:), ...
    'VariableNames',{'H_over_a','strict_BIC','nonzero_RA_control'}), ...
    fullfile(dataDir,'fig3_exterior_norm.csv'));

%% Visual system
ink=[31,36,43]/255; mid=[101,113,128]/255;
gridColor=[218,223,228]/255; blue=[33,102,172]/255;
orange=[217,119,6]/255; red=[199,53,43]/255;
cyan=[53,154,178]/255; purple=[139,103,184]/255;

fig=figure('Visible','off','Color','w','Units','inches', ...
    'Position',[.5,.5,7.15,5.20]);
set(fig,'DefaultAxesFontName','Arial','DefaultTextFontName','Arial', ...
    'DefaultAxesFontSize',7.2,'DefaultTextFontSize',7.2, ...
    'DefaultAxesLineWidth',.72,'DefaultLineLineWidth',1.15);

%% (a) Isolated strict null on the Rayleigh line
axA=axes(fig,'Position',[.065,.630,.390,.285]);
dk=1e4*Tnull.deltaKappa;
semilogy(axA,dk,Tnull.sigmaRelative,'-','Color',blue,'LineWidth',1.25); hold(axA,'on');
plot(axA,dk(1:10:end),Tnull.sigmaRelative(1:10:end),'o', ...
    'Color',blue,'MarkerFaceColor','w','MarkerSize',3.0);
[~,i0]=min(abs(Tnull.deltaKappa));
plot(axA,dk(i0),Tnull.sigmaRelative(i0),'o','MarkerSize',5.8, ...
    'MarkerFaceColor',red,'MarkerEdgeColor','w','LineWidth',.8);
text(axA,.55,2.2e-12,'strict BIC','Color',red);
xlabel(axA,'(\kappa-\kappa_{BIC})\times10^4');
ylabel(axA,'\sigma_{min}/\sigma_{max}');
xlim(axA,[-4,4]); ylim(axA,[5e-16,8e-3]);
yticks(axA,[1e-15,1e-11,1e-7,1e-3]);
style_axes(axA,ink,gridColor,true); panel_label(axA,'a',ink);

%% (b,c) Exact two-groove radiation closures of the strict eigenmode
axB=axes(fig,'Position',[.525,.640,.190,.270]);
draw_phasor(axB,Rdata.phasors(orders_index(Rdata.phasorOrders,0),:), ...
    Rdata.totalPhasors(orders_index(Rdata.phasorOrders,0)), ...
    'open, n=0',cyan,purple,red,ink,mid,true);
panel_label(axB,'b',ink);

axC=axes(fig,'Position',[.775,.640,.190,.270]);
draw_phasor(axC,Rdata.phasors(orders_index(Rdata.phasorOrders,-1),:), ...
    Rdata.totalPhasors(orders_index(Rdata.phasorOrders,-1)), ...
    'Rayleigh, n=-1',cyan,purple,red,ink,mid,false);
panel_label(axC,'c',ink);

%% (d) Homogeneous strict eigenfield (no incident wave)
axD=axes(fig,'Position',[.065,.090,.555,.420]);
state=Rdata.fieldStates(1);
h=imagesc(axD,Rdata.xGrid,Rdata.yGrid,state.pmag); hold(axD,'on');
set(axD,'YDir','normal'); set(h,'AlphaData',isfinite(state.pmag));
colormap(axD,parula(256)); clim(axD,[0,1]);
axis(axD,'image'); xlim(axD,[0,1]); ylim(axD,[-.75,1.05]);
draw_groove_outlines(axD,x,ink);
xlabel(axD,'x/a'); ylabel(axD,'y/a');
text(axD,.03,.95,'homogeneous eigenfield','Units','normalized', ...
    'Color','w','FontWeight','bold','VerticalAlignment','top');
cb=colorbar(axD,'eastoutside'); cb.Label.String='|p| / max_{groove}|p|';
cb.FontName='Arial'; cb.FontSize=7;
style_axes(axD,ink,gridColor,false); panel_label(axD,'d',ink);

%% (e) Square-integrability test
axE=axes(fig,'Position',[.720,.120,.245,.345]);
plot(axE,H,Estrict,'-','Color',blue,'LineWidth',1.45); hold(axE,'on');
plot(axE,H,Econtrol,'--','Color',orange,'LineWidth',1.35);
plot(axE,[0,20],[Estrict(end),Estrict(end)],':','Color',mid,'LineWidth',.8);
xlabel(axE,'Exterior height, H/a');
ylabel(axE,'Exterior norm, N(H)');
legend(axE,{'strict BIC','A_{RA}=0.05 control'}, ...
    'Location','northwest','Box','off');
text(axE,11,Estrict(end)*1.16,'saturates','Color',blue);
text(axE,10,Econtrol(round(end*.72))+.007,'\propto H','Color',orange);
xlim(axE,[0,20]); ylim(axE,[0,max(Econtrol)*1.08]);
style_axes(axE,ink,gridColor,true); panel_label(axE,'e',ink);

%% Export
drawnow;
stem=fullfile(figureDir,'fig3_strict_eigenmode_proof_matlab');
exportgraphics(fig,[stem,'.pdf'],'ContentType','vector','BackgroundColor','white');
exportgraphics(fig,[stem,'.svg'],'ContentType','vector','BackgroundColor','white');
exportgraphics(fig,[stem,'.png'],'Resolution',400,'BackgroundColor','white');
exportgraphics(fig,[stem,'.tiff'],'Resolution',600,'BackgroundColor','white');
close(fig);
fprintf('Exported %s.[pdf|svg|png|tiff]\n',stem);

function id=orders_index(orders,target)
id=find(orders==target,1); assert(~isempty(id));
end

function draw_phasor(ax,z,total,titleText,c1,c2,red,ink,mid,showLegend)
axis(ax,'equal'); axis(ax,'off'); hold(ax,'on');
scale=max(abs(z)); zn=z/max(scale,eps);
plot(ax,[-1.1,1.1],[0,0],'-','Color',[.86,.88,.90],'LineWidth',.6, ...
    'HandleVisibility','off');
plot(ax,[0,0],[-1.1,1.1],'-','Color',[.86,.88,.90],'LineWidth',.6, ...
    'HandleVisibility','off');
quiver(ax,0,0,real(zn(1)),imag(zn(1)),0,'Color',c1,'LineWidth',1.8, ...
    'MaxHeadSize',.28,'HandleVisibility','off');
quiver(ax,real(zn(1)),imag(zn(1)),real(zn(2)),imag(zn(2)),0, ...
    'Color',c2,'LineWidth',1.8,'MaxHeadSize',.28,'HandleVisibility','off');
plot(ax,0,0,'o','MarkerSize',5.2,'MarkerFaceColor','w', ...
    'MarkerEdgeColor',red,'LineWidth',1.0,'HandleVisibility','off');
text(ax,0,1.23,titleText,'HorizontalAlignment','center','Color',ink, ...
    'FontSize',7);
rel=abs(total)/max(sum(abs(z)),eps);
text(ax,0,-1.28,sprintf('|\\Sigma|/\\Sigma|s_j| = %.1e',rel), ...
    'HorizontalAlignment','center','Color',red,'FontSize',7);
text(ax,1.03,.08,'Re','HorizontalAlignment','right','Color',mid);
text(ax,.06,1.02,'Im','VerticalAlignment','top','Color',mid);
h1=plot(ax,nan,nan,'-','Color',c1,'LineWidth',1.8,'DisplayName','Groove 1');
h2=plot(ax,nan,nan,'-','Color',c2,'LineWidth',1.8,'DisplayName','Groove 2');
if showLegend
    legend(ax,[h1,h2],{'Groove 1','Groove 2'},'Location','southoutside', ...
        'Orientation','horizontal','Box','off');
else
    set([h1,h2],'HandleVisibility','off');
end
xlim(ax,[-1.15,1.15]); ylim(ax,[-1.45,1.35]);
end

function draw_groove_outlines(ax,x,ink)
margin=(1-sum(x(4:6)))/2;
xl=[margin,margin+x(4)+x(6)];
for ell=1:2
    rectangle(ax,'Position',[xl(ell),-x(ell+1),x(ell+3),x(ell+1)], ...
        'EdgeColor',ink,'LineWidth',.75,'FaceColor','none');
end
plot(ax,[0,1],[0,0],'-','Color',ink,'LineWidth',.75);
end

function panel_label(ax,letter,ink)
text(ax,0,1.08,['(',letter,')'],'Units','normalized','FontWeight','bold', ...
    'FontSize',9,'Color',ink,'HorizontalAlignment','left', ...
    'VerticalAlignment','top','Clipping','off');
end

function style_axes(ax,ink,gridColor,useGrid)
set(ax,'Box','off','Layer','top','TickDir','out','TickLength',[.018,.018], ...
    'XColor',ink,'YColor',ink,'LineWidth',.72,'FontName','Arial', ...
    'FontSize',7.2,'XMinorGrid','off','YMinorGrid','off');
if useGrid
    grid(ax,'on'); set(ax,'GridColor',gridColor,'GridAlpha',.52);
else
    grid(ax,'off');
end
end
