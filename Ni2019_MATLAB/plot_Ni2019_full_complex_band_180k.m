%% PRL-style full complex-band diagram for the 180-kHz Rayleigh BIC
clear; clc;
rootDir=fileparts(mfilename('fullpath'));
dataDir=fullfile(rootDir,'results','full_complex_band_180k');
figDir=fullfile(rootDir,'..','arxiv_theory_paper','figures');
if ~exist(figDir,'dir'), mkdir(figDir); end
D=load(fullfile(dataDir,'positive_pole_census.mat'));
T=load(fullfile(dataDir,'target_rayleigh_band.mat'));
S=load(fullfile(rootDir,'results','StrictRayleighBIC_180kHz_7p10deg_final.mat'));

navy=[20,73,122]/255; orange=[220,119,20]/255;
red=[211,62,50]/255; gray=[117,127,136]/255; light=[220,225,229]/255;
cmap=interp1([0,.48,1],[16,30,56;20,132,128;246,211,99]/255,linspace(0,1,256));

fig=figure('Color','w','Units','centimeters','Position',[2,2,18.3,11.7]);
ax=axes(fig,'Position',[.105,.135,.78,.80]); hold(ax,'on');
set(fig,'InvertHardcopy','off');

% Rayleigh/light lines in the first Brillouin zone.
k=linspace(-.5,.5,1001);
plot(ax,k,abs(k),'Color',gray,'LineWidth',.9,'LineStyle','--');
plot(ax,k,abs(k-1),'Color',orange,'LineWidth',1.0,'LineStyle','--');
plot(ax,k,abs(k+1),'Color',orange,'LineWidth',1.0,'LineStyle','--');

% Physical-sheet groove-pole census.  Thin lines show only locally connected
% segments; colored points carry the independently refined pole data.
rows=[];
for j=1:numel(D.tracks)
    tr=D.tracks(j);
    if numel(tr.kappa)<3, continue; end
    reAll=real(tr.Omega(:));
    % The target Rayleigh-sheet branch is plotted below from its dedicated
    % square-root continuation.  Suppress the duplicate physical-sheet
    % segment that the independent census rediscovers near the zone edge.
    isTargetDuplicate=min(tr.kappa)>.30 && ...
        median(reAll,'omitnan')>.88 && median(reAll,'omitnan')<.91;
    if isTargetDuplicate, continue; end
    for sign=[-1,1]
        kap=sign*tr.kappa(:); om=tr.Omega(:); q=abs(tr.Q(:));
        [kap,id]=sort(kap); om=om(id); q=q(id);
        visible=real(om)<=1.0005;
        kap=kap(visible); om=om(visible); q=q(visible);
        if isempty(kap), continue; end
        plot(ax,kap,real(om),'Color',[.68,.71,.73],'LineWidth',.55);
        bound=real(om)<abs(kap)-4e-4;
        if any(bound)
            scatter(ax,kap(bound),real(om(bound)),13,[.43,.48,.52],'filled', ...
                'MarkerEdgeColor','none');
        end
        radiative=~bound;
        if any(radiative)
            scatter(ax,kap(radiative),real(om(radiative)),13, ...
                clip_logq(q(radiative)),'filled','MarkerEdgeColor','none', ...
                'MarkerFaceAlpha',.92);
        end
        rows=[rows;repmat(j,numel(kap),1),kap,real(om),imag(om),q]; %#ok<AGROW>
    end
end

% Open symbols identify square-root endpoints; no line is drawn across them.
branchPoints=[.27,.7293;.42,.4215;.43,.5775];
for b=1:size(branchPoints,1)
    plot(ax,[branchPoints(b,1),-branchPoints(b,1)], ...
        [branchPoints(b,2),branchPoints(b,2)],'o','MarkerSize',4.3, ...
        'MarkerFaceColor','w','MarkerEdgeColor',gray,'LineWidth',.65);
end
text(ax,.285,.746,'branch point','Color',gray,'FontSize',6.8);

% Target n=-1 Rayleigh sheet and its reciprocal n=+1 partner.
valid=isfinite(T.kappa)&isfinite(T.Omega);
kt=T.kappa(valid); ot=T.Omega(valid); qt=T.Q(valid);
for sign=[-1,1]
    kap=sign*kt(:); om=ot(:); q=qt(:);
    [kap,id]=sort(kap); om=om(id); q=q(id);
    plot(ax,kap,real(om),'Color',navy,'LineWidth',1.45);
    scatter(ax,kap,real(om),16,clip_logq(q),'filled', ...
        'MarkerEdgeColor',navy,'LineWidth',.2);
    rows=[rows;repmat(100,numel(kap),1),kap,real(om),imag(om),q]; %#ok<AGROW>
end

kB=S.x(1); oB=1-kB;
plot(ax,[kB,-kB],[oB,oB],'o','MarkerSize',6.5,'MarkerFaceColor',red, ...
    'MarkerEdgeColor','w','LineWidth',.7);
text(ax,kB+.025,oB+.025,'Rayleigh BIC','Color',red,'FontSize',7.3, ...
    'HorizontalAlignment','left','VerticalAlignment','bottom');
text(ax,.29,.345,'guided mode','Color',[.37,.42,.46],'FontSize',7.0, ...
    'HorizontalAlignment','left');

% Minimal direct labels for channel thresholds.
text(ax,.39,.405,'n=0','Color',gray,'FontSize',7.1,'Rotation',42);
text(ax,.31,.705,'n=-1','Color',orange,'FontSize',7.1,'Rotation',-42);
text(ax,-.31,.705,'n=+1','Color',orange,'FontSize',7.1,'Rotation',42);

xlim(ax,[-.5,.5]); ylim(ax,[0,1.0]);
xticks(ax,-.5:.25:.5); yticks(ax,0:.2:1.0);
xlabel(ax,'normalized Bloch wave number, \kappa=k_xa/2\pi');
ylabel(ax,'normalized eigenfrequency, Re \Omega');
set(ax,'FontName','Arial','FontSize',7.5,'LineWidth',.75, ...
    'TickDir','out','Box','off','Layer','top','XGrid','on','YGrid','on', ...
    'GridColor',light,'GridAlpha',.62,'Color','w', ...
    'XColor',[.12,.14,.16],'YColor',[.12,.14,.16]);
colormap(ax,cmap); caxis(ax,[0,6]);
cb=colorbar(ax,'eastoutside'); cb.Ticks=[0,2,4,6];
cb.Label.String='log_{10} Q  (radiative modes)'; cb.FontName='Arial'; cb.FontSize=7.3;
cb.Color=[.12,.14,.16];

% Local inset: branch-point crossing and diverging lifetime.
inset=axes(fig,'Position',[.185,.645,.215,.215]); hold(inset,'on');
plot(inset,kt,real(ot),'Color',navy,'LineWidth',1.25);
plot(inset,kt,1-kt,'Color',orange,'LineWidth',1.0,'LineStyle','--');
scatter(inset,kt,real(ot),10,clip_logq(qt),'filled','MarkerEdgeColor','none');
plot(inset,kB,oB,'o','MarkerSize',5.2,'MarkerFaceColor',red, ...
    'MarkerEdgeColor','w','LineWidth',.6);
xlim(inset,[.055,.165]); ylim(inset,[.875,.905]);
xticks(inset,[.06,.11,.16]); yticks(inset,[.88,.89,.90]);
xlabel(inset,'\kappa'); ylabel(inset,'Re \Omega');
set(inset,'FontName','Arial','FontSize',6.7,'LineWidth',.65,'TickDir','out', ...
    'Box','on','Color','w','XColor',[.12,.14,.16],'YColor',[.12,.14,.16]);
text(inset,.065,.902,'target Rayleigh sheet','FontSize',6.8,'Color',navy);

% Export numerical points and vector/raster artwork.
bandTable=array2table(rows,'VariableNames', ...
    {'track_id','kappa','Omega_real','Omega_imag','Q'});
writetable(bandTable,fullfile(dataDir,'full_complex_band_points.csv'));
exportgraphics(fig,fullfile(figDir,'fig2b_full_complex_band.pdf'), ...
    'ContentType','vector');
exportgraphics(fig,fullfile(figDir,'fig2b_full_complex_band.png'), ...
    'Resolution',600);
savefig(fig,fullfile(figDir,'fig2b_full_complex_band.fig'));
close(fig);

function value=clip_logq(q)
value=log10(max(abs(q),1)); value=min(max(value,0),6);
end
