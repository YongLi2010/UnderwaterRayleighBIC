%% PRL Figure 4 -- driven scattering and observable acoustic fields
clear; close all; clc;

thisFile=mfilename('fullpath'); solverDir=fileparts(thisFile); repoDir=fileparts(solverDir);
resultDir=fullfile(solverDir,'results');
figureDir=fullfile(repoDir,'arxiv_theory_paper','figures');
dataDir=fullfile(repoDir,'arxiv_theory_paper','figure_data');
if ~exist(figureDir,'dir'), mkdir(figureDir); end
addpath(solverDir);

S=load(fullfile(resultDir,'StrictRayleighBIC_200kHz_min1mm.mat'));
x=S.xFinal(:).'; aPhysical=S.aPhysical; cWater=S.cWater;
kappaBIC=x(1); thetaBIC=asind(kappaBIC/(1-kappaBIC));

% Fully converged angle-frequency atlas.  The earlier lower-order maps are
% retained for discovery only; every quantitative scattering panel here
% uses the same N=313,K=39 truncation as the strict BIC verification.
mapFile=fullfile(resultDir,'ConvergedScatteringMap_N313K39_73x145.mat');
thetaMap=linspace(4.7,6.3,73);
frequencyMap=linspace(199.4e3,200.6e3,145);
if exist(mapFile,'file')
    HM=load(mapFile);
else
    nT=numel(thetaMap); nF=numel(frequencyMap); nP=nT*nF;
    etaM1Flat=nan(nP,1); eta0Flat=nan(nP,1);
    fprintf('Computing converged scattering atlas: %d x %d, N=313,K=39.\n',nT,nF);
    tic;
    parfor idx=1:nP
        [it,jf]=ind2sub([nT,nF],idx);
        R=solve_point(thetaMap(it),frequencyMap(jf),x,aPhysical,cWater,313,39);
        etaM1Flat(idx)=R.eta(R.orders==-1);
        eta0Flat(idx)=R.eta(R.orders==0);
    end
    elapsedSeconds=toc;
    etaM1=reshape(etaM1Flat,[nT,nF]); eta0=reshape(eta0Flat,[nT,nF]);
    rayleighFrequency=cWater./(aPhysical*(1+sind(thetaMap)));
    save(mapFile,'thetaMap','frequencyMap','etaM1','eta0', ...
        'rayleighFrequency','elapsedSeconds','-v7.3');
    HM=load(mapFile);
    fprintf('Converged atlas completed in %.2f s.\n',elapsedSeconds);
end

% Representative states: the strongest resolved local response in the
% high-order atlas and the audited off-Rayleigh anomalous-reflection point.
[~,imax]=max(HM.etaM1(:));
[itMax,jfMax]=ind2sub(size(HM.etaM1),imax);
points(1)=struct('name','resolved near-BIC response', ...
    'theta',HM.thetaMap(itMax),'frequency',HM.frequencyMap(jfMax));
points(2)=struct('name','anomalous reflection','theta',32.6328, ...
    'frequency',202.430e3);

for j=1:2
    Rh=solve_point(points(j).theta,points(j).frequency,x,aPhysical,cWater,313,39);
    ids=[find(Rh.orders==-1,1),find(Rh.orders==0,1),find(Rh.orders==1,1)];
    points(j).eta=Rh.eta(ids);
    points(j).A=Rh.A(ids);
end

% Exact homogeneous BIC field.  At the BIC itself a driven total field is not
% uniquely defined; the physically meaningful field is the null vector of the
% homogeneous strict operator, with both n=0 and n=-1 radiation amplitudes zero.
cfgBIC=struct('a',1,'lambda',1/(1-kappaBIC), ...
    'theta_i_deg',thetaBIC,'depths',x(2:3),'widths',x(4:5), ...
    'gaps',x(6),'N',313,'K',39,'solve_scattering',false);
RBIC=ni2019_strict_rayleigh_operator(cfgBIC,'TargetOrder',-1, ...
    'KyTolerance',1e-6);
assert(RBIC.sigma_ratio<2e-15,'Strict BIC eigenvector failed validation.');

xField=linspace(-3,3,721); yField=linspace(-.72,2.15,401);
[X,Y]=meshgrid(xField,yField);
pBIC=homogeneous_bic_field(RBIC,xField,yField,kappaBIC);
pScale=max(abs(pBIC(:)),[],'omitnan');
pBIC=pBIC/max(pScale,eps);
reBIC=real(pBIC);
logAbsBIC=log10(max(abs(pBIC),1e-5));

% Source-data exports for the quantitative atlas.
writematrix(HM.thetaMap(:),fullfile(dataDir,'fig4_theta_deg.csv'));
writematrix(HM.frequencyMap(:),fullfile(dataDir,'fig4_frequency_hz.csv'));
writematrix(HM.etaM1,fullfile(dataDir,'fig4_eta_m1.csv'));
writematrix(HM.eta0,fullfile(dataDir,'fig4_eta_0.csv'));

% Machine-readable selected-state table.
etaSelected=reshape([points.eta],3,[]).';
stateTable=table(string({points.name}).', [points.theta].', ...
    [points.frequency].'/1e3,etaSelected(:,1),etaSelected(:,2),etaSelected(:,3), ...
    'VariableNames',{'state','theta_deg','frequency_kHz','eta_m1','eta_0','eta_p1'});
writetable(stateTable,fullfile(dataDir,'fig4_selected_states.csv'));

%% Visual system
ink=[31,36,43]/255; gridColor=[218,223,228]/255;
blue=[33,102,172]/255; orange=[217,119,6]/255; red=[199,53,43]/255;
teal=[38,148,139]/255; purple=[123,92,159]/255;

fig=figure('Visible','off','Color','w','Units','inches','Position',[.5,.5,7.15,6.25]);
set(fig,'DefaultAxesFontName','Arial','DefaultTextFontName','Arial', ...
    'DefaultAxesFontSize',7.2,'DefaultTextFontSize',7.2, ...
    'DefaultAxesLineWidth',.72,'DefaultLineLineWidth',1.15);

%% (a,b) Angle-frequency diffraction efficiencies
axA=axes(fig,'Position',[.072,.665,.348,.265]);
imagesc(axA,HM.thetaMap,HM.frequencyMap/1e3,HM.etaM1.');
set(axA,'YDir','normal'); colormap(axA,parula(256));
clim(axA,[0,max(HM.etaM1(:))]); hold(axA,'on');
plot(axA,HM.thetaMap,HM.rayleighFrequency/1e3,'w--','LineWidth',1.15);
plot(axA,thetaBIC,200,'o','MarkerSize',5.3,'MarkerFaceColor',red,'MarkerEdgeColor','w');
xlabel(axA,'Angle, \theta (deg)'); ylabel(axA,'Frequency (kHz)');
title(axA,'\eta_{-1}','FontWeight','normal');
cb=colorbar(axA,'eastoutside'); cb.FontSize=7;
style_axes(axA,ink,gridColor,false); panel_label(axA,'a',ink);

axB=axes(fig,'Position',[.555,.665,.360,.265]);
imagesc(axB,HM.thetaMap,HM.frequencyMap/1e3,HM.eta0.');
set(axB,'YDir','normal'); colormap(axB,parula(256));
clim(axB,[min(HM.eta0(:)),1]); hold(axB,'on');
plot(axB,HM.thetaMap,HM.rayleighFrequency/1e3,'w--','LineWidth',1.15);
plot(axB,thetaBIC,200,'o','MarkerSize',5.3,'MarkerFaceColor',red,'MarkerEdgeColor','w');
xlabel(axB,'Angle, \theta (deg)'); ylabel(axB,'Frequency (kHz)');
title(axB,'\eta_0','FontWeight','normal');
cb=colorbar(axB,'eastoutside'); cb.FontSize=7;
style_axes(axB,ink,gridColor,false); panel_label(axB,'b',ink);

%% (c) Energy conversion at a representative angle
axC=axes(fig,'Position',[.070,.405,.390,.175]);
[~,ith]=min(abs(HM.thetaMap-points(1).theta));
plot(axC,HM.frequencyMap/1e3,HM.etaM1(ith,:),'-', ...
    'Color',orange,'LineWidth',1.35); hold(axC,'on');
plot(axC,HM.frequencyMap/1e3,HM.eta0(ith,:),'-', ...
    'Color',blue,'LineWidth',1.35);
plot(axC,HM.frequencyMap/1e3,HM.etaM1(ith,:)+HM.eta0(ith,:), ...
    '--','Color',[.35,.38,.42], ...
    'LineWidth',.9);
xline(axC,points(1).frequency/1e3,':','Color',red,'LineWidth',1.0);
xlabel(axC,'Frequency (kHz)'); ylabel(axC,'Power fraction');
legend(axC,{'\eta_{-1}','\eta_0','closure'},'Location','best','Box','off');
xlim(axC,[199.4,200.6]); ylim(axC,[0,1.04]);
style_axes(axC,ink,gridColor,true); panel_label(axC,'c',ink);

%% (d) Floquet power partition at the two selected states
axD=axes(fig,'Position',[.585,.405,.315,.175]);
Ybar=etaSelected;
b=bar(axD,Ybar,'grouped','LineWidth',.55);
b(1).FaceColor=orange; b(2).FaceColor=blue; b(3).FaceColor=purple;
set(axD,'XTick',1:2,'XTickLabel',{'Rayleigh side','off-Rayleigh'}, ...
    'XTickLabelRotation',0);
ylabel(axD,'Diffraction efficiency'); ylim(axD,[0,1.04]);
legend(axD,{'n=-1','n=0','n=+1'},'Location','northoutside', ...
    'Orientation','horizontal','Box','off');
style_axes(axD,ink,gridColor,true); panel_label(axD,'d',ink);

%% (e,f) Exact homogeneous pressure field at the strict BIC
axE=axes(fig,'Position',[.055,.065,.425,.250]);
draw_bic_real_field(axE,reBIC,xField,yField,x,ink);
title(axE,sprintf('strict BIC: %.3f^\\circ, 200.000 kHz',thetaBIC), ...
    'FontWeight','normal');
style_axes(axE,ink,gridColor,false); panel_label(axE,'e',ink);

axF=axes(fig,'Position',[.550,.065,.425,.250]);
draw_bic_log_field(axF,logAbsBIC,xField,yField,x,ink);
title(axF,'evanescent localization','FontWeight','normal');
style_axes(axF,ink,gridColor,false); panel_label(axF,'f',ink);

%% Export
drawnow;
stem=fullfile(figureDir,'fig4_scattering_fields_matlab');
exportgraphics(fig,[stem,'.pdf'],'ContentType','vector','BackgroundColor','white');
exportgraphics(fig,[stem,'.svg'],'ContentType','vector','BackgroundColor','white');
exportgraphics(fig,[stem,'.png'],'Resolution',400,'BackgroundColor','white');
exportgraphics(fig,[stem,'.tiff'],'Resolution',600,'BackgroundColor','white');
close(fig);
fprintf('Selected state 1 eta[-1,0,+1] = %.6g %.6g %.6g\n',points(1).eta);
fprintf('Selected state 2 eta[-1,0,+1] = %.6g %.6g %.6g\n',points(2).eta);
fprintf('Exported %s.[pdf|svg|png|tiff]\n',stem);

function R=solve_point(thetaDeg,frequencyHz,x,aPhysical,cWater,N,K)
Omega=frequencyHz*aPhysical/cWater;
cfg=struct('a',1,'lambda',1/Omega,'theta_i_deg',thetaDeg, ...
    'depths',x(2:3),'widths',x(4:5),'gaps',x(6), ...
    'N',N,'K',K,'solve_scattering',true);
R=ni2019_modal_solver(cfg);
end

function draw_bic_real_field(ax,P,xg,yg,x,ink)
h=imagesc(ax,xg,yg,P); set(ax,'YDir','normal'); hold(ax,'on');
set(h,'AlphaData',isfinite(P)); colormap(ax,diverging_map(256)); clim(ax,[-1,1]);
draw_periodic_grooves(ax,x,-3,3,ink);
xlabel(ax,'x/a'); ylabel(ax,'y/a'); xlim(ax,[-3,3]); ylim(ax,[-.72,2.15]);
cb=colorbar(ax,'eastoutside'); cb.Label.String='Re(p)/max|p|'; cb.FontSize=7;
end

function draw_bic_log_field(ax,P,xg,yg,x,ink)
h=imagesc(ax,xg,yg,P); set(ax,'YDir','normal'); hold(ax,'on');
set(h,'AlphaData',isfinite(P)); colormap(ax,parula(256)); clim(ax,[-5,0]);
draw_periodic_grooves(ax,x,-3,3,ink);
xlabel(ax,'x/a'); ylabel(ax,'y/a'); xlim(ax,[-3,3]); ylim(ax,[-.72,2.15]);
cb=colorbar(ax,'eastoutside'); cb.Label.String='log_{10}|p|/max|p|'; cb.FontSize=7;
end

function p=homogeneous_bic_field(R,xg,yg,kappa)
[X,Y]=meshgrid(xg,yg); p=nan(size(X)); above=Y>=0;
pa=complex(zeros(nnz(above),1));
A=R.mode.A; op=R.full_operator;
for n=1:numel(op.orders)
    pa=pa+A(n).*exp(-1i*op.kx(n).*X(above)-1i*op.ky(n).*Y(above));
end
p(above)=pa;
C=reshape(R.mode.surface_coefficients,op.K,op.L);
for cellId=-3:2
    bloch=exp(-1i*2*pi*kappa*cellId);
    for ell=1:op.L
        xl=cellId+op.xleft(ell);
        inside=Y<0 & Y>=-op.depths(ell) & X>=xl & X<=xl+op.widths(ell);
        if ~any(inside,'all'), continue; end
        u=X(inside)-xl; yy=Y(inside); pg=complex(zeros(size(u)));
        for q=0:op.K-1
            alpha=q*pi/op.widths(ell);
            beta=groove_beta_local(op.k0^2-alpha^2);
            vertical=cos(beta*(yy+op.depths(ell)))/cos(beta*op.depths(ell));
            pg=pg+bloch*C(q+1,ell).*cos(alpha*u).*vertical;
        end
        p(inside)=pg;
    end
end
end

function b=groove_beta_local(z)
if real(z)>=0
    b=sqrt(real(z));
else
    b=1i*sqrt(-real(z));
end
end

function draw_periodic_grooves(ax,x,cmin,cmax,ink)
margin=(1-sum(x(4:6)))/2;
for cellId=cmin:cmax-1
    xl1=cellId+margin; xl2=xl1+x(4)+x(6);
    rectangle(ax,'Position',[xl1,-x(2),x(4),x(2)], ...
        'EdgeColor',ink,'FaceColor','none','LineWidth',.55);
    rectangle(ax,'Position',[xl2,-x(3),x(5),x(3)], ...
        'EdgeColor',ink,'FaceColor','none','LineWidth',.55);
end
plot(ax,[cmin,cmax],[0,0],'-','Color',ink,'LineWidth',.65);
end

function map=diverging_map(n)
anchors=[30,70,140; 91,150,196; 245,247,248; 220,112,83; 156,35,49]/255;
t=linspace(0,1,size(anchors,1)); tq=linspace(0,1,n);
map=[interp1(t,anchors(:,1),tq).',interp1(t,anchors(:,2),tq).', ...
    interp1(t,anchors(:,3),tq).'];
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
