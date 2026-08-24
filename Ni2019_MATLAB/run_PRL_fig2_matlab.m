%% PRL Figure 2 -- spectral formation of the underwater acoustic Rayleigh BIC
% Every quantitative panel is generated from the audited final geometry.
% The figure separates global spectroscopy, pole motion, linewidth collapse,
% and fixed-angle spectral evolution.  Strict channel cancellation is kept
% for Fig. 3 so that this figure tells one formation story.
clear; close all; clc;

thisFile = mfilename('fullpath');
solverDir = fileparts(thisFile);
repoDir = fileparts(solverDir);
resultDir = fullfile(solverDir,'results');
paperDir = fullfile(repoDir,'arxiv_theory_paper');
figureDir = fullfile(paperDir,'figures');
dataDir = fullfile(paperDir,'figure_data');
if ~exist(figureDir,'dir'), mkdir(figureDir); end
if ~exist(dataDir,'dir'), mkdir(dataDir); end
addpath(solverDir);

S = load(fullfile(resultDir,'StrictRayleighBIC_200kHz_min1mm.mat'));
M = load(fullfile(resultDir,'DenseFinalRayleighMap_241x321.mat'));
P = readtable(fullfile(repoDir,'advanced_analysis','final_root_pole', ...
    'final_root_pole_track.csv'));
x = S.xFinal(:).';
kappaBIC = x(1);
OmegaBIC = 1-kappaBIC;
thetaBIC = asind(kappaBIC/OmegaBIC);
fBIC = S.fTarget;
aPhysical = S.aPhysical;
cWater = S.cWater;

%% High-resolution fixed-angle cuts (MATLAB modal matching)
cutFile = fullfile(resultDir,'Fig2FixedAngleCuts_1201.mat');
angleOffsetDeg = [-0.30,-0.15,-0.05,0.05,0.15,0.30];
cutAnglesDeg = thetaBIC + angleOffsetDeg;
cutFrequencyHz = linspace(198.8e3,201.2e3,1201);
if exist(cutFile,'file')
    C = load(cutFile);
    sameGrid = isequal(size(C.cutAnglesDeg),size(cutAnglesDeg)) && ...
        max(abs(C.cutAnglesDeg-cutAnglesDeg))<1e-12 && ...
        isequal(size(C.cutFrequencyHz),size(cutFrequencyHz)) && ...
        max(abs(C.cutFrequencyHz-cutFrequencyHz))<1e-9;
else
    sameGrid = false;
end

if ~sameGrid
    nA = numel(cutAnglesDeg);
    nF = numel(cutFrequencyHz);
    nPoint = nA*nF;
    cutAbsAm1Flat = nan(nPoint,1);
    cutGrooveEnergyFlat = nan(nPoint,1);
    cutConditionFlat = nan(nPoint,1);
    fprintf('Computing Fig. 2 fixed-angle cuts: %d x %d points.\n',nA,nF);
    tic;
    parfor idx = 1:nPoint
        [ia,jf] = ind2sub([nA,nF],idx);
        R = solve_driven_point(cutAnglesDeg(ia),cutFrequencyHz(jf),x, ...
            aPhysical,cWater,61,9);
        idm1 = find(R.orders==-1,1);
        cutAbsAm1Flat(idx) = abs(R.A(idm1));
        cutGrooveEnergyFlat(idx) = norm(R.groove_surface_coefficients(:))^2;
        cutConditionFlat(idx) = R.condition_number;
    end
    elapsedSeconds = toc;
    cutAbsAm1 = reshape(cutAbsAm1Flat,[nA,nF]);
    cutGrooveEnergy = reshape(cutGrooveEnergyFlat,[nA,nF]);
    cutCondition = reshape(cutConditionFlat,[nA,nF]);
    cutRayleighHz = cWater./(aPhysical*(1+sind(cutAnglesDeg)));
    save(cutFile,'cutAnglesDeg','angleOffsetDeg','cutFrequencyHz', ...
        'cutAbsAm1','cutGrooveEnergy','cutCondition','cutRayleighHz', ...
        'elapsedSeconds','-v7.3');
    writematrix(cutFrequencyHz(:),fullfile(dataDir,'fig2e_frequency_hz.csv'));
    writematrix(cutAnglesDeg(:),fullfile(dataDir,'fig2e_angle_deg.csv'));
    writematrix(cutGrooveEnergy,fullfile(dataDir,'fig2e_groove_energy.csv'));
    writematrix(cutAbsAm1,fullfile(dataDir,'fig2e_abs_am1.csv'));
    fprintf('Fixed-angle cuts completed in %.2f s.\n',elapsedSeconds);
else
    cutAbsAm1 = C.cutAbsAm1;
    cutGrooveEnergy = C.cutGrooveEnergy;
    cutCondition = C.cutCondition;
    cutRayleighHz = C.cutRayleighHz;
end

assert(all(isfinite(M.absA(:))) && all(M.absA(:)>=0), ...
    'The dense driven map contains invalid amplitudes.');
assert(all(isfinite(cutGrooveEnergy(:))) && all(cutGrooveEnergy(:)>0), ...
    'The fixed-angle groove-response data are invalid.');

%% PRL visual system
ink = [31,36,43]/255;
mid = [101,113,128]/255;
gridColor = [218,223,228]/255;
blue = [33,102,172]/255;
orange = [217,119,6]/255;
red = [199,53,43]/255;
cutColors = [49,90,155; 60,133,170; 83,168,163; ...
    202,150,61; 215,99,73; 151,73,128]/255;

fig = figure('Visible','off','Color','w','Units','inches', ...
    'Position',[0.5,0.5,7.15,6.15]);
set(fig,'DefaultAxesFontName','Arial','DefaultTextFontName','Arial', ...
    'DefaultAxesFontSize',7.2,'DefaultTextFontSize',7.2, ...
    'DefaultAxesLineWidth',0.72,'DefaultLineLineWidth',1.15);

%% (a) Periodic two-groove acoustic metasurface
axA = axes(fig,'Position',[.028,.705,.292,.225]);
draw_structure(axA,x,ink,mid);
panel_label(axA,'a',ink);

%% (b) Global angle-frequency spectroscopy
axB = axes(fig,'Position',[.355,.720,.235,.195]);
logAm1 = log10(max(squeeze(M.absA(:,:,1)),1e-6));
imagesc(axB,M.thetaValues,M.frequencyValues/1e3,logAm1.');
set(axB,'YDir','normal');
colormap(axB,parula(256));
clim(axB,[-4.3,0.4]);
hold(axB,'on');
plot(axB,M.thetaValues,M.rayleighFrequency/1e3,'w--','LineWidth',1.25);

finitePole = P.delta_kappa<0 & isfinite(P.Re_Omega);
thetaPole = asind(P.kappa(finitePole)./P.Re_Omega(finitePole));
fPoleKHz = P.Re_Omega(finitePole)*cWater/aPhysical/1e3;
plot(axB,thetaPole,fPoleKHz,'-','Color',[1,.55,.20], ...
    'LineWidth',1.15);
plot(axB,thetaBIC,fBIC/1e3,'o','MarkerSize',5.8, ...
    'MarkerFaceColor',red,'MarkerEdgeColor','w','LineWidth',0.8);
text(axB,thetaBIC+0.10,fBIC/1e3+0.15,'BIC','Color','w', ...
    'FontSize',7.2,'HorizontalAlignment','left');
xlim(axB,[4.4,6.2]); ylim(axB,[198.5,201.5]);
xlabel(axB,'Angle, \theta (deg)'); ylabel(axB,'Frequency (kHz)');
cb = colorbar(axB,'eastoutside');
cb.Label.String = 'log_{10}|A_{-1}|';
cb.FontName = 'Arial'; cb.FontSize = 7.0;
style_axes(axB,ink,gridColor,false);
panel_label(axB,'b',ink);

%% (c) Specular background channel in the same driven scan
axC = axes(fig,'Position',[.685,.720,.235,.195]);
logA0 = log10(max(squeeze(M.absA(:,:,2)),1e-6));
imagesc(axC,M.thetaValues,M.frequencyValues/1e3,logA0.');
set(axC,'YDir','normal');
colormap(axC,parula(256));
clim(axC,[-0.18,0]);
hold(axC,'on');
plot(axC,M.thetaValues,M.rayleighFrequency/1e3,'w--','LineWidth',1.15);
plot(axC,thetaPole,fPoleKHz,'-','Color',[1,.55,.20],'LineWidth',1.05);
plot(axC,thetaBIC,fBIC/1e3,'o','MarkerSize',5.3, ...
    'MarkerFaceColor',red,'MarkerEdgeColor','w','LineWidth',0.75);
xlim(axC,[4.4,6.2]); ylim(axC,[198.5,201.5]);
xlabel(axC,'Angle, \theta (deg)'); ylabel(axC,'Frequency (kHz)');
cb0 = colorbar(axC,'eastoutside');
cb0.Label.String = 'log_{10}|A_0|';
cb0.FontName = 'Arial'; cb0.FontSize = 7.0;
style_axes(axC,ink,gridColor,false);
panel_label(axC,'c',ink);

%% (d) Real part of the same-geometry pole
axD = axes(fig,'Position',[.065,.395,.405,.225]);
validPole = P.delta_kappa<0 & isfinite(P.Re_Omega);
kPlot = linspace(min(P.kappa(validPole)),kappaBIC,400);
plot(axD,kPlot,1-kPlot,'--','Color',orange,'LineWidth',1.25); hold(axD,'on');
plot(axD,P.kappa(validPole),P.Re_Omega(validPole),'o-', ...
    'Color',blue,'MarkerFaceColor','w','MarkerEdgeColor',blue, ...
    'MarkerSize',3.3,'LineWidth',1.15);
plot(axD,kappaBIC,OmegaBIC,'o','MarkerSize',5.8, ...
    'MarkerFaceColor',red,'MarkerEdgeColor','w','LineWidth',0.8);
text(axD,kappaBIC-8e-5,OmegaBIC+7.5e-5,'BIC','Color',red, ...
    'HorizontalAlignment','right');
text(axD,kappaBIC-8.8e-4,OmegaBIC+8.7e-4,'Rayleigh line', ...
    'Color',orange,'HorizontalAlignment','left');
xlabel(axD,'Bloch wavenumber, \kappa');
ylabel(axD,'Re \Omega_p');
xlim(axD,[kappaBIC-1.03e-3,kappaBIC+3e-5]);
ylim(axD,[OmegaBIC-1.2e-4,OmegaBIC+1.14e-3]);
style_axes(axD,ink,gridColor,true);
panel_label(axD,'d',ink);

%% (e) Linewidth collapse and Q divergence
axE = axes(fig,'Position',[.555,.395,.355,.225]);
dk = abs(P.delta_kappa(validPole));
imOmega = abs(P.Im_Omega(validPole));
Q = P.Q(validPole);
[dk,order] = sort(dk);
imOmega = imOmega(order); Q = Q(order);
yyaxis(axE,'left');
plot(axE,dk,imOmega,'o-','Color',blue,'MarkerFaceColor','w', ...
    'MarkerSize',3.2,'LineWidth',1.15);
set(axE,'XScale','log','YScale','log');
ylabel(axE,'|Im \Omega_p|','Color',blue);
axE.YAxis(1).Color = blue;
yyaxis(axE,'right');
plot(axE,dk,Q,'s-','Color',orange,'MarkerFaceColor','w', ...
    'MarkerSize',3.1,'LineWidth',1.15);
set(axE,'XScale','log','YScale','log');
ylabel(axE,'Q','Color',orange);
axE.YAxis(2).Color = orange;
xlabel(axE,'|\kappa-\kappa_{BIC}|');
xlim(axE,[8e-6,1.2e-3]);
text(axE,1.25e-5,4.5e11,'Q \rightarrow \infty','Color',orange, ...
    'HorizontalAlignment','left');
style_axes(axE,ink,gridColor,true);
panel_label(axE,'e',ink);

%% (f) Fixed-angle spectral evolution toward the BIC
axF = axes(fig,'Position',[.065,.080,.845,.225]);
hold(axF,'on');
for ia = 1:numel(cutAnglesDeg)
    y = log10(cutGrooveEnergy(ia,:));
    plot(axF,cutFrequencyHz/1e3,y,'Color',cutColors(ia,:), ...
        'LineWidth',1.15,'DisplayName',sprintf('%+.2f^\\circ',angleOffsetDeg(ia)));
    xline(axF,cutRayleighHz(ia)/1e3,':','Color', ...
        0.72*cutColors(ia,:)+0.28,'LineWidth',0.72,'HandleVisibility','off');
end
xline(axF,fBIC/1e3,'--','Color',red,'LineWidth',1.05, ...
    'HandleVisibility','off');
text(axF,fBIC/1e3+0.012,max(log10(cutGrooveEnergy(:)))-0.12, ...
    'BIC limit','Color',red,'VerticalAlignment','top');
xlabel(axF,'Frequency (kHz)');
ylabel(axF,'log_{10} ||C_{surf}||_2^2');
xlim(axF,[198.8,201.2]);
xticks(axF,[199,199.5,200,200.5,201]);
lgd = legend(axF,'Location','northoutside','Orientation','horizontal', ...
    'NumColumns',6,'Box','off');
lgd.Title.String = '\theta-\theta_{BIC}';
style_axes(axF,ink,gridColor,true);
panel_label(axF,'f',ink);

%% Export
drawnow;
stem = fullfile(figureDir,'fig2_rayleigh_bic_formation_matlab');
exportgraphics(fig,[stem,'.pdf'],'ContentType','vector', ...
    'BackgroundColor','white');
exportgraphics(fig,[stem,'.svg'],'ContentType','vector', ...
    'BackgroundColor','white');
exportgraphics(fig,[stem,'.png'],'Resolution',400, ...
    'BackgroundColor','white');
exportgraphics(fig,[stem,'.tiff'],'Resolution',600, ...
    'BackgroundColor','white');
close(fig);
fprintf('Exported %s.[pdf|svg|png|tiff]\n',stem);

function R = solve_driven_point(thetaDeg,frequencyHz,x,aPhysical,cWater,N,K)
Omega = frequencyHz*aPhysical/cWater;
cfg = struct('a',1,'lambda',1/Omega,'theta_i_deg',thetaDeg, ...
    'depths',x(2:3),'widths',x(4:5),'gaps',x(6), ...
    'N',N,'K',K,'solve_scattering',true);
R = ni2019_modal_solver(cfg);
end

function panel_label(ax,letter,ink)
text(ax,0,1.085,['(',letter,')'],'Units','normalized', ...
    'FontWeight','bold','FontSize',9,'Color',ink, ...
    'HorizontalAlignment','left','VerticalAlignment','top', ...
    'Clipping','off');
end

function style_axes(ax,ink,gridColor,useGrid)
set(ax,'Box','off','Layer','top','TickDir','out','TickLength',[.018,.018], ...
    'XColor',ink,'YColor',ink,'LineWidth',.72,'FontName','Arial', ...
    'FontSize',7.2,'XMinorGrid','off','YMinorGrid','off');
if useGrid
    grid(ax,'on');
    set(ax,'GridColor',gridColor,'GridAlpha',.52,'GridLineStyle','-');
else
    grid(ax,'off');
end
end

function draw_structure(ax,x,ink,mid)
axis(ax,'equal'); axis(ax,'off'); hold(ax,'on');
nCell = 8;
zMax = 1.35;
thick = 0.74;
frontLeft = [77,126,145]/255;
frontRight = [111,94,146]/255;
topFront = [177,218,220]/255;
topBack = [210,200,229]/255;

% Shadow.
p = proj_points([0.10,0,-thick-0.08; nCell+0.15,0,-thick-0.08; ...
    nCell+0.33,0,-thick-0.14; 0.30,0,-thick-0.14]);
patch(ax,p(:,1),p(:,2),[.72,.77,.80],'EdgeColor','none','FaceAlpha',.20);

% Front material gradient.
for j = 1:32
    xa = nCell*(j-1)/32; xb = nCell*j/32;
    t = (j-.5)/32; col = (1-t)*frontLeft+t*frontRight;
    p = proj_points([xa,0,0; xb,0,0; xb,0,-thick; xa,0,-thick]);
    patch(ax,p(:,1),p(:,2),col,'EdgeColor','none');
end
% Top material gradient.
for j = 1:18
    za = zMax*(j-1)/18; zb = zMax*j/18;
    t = (j-.5)/18; col = (1-t)*topFront+t*topBack;
    p = proj_points([0,za,0; nCell,za,0; nCell,zb,0; 0,zb,0]);
    patch(ax,p(:,1),p(:,2),col,'EdgeColor','none');
end
% Side and outline.
p = proj_points([nCell,0,0; nCell,zMax,0; nCell,zMax,-thick; nCell,0,-thick]);
patch(ax,p(:,1),p(:,2),[.43,.41,.61],'EdgeColor',ink,'LineWidth',.7);
line_outline(ax,[0,0,0; nCell,0,0; nCell,zMax,0; 0,zMax,0; 0,0,0],ink,.72);
line_outline(ax,[0,0,0; 0,0,-thick; nCell,0,-thick; nCell,0,0],ink,.72);

w1 = x(4); w2 = x(5); gap = x(6); d1 = x(2); d2 = x(3);
margin = (1-w1-w2-gap)/2;
for cellId = 0:nCell-1
    x1a = cellId+margin; x1b = x1a+w1;
    x2a = x1b+gap; x2b = x2a+w2;
    draw_groove(ax,x1a,x1b,zMax,d1,ink);
    draw_groove(ax,x2a,x2b,zMax,d2,ink);
end
for xb = 1:nCell-1
    line_outline(ax,[xb,0,.018; xb,zMax,.018],[1,1,1],.55);
end

% Period marker.
u0 = proj_points([3.5,0,-.90]);
u1 = proj_points([4.5,0,-.90]);
plot(ax,[u0(1),u1(1)],[u0(2),u1(2)],'-','Color',mid,'LineWidth',.8);
plot(ax,u0(1),u0(2),'<','Color',mid,'MarkerFaceColor',mid,'MarkerSize',4);
plot(ax,u1(1),u1(2),'>','Color',mid,'MarkerFaceColor',mid,'MarkerSize',4);
text(ax,mean([u0(1),u1(1)]),u0(2),'a','FontName','Times New Roman', ...
    'FontAngle','italic','FontSize',8,'Color',ink,'BackgroundColor','w', ...
    'Margin',.2,'HorizontalAlignment','center','VerticalAlignment','middle');
xlim(ax,[-.15,nCell+.40*zMax+.18]); ylim(ax,[-.94,.24*zMax+.15]);
end

function draw_groove(ax,xa,xb,zMax,depth,ink)
% White aperture and front cut: the groove is a void, not another material.
p = proj_points([xa,0,.014; xb,0,.014; xb,zMax,.014; xa,zMax,.014]);
patch(ax,p(:,1),p(:,2),[.97,.985,.99],'EdgeColor',ink,'LineWidth',.42);
p = proj_points([xa,0,0; xb,0,0; xb,0,-depth; xa,0,-depth]);
patch(ax,p(:,1),p(:,2),'w','EdgeColor','none');
line_outline(ax,[xa,0,0; xa,0,-depth; xb,0,-depth; xb,0,0],ink,.62);
end

function line_outline(ax,xyz,color,lw)
p = proj_points(xyz);
plot(ax,p(:,1),p(:,2),'-','Color',color,'LineWidth',lw);
end

function uv = proj_points(xzy)
% Input columns are x, z, y.
uv = [xzy(:,1)+.40*xzy(:,2), xzy(:,3)+.24*xzy(:,2)];
end
