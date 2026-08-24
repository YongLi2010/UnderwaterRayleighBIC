%% Scattering characterization around the strict single-Rayleigh BIC
% Produces experimentally interpretable angle-frequency maps, complex
% diffraction-amplitude cuts, and driven near-field maps.  The exact BIC is
% a dark homogeneous state, so the driven field is evaluated at nearby
% regular scattering points rather than at the singular BIC coordinate.
clear; close all; clc;

rootDir=fileparts(mfilename('fullpath'));
outputDir=fullfile(rootDir,'results');
if ~exist(outputDir,'dir'), mkdir(outputDir); end

S=load(fullfile(outputDir,'StrictRayleighBIC_optimized_practical.mat'));
% Use the (121,15) root because it is also the geometry used for the
% independent scattering-pole continuation.  Moderate truncation below is
% chosen for dense maps; representative points are re-evaluated at N=101.
x=S.xSequence(1,:);
kappaBIC=x(1);
OmegaBIC=1-kappaBIC;
thetaBIC=asind(kappaBIC/OmegaBIC);

mapN=81; mapK=11;
thetaValues=linspace(thetaBIC-.60,thetaBIC+.60,61);
OmegaValues=linspace(OmegaBIC-.0045,OmegaBIC+.0045,81);
nTheta=numel(thetaValues); nOmega=numel(OmegaValues);

absA=nan(nTheta,nOmega,3); phaseA=absA;
eta=nan(nTheta,nOmega,3); grooveNorm=nan(nTheta,nOmega);
conditionNumber=nan(nTheta,nOmega);

fprintf('Dense scattering map: %d angles x %d frequencies\n',nTheta,nOmega);
for it=1:nTheta
    for io=1:nOmega
        R=solve_point(thetaValues(it),OmegaValues(io),x,mapN,mapK);
        ids=[find(R.orders==-1,1),find(R.orders==0,1),find(R.orders==1,1)];
        absA(it,io,:)=abs(R.A(ids));
        phaseA(it,io,:)=angle(R.A(ids));
        eta(it,io,:)=R.eta(ids);
        grooveNorm(it,io)=norm(R.groove_surface_coefficients(:));
        conditionNumber(it,io)=R.condition_number;
    end
end

% Select a regular high-near-field point below the -1 opening.  A condition
% cap prevents a nearly singular numerical solve from being mistaken for a
% physically excitable field maximum.
[TH,OM]=ndgrid(thetaValues,OmegaValues);
rayleighOmega=1./(1+sind(TH));
regular=conditionNumber<1e9 & isfinite(grooveNorm);
below=OM<rayleighOmega-2e-6;
hotMetric=grooveNorm;
hotMetric(~(regular & below))=-inf;
[~,hotId]=max(hotMetric(:));
[hotIt,hotIo]=ind2sub(size(hotMetric),hotId);
hotSeed=[thetaValues(hotIt),OmegaValues(hotIo)];

% Select the strongest -1 diffraction point above its Rayleigh opening.
etaM=eta(:,:,1);
above=OM>rayleighOmega+2e-6;
routeMetric=etaM;
routeMetric(~(regular & above))=-inf;
[~,routeId]=max(routeMetric(:));
[routeIt,routeIo]=ind2sub(size(routeMetric),routeId);
routeSeed=[thetaValues(routeIt),OmegaValues(routeIo)];

% Adaptive local grids resolve narrow quasi-BIC features without requiring
% an impractically fine uniform two-dimensional map.
% The final point coordinates use the same (121,15) truncation as the root
% and pole track, so the quoted experimental coordinates are not shifted by
% the cheaper overview-map discretization.
hotPoint=refine_point(hotSeed,x,'groove',121,15,regular_bounds(thetaBIC,OmegaBIC));
routePoint=refine_point(routeSeed,x,'eta_minus1',121,15,regular_bounds(thetaBIC,OmegaBIC));

% An off-resonant reference at the same angle as the hot point makes the
% near-field comparison directly interpretable in a frequency sweep.
referencePoint=[hotPoint(1),OmegaBIC-.0040];
points=[referencePoint;hotPoint;routePoint];
pointNames={'off-resonant reference','near-field quasi-BIC','strongest -1 routing'};

%% Figure 1: angle-frequency amplitude, efficiency, and near-field maps
navy=[.035 .10 .30]; orange=[1 .36 .08]; red=[.90 .05 .06];
fig1=figure('Color','w','Position',[40 40 1420 780]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
mapTitles={'log_{10}|A_{-1}|','log_{10}|A_0|','log_{10}|A_{+1}|', ...
    '\eta_{-1}','\eta_0','log_{10}||C_{groove}||_2'};
mapData={log10(max(absA(:,:,1),1e-8)),log10(max(absA(:,:,2),1e-8)), ...
    log10(max(absA(:,:,3),1e-8)),eta(:,:,1),eta(:,:,2), ...
    log10(max(grooveNorm,1e-8))};
for j=1:6
    ax=nexttile;
    imagesc(thetaValues,OmegaValues,mapData{j}.'); axis xy; hold on;
    plot(thetaValues,1./(1+sind(thetaValues)),'w--','LineWidth',1.3);
    plot(thetaBIC,OmegaBIC,'p','MarkerSize',10,'MarkerFaceColor',red, ...
        'MarkerEdgeColor','w');
    plot(hotPoint(1),hotPoint(2),'o','MarkerSize',7,'MarkerFaceColor',orange, ...
        'MarkerEdgeColor','w');
    plot(routePoint(1),routePoint(2),'s','MarkerSize',7,'MarkerFaceColor',navy, ...
        'MarkerEdgeColor','w');
    xlabel('incidence angle \theta (deg)'); ylabel('\Omega=fa/c');
    title(mapTitles{j}); colorbar; box on;
    if j==4 || j==5, clim([0 1]); end
end
colormap(fig1,ni2019_viridis(256));
sgtitle(sprintf('Driven scattering around the single-Rayleigh BIC:  theta_B=%.6f deg, Omega_B=%.9f', ...
    thetaBIC,OmegaBIC),'Color','k');
style_figure(fig1);
mapFile=fullfile(outputDir,'SingleRayleighBIC_scattering_maps.png');
exportgraphics(fig1,mapFile,'Resolution',220);

%% Figure 2: experimental frequency and angle cuts of complex A_n
freqTheta=hotPoint(1);
freqCenter=hotPoint(2);
freqValues=linspace(freqCenter-5e-4,freqCenter+5e-4,801);
freqCut=evaluate_cut(repmat(freqTheta,size(freqValues)),freqValues,x,121,15);

angleOmega=routePoint(2);
angleCenter=routePoint(1);
angleValues=linspace(angleCenter-.12,angleCenter+.12,801);
angleCut=evaluate_cut(angleValues,repmat(angleOmega,size(angleValues)),x,121,15);

fig2=figure('Color','w','Position',[50 45 1450 800]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
plot_cut_row(freqValues,freqCut,'Omega',sprintf('theta=%.6f deg',freqTheta), ...
    [OmegaBIC OmegaBIC],navy,orange,red);
plot_cut_row(angleValues,angleCut,'theta (deg)',sprintf('Omega=%.9f',angleOmega), ...
    [thetaBIC thetaBIC],navy,orange,red);
sgtitle('Complex diffraction amplitudes and power efficiencies on experimental cuts','Color','k');
style_figure(fig2);
cutFile=fullfile(outputDir,'SingleRayleighBIC_A_order_cuts.png');
exportgraphics(fig2,cutFile,'Resolution',220);

%% Figure 3: driven pressure and intensity fields at representative points
fieldN=121; fieldK=15;
xField=linspace(0,1,361);
yField=linspace(-max(x(2:3)),1.0,401);
fieldResults=cell(3,1); fields=cell(3,1);
fig3=figure('Color','w','Position',[35 80 1540 520]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for j=1:3
    fieldResults{j}=solve_point(points(j,1),points(j,2),x,fieldN,fieldK);
    fields{j}=ni2019_reconstruct_field(fieldResults{j},xField,yField,'total');
    F=fields{j};
    level=log10(max(abs(F.p),1e-4));
    ax=nexttile;
    imagesc(xField,yField,level); axis xy; hold on;
    clim([-3 2.5]);
    draw_grooves(fieldResults{j});
    stepX=22; stepY=24;
    qx=F.X(1:stepY:end,1:stepX:end);
    qy=F.Y(1:stepY:end,1:stepX:end);
    qIx=F.Ix(1:stepY:end,1:stepX:end);
    qIy=F.Iy(1:stepY:end,1:stepX:end);
    qScale=sqrt(qIx.^2+qIy.^2);
    qIx=qIx./max(qScale,1e-12); qIy=qIy./max(qScale,1e-12);
    quiver(qx,qy,qIx,qIy,.45,'w','LineWidth',.65);
    ids=[find(fieldResults{j}.orders==-1,1),find(fieldResults{j}.orders==0,1), ...
        find(fieldResults{j}.orders==1,1)];
    title(sprintf('%s\n\\theta=%.6f^\\circ, \\Omega=%.9f\n|A|=[%.3g, %.3g, %.3g]', ...
        pointNames{j},points(j,1),points(j,2),abs(fieldResults{j}.A(ids))));
    xlabel('x/a'); ylabel('y/a'); box on;
end
colormap(fig3,ni2019_viridis(256));
cb=colorbar; cb.Layout.Tile='east'; cb.Label.String='log_{10}|p/p_{inc}|';
sgtitle('Driven total-pressure fields; arrows show normalized time-averaged intensity','Color','k');
style_figure(fig3);
fieldFile=fullfile(outputDir,'SingleRayleighBIC_driven_fields.png');
exportgraphics(fig3,fieldFile,'Resolution',220);

%% Machine-readable data and representative-point table
summary=point_table(points,pointNames,fieldResults);
csvFile=fullfile(outputDir,'SingleRayleighBIC_scattering_points.csv');
writetable(summary,csvFile);
matFile=fullfile(outputDir,'SingleRayleighBIC_scattering_characterization.mat');
save(matFile,'x','kappaBIC','OmegaBIC','thetaBIC','thetaValues','OmegaValues', ...
    'absA','phaseA','eta','grooveNorm','conditionNumber','hotPoint','routePoint', ...
    'points','pointNames','freqValues','freqCut','angleValues','angleCut','summary');

fprintf('\nSingle-Rayleigh scattering characterization complete\n');
fprintf('  BIC: theta=%.9f deg, Omega=%.12f\n',thetaBIC,OmegaBIC);
disp(summary(:,{'point','theta_deg','Omega','A_m1_abs','A_0_abs','A_p1_abs', ...
    'eta_m1','eta_0','eta_p1','groove_surface_norm'}));
fprintf('  %s\n  %s\n  %s\n  %s\n  %s\n',mapFile,cutFile,fieldFile,csvFile,matFile);

function R=solve_point(theta,Omega,x,N,K)
cfg=struct('a',1,'lambda',1/Omega,'theta_i_deg',theta, ...
    'depths',x(2:3),'widths',x(4:5),'gaps',x(6), ...
    'N',N,'K',K,'solve_scattering',true);
R=ni2019_modal_solver(cfg);
end

function bounds=regular_bounds(thetaBIC,OmegaBIC)
bounds=[thetaBIC-.70 thetaBIC+.70 OmegaBIC-.006 OmegaBIC+.006];
end

function point=refine_point(seed,x,metricName,N,K,bounds)
point=seed;
dTheta=.06; dOmega=6e-4;
for pass=1:4
    tv=linspace(max(bounds(1),point(1)-dTheta),min(bounds(2),point(1)+dTheta),25);
    ov=linspace(max(bounds(3),point(2)-dOmega),min(bounds(4),point(2)+dOmega),25);
    best=-inf; bestPoint=point;
    for it=1:numel(tv)
        for io=1:numel(ov)
            R=solve_point(tv(it),ov(io),x,N,K);
            if ~isfinite(R.condition_number) || R.condition_number>1e9, continue; end
            iM=find(R.orders==-1,1);
            rayleigh=1/(1+sind(tv(it)));
            switch metricName
                case 'groove'
                    if ov(io)>=rayleigh-2e-7, continue; end
                    value=norm(R.groove_surface_coefficients(:));
                case 'eta_minus1'
                    if ov(io)<=rayleigh+2e-7, continue; end
                    value=R.eta(iM);
                otherwise
                    error('Unknown refinement metric.');
            end
            if value>best, best=value; bestPoint=[tv(it),ov(io)]; end
        end
    end
    point=bestPoint;
    dTheta=dTheta/5; dOmega=dOmega/5;
end
end

function cut=evaluate_cut(theta,Omega,x,N,K)
n=numel(theta);
cut.A=complex(nan(n,3)); cut.eta=nan(n,3); cut.groove_norm=nan(n,1);
cut.condition_number=nan(n,1);
for j=1:n
    R=solve_point(theta(j),Omega(j),x,N,K);
    ids=[find(R.orders==-1,1),find(R.orders==0,1),find(R.orders==1,1)];
    cut.A(j,:)=R.A(ids); cut.eta(j,:)=R.eta(ids);
    cut.groove_norm(j)=norm(R.groove_surface_coefficients(:));
    cut.condition_number(j)=R.condition_number;
end
end

function plot_cut_row(xv,cut,xlabelText,subtitle,reference,navy,orange,red)
colors={navy,orange,red}; labels={'n=-1','n=0','n=+1'};
nexttile;
for n=1:3
    semilogy(xv,max(abs(cut.A(:,n)),1e-8),'Color',colors{n},'LineWidth',1.35); hold on;
end
xline(reference(1),'k--','BIC','LabelVerticalAlignment','bottom');
grid on; box on; xlabel(xlabelText); ylabel('|A_n|'); title(['amplitude, ' subtitle]);
legend(labels,'Location','best');

nexttile;
for n=1:3
    plot(xv,unwrap(angle(cut.A(:,n)))*180/pi,'Color',colors{n},'LineWidth',1.2); hold on;
end
xline(reference(1),'k--'); grid on; box on; xlabel(xlabelText);
ylabel('unwrapped phase (deg)'); title(['phase, ' subtitle]);

nexttile;
for n=1:3
    plot(xv,cut.eta(:,n),'Color',colors{n},'LineWidth',1.35); hold on;
end
xline(reference(1),'k--'); grid on; box on; xlabel(xlabelText);
ylabel('\eta_n'); ylim([0 1.02]); title(['normal-power efficiency, ' subtitle]);
end

function draw_grooves(R)
yline(0,'w-','LineWidth',1.1);
for ell=1:numel(R.widths)
    rectangle('Position',[R.xleft(ell),-R.depths(ell),R.widths(ell),R.depths(ell)], ...
        'EdgeColor','w','LineWidth',1.1);
end
end

function style_figure(fig)
axesList=findall(fig,'Type','axes');
for j=1:numel(axesList)
    ax=axesList(j);
    ax.Color='w'; ax.XColor='k'; ax.YColor='k';
    ax.GridColor=[.72 .72 .72]; ax.MinorGridColor=[.84 .84 .84];
    ax.Title.Color='k'; ax.XLabel.Color='k'; ax.YLabel.Color='k';
end
legends=findall(fig,'Type','legend');
for j=1:numel(legends)
    legends(j).Color='w'; legends(j).TextColor='k'; legends(j).EdgeColor=[.3 .3 .3];
end
end

function T=point_table(points,names,results)
n=size(points,1);
point=string(names(:)); theta_deg=points(:,1); Omega=points(:,2);
kappa=Omega.*sind(theta_deg);
A_m1=complex(zeros(n,1)); A_0=A_m1; A_p1=A_m1;
eta_m1=zeros(n,1); eta_0=eta_m1; eta_p1=eta_m1;
groove_surface_norm=zeros(n,1); condition_number=zeros(n,1);
for j=1:n
    R=results{j}; ids=[find(R.orders==-1,1),find(R.orders==0,1),find(R.orders==1,1)];
    A_m1(j)=R.A(ids(1)); A_0(j)=R.A(ids(2)); A_p1(j)=R.A(ids(3));
    eta_m1(j)=R.eta(ids(1)); eta_0(j)=R.eta(ids(2)); eta_p1(j)=R.eta(ids(3));
    groove_surface_norm(j)=norm(R.groove_surface_coefficients(:));
    condition_number(j)=R.condition_number;
end
T=table(point,theta_deg,Omega,kappa,real(A_m1),imag(A_m1),abs(A_m1),angle(A_m1)*180/pi, ...
    real(A_0),imag(A_0),abs(A_0),angle(A_0)*180/pi, ...
    real(A_p1),imag(A_p1),abs(A_p1),angle(A_p1)*180/pi, ...
    eta_m1,eta_0,eta_p1,groove_surface_norm,condition_number, ...
    'VariableNames',{'point','theta_deg','Omega','kappa', ...
    'A_m1_real','A_m1_imag','A_m1_abs','A_m1_phase_deg', ...
    'A_0_real','A_0_imag','A_0_abs','A_0_phase_deg', ...
    'A_p1_real','A_p1_imag','A_p1_abs','A_p1_phase_deg', ...
    'eta_m1','eta_0','eta_p1','groove_surface_norm','condition_number'});
end
