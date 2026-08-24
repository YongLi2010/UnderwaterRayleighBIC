%% PRL-level article figure plates for the 200 kHz acoustic Rayleigh BIC
% Ten self-contained, multi-panel article figures are generated.  Every
% driven-scattering and field panels use the final manufacturable w>=1 mm
% geometry.  Strict-mode panels use its polished N=313,K=39 root.  The
% moderate N=121,K=15 pole track is explicitly a convergence diagnostic.
% The figures distinguish reciprocal signed-order routing from a genuinely
% nonreciprocal one-way device.
clear; close all; clc;

rootDir=fileparts(mfilename('fullpath'));
resultDir=fullfile(rootDir,'results');
galleryDir=fullfile(resultDir,'PRL_gallery');
if ~exist(galleryDir,'dir'), mkdir(galleryDir); end

design=load(fullfile(resultDir,'StrictRayleighBIC_200kHz_min1mm.mat'));
poleData=load(fullfile(resultDir,'StrictRayleighBIC_200kHz_min1mm_poles.mat'));
singleData=load(fullfile(resultDir,'SingleGroove_strict_Rayleigh_search.mat'));

xPhysical=design.xFinal;
xPole=design.xSequence(1,:); % exact root used by the independent pole track
f0=design.fTarget; c0=design.cWater; aPhysical=design.aPhysical;
kappa0=xPhysical(1); Omega0=1-kappa0;
theta0=asind(kappa0/Omega0);
designPlot=design;

navy=[.035 .10 .30]; orange=[1 .36 .08]; red=[.90 .05 .06];
green=[.10 .55 .38]; gray=[.50 .50 .50]; cyan=[.12 .62 .72];
lightBlue=[.89 .96 .99]; solidGray=[.70 .72 .74];

cacheFile=fullfile(resultDir,'PRL_article_figure_data_v5.mat');
recalculate=strcmp(getenv('NI2019_RECALC_PRL_GALLERY'),'1') || ~exist(cacheFile,'file');
if recalculate
    fprintf('Computing PRL gallery data for the 200 kHz design...\n');

    %% Angle-frequency scattering atlas
    thetaValues=linspace(3.8,7.0,57);
    frequencyValues=linspace(196e3,205e3,73);
    nTheta=numel(thetaValues); nFrequency=numel(frequencyValues);
    absA=nan(nTheta,nFrequency,3); phaseA=absA; eta=absA;
    grooveNorm=nan(nTheta,nFrequency); conditionNumber=grooveNorm;
    for it=1:nTheta
        for jf=1:nFrequency
            R=solve_physical(thetaValues(it),frequencyValues(jf),xPhysical,aPhysical,c0,61,9);
            ids=order_ids(R);
            absA(it,jf,:)=abs(R.A(ids)); phaseA(it,jf,:)=angle(R.A(ids));
            eta(it,jf,:)=R.eta(ids);
            grooveNorm(it,jf)=norm(R.groove_surface_coefficients(:));
            conditionNumber(it,jf)=R.condition_number;
        end
    end
    [TH,FF]=ndgrid(thetaValues,frequencyValues);
    fRayleigh=c0./(aPhysical*(1+sind(TH)));
    regular=conditionNumber<1e9 & isfinite(grooveNorm);
    below=FF<fRayleigh-20; above=FF>fRayleigh+20;
    % High-truncation operating points.  The practical route point was
    % optimized at N=313,K=39 and independently retuned at N=345,K=43,
    % giving 30.6% and 30.8% conversion, respectively.  The closer point
    % illustrates the ideal high-Q/vanishing-bandwidth limit.
    hotPoint=[5.270442546972,200.005779468e3];
    routePoint=[6.000000000000,200.216874212e3];
    referencePoint=[hotPoint(1),199.600e3];
    routeConvergenceK=[15 23 35 39];
    routeConvergenceFrequency=[200.214750483 200.216070706 ...
        200.216755666 200.216874212];
    routeConvergenceEfficiency=[.307369214748 .306622452885 ...
        .306203460852 .306128029770];
    fieldPoints=[referencePoint;hotPoint;routePoint];

    %% Driven fields and strict homogeneous BIC field
    xField=linspace(0,1,321); yField=linspace(-max(xPhysical(2:3)),1,361);
    fieldResults=cell(3,1); drivenFields=cell(3,1);
    for j=1:3
        fieldResults{j}=solve_physical(fieldPoints(j,1),fieldPoints(j,2), ...
            xPhysical,aPhysical,c0,313,39);
        drivenFields{j}=ni2019_reconstruct_field(fieldResults{j},xField,yField,'total');
    end
    strictCfg=make_cfg(xPhysical,313,39);
    strictMode=ni2019_strict_rayleigh_operator(strictCfg,'TargetOrder',-1);
    strictResult=strict_to_field_result(strictMode);
    bicField=ni2019_reconstruct_field(strictResult,xField,yField,'scattered');

    %% Exact two-groove radiation phasors for A0 and the Rayleigh order
    % The full operator is [pressure-continuity rows; velocity-continuity
    % rows] acting on [A;C].  Decompose the velocity equation of each
    % outgoing order into the two groove contributions.
    nGroove=strictMode.full_operator.n_groove;
    nFloquet=strictMode.full_operator.n_floquet;
    velocityBlock=strictMode.F_full_raw(nGroove+1:end,nFloquet+1:end);
    Cscaled=strictMode.mode.C_scaled;
    phasorOrders=[0 -1]; phasors=complex(zeros(2,2));
    for n=1:2
        row=find(strictMode.full_operator.orders==phasorOrders(n),1);
        for ell=1:2
            cols=(ell-1)*strictMode.full_operator.K+(1:strictMode.full_operator.K);
            phasors(n,ell)=velocityBlock(row,cols)*Cscaled(cols);
        end
    end

    %% Strict-mode exterior energy and nonzero-grazing control
    heightValues=linspace(0,20,161);
    A=strictMode.mode.A; ky=strictMode.full_operator.ky;
    evanescent=imag(ky)<-1e-10;
    exteriorEnergy=zeros(size(heightValues));
    for n=find(evanescent).'
        alpha=-imag(ky(n));
        exteriorEnergy=exteriorEnergy+abs(A(n))^2*(1-exp(-2*alpha*heightValues))/(2*alpha);
    end
    grazingControl=exteriorEnergy+.05^2*heightValues;

    %% Fixed-angle spectral cuts, following the arXiv Fig. 2(c) layout
    cutTheta=theta0+[-.60 -.30 0 .30 .60];
    cutFrequency=linspace(198e3,202e3,361);
    cutA0Phase=nan(numel(cutTheta),numel(cutFrequency));
    cutEtaAnom=cutA0Phase;
    for it=1:numel(cutTheta)
        for jf=1:numel(cutFrequency)
            R=solve_physical(cutTheta(it),cutFrequency(jf),xPhysical, ...
                aPhysical,c0,61,9);
            i0=find(R.orders==0,1); im=find(R.orders==-1,1);
            cutA0Phase(it,jf)=angle(R.A(i0));
            cutEtaAnom(it,jf)=R.eta(im);
        end
        cutA0Phase(it,:)=unwrap(cutA0Phase(it,:));
    end
    cutRayleigh=c0./(aPhysical*(1+sind(cutTheta)));

    %% Direction-paired scattering at 1% above the BIC frequency
    directionFrequency=routePoint(2);
    directionTheta=linspace(routePoint(1)-.12,routePoint(1)+.12,481);
    etaPlus=zeros(size(directionTheta)); etaMinus=etaPlus;
    eta0Plus=etaPlus; eta0Minus=etaPlus; phasePlus=etaPlus; phaseMinus=etaPlus;
    ampPlus=etaPlus; ampMinus=etaPlus;
    for j=1:numel(directionTheta)
        Rp=solve_physical(directionTheta(j),directionFrequency,xPhysical,aPhysical,c0,313,39);
        Rm=solve_physical(-directionTheta(j),directionFrequency,xPhysical,aPhysical,c0,313,39);
        im=find(Rp.orders==-1,1); ip=find(Rm.orders==1,1);
        i0p=find(Rp.orders==0,1); i0m=find(Rm.orders==0,1);
        etaPlus(j)=Rp.eta(im); etaMinus(j)=Rm.eta(ip);
        eta0Plus(j)=Rp.eta(i0p); eta0Minus(j)=Rm.eta(i0m);
        ampPlus(j)=abs(Rp.A(im)); ampMinus(j)=abs(Rm.A(ip));
        phasePlus(j)=angle(Rp.A(im)); phaseMinus(j)=angle(Rm.A(ip));
    end
    reciprocityError=abs(etaPlus-etaMinus);

    % Full two-port power-normalized scattering matrices at the device point.
    RpRoute=solve_physical(routePoint(1),routePoint(2),xPhysical,aPhysical,c0,313,39);
    RmRoute=solve_physical(-routePoint(1),routePoint(2),xPhysical,aPhysical,c0,313,39);
    Splus=open_scattering_matrix(RpRoute,[0 -1]);
    Sminus=open_scattering_matrix(RmRoute,[0 1]);
    matrixReciprocityError=norm(Splus-Sminus.','fro');
    matrixUnitarityError=[norm(Splus'*Splus-eye(2),'fro'), ...
        norm(Sminus'*Sminus-eye(2),'fro')];

    %% Fabrication tolerance at the exact (121,15) root
    errorValuesMm=linspace(-.08,.08,31);
    toleranceSigma=nan(numel(errorValuesMm));
    for iw=1:numel(errorValuesMm)
        for id=1:numel(errorValuesMm)
            xt=xPole;
            aPole=(1-xPole(1))*c0/f0;
            xt(5)=xPole(5)+1e-3*errorValuesMm(iw)/aPole;
            xt(3)=xPole(3)+1e-3*errorValuesMm(id)/aPole;
            Rt=ni2019_strict_rayleigh_operator(make_cfg(xt,121,15),'TargetOrder',-1);
            toleranceSigma(id,iw)=Rt.sigma_ratio;
        end
    end

    %% Finite-array far-field envelope and output angle
    frequencyRoute=linspace(f0,210e3,201);
    sinOut= sind(theta0)-c0./(frequencyRoute*aPhysical);
    thetaOut=asind(max(-1,min(1,sinOut)));
    thetaFar=linspace(-90,-55,701);
    arrayPeriods=[10 20 40]; arrayFactor=nan(numel(arrayPeriods),numel(thetaFar));
    fArray=routePoint(2); kArray=2*pi*fArray/c0;
    targetSin=sind(routePoint(1))-c0/(fArray*aPhysical);
    for m=1:numel(arrayPeriods)
        psi=kArray*aPhysical*(sind(thetaFar)-targetSin);
        den=sin(psi/2); num=sin(arrayPeriods(m)*psi/2);
        af=ones(size(psi)); nz=abs(den)>1e-12; af(nz)=abs(num(nz)./(arrayPeriods(m)*den(nz))).^2;
        arrayFactor(m,:)=af;
    end

    %% Independent off-Rayleigh device point on the same unit cell
    % This operating point is deliberately kept separate from the
    % Rayleigh-BIC story: it is a high-efficiency two-channel reflection
    % point well above the n=-1 Rayleigh threshold, not a BIC-enhanced
    % state.  The N/K sweep is retained to expose modal-truncation error.
    offRayleighPoint=[32.6328,202.430e3];
    offRayleighNK=[81 11;121 15;161 20;201 25;313 39;401 50];
    offRayleighEta=nan(size(offRayleighNK,1),1);
    offRayleighEta0=offRayleighEta;
    offRayleighCondition=offRayleighEta;
    offRayleighEnergyError=offRayleighEta;
    for j=1:size(offRayleighNK,1)
        R=solve_physical(offRayleighPoint(1),offRayleighPoint(2),xPhysical, ...
            aPhysical,c0,offRayleighNK(j,1),offRayleighNK(j,2));
        im=find(R.orders==-1,1); i0=find(R.orders==0,1);
        offRayleighEta(j)=R.eta(im); offRayleighEta0(j)=R.eta(i0);
        offRayleighCondition(j)=R.condition_number;
        offRayleighEnergyError(j)=R.energy_error;
    end

    % Publication-resolution cuts use N=201,K=25; the convergence panel
    % reports the higher-truncation values independently.
    offRayleighFrequency=linspace(202.30e3,202.56e3,261);
    offRayleighEtaFrequency=nan(size(offRayleighFrequency));
    offRayleighEta0Frequency=offRayleighEtaFrequency;
    offRayleighOutputFrequency=offRayleighEtaFrequency;
    for j=1:numel(offRayleighFrequency)
        R=solve_physical(offRayleighPoint(1),offRayleighFrequency(j),xPhysical, ...
            aPhysical,c0,201,25);
        im=find(R.orders==-1,1); i0=find(R.orders==0,1);
        offRayleighEtaFrequency(j)=R.eta(im);
        offRayleighEta0Frequency(j)=R.eta(i0);
        offRayleighOutputFrequency(j)=R.theta_deg(im);
    end

    offRayleighAngle=linspace(25,40,301);
    offRayleighEtaAnglePlus=nan(size(offRayleighAngle));
    offRayleighEtaAngleMinus=offRayleighEtaAnglePlus;
    offRayleighEta0AnglePlus=offRayleighEtaAnglePlus;
    offRayleighEta0AngleMinus=offRayleighEtaAnglePlus;
    offRayleighOutputAnglePlus=offRayleighEtaAnglePlus;
    offRayleighOutputAngleMinus=offRayleighEtaAnglePlus;
    for j=1:numel(offRayleighAngle)
        Rp=solve_physical(offRayleighAngle(j),offRayleighPoint(2),xPhysical, ...
            aPhysical,c0,201,25);
        Rm=solve_physical(-offRayleighAngle(j),offRayleighPoint(2),xPhysical, ...
            aPhysical,c0,201,25);
        im=find(Rp.orders==-1,1); i0p=find(Rp.orders==0,1);
        ip=find(Rm.orders==1,1); i0m=find(Rm.orders==0,1);
        offRayleighEtaAnglePlus(j)=Rp.eta(im);
        offRayleighEtaAngleMinus(j)=Rm.eta(ip);
        offRayleighEta0AnglePlus(j)=Rp.eta(i0p);
        offRayleighEta0AngleMinus(j)=Rm.eta(i0m);
        offRayleighOutputAnglePlus(j)=Rp.theta_deg(im);
        offRayleighOutputAngleMinus(j)=Rm.theta_deg(ip);
    end

    save(cacheFile,'thetaValues','frequencyValues','absA','phaseA','eta', ...
        'grooveNorm','conditionNumber','hotPoint','routePoint','referencePoint', ...
        'routeConvergenceK','routeConvergenceFrequency','routeConvergenceEfficiency', ...
        'fieldPoints','fieldResults','drivenFields','strictMode','bicField', ...
        'phasorOrders','phasors','heightValues','exteriorEnergy','grazingControl', ...
        'cutTheta','cutFrequency','cutA0Phase','cutEtaAnom','cutRayleigh', ...
        'directionFrequency','directionTheta','etaPlus','etaMinus','eta0Plus', ...
        'eta0Minus','phasePlus','phaseMinus','ampPlus','ampMinus','reciprocityError', ...
        'Splus','Sminus','matrixReciprocityError','matrixUnitarityError', ...
        'errorValuesMm','toleranceSigma','frequencyRoute','thetaOut','thetaFar', ...
        'arrayPeriods','arrayFactor','fArray','offRayleighPoint','offRayleighNK', ...
        'offRayleighEta','offRayleighEta0','offRayleighCondition', ...
        'offRayleighEnergyError','offRayleighFrequency','offRayleighEtaFrequency', ...
        'offRayleighEta0Frequency','offRayleighOutputFrequency','offRayleighAngle', ...
        'offRayleighEtaAnglePlus','offRayleighEtaAngleMinus', ...
        'offRayleighEta0AnglePlus','offRayleighEta0AngleMinus', ...
        'offRayleighOutputAnglePlus','offRayleighOutputAngleMinus');
else
    load(cacheFile);
end

% Reconstruct this lightweight helper after either compute or cache load.
strictResult=strict_to_field_result(strictMode);

P=poleData.P;
figureFiles=strings(10,1);

%% FIGURE 1 - Geometry and single-Rayleigh channel topology
fig=paper_figure([40 40 1420 800]); tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; draw_unit_cell_mm(designPlot,lightBlue,solidGray,navy,orange,red); title('(a) Manufacturable acoustic unit cell');
nexttile; draw_rayleigh_topology(kappa0,Omega0,navy,orange,red,gray); title('(b) Rayleigh-channel topology');
nexttile; draw_channel_bars(kappa0,Omega0,navy,orange,red); title('(c) Channel status at the BIC pair');
nexttile; draw_signed_channel_cartoon(theta0,navy,orange,red); title('(d) Signed-order switching');
sgtitle('Figure 1 | Acoustic structure and the single-Rayleigh condition','Color','k');
style_figure(fig); figureFiles(1)=fullfile(galleryDir,'PRL_Fig01_geometry_channels.png'); exportgraphics(fig,figureFiles(1),'Resolution',220);

%% FIGURE 2 - Strict BIC evidence
fig=paper_figure([40 40 1420 800]); tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile;
plot(design.rootTruncations(:,2),design.xSequence(:,2),'-o','Color',navy,'LineWidth',1.6); hold on;
plot(design.rootTruncations(:,2),design.xSequence(:,3),'-s','Color',orange,'LineWidth',1.6);
plot(design.rootTruncations(:,2),design.xSequence(:,4),'-^','Color',green,'LineWidth',1.5);
plot(design.rootTruncations(:,2),design.xSequence(:,5),'-v','Color',red,'LineWidth',1.5);
plot(design.rootTruncations(:,2),design.xSequence(:,6),'-d','Color',gray,'LineWidth',1.4);
xlabel('groove truncation K'); ylabel('normalized geometry'); grid on; box on;
legend('d_1/a','d_2/a','w_1/a','w_2/a','g/a','Location','eastoutside'); title('(a) Root-parameter convergence');
nexttile;
semilogy(design.crossTruncations(:,2),design.crossSigma,'-s','Color',navy,'LineWidth',1.6); hold on;
semilogy(design.rootTruncations(:,2),design.rootSigma,'o','Color',red,'MarkerFaceColor',red);
yline(1e-8,'--','Color',gray); grid on; box on; xlabel('K'); ylabel('\sigma_{min}/\sigma_{max}');
legend('fixed geometry','reoptimized root','10^{-8}','Location','best'); title('(b) Strict compatibility');
nexttile;
orders=strictMode.full_operator.orders; id=abs(orders)<=5;
bar(orders(id),log10(max(abs(strictMode.mode.A(id)),1e-16)),'FaceColor',navy); hold on;
plot([-1 0],[-16 -16],'o','Color',red,'MarkerFaceColor',red,'MarkerSize',7);
xlabel('Floquet order n'); ylabel('log_{10}|A_n^{BIC}|'); ylim([-16 1]); grid on; box on;
title('(c) Homogeneous BIC radiation amplitudes');
nexttile;
plot(heightValues,exteriorEnergy,'Color',navy,'LineWidth',1.8); hold on;
plot(heightValues,grazingControl,'--','Color',orange,'LineWidth',1.8);
xlabel('exterior height H/a'); ylabel('modal exterior-energy proxy'); grid on; box on;
legend('strict BIC','nonzero-grazing control','Location','northwest'); title('(d) Square-integrability test');
sgtitle('Figure 2 | A_{-1}=A_0=0, not merely zero grazing flux','Color','k');
style_figure(fig); figureFiles(2)=fullfile(galleryDir,'PRL_Fig02_strict_bic_evidence.png'); exportgraphics(fig,figureFiles(2),'Resolution',220);

%% FIGURE 3 - Leaky pole, Rayleigh branch point, and Q scaling
fig=paper_figure([40 40 1450 780]); tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
valid=1:numel(P.kappa)-1;
nexttile;
plot(P.kappa,real(P.Omega),'-o','Color',navy,'LineWidth',1.6,'MarkerSize',3); hold on;
plot(P.kappa,1-P.kappa,'--','Color',orange,'LineWidth',1.5);
plot(P.kappa(end),real(P.Omega(end)),'p','Color',red,'MarkerFaceColor',red,'MarkerSize',10);
xlabel('\kappa'); ylabel('Re \Omega'); grid on; box on; legend('leaky pole','n=-1 Rayleigh line','BIC');
title('(a) Pole reaches the branch point');
nexttile;
yyaxis left; loglog(abs(P.delta_kappa(valid)),imag(P.Omega(valid)),'o-','Color',navy,'LineWidth',1.5); ylabel('Im \Omega');
yyaxis right; loglog(abs(P.delta_kappa(valid)),P.Q(valid),'s-','Color',orange,'LineWidth',1.5); ylabel('Q');
xlabel('|\Delta\kappa|'); grid on; box on; title('(b) Linewidth collapse and Q divergence');
nexttile;
cutColors=ni2019_viridis(numel(cutTheta));
for j=1:numel(cutTheta)
    plot(cutFrequency/1e3,cutA0Phase(j,:)*180/pi,'Color',cutColors(j,:),'LineWidth',1.35); hold on;
    xline(cutRayleigh(j)/1e3,':','Color',cutColors(j,:),'LineWidth',.9);
end
xlabel('frequency (kHz)'); ylabel('arg A_0 (deg)'); grid on; box on;
legend(compose('theta = %.2f deg',cutTheta),'Location','best');
title('(c) Fixed-angle specular-phase cuts');
nexttile;
for j=1:numel(cutTheta)
    plot(cutFrequency/1e3,cutEtaAnom(j,:),'Color',cutColors(j,:),'LineWidth',1.35); hold on;
    xline(cutRayleigh(j)/1e3,':','Color',cutColors(j,:),'LineWidth',.9);
end
xlabel('frequency (kHz)'); ylabel('\eta_{-1}'); grid on; box on; ylim([0 .36]);
title('(d) Fixed-angle anomalous-order cuts');
sgtitle('Figure 3 | Rayleigh branch point, pole collapse, and acoustic scattering cuts','Color','k');
style_figure(fig); figureFiles(3)=fullfile(galleryDir,'PRL_Fig03_pole_q_scaling.png'); exportgraphics(fig,figureFiles(3),'Resolution',220);

%% FIGURE 4 - Driven scattering atlas
fig=paper_figure([30 30 1550 850]); tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
data={log10(max(absA(:,:,1),1e-8)),log10(max(absA(:,:,2),1e-8)), ...
    log10(max(absA(:,:,3),1e-8)),eta(:,:,1),eta(:,:,2),log10(max(grooveNorm,1e-8))};
titles={'log_{10}|A_{-1}|','log_{10}|A_0|','log_{10}|A_{+1}|', ...
    '\eta_{-1}','\eta_0','log_{10}||C_{groove}||'};
fRA=c0./(aPhysical*(1+sind(thetaValues)));
for j=1:6
    nexttile; imagesc(thetaValues,frequencyValues/1e3,data{j}.'); axis xy; hold on;
    plot(thetaValues,fRA/1e3,'w--','LineWidth',1.2); plot(theta0,f0/1e3,'p','Color',red,'MarkerFaceColor',red,'MarkerSize',9);
    plot(hotPoint(1),hotPoint(2)/1e3,'o','MarkerFaceColor',orange,'MarkerEdgeColor','w');
    plot(routePoint(1),routePoint(2)/1e3,'s','MarkerFaceColor',navy,'MarkerEdgeColor','w');
    xlabel('\theta (deg)'); ylabel('f (kHz)'); title(titles{j}); colorbar; box on;
    if j==4 || j==5, clim([0 1]); end
end
colormap(fig,ni2019_viridis(256)); sgtitle('Figure 4 | A_n, power efficiencies, and near-field response','Color','k');
style_figure(fig); figureFiles(4)=fullfile(galleryDir,'PRL_Fig04_scattering_atlas.png'); exportgraphics(fig,figureFiles(4),'Resolution',220);

%% FIGURE 5 - Representative driven fields and the dark BIC eigenfield
fig=paper_figure([30 35 1500 820]); tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
fieldNames={'off resonance','ideal high-Q quasi-BIC','practical 30% Rayleigh router'};
for j=1:3
    nexttile; F=drivenFields{j}; level=log10(max(abs(F.p),1e-5)); imagesc(F.x,F.y,level); axis xy; hold on;
    draw_grooves(fieldResults{j},'w'); overlay_intensity(F,'w');
    ids=order_ids(fieldResults{j}); maxp=max(abs(F.p),[],'all','omitnan');
    title(sprintf('(%c) %s: theta=%.4f deg, f=%.3f kHz\nmax|p|=%.2g, |A|=[%.2g %.2g %.2g]', ...
        char('a'+j-1),fieldNames{j},fieldPoints(j,1),fieldPoints(j,2)/1e3,maxp,abs(fieldResults{j}.A(ids))));
    xlabel('x/a'); ylabel('y/a'); colorbar; box on;
end
nexttile; level=log10(max(abs(bicField.p)/max(abs(bicField.p),[],'all'),1e-6));
imagesc(bicField.x,bicField.y,level); axis xy; hold on; draw_grooves(strictResult,'w');
xlabel('x/a'); ylabel('y/a'); colorbar; box on; title('(d) dark homogeneous BIC eigenfield (normalized)');
colormap(fig,ni2019_viridis(256)); sgtitle('Figure 5 | From a driven quasi-BIC to the nonradiating eigenstate','Color','k');
style_figure(fig); figureFiles(5)=fullfile(galleryDir,'PRL_Fig05_driven_and_bic_fields.png'); exportgraphics(fig,figureFiles(5),'Resolution',220);

%% FIGURE 6 - Reciprocal momentum-channel locking
fig=paper_figure([40 40 1420 800]); tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; draw_signed_channel_cartoon(theta0,navy,orange,red); title('(a) Momentum sign selects the anomalous order');
nexttile;
plot(directionTheta,etaPlus,'-','Color',navy,'LineWidth',1.8); hold on;
plot(directionTheta,etaMinus,'--','Color',orange,'LineWidth',1.8);
xlabel('|\theta| (deg)'); ylabel('anomalous-order efficiency'); grid on; box on;
legend('+\theta: \eta_{-1}','-\theta: \eta_{+1}','Location','best'); title(sprintf('(b) Reciprocal pair at %.1f kHz',directionFrequency/1e3));
nexttile;
draw_scattering_matrix_pair(Splus,Sminus); title('(c) Power-normalized two-port matrices');
nexttile;
draw_reciprocity_residual(Splus,Sminus,matrixReciprocityError,matrixUnitarityError);
title('(d) Reciprocity and energy-conservation audit');
colormap(fig,ni2019_viridis(256));
sgtitle('Figure 6 | Reciprocal signed-order routing, not acoustic isolation','Color','k');
style_figure(fig); figureFiles(6)=fullfile(galleryDir,'PRL_Fig06_reciprocal_channel_locking.png'); exportgraphics(fig,figureFiles(6),'Resolution',220);

%% FIGURE 7 - Why the second cavity is essential
fig=paper_figure([40 40 1420 800]); tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; draw_ablation_schematic(designPlot,lightBlue,solidGray,navy,red); title('(a) One storage cavity + one phase-trim cavity');
nexttile;
semilogy(singleData.rootTruncations(:,2),singleData.rootSigma,'-o','Color',navy,'LineWidth',1.7); hold on;
semilogy(design.rootTruncations(:,2),design.rootSigma,'-s','Color',red,'LineWidth',1.7);
yline(1e-8,'--','Color',gray); xlabel('K'); ylabel('\sigma_{min}/\sigma_{max}'); grid on; box on;
legend('single wide cavity','two cavities','strict criterion','Location','best'); title('(b) High transverse modes do not replace cavity 2');
nexttile; draw_phasors(phasors(1,:),navy,orange,red); title('(c) A_0 radiation cancellation');
nexttile; draw_phasors(phasors(2,:),navy,orange,red); title('(d) A_{-1} grazing-amplitude cancellation');
sgtitle('Figure 7 | The small cavity supplies an independent radiation phase','Color','k');
style_figure(fig); figureFiles(7)=fullfile(galleryDir,'PRL_Fig07_second_cavity_mechanism.png'); exportgraphics(fig,figureFiles(7),'Resolution',220);

%% FIGURE 8 - Manufacturability, tolerance, and environmental tuning
fig=paper_figure([40 40 1420 800]); tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile;
bar([designPlot.depthsMm(:),designPlot.widthsMm(:)]); hold on; yline(1,'--','1 mm','Color',red);
set(gca,'XTickLabel',{'groove 1','groove 2'}); ylabel('dimension (mm)'); grid on; box on;
legend('depth','width','Location','best'); title('(a) 200 kHz physical dimensions');
nexttile;
imagesc(errorValuesMm,errorValuesMm,log10(max(toleranceSigma,1e-16))); axis xy; hold on; plot(0,0,'p','Color',red,'MarkerFaceColor',red,'MarkerSize',9);
xlabel('\Delta w_2 (mm)'); ylabel('\Delta d_2 (mm)'); colorbar; title('(b) Fixed-operating-point tolerance'); box on;
nexttile;
cValues=linspace(1450,1530,161); fShift=f0*cValues/c0/1e3;
plot(cValues,fShift,'Color',navy,'LineWidth',1.8); hold on; plot(c0,f0/1e3,'o','Color',red,'MarkerFaceColor',red);
xlabel('water sound speed (m/s)'); ylabel('shifted BIC frequency (kHz)'); grid on; box on; title('(c) Environmental frequency tuning');
nexttile;
yyaxis left; plot(routeConvergenceK,routeConvergenceEfficiency,'-o','Color',navy,'LineWidth',1.6,'MarkerFaceColor',navy);
ylabel('optimized \eta_{-1} at theta=6 deg'); ylim([.29 .32]);
yyaxis right; plot(routeConvergenceK,routeConvergenceFrequency,'-s','Color',red,'LineWidth',1.5,'MarkerFaceColor',red);
ylabel('optimized frequency (kHz)'); xlabel('groove truncation K'); grid on; box on;
title('(d) Retuned device-point convergence');
colormap(fig,ni2019_viridis(256)); sgtitle('Figure 8 | A manufacturable underwater design and its tolerance budget','Color','k');
style_figure(fig); figureFiles(8)=fullfile(galleryDir,'PRL_Fig08_manufacturability_tolerance.png'); exportgraphics(fig,figureFiles(8),'Resolution',220);

%% FIGURE 9 - Underwater experiment and modal extraction
fig=paper_figure([40 40 1480 800]); tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; draw_experiment_schematic(designPlot,lightBlue,solidGray,navy,orange,red); title('(a) Water-tank measurement geometry');
nexttile; draw_scan_plan(theta0,f0,navy,orange,red); title('(b) Coarse-to-fine scan plan');
nexttile; draw_fourier_extraction(fieldResults{3},drivenFields{3},navy,orange,red); title('(c) Complex near field to A_n and \eta_n');
nexttile; draw_experiment_workflow(navy,orange,red,green); title('(d) PRL-grade evidence chain');
sgtitle('Figure 9 | Water-tank experiment and data-reduction protocol','Color','k');
style_figure(fig); figureFiles(9)=fullfile(galleryDir,'PRL_Fig09_underwater_experiment.png'); exportgraphics(fig,figureFiles(9),'Resolution',220);

%% FIGURE 10 - Near-unity off-Rayleigh operation on the same unit cell
fig=paper_figure([40 40 1450 800]); tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile;
thetaRegime=linspace(0,40,401);
fRayleighRegime=c0./(aPhysical*(1+sind(thetaRegime)))/1e3;
plot(thetaRegime,fRayleighRegime,'--','Color',orange,'LineWidth',1.6); hold on;
plot(theta0,f0/1e3,'p','Color',red,'MarkerFaceColor',red,'MarkerSize',10);
plot(routePoint(1),routePoint(2)/1e3,'s','Color',navy,'MarkerFaceColor',navy,'MarkerSize',8);
plot(offRayleighPoint(1),offRayleighPoint(2)/1e3,'d','Color',green,'MarkerFaceColor',green,'MarkerSize',9);
text(theta0-.7,f0/1e3-4,'strict BIC','Color',red,'FontSize',9, ...
    'HorizontalAlignment','right');
text(routePoint(1)+1,routePoint(2)/1e3+5,'near-RA: 30.6%','Color',navy,'FontSize',9);
offRayleighRA=c0/(aPhysical*(1+sind(offRayleighPoint(1))))/1e3;
text(offRayleighPoint(1)-8,offRayleighPoint(2)/1e3-8, ...
    sprintf('off-RA: 99.96%%\\newline\\Delta f=%.1f kHz',offRayleighPoint(2)/1e3-offRayleighRA), ...
    'Color',green,'FontSize',9);
xlabel('|\theta_i| (deg)'); ylabel('frequency (kHz)'); grid on; box on;
legend('n=-1 Rayleigh threshold','strict BIC','near-RA router','off-RA device','Location','southwest');
title('(a) Distinct operating regimes of the same unit cell'); xlim([0 40]); ylim([130 220]);
nexttile;
convX=1:size(offRayleighNK,1);
plot(convX,offRayleighEta,'-o','Color',green,'LineWidth',1.8,'MarkerFaceColor',green); hold on;
plot(convX,offRayleighEta0,'-s','Color',gray,'LineWidth',1.6,'MarkerFaceColor',gray);
yline(.99,'--','Color',orange,'LineWidth',1.1);
set(gca,'XTick',convX,'XTickLabel',compose('%d/%d',offRayleighNK(:,1),offRayleighNK(:,2)));
xlabel('modal truncation N/K'); ylabel('normal-power efficiency'); grid on; box on;
legend('target anomalous order','specular n=0','99% reference','Location','best'); ylim([0 1.05]);
title(sprintf('(b) Convergence to \\eta_{-1}=%.4f (N=%d,K=%d)',offRayleighEta(end),offRayleighNK(end,1),offRayleighNK(end,2)));
nexttile;
plot(offRayleighFrequency/1e3,offRayleighEtaFrequency,'Color',green,'LineWidth',1.8); hold on;
plot(offRayleighFrequency/1e3,offRayleighEta0Frequency,'Color',gray,'LineWidth',1.6);
yline(.90,':','Color',orange,'HandleVisibility','off'); yline(.99,'--','Color',orange,'HandleVisibility','off');
i90=find(offRayleighEtaFrequency>=.90); i99=find(offRayleighEtaFrequency>=.99);
f90=offRayleighFrequency([i90(1) i90(end)])/1e3;
f99=offRayleighFrequency([i99(1) i99(end)])/1e3;
xline(f90(1),':','Color',gray,'HandleVisibility','off'); xline(f90(2),':','Color',gray,'HandleVisibility','off');
xline(f99(1),'--','Color',gray,'HandleVisibility','off'); xline(f99(2),'--','Color',gray,'HandleVisibility','off');
[~,idFreq]=min(abs(offRayleighFrequency-offRayleighPoint(2)));
plot(offRayleighPoint(2)/1e3,offRayleighEtaFrequency(idFreq),'d','Color',red,'MarkerFaceColor',red,'MarkerSize',9);
xlabel('frequency (kHz)'); ylabel('normal-power efficiency'); grid on; box on; ylim([0 1.05]);
legend('target n=-1','specular n=0','operating point','Location','best');
text(202.385,.54,sprintf('\\eta>0.90: %.3f kHz',diff(f90)),'Color',orange,'FontSize',9);
text(202.425,.42,sprintf('\\eta>0.99: %.3f kHz',diff(f99)),'Color',orange,'FontSize',9);
title('(c) Narrowband high-efficiency reflection (N=201,K=25)');
nexttile;
plot(offRayleighAngle,offRayleighEtaAnglePlus,'Color',navy,'LineWidth',1.8); hold on;
plot(offRayleighAngle,offRayleighEtaAngleMinus,'--','Color',orange,'LineWidth',1.8);
plot(offRayleighAngle,offRayleighEta0AnglePlus,':','Color',gray,'LineWidth',1.4);
plot(offRayleighAngle,offRayleighEta0AngleMinus,':','Color',gray,'LineWidth',1.4,'HandleVisibility','off');
yline(.90,':','Color',orange,'HandleVisibility','off'); yline(.99,'--','Color',orange,'HandleVisibility','off');
xline(offRayleighPoint(1),'--','Color',red,'HandleVisibility','off');
[~,idAngle]=min(abs(offRayleighAngle-offRayleighPoint(1)));
plot(offRayleighPoint(1),offRayleighEtaAnglePlus(idAngle),'d','Color',red,'MarkerFaceColor',red,'MarkerSize',9);
xlabel('|\theta_i| (deg)'); ylabel('normal-power efficiency'); grid on; box on; ylim([0 1.05]);
legend('+\theta_i: target n=-1','-\theta_i: target n=+1','specular n=0','Location','best');
text(25.4,.77,sprintf('+\\theta_i -> n=-1, \\theta_o=%.2f deg',offRayleighOutputAnglePlus(idAngle)), ...
    'Color',navy,'FontSize',9);
text(25.4,.67,sprintf('-\\theta_i -> n=+1, \\theta_o=+%.2f deg',abs(offRayleighOutputAngleMinus(idAngle))), ...
    'Color',orange,'FontSize',9);
title(sprintf('(d) Reciprocal negative reflection at %.3f kHz',offRayleighPoint(2)/1e3));
sgtitle('Figure 10 | Same unit cell: near-unity off-Rayleigh anomalous reflection','Color','k');
style_figure(fig); figureFiles(10)=fullfile(galleryDir,'PRL_Fig10_device_operation.png'); exportgraphics(fig,figureFiles(10),'Resolution',220);

%% Contact sheet and figure index
fig=figure('Color','w','Position',[20 20 1800 2200]); tiledlayout(5,2,'TileSpacing','compact','Padding','compact');
for j=1:10
    nexttile; im=imread(figureFiles(j)); image(im); axis image off; title(sprintf('Article Figure %d',j),'Color','k');
end
contactFile=fullfile(galleryDir,'PRL_Figure_gallery_contact_sheet.png'); exportgraphics(fig,contactFile,'Resolution',150);

story={ ...
    'Geometry and single-Rayleigh channel topology'; ...
    'Strict zero-radiation and square-integrability evidence'; ...
    'Leaky pole, Rayleigh branch point, Q scaling, and spectral cuts'; ...
    'Theory atlas for water-tank A_n, efficiencies, and near field'; ...
    'Off-resonant, ideal-qBIC, practical-router, and dark eigenfields'; ...
    'Reciprocal signed-order locking and full scattering-matrix audit'; ...
    'Independent phase role of the second cavity'; ...
    'Manufacturability, tolerance, sound-speed drift, and route convergence'; ...
    'Underwater experiment and modal extraction'; ...
    'Near-unity off-Rayleigh anomalous reflection and reciprocal angle window'};
figureIndex=table((1:10).',story,figureFiles,'VariableNames',{'figure','story','file'});
writetable(figureIndex,fullfile(galleryDir,'PRL_figure_index.csv'));
routeTable=table(routeConvergenceK(:),routeConvergenceFrequency(:), ...
    routeConvergenceEfficiency(:),'VariableNames', ...
    {'K','optimized_frequency_kHz','eta_anomalous'});
writetable(routeTable,fullfile(galleryDir,'PRL_route_convergence.csv'));
save(fullfile(galleryDir,'PRL_figure_gallery.mat'),'figureIndex','contactFile','cacheFile');
fprintf('\nPRL article-figure gallery complete: %s\n',galleryDir);
disp(figureIndex);

%% Helpers
function cfg=make_cfg(x,N,K)
Omega=1-x(1);
cfg=struct('a',1,'lambda',1/Omega,'theta_i_deg',asind(x(1)/Omega), ...
    'depths',x(2:3),'widths',x(4:5),'gaps',x(6),'N',N,'K',K,'solve_scattering',false);
end

function R=solve_physical(theta,f,x,aPhysical,c,N,K)
Omega=f*aPhysical/c;
cfg=struct('a',1,'lambda',1/Omega,'theta_i_deg',theta, ...
    'depths',x(2:3),'widths',x(4:5),'gaps',x(6),'N',N,'K',K,'solve_scattering',true);
R=ni2019_modal_solver(cfg);
end

function ids=order_ids(R)
ids=[find(R.orders==-1,1),find(R.orders==0,1),find(R.orders==1,1)];
end

function S=open_scattering_matrix(R,orderList)
% Rebuild the pressure-to-pressure reflection matrix for arbitrary incoming
% open Floquet orders, then convert it to power-normalized port amplitudes.
n=numel(R.orders); ids=zeros(1,numel(orderList));
for j=1:numel(orderList), ids(j)=find(R.orders==orderList(j),1); end
E=complex(zeros(n,numel(ids)));
for j=1:numel(ids), E(ids(j),j)=1; end
ky=real(R.ky(ids));
A=R.system_matrix\(E*diag(2*ky))-E;
S=diag(sqrt(ky))*A(ids,:)*diag(1./sqrt(ky));
end

function point=refine_physical(seed,x,aPhysical,c,metric)
point=seed; dt=.10; df=300;
for pass=1:3
    tv=linspace(point(1)-dt,point(1)+dt,21); fv=linspace(point(2)-df,point(2)+df,21);
    best=-inf; bestPoint=point;
    for it=1:numel(tv)
        fRA=c/(aPhysical*(1+sind(tv(it))));
        for jf=1:numel(fv)
            R=solve_physical(tv(it),fv(jf),x,aPhysical,c,81,11);
            if ~isfinite(R.condition_number) || R.condition_number>1e9, continue; end
            switch metric
                case 'groove'
                    if fv(jf)>=fRA-5, continue; end
                    value=norm(R.groove_surface_coefficients(:));
                otherwise
                    if fv(jf)<=fRA+5, continue; end
                    value=R.eta(R.orders==-1);
            end
            if value>best, best=value; bestPoint=[tv(it),fv(jf)]; end
        end
    end
    point=bestPoint; dt=dt/5; df=df/5;
end
end

function out=strict_to_field_result(S)
op=S.full_operator;
out=struct('orders',op.orders,'A',S.mode.A,'eta',zeros(op.N,1), ...
    'theta_deg',nan(op.N,1),'is_propagating',false(op.N,1), ...
    'total_efficiency',0,'energy_error',0,'k0',op.k0,'kx',op.kx,'ky',op.ky, ...
    'ky_incident',op.k0*cosd(S.cfg.theta_i_deg),'a',op.a,'lambda',op.lambda, ...
    'theta_i_deg',S.cfg.theta_i_deg,'widths',op.widths,'depths',op.depths, ...
    'gaps',op.gaps,'xleft',op.xleft,'N',op.N,'K',op.K, ...
    'condition_number',nan,'system_matrix',[],'forcing',[], ...
    'groove_surface_coefficients',reshape(S.groove.surface_coefficients,op.K,op.L), ...
    'groove_mode_indices',(0:op.K-1).');
end

function fig=paper_figure(position)
fig=figure('Color','w','Position',position);
end

function style_figure(fig)
axesList=findall(fig,'Type','axes');
for j=1:numel(axesList)
    ax=axesList(j); ax.Color='w'; ax.XColor='k'; ax.YColor='k';
    ax.GridColor=[.72 .72 .72]; ax.MinorGridColor=[.84 .84 .84];
    ax.Title.Color='k'; ax.XLabel.Color='k'; ax.YLabel.Color='k';
end
legends=findall(fig,'Type','legend');
for j=1:numel(legends)
    legends(j).Color='w'; legends(j).TextColor='k'; legends(j).EdgeColor=[.3 .3 .3];
end
texts=findall(fig,'Type','text');
for j=1:numel(texts)
    color=texts(j).Color;
    if isnumeric(color) && numel(color)==3 && max(color)-min(color)<.04
        texts(j).Color='k';
    end
end
end

function draw_unit_cell_mm(D,water,solid,navy,orange,red)
a=1e3*D.aPhysical; d=D.depthsMm; w=D.widthsMm; g=D.gapMm;
occupied=sum(w)+g; x0=(a-occupied)/2; x2=x0+w(1)+g;
rectangle('Position',[0,0,a,3.2],'FaceColor',water,'EdgeColor','none'); hold on;
rectangle('Position',[0,-max(d)-.6,a,max(d)+.6],'FaceColor',solid,'EdgeColor','k');
rectangle('Position',[x0,-d(1),w(1),d(1)],'FaceColor','w','EdgeColor',navy,'LineWidth',1.5);
rectangle('Position',[x2,-d(2),w(2),d(2)],'FaceColor','w','EdgeColor',red,'LineWidth',1.5);
quiver(a*.52,2.7,-.55,-2.0,0,'Color',navy,'LineWidth',2,'MaxHeadSize',.5);
quiver(a*.50,.1,-1.5,.15,0,'Color',orange,'LineWidth',2,'MaxHeadSize',.5);
text(a*.58,2.55,'incident n=0','Color',navy); text(.15,.35,'n=-1 grazing','Color',orange);
text(x0+w(1)/2,-d(1)/2,'storage','HorizontalAlignment','center','Color',navy);
text(x2+w(2)/2,-d(2)/2,'phase trim','HorizontalAlignment','center','Rotation',90,'Color',red);
xlim([0 a]); ylim([-max(d)-.6 3.2]); axis equal; xlabel('x (mm)'); ylabel('y (mm)'); box on;
end

function draw_rayleigh_topology(k0,O0,navy,orange,red,gray)
k=linspace(-.32,.32,401); hold on;
plot(k,abs(k-1),'Color',orange,'LineWidth',1.8);
plot(k,abs(k),'Color',gray,'LineWidth',1.4);
plot(k,abs(k+1),'Color',navy,'LineWidth',1.8);
plot([k0 -k0],[O0 O0],'p','Color',red,'MarkerFaceColor',red,'MarkerSize',9);
fill([-.32 .32 .32 -.32],[.86 .86 1.08 1.08],[.94 .94 .94],'EdgeColor','none','FaceAlpha',.25);
xlabel('\kappa'); ylabel('\Omega=fa/c'); xlim([-.32 .32]); ylim([.70 1.15]); grid on; box on;
legend('n=-1 Rayleigh','n=0 light line','n=+1 Rayleigh','BIC pair','Location','south');
end

function draw_channel_bars(k,O,navy,orange,red)
orders=-1:1; kvalsPlus=k+orders; kvalsMinus=-k+orders;
kyPlus=sqrt(complex(O^2-kvalsPlus.^2)); kyMinus=sqrt(complex(O^2-kvalsMinus.^2));
vals=[real(kyPlus(:)),real(kyMinus(:))]; bar(orders,vals); hold on;
yline(0,'k-'); xlabel('order n'); ylabel('Re(k_{y,n})a/2\pi'); grid on; box on;
legend('+\theta','-\theta','Location','north');
text(-1,.04,'RA at +\theta','HorizontalAlignment','center','Color',orange);
text(1,.04,'RA at -\theta','HorizontalAlignment','center','Color',navy);
end

function draw_signed_channel_cartoon(theta,navy,orange,red)
axis off; hold on; plot([.05 .95],[.42 .42],'k-','LineWidth',2);
quiver(.28,.88,.12,-.40,0,'Color',navy,'LineWidth',2,'MaxHeadSize',.5);
quiver(.72,.88,-.12,-.40,0,'Color',orange,'LineWidth',2,'MaxHeadSize',.5);
quiver(.40,.44,-.33,.02,0,'Color',red,'LineWidth',2,'MaxHeadSize',.5);
quiver(.60,.44,.33,.02,0,'Color',red,'LineWidth',2,'MaxHeadSize',.5);
text(.16,.90,sprintf('+theta = %.3f deg',theta),'Color',navy);
text(.66,.90,sprintf('-theta = %.3f deg',theta),'Color',orange);
text(.05,.32,'selected n=-1 port','Color',red); text(.67,.32,'selected n=+1 port','Color',red);
text(.5,.12,'reciprocal pair: direction swaps, isolation does not occur','HorizontalAlignment','center');
xlim([0 1]); ylim([0 1]);
end

function draw_grooves(R,color)
yline(0,[color '-'],'LineWidth',1);
for ell=1:numel(R.widths)
    rectangle('Position',[R.xleft(ell),-R.depths(ell),R.widths(ell),R.depths(ell)], ...
        'EdgeColor',color,'LineWidth',1.1);
end
end

function overlay_intensity(F,color)
sy=28; sx=25; X=F.X(1:sy:end,1:sx:end); Y=F.Y(1:sy:end,1:sx:end);
Ix=F.Ix(1:sy:end,1:sx:end); Iy=F.Iy(1:sy:end,1:sx:end);
mag=sqrt(Ix.^2+Iy.^2); Ix=Ix./max(mag,1e-12); Iy=Iy./max(mag,1e-12);
quiver(X,Y,Ix,Iy,.42,'Color',color,'LineWidth',.55);
end

function draw_ablation_schematic(D,water,solid,navy,red)
axis off; hold on;
for row=1:2
    y=.72-(row-1)*.48; rectangle('Position',[.06,y-.18,.88,.18],'FaceColor',solid,'EdgeColor','k');
    rectangle('Position',[.18,y-.15,.47,.15],'FaceColor','w','EdgeColor',navy,'LineWidth',1.5);
    if row==1, rectangle('Position',[.76,y-.07,.12,.07],'FaceColor','w','EdgeColor',red,'LineWidth',1.5); end
end
text(.05,.93,'two cavities: exact A_0 and A_{RA} cancellation','Color',red);
text(.05,.45,'single wide cavity: transverse modes retained, residual remains finite','Color',navy);
xlim([0 1]); ylim([0 1]);
end

function draw_phasors(z,navy,orange,red)
scale=max(abs(z)); if scale==0, scale=1; end; z=z/scale;
quiver(0,0,real(z(1)),imag(z(1)),0,'Color',navy,'LineWidth',2,'MaxHeadSize',.35); hold on;
quiver(real(z(1)),imag(z(1)),real(z(2)),imag(z(2)),0,'Color',orange,'LineWidth',2,'MaxHeadSize',.35);
plot([0 real(sum(z))],[0 imag(sum(z))],'--','Color',red,'LineWidth',1.5);
plot(0,0,'ko','MarkerFaceColor','k'); plot(real(sum(z)),imag(sum(z)),'p','Color',red,'MarkerFaceColor',red);
axis equal; lim=1.25; xlim([-lim lim]); ylim([-lim lim]); grid on; box on;
xlabel('Re contribution'); ylabel('Im contribution');
legend('groove 1','groove 2','residual','Location','best');
end

function draw_scattering_matrix_pair(Splus,Sminus)
M=[abs(Splus),nan(2,1),abs(Sminus.')];
imagesc(M); axis image; box on; clim([0 1]); colorbar;
set(gca,'XTick',[1 2 4 5],'XTickLabel',{'S_+:0','S_+:-1','S_-^T:0','S_-^T:+1'}, ...
    'YTick',[1 2],'YTickLabel',{'port 1','port 2'});
for r=1:2
    for c=[1 2 4 5]
        if c==3 || isnan(M(r,c)), continue; end
        value=M(r,c); text(c,r,sprintf('%.3f',value),'HorizontalAlignment','center', ...
            'Color',conditional_text_color(value),'FontWeight','bold');
    end
end
xlabel('input port / reciprocal representation'); ylabel('output port');
end

function draw_reciprocity_residual(Splus,Sminus,reciprocityError,unitarityError)
E=abs(Splus-Sminus.'); imagesc(log10(max(E,1e-18))); axis image; box on;
clim([-18 -10]); colorbar; set(gca,'XTick',1:2,'YTick',1:2);
xlabel('input port'); ylabel('output port');
for r=1:2
    for c=1:2
        text(c,r,sprintf('%.1e',E(r,c)),'HorizontalAlignment','center','Color','w','FontWeight','bold');
    end
end
text(2.72,.72,sprintf('||S_+-S_-^T||_F = %.1e',reciprocityError),'Color','k');
text(2.72,1.18,sprintf('unitarity: %.1e / %.1e',unitarityError(1),unitarityError(2)),'Color','k');
xlim([.5 4.8]);
end

function color=conditional_text_color(value)
if value>.55, color='k'; else, color='w'; end
end

function draw_experiment_schematic(D,water,solid,navy,orange,red)
axis off; hold on; rectangle('Position',[.02,.08,.96,.82],'FaceColor',water,'EdgeColor',navy,'LineWidth',1.4);
rectangle('Position',[.34,.30,.44,.08],'FaceColor',solid,'EdgeColor','k');
for j=0:7, rectangle('Position',[.36+.05*j,.30,.025,.06],'FaceColor','w','EdgeColor','none'); end
rectangle('Position',[.08,.52,.10,.20],'Curvature',[.2 .2],'FaceColor',orange,'EdgeColor','k');
quiver(.19,.62,.20,-.18,0,'Color',navy,'LineWidth',2,'MaxHeadSize',.5);
plot(.62+.22*cosd(linspace(10,170,80)),.42+.22*sind(linspace(10,170,80)),'--','Color',red,'LineWidth',1.5);
plot(.62,.66,'o','Color',red,'MarkerFaceColor',red); plot([.38 .76],[.45 .45],'Color',navy,'LineWidth',1.2);
text(.07,.76,'200 kHz transmitter'); text(.50,.24,'M=20,40,60-period samples');
text(.66,.78,'hydrophone arc'); text(.40,.47,'near-field raster scan');
xlim([0 1]); ylim([0 1]);
end

function draw_scan_plan(theta,f0,navy,orange,red)
axis off; hold on;
steps={'1  ToF sound-speed calibration','2  190-205 kHz x -7 to +7 deg map', ...
    '3  10-50 Hz fine scan near RA','4  complex p(x,y) raster + ring-down'};
for j=1:4
    y=.88-(j-1)*.22; rectangle('Position',[.08,y-.10,.84,.13],'Curvature',.04,'EdgeColor',[.6 .6 .6],'FaceColor',[.97 .97 .97]);
    text(.12,y-.035,steps{j},'Color','k');
    if j<4, quiver(.50,y-.11,0,-.08,0,'Color',red,'LineWidth',1.5,'MaxHeadSize',.8); end
end
text(.5,.03,sprintf('target: f=%.1f kHz, theta=+/-%.3f deg',f0/1e3,theta),'HorizontalAlignment','center','Color',navy);
xlim([0 1]); ylim([0 1]);
end

function draw_fourier_extraction(R,F,navy,orange,red)
axis off; hold on; x=F.x;
[~,iy]=min(abs(F.y-.25)); trace=real(F.p(iy,:));
trace=trace/max(max(abs(trace)),eps); xp=.08+.84*x; yp=.76+.13*trace;
plot(xp,yp,'Color',navy,'LineWidth',1.5); plot([.08 .92],[.76 .76],'k-','LineWidth',.8);
text(.08,.94,'measured complex trace p(x,y_0)','Color',navy);
quiver(.50,.61,0,-.12,0,'Color',red,'LineWidth',1.5,'MaxHeadSize',.7);
text(.50,.64,'least-squares Floquet fit','HorizontalAlignment','center','Color',red);
ids=order_ids(R); orderLabels={'n=-1','n=0','n=+1'}; colors={orange,navy,red};
for j=1:3
    left=.05+(j-1)*.32;
    rectangle('Position',[left,.14,.27,.27],'FaceColor',[.97 .97 .97], ...
        'EdgeColor',colors{j},'LineWidth',1.4);
    text(left+.135,.34,orderLabels{j},'HorizontalAlignment','center','Color',colors{j});
    text(left+.135,.25,sprintf('|A|=%.3g',abs(R.A(ids(j)))),'HorizontalAlignment','center');
    text(left+.135,.18,sprintf('eta=%.3f',R.eta(ids(j))),'HorizontalAlignment','center');
end
text(.50,.04,'propagating eta_n = Re(k_{y,n}) |A_n|^2 / k_{y,inc}', ...
    'HorizontalAlignment','center');
xlim([0 1]); ylim([0 1]);
end

function draw_experiment_workflow(navy,orange,red,green)
axis off; hold on; labels={'complex far field','A_n and eta_n','near-field hotspot','ring-down Q','finite-M scaling'};
colors={navy,orange,red,green,[.4 .2 .65]};
for j=1:5
    x=.04+(j-1)*.19; rectangle('Position',[x,.40,.15,.20],'Curvature',.08,'FaceColor',[.97 .97 .97],'EdgeColor',colors{j},'LineWidth',1.5);
    text(x+.075,.50,labels{j},'HorizontalAlignment','center','Color',colors{j});
    if j<5, quiver(x+.15,.50,.035,0,0,'Color','k','LineWidth',1.2,'MaxHeadSize',1); end
end
text(.5,.75,'one evidence chain, not one spectacular but ambiguous spectrum','HorizontalAlignment','center');
text(.5,.20,'controls: flat plate | single cavity | blocked small cavity | M=20/40/60','HorizontalAlignment','center');
xlim([0 1]); ylim([0 1]);
end

function draw_router_truth_table(theta,navy,orange,red)
axis off; hold on;
rectangle('Position',[.08,.55,.84,.24],'FaceColor',[.97 .97 .97],'EdgeColor',navy);
rectangle('Position',[.08,.20,.84,.24],'FaceColor',[.97 .97 .97],'EdgeColor',orange);
text(.13,.70,sprintf('input +%.3f deg',theta),'Color',navy); text(.58,.70,'selected anomalous port: n=-1 (left)','Color',red);
text(.13,.35,sprintf('input -%.3f deg',theta),'Color',orange); text(.58,.35,'selected anomalous port: n=+1 (right)','Color',red);
text(.50,.06,'time-reversal partners: same efficiency, opposite signed order','HorizontalAlignment','center');
xlim([0 1]); ylim([0 1]);
end
