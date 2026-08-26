%% Isolated-groove complex eigenfrequencies in an infinite rigid baffle
% The periodic Floquet sum is taken to its continuous a->infinity limit.
% No effective end correction is used; transverse groove modes q=0...K-1
% and the full outgoing angular spectrum are matched directly.
clear; close all; clc;

thisFile=mfilename('fullpath'); solverDir=fileparts(thisFile);
resultDir=fullfile(solverDir,'results','isolated_groove_eigenmodes');
if ~exist(resultDir,'dir'), mkdir(resultDir); end
addpath(solverDir);
D=load(fullfile(solverDir,'results', ...
    'StrictRayleighBIC_200kHz_min1mm.mat'));
x=D.xFinal(:).'; a=D.aPhysical; cWater=D.cWater;
names={'wide_q1','narrow_q0'};
width=[a*x(4),a*x(5)]; depth=[a*x(2),a*x(3)];
seed=[200.8+2.45i,213.1+66.1i]*1e3;
Klist=[10 14 18 22 30 38 46 62].';
tailMax=300; panelsPerPi=4;

rows=[]; poles=cell(2,numel(Klist));
for cavity=1:2
    fSeed=seed(cavity);
    for j=1:numel(Klist)
        K=Klist(j);
        pole=ni2019_find_isolated_groove_pole(width(cavity), ...
            depth(cavity),K,cWater,fSeed,'OuterIterations',7, ...
            'SMax',tailMax,'PanelsPerPi',panelsPerPi);
        fSeed=pole.frequency_hz; poles{cavity,j}=pole;
        [mainFraction,mainIndex]=max(pole.transverse_fraction);
        surfaceMain=max(pole.surface_pressure_fraction);
        rows=[rows;{names{cavity},K,real(fSeed),imag(fSeed),pole.Q, ...
            pole.sigma_ratio,pole.raw_residual,mainIndex-1,mainFraction, ...
            surfaceMain}]; %#ok<AGROW>
        fprintf('%s K=%d: f=(%.9f%+.9fi) kHz, Q=%.6f, sigma=%.3e\n', ...
            names{cavity},K,real(fSeed)/1e3,imag(fSeed)/1e3, ...
            pole.Q,pole.sigma_ratio);
    end
end
convergence=cell2table(rows,'VariableNames',{'cavity','K','Re_f_Hz', ...
    'Im_f_Hz','Q','sigma_ratio','raw_residual','dominant_q', ...
    'bottom_modal_fraction','surface_modal_fraction'});
writetable(convergence,fullfile(resultDir,'isolated_groove_K_convergence.csv'));

% Edge singularities make modal convergence algebraic.  Report a transparent
% 1/K extrapolation of the four highest truncations and retain the K=62 value.
summaryRows=[];
for cavity=1:2
    select=strcmp(convergence.cavity,names{cavity});
    tableCase=convergence(select,:); fitRows=height(tableCase)-3:height(tableCase);
    fitReal=polyfit(1./tableCase.K(fitRows),tableCase.Re_f_Hz(fitRows),1);
    fitImag=polyfit(1./tableCase.K(fitRows),tableCase.Im_f_Hz(fitRows),1);
    fInf=fitReal(2)+1i*fitImag(2); QInf=real(fInf)/(2*imag(fInf));
    fLast=tableCase.Re_f_Hz(end)+1i*tableCase.Im_f_Hz(end);
    summaryRows=[summaryRows;{names{cavity},width(cavity)*1e3, ...
        depth(cavity)*1e3,real(fLast),imag(fLast),tableCase.Q(end), ...
        real(fInf),imag(fInf),QInf,100*(real(fInf)-200e3)/200e3}]; %#ok<AGROW>
end
summary=cell2table(summaryRows,'VariableNames',{'cavity','width_mm', ...
    'depth_mm','Re_f_K62_Hz','Im_f_K62_Hz','Q_K62', ...
    'Re_f_extrapolated_Hz','Im_f_extrapolated_Hz','Q_extrapolated', ...
    'detuning_from_200k_percent'});
writetable(summary,fullfile(resultDir,'isolated_groove_summary.csv'));
save(fullfile(resultDir,'isolated_groove_eigenmodes.mat'),'D','width', ...
    'depth','Klist','tailMax','panelsPerPi','poles','convergence','summary');

fig=figure('Visible','off','Color','w','Units','inches', ...
    'Position',[.5,.5,6.6,2.65]);
set(fig,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica', ...
    'DefaultAxesFontSize',8,'DefaultTextFontSize',8, ...
    'DefaultAxesLineWidth',.65,'DefaultLineLineWidth',1.1);
tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
colors=[.08,.36,.58;.82,.32,.12];
nexttile; hold on;
for cavity=1:2
    select=strcmp(convergence.cavity,names{cavity});
    plot(1./convergence.K(select),convergence.Re_f_Hz(select)/1e3, ...
        'o-','Color',colors(cavity,:),'MarkerFaceColor',colors(cavity,:));
end
yline(200,'--','Color',[.35,.35,.35]);
xlabel('1/K'); ylabel('Re f_p (kHz)'); box on;
legend('wide groove','narrow groove','BIC','Location','best','Box','off');
nexttile; hold on;
for cavity=1:2
    select=strcmp(convergence.cavity,names{cavity});
    plot(convergence.Re_f_Hz(select)/1e3, ...
        convergence.Im_f_Hz(select)/1e3,'o-','Color',colors(cavity,:), ...
        'MarkerFaceColor',colors(cavity,:));
end
xlabel('Re f_p (kHz)'); ylabel('Im f_p (kHz)'); box on;
exportgraphics(fig,fullfile(resultDir,'isolated_groove_convergence.pdf'), ...
    'ContentType','vector');
exportgraphics(fig,fullfile(resultDir,'isolated_groove_convergence.png'), ...
    'Resolution',300);

disp(summary);
