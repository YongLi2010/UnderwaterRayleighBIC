%% Ideal infinite-period plane-wave scattering fields for Fig. 5
% Data only. Figure rendering is performed by the selected Python backend.
% At the exact dark states the forced system is singular, so the first two
% columns use the stable +0.1-Hz limiting background scattering solution.
clear; clc;
rootDir=fileparts(mfilename('fullpath'));
outDir=fullfile(rootDir,'results','fig5_infinite_period_fields_180k');
if ~exist(outDir,'dir'), mkdir(outDir); end
addpath(rootDir);

S=load(fullfile(rootDir,'results', ...
    'StrictRayleighBIC_180kHz_7p10deg_final.mat'));
cWater=1500;
stateLabel={'Gamma_limit','Rayleigh_BIC_limit','above_BIC'};
thetaDeg=[0,7.09966291282572,10.05404];
frequencyHz=[178780.3005,180000.1,180498.322];

xOverA=linspace(-2.5,2.5,1001);
yOverA=linspace(0.12,6.0,601);
[X,Y]=meshgrid(xOverA,yOverA);
writematrix(xOverA(:),fullfile(outDir,'x_over_a.csv'));
writematrix(yOverA(:),fullfile(outDir,'y_over_a.csv'));

orderRows=cell(3,1);
summary=table('Size',[3,8],'VariableTypes',repmat({'double'},1,8), ...
    'VariableNames',{'state','theta_deg','frequency_hz','condition_number', ...
    'energy_error','eta_m1','eta_0','eta_other'});

for state=1:3
    cfg=S.cfg;
    cfg.solve_scattering=true;
    cfg.lambda=1/(frequencyHz(state)*S.a/cWater);
    cfg.theta_i_deg=thetaDeg(state);
    R=ni2019_modal_solver(cfg);

    pIncident=exp(-1i*R.kx(R.orders==0).*X + 1i*R.ky_incident.*Y);
    pReflected=complex(zeros(size(X)));
    for idx=1:numel(R.orders)
        pReflected=pReflected+R.A(idx).* ...
            exp(-1i*R.kx(idx).*X-1i*R.ky(idx).*Y);
    end
    prefix=lower(stateLabel{state});
    writematrix(real(pIncident),fullfile(outDir,[prefix '_incident_real.csv']));
    writematrix(imag(pIncident),fullfile(outDir,[prefix '_incident_imag.csv']));
    writematrix(real(pReflected),fullfile(outDir,[prefix '_reflected_real.csv']));
    writematrix(imag(pReflected),fullfile(outDir,[prefix '_reflected_imag.csv']));

    keep=R.orders>=-3 & R.orders<=3;
    isPropagating=abs(imag(R.ky(keep)))<1e-9 & real(R.ky(keep))>1e-9;
    orderRows{state}=table(repmat(state,nnz(keep),1),R.orders(keep), ...
        real(R.kx(keep)./R.k0),real(R.ky(keep)./R.k0), ...
        imag(R.ky(keep)./R.k0),abs(R.A(keep)),angle(R.A(keep)), ...
        R.eta(keep),isPropagating, ...
        'VariableNames',{'state','order','kx_over_k0','ky_real_over_k0', ...
        'ky_imag_over_k0','amplitude_abs','amplitude_phase_rad', ...
        'power_fraction','is_propagating'});

    id0=find(R.orders==0,1); idm=find(R.orders==-1,1);
    etaOther=sum(R.eta)-R.eta(id0)-R.eta(idm);
    summary{state,:}=[state,thetaDeg(state),frequencyHz(state), ...
        R.condition_number,R.energy_error,R.eta(idm),R.eta(id0),etaOther];
end

orderTable=vertcat(orderRows{:});
writetable(orderTable,fullfile(outDir,'floquet_orders.csv'));
writetable(summary,fullfile(outDir,'state_summary.csv'));
save(fullfile(outDir,'infinite_period_scattering_fields.mat'), ...
    'S','thetaDeg','frequencyHz','stateLabel','xOverA','yOverA', ...
    'orderTable','summary','-v7.3');
disp(summary);
