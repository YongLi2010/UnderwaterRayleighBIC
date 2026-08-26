%% Export homogeneous eigenfields along the 180-kHz Rayleigh-BIC branch
% No incident field is used anywhere in this script.  The two off-BIC states
% are outgoing complex-frequency eigenmodes from the stored pole branch; the
% center state is the strict zero-radiation homogeneous solution.
clear; clc;

rootDir = fileparts(mfilename('fullpath'));
outDir = fullfile(rootDir,'results','fig2_eigenmode_fields_180k');
if ~exist(outDir,'dir'), mkdir(outDir); end

S = load(fullfile(rootDir,'results', ...
    'StrictRayleighBIC_180kHz_7p10deg_final.mat'));
E = readtable(fullfile(rootDir,'results','kx_channel_evolution_180k', ...
    'kx_channel_evolution.csv'));
cfg = S.cfg;
kappaBIC = S.x(1);
OmegaBIC = 1-kappaBIC;
thetaData = asind(E.kappa./E.Omega_real);
thetaBIC = asind(kappaBIC/OmegaBIC);

% Existing continued-pole samples nearest to symmetric +/-0.5-degree views.
targetTheta = [thetaBIC-.5,thetaBIC+.5];
[~,idLeft] = min(abs(thetaData-targetTheta(1)));
[~,idRight] = min(abs(thetaData-targetTheta(2)));
ids = [idLeft,idRight];

x = linspace(0,1,321);
y = linspace(-0.74,1.05,361);
[X,Y] = meshgrid(x,y);
% Left outgoing eigenmode.
states = continued_mode(E(ids(1),:),cfg,X,Y,'below Rayleigh line');

% Strict endpoint: both A_0 and A_-1 are removed from the homogeneous
% operator before the null vector is reconstructed.
R = ni2019_strict_rayleigh_operator(cfg,'TargetOrder',-1, ...
    'EnforceAllFiniteOpen',true,'KyTolerance',1e-6,'Verbose',true);
modeCenter = struct('A',R.mode.A,'surface', ...
    reshape(R.mode.surface_coefficients,cfg.K,numel(cfg.widths)), ...
    'kx',R.full_operator.kx,'ky',R.full_operator.ky, ...
    'k0',R.full_operator.k0);
pCenter = homogeneous_pressure(modeCenter,cfg,X,Y);
states(2) = pack_state('Rayleigh BIC',thetaBIC,kappaBIC, ...
    complex(OmegaBIC,0),Inf,R.sigma_ratio,R.residual.full_raw, ...
    R.mode.A,R.full_operator.orders,R.full_operator.ky,pCenter,cfg,X,Y);

% Right outgoing eigenmode.
states(3) = continued_mode(E(ids(2),:),cfg,X,Y,'above Rayleigh line');

% Fix the arbitrary eigenvector scale by one common physical rule:
% integral_grooves |p|^2 dS = 1.  A single global reference then sets the
% plotted logarithmic color scale for all three states.
for j = 1:3
    groove = groove_mask(cfg,X,Y);
    dx = mean(diff(x)); dy = mean(diff(y));
    Eg = sum(abs(states(j).p(groove)).^2,'omitnan')*dx*dy;
    states(j).groove_energy_before_normalization = Eg;
    states(j).p = states(j).p/sqrt(max(Eg,eps));
end
globalMax = max(arrayfun(@(s)max(abs(s.p),[],'all','omitnan'),states));
for j = 1:3
    states(j).p_normalized = states(j).p/globalMax;
    states(j).log_magnitude = log10(max(abs(states(j).p_normalized),1e-7));
end

meta = struct2table(rmfield(states,{'p','p_normalized','log_magnitude'}));
writetable(meta,fullfile(outDir,'field_metadata.csv'));
save(fullfile(outDir,'eigenmode_fields.mat'),'states','x','y','X','Y', ...
    'cfg','thetaBIC','kappaBIC','OmegaBIC','globalMax');

fprintf('Exported homogeneous eigenfields:\n');
for j = 1:3
    fprintf('  theta=%8.5f deg, ReOmega=%.10f, ImOmega=%+.3e, ', ...
        states(j).theta_deg,real(states(j).Omega),imag(states(j).Omega));
    fprintf('|A0|/||C||=%.3e, |A-1|/||C||=%.3e\n', ...
        states(j).A0_over_C,states(j).Am1_over_C);
end

function state = continued_mode(row,cfg,X,Y,label)
kappa = row.kappa;
Omega = complex(row.Omega_real,row.Omega_imag);
qbar = complex(row.q_real,row.q_imag);
op = ni2019_full_eigen_operator_complex(cfg,2*pi*Omega,2*pi*kappa, ...
    -1,2*pi*qbar);
[~,Sv,V] = svd(op.Fscaled,'econ');
singular = diag(Sv);
z = V(:,end)./transpose(op.column_scale);
C = z(op.N+1:end);
z = z/max(norm(C),eps);
A = z(1:op.N);
Cscaled = z(op.N+1:end);
surface = reshape(Cscaled.*op.cos_depth_normalized,cfg.K,numel(cfg.widths));
mode = struct('A',A,'surface',surface,'kx',op.kx,'ky',op.ky,'k0',op.k0);
p = homogeneous_pressure(mode,cfg,X,Y);
theta = asind(kappa/real(Omega));
Q = real(Omega)/(2*abs(imag(Omega)));
rawResidual = norm(op.F*z)/max(norm(z),eps);
state = pack_state(label,theta,kappa,Omega,Q, ...
    singular(end)/singular(1),rawResidual,A,op.orders,op.ky,p,cfg,X,Y);
end

function state = pack_state(label,theta,kappa,Omega,Q,sigma,raw,A,orders,ky,p,cfg,X,Y)
groove = groove_mask(cfg,X,Y);
Cproxy = sqrt(sum(abs(p(groove)).^2,'omitnan'));
id0 = orders==0; idm1 = orders==-1;
tol = 1e-7*max(abs(2*pi*Omega),1);
finiteOpen = abs(imag(ky))<=tol & real(ky)>tol;
threshold = abs(ky)<=tol;
state = struct('label',label,'theta_deg',theta,'kappa',kappa,'Omega',Omega, ...
    'Q',Q,'sigma_ratio',sigma,'raw_residual',raw, ...
    'A0_over_C',abs(A(id0))/max(Cproxy,eps), ...
    'Am1_over_C',abs(A(idm1))/max(Cproxy,eps), ...
    'finite_open_orders',strtrim(sprintf('%d ',orders(finiteOpen))), ...
    'threshold_orders',strtrim(sprintf('%d ',orders(threshold))), ...
    'groove_energy_before_normalization',NaN,'p',p, ...
    'p_normalized',[],'log_magnitude',[]);
end

function p = homogeneous_pressure(mode,cfg,X,Y)
p = complex(nan(size(X)));
above = Y>=0;
pa = complex(zeros(nnz(above),1));
for n = 1:numel(mode.A)
    pa = pa + mode.A(n).*exp(-1i*mode.kx(n).*X(above) - ...
        1i*mode.ky(n).*Y(above));
end
p(above) = pa;

xleft = groove_left_edges(cfg);
for ell = 1:numel(cfg.widths)
    inside = Y<0 & Y>=-cfg.depths(ell) & ...
        X>=xleft(ell) & X<=xleft(ell)+cfg.widths(ell);
    if ~any(inside,'all'), continue; end
    u = X(inside)-xleft(ell);
    yy = Y(inside);
    pg = complex(zeros(size(u)));
    for q = 0:cfg.K-1
        alpha = q*pi/cfg.widths(ell);
        beta = sqrt(mode.k0^2-alpha^2);
        vertical = cos(beta*(yy+cfg.depths(ell)))/cos(beta*cfg.depths(ell));
        pg = pg + mode.surface(q+1,ell).*cos(alpha*u).*vertical;
    end
    p(inside) = pg;
end
end

function mask = groove_mask(cfg,X,Y)
xleft = groove_left_edges(cfg);
mask = false(size(X));
for ell = 1:numel(cfg.widths)
    mask = mask | (Y<0 & Y>=-cfg.depths(ell) & ...
        X>=xleft(ell) & X<=xleft(ell)+cfg.widths(ell));
end
end

function xleft = groove_left_edges(cfg)
occupied = sum(cfg.widths)+sum(cfg.gaps);
if isfield(cfg,'x0'), x0=cfg.x0; else, x0=.5*(cfg.a-occupied); end
xleft = zeros(1,numel(cfg.widths)); xleft(1)=x0;
for ell=2:numel(cfg.widths)
    xleft(ell)=xleft(ell-1)+cfg.widths(ell-1)+cfg.gaps(ell-1);
end
end
