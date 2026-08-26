%% Export homogeneous eigenfields at the Gamma BIC, Rayleigh BIC, and beyond it
% Numerical source for Fig. 2e.  These are homogeneous complex-frequency
% eigenmodes; no incident field is introduced.  The Gamma sample is the
% symmetry-protected BIC on the same physical target branch, the
% central state is reconstructed from the strict zero-radiation operator,
% and the post-BIC state lies on the physical sheet with n=-1 open.
clear; clc;

rootDir = fileparts(mfilename('fullpath'));
outDir = fullfile(rootDir,'results','fig2_three_state_fields_180k');
if ~exist(outDir,'dir'), mkdir(outDir); end

S = load(fullfile(rootDir,'results', ...
    'StrictRayleighBIC_180kHz_7p10deg_final.mat'));
P = load(fullfile(rootDir,'results','fig2_three_state_fields_180k', ...
    'high_truncation_pole_samples.mat'));
T = load(fullfile(rootDir,'results','target_branch_physical_bz_180k', ...
    'target_branch_physical_bz.mat'));
cfg = S.cfg;
cfg.solve_scattering = false;
kappaBIC = S.x(1);
OmegaBIC = 1-kappaBIC;

% The Gamma comparison is obtained from the all-channel outgoing/decaying
% operator on the physical target branch.
kappaBeyondRequested = P.beyond.kappa;

x = linspace(0,5,1001);
y = linspace(-0.74,3.2,641);
[X,Y] = meshgrid(x,y);

gammaSeed = complex(T.data.Omega_real(1),T.data.Omega_imag(1));
gammaPole = ni2019_refine_outgoing_pole_kappa(cfg,0,gammaSeed, ...
    'OuterIterations',12,'Display','off');
states = physical_mode(gammaPole,cfg,X,Y,'Gamma BIC');

% Strict endpoint: every finite-flux channel and the n=-1 threshold
% pressure amplitude are removed before obtaining the null vector.
R = ni2019_strict_rayleigh_operator(cfg,'TargetOrder',-1, ...
    'EnforceAllFiniteOpen',true,'KyTolerance',1e-6,'Verbose',true);
modeBIC = struct('A',R.mode.A,'surface', ...
    reshape(R.mode.surface_coefficients,cfg.K,numel(cfg.widths)), ...
    'kx',R.full_operator.kx,'ky',R.full_operator.ky, ...
    'k0',R.full_operator.k0,'kappa',kappaBIC);
pBIC = homogeneous_pressure(modeBIC,cfg,X,Y);
states(2) = pack_state('Rayleigh BIC',kappaBIC,complex(OmegaBIC,0), ...
    Inf,R.sigma_ratio,R.residual.full_raw,R.mode.A, ...
    R.full_operator.orders,R.full_operator.kx,R.full_operator.ky,pBIC);

states(3) = continued_mode(P.beyond.kappa,P.beyond.Omega, ...
    P.beyond.qbar,cfg,X,Y,'beyond BIC');

% Remove the arbitrary eigenvector scale with the same physical rule in all
% panels, then use one global pressure reference and one common color scale.
groove = groove_mask(cfg,X,Y);
dx = mean(diff(x)); dy = mean(diff(y));
nPeriods = max(x)-min(x);
for j = 1:numel(states)
    Eg = sum(abs(states(j).p(groove)).^2,'omitnan')*dx*dy/nPeriods;
    states(j).groove_energy_before_normalization = Eg;
    states(j).p = states(j).p/sqrt(max(Eg,eps));
end
globalMax = max(arrayfun(@(s)max(abs(s.p),[],'all','omitnan'),states));
for j = 1:numel(states)
    states(j).p_normalized = states(j).p/globalMax;
    states(j).log_magnitude = log10(max(abs(states(j).p_normalized),1e-7));
end

meta = struct2table(rmfield(states,{'p','p_normalized','log_magnitude'}));
writetable(meta,fullfile(outDir,'field_metadata.csv'));
save(fullfile(outDir,'eigenmode_fields_gamma_bic_beyond.mat'), ...
    'states','x','y','X','Y','cfg','kappaBIC','OmegaBIC','globalMax', ...
    'kappaBeyondRequested');

fprintf('Exported homogeneous Fig. 2e eigenfields:\n');
for j = 1:numel(states)
    fprintf('  %-13s kappa=%+.10f, Omega=%.10f%+.3ei, Q=%.4g, ', ...
        states(j).label,states(j).kappa,real(states(j).Omega), ...
        imag(states(j).Omega),states(j).Q);
    fprintf('|A0|/||C||=%.3e, |A-1|/||C||=%.3e\n', ...
        states(j).A0_over_C,states(j).Am1_over_C);
end

function state = physical_mode(pole,cfg,X,Y,label)
op = pole.operator;
z = pole.mode_vector;
A = z(1:op.N);
Cscaled = z(op.N+1:end);
surface = reshape(Cscaled.*op.cos_depth_normalized,cfg.K,numel(cfg.widths));
mode = struct('A',A,'surface',surface,'kx',op.kx,'ky',op.ky, ...
    'k0',op.k0,'kappa',pole.kappa);
p = homogeneous_pressure(mode,cfg,X,Y);
state = pack_state(label,pole.kappa,pole.Omega,pole.Q, ...
    pole.sigma_ratio,pole.raw_residual,A,op.orders,op.kx,op.ky,p);
end

function state = continued_mode(kappa,Omega,qbar,cfg,X,Y,label)
op = ni2019_full_eigen_operator_complex(cfg,2*pi*Omega,2*pi*kappa, ...
    -1,2*pi*qbar);
[~,Sv,V] = svd(op.Fscaled,'econ');
singular = diag(Sv);
yScaled = V(:,end);
z = yScaled./transpose(op.column_scale);
C = z(op.N+1:end);
z = z/max(norm(C),eps);
A = z(1:op.N);
Cscaled = z(op.N+1:end);
surface = reshape(Cscaled.*op.cos_depth_normalized,cfg.K,numel(cfg.widths));
mode = struct('A',A,'surface',surface,'kx',op.kx,'ky',op.ky, ...
    'k0',op.k0,'kappa',kappa);
p = homogeneous_pressure(mode,cfg,X,Y);
Q = real(Omega)/(2*abs(imag(Omega)));
rawResidual = norm(op.F*z)/max(norm(z),eps);
state = pack_state(label,kappa,Omega,Q,singular(end)/singular(1), ...
    rawResidual,A,op.orders,op.kx,op.ky,p);
end

function state = pack_state(label,kappa,Omega,Q,sigma,raw,A,orders,kx,ky,p)
idp1 = orders==1; id0 = orders==0; idm1 = orders==-1;
state = struct('label',label,'kappa',kappa, ...
    'theta_deg',asind(real(kappa/Omega)),'Omega',Omega,'Q',Q, ...
    'sigma_ratio',sigma,'raw_residual',raw, ...
    'Ap1_over_C',abs(A(idp1)),'A0_over_C',abs(A(id0)), ...
    'Am1_over_C',abs(A(idm1)), ...
    'Ap1_complex',A(idp1),'A0_complex',A(id0),'Am1_complex',A(idm1), ...
    'kxp1',kx(idp1),'kx0',kx(id0),'kxm1',kx(idm1), ...
    'kyp1',ky(idp1),'ky0',ky(id0),'kym1',ky(idm1), ...
    'A_all',A,'orders_all',orders,'kx_all',kx,'ky_all',ky, ...
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
xCell = floor(X+1e-12);
xLocal = X-xCell;
phase = exp(-1i*2*pi*mode.kappa.*xCell);
for ell = 1:numel(cfg.widths)
    inside = Y<0 & Y>=-cfg.depths(ell) & ...
        xLocal>=xleft(ell) & xLocal<=xleft(ell)+cfg.widths(ell);
    if ~any(inside,'all'), continue; end
    u = xLocal(inside)-xleft(ell);
    yy = Y(inside);
    pg = complex(zeros(size(u)));
    for q = 0:cfg.K-1
        alpha = q*pi/cfg.widths(ell);
        beta = sqrt(mode.k0^2-alpha^2);
        vertical = cos(beta*(yy+cfg.depths(ell)))/cos(beta*cfg.depths(ell));
        pg = pg + mode.surface(q+1,ell).*cos(alpha*u).*vertical;
    end
    p(inside) = phase(inside).*pg;
end
end

function mask = groove_mask(cfg,X,Y)
xleft = groove_left_edges(cfg);
xLocal = X-floor(X+1e-12);
mask = false(size(X));
for ell = 1:numel(cfg.widths)
    mask = mask | (Y<0 & Y>=-cfg.depths(ell) & ...
        xLocal>=xleft(ell) & xLocal<=xleft(ell)+cfg.widths(ell));
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
