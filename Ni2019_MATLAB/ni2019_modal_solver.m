function out = ni2019_modal_solver(cfg)
%NI2019_MODAL_SOLVER Modal-matching solver for Ni et al., PRB 100, 104104 (2019).
%
% The grating is periodic in x, rigid at y=0, and contains rectangular
% grooves extending into y<0. Pressure uses the exp(+j*omega*t)
% convention used in the paper. All lengths may be dimensional or
% normalized, provided that lambda, a, widths, depths, and positions use
% the same unit.
%
% Required cfg fields:
%   lambda, a, theta_i_deg, widths, depths, gaps
% Optional:
%   N (odd Floquet truncation, default 101)
%   K (groove modes, default 10)
%   x0 (left edge of first groove; default centers the groove group)

arguments
    cfg struct
end

lambda = cfg.lambda;
a = cfg.a;
theta_i = deg2rad(cfg.theta_i_deg);
widths = cfg.widths(:).';
depths = cfg.depths(:).';
L = numel(widths);

if numel(depths) ~= L
    error('widths and depths must contain the same number of grooves.');
end
if ~isfield(cfg, 'gaps'), cfg.gaps = zeros(1, max(0,L-1)); end
gaps = cfg.gaps(:).';
if numel(gaps) ~= max(0,L-1)
    error('gaps must contain L-1 edge-to-edge distances.');
end
if ~isfield(cfg, 'N'), cfg.N = 101; end
if ~isfield(cfg, 'K'), cfg.K = 10; end
if ~isfield(cfg, 'solve_scattering'), cfg.solve_scattering = true; end
if mod(cfg.N,2) ~= 1
    error('N must be odd so that orders -(N-1)/2:(N-1)/2 are retained.');
end

occupied = sum(widths) + sum(gaps);
if occupied >= a
    error('Grooves and gaps occupy %.4g, which is not smaller than a=%.4g.', occupied, a);
end
if ~isfield(cfg, 'x0')
    x0 = 0.5*(a-occupied); % Translation changes phases, not efficiencies.
else
    x0 = cfg.x0;
end
xleft = zeros(1,L);
xleft(1) = x0;
for ell = 2:L
    xleft(ell) = xleft(ell-1) + widths(ell-1) + gaps(ell-1);
end

k0 = 2*pi/lambda;
kx_i = k0*sin(theta_i);
ky_i = k0*cos(theta_i);
P = (cfg.N-1)/2;
orders = (-P:P).';
G = 2*pi*orders/a;
kx = kx_i + G;
ky = radiation_sqrt(k0^2-kx.^2); % reflected exp(-j*ky*y)

% B maps reflected Floquet amplitudes to cosine pressure coefficients at
% each groove. D maps groove cosine velocities to global Floquet tests.
K = cfg.K;
nUnknownGroove = L*K;
B = complex(zeros(nUnknownGroove,cfg.N));
binc = complex(zeros(nUnknownGroove,1));
D = complex(zeros(cfg.N,nUnknownGroove));
yadm = complex(zeros(nUnknownGroove,1)); % rho*omega*v/p
normCos = [1, 0.5*ones(1,K-1)];

row = 0;
for ell = 1:L
    t = widths(ell);
    d = depths(ell);
    xl = xleft(ell);
    for q = 0:K-1
        row = row+1;
        alpha = q*pi/t;
        beta = groove_sqrt(k0^2-alpha^2);
        binc(row) = cosine_projection(kx_i, q, xl, t);
        for n = 1:cfg.N
            B(row,n) = cosine_projection(kx(n), q, xl, t);
        end
        for j = 1:cfg.N
            % (1/a) integral exp(+j*G_j*x)*cos(alpha*(x-xl)) dx
            D(j,row) = (t/a)*cosine_projection(-kx(j), q, xl, t);
        end
        yadm(row) = -1i*beta*tan(beta*d);
        B(row,:) = B(row,:)/normCos(q+1);
        binc(row) = binc(row)/normCos(q+1);
    end
end

S = D.*transpose(yadm); % implicit expansion: each column times admittance
systemMatrix = diag(ky) - S*B;
rhs = ky_i*(orders==0) + S*binc;
if cfg.solve_scattering
    A = systemMatrix\rhs; % incident amplitude A0+ = 1
else
    A = complex(nan(cfg.N,1));
end

isProp = abs(imag(ky)) < 1e-10*k0 & real(ky) > 0;
eta = zeros(cfg.N,1);
eta(isProp) = real(ky(isProp))/ky_i .* abs(A(isProp)).^2;
theta = nan(cfg.N,1);
theta(isProp) = asind(real(kx(isProp))/k0);

out = struct();
out.orders = orders;
out.A = A;
out.eta = eta;
out.theta_deg = theta;
out.is_propagating = isProp;
out.total_efficiency = sum(eta);
out.energy_error = abs(1-out.total_efficiency);
out.k0 = k0;
out.kx = kx;
out.ky = ky;
out.ky_incident = ky_i;
out.a = a;
out.lambda = lambda;
out.theta_i_deg = cfg.theta_i_deg;
out.widths = widths;
out.depths = depths;
out.gaps = gaps;
out.xleft = xleft;
out.N = cfg.N;
out.K = cfg.K;
out.condition_number = cond(systemMatrix);
out.system_matrix = systemMatrix;
out.forcing = rhs;
out.groove_surface_coefficients = reshape(binc+B*A,K,L);
out.groove_mode_indices = (0:K-1).';
end

function value = cosine_projection(kx, q, xl, t)
% (1/t) integral_xl^(xl+t) exp(-j*kx*x)*cos(q*pi*(x-xl)/t) dx
alpha = q*pi/t;
value = exp(-1i*kx*xl)*0.5*(expint_stable(kx-alpha,t) + ...
                                  expint_stable(kx+alpha,t));
end

function value = expint_stable(kappa,t)
% (1/t) integral_0^t exp(-j*kappa*u) du, stable at kappa=0.
z = kappa*t;
if abs(z) < 1e-8
    value = 1 - 1i*z/2 - z^2/6;
else
    value = exp(-1i*z/2)*sin(z/2)/(z/2);
end
end

function value = radiation_sqrt(z)
% Branch needed by exp(-j*ky*y): evanescent orders decay for y>0.
value = complex(zeros(size(z)));
prop = real(z) >= 0;
value(prop) = sqrt(real(z(prop)));
value(~prop) = -1i*sqrt(-real(z(~prop)));
end

function value = groove_sqrt(z)
% Positive-imaginary branch is convenient in cos(beta*(y+d)).
if real(z) >= 0
    value = sqrt(real(z));
else
    value = 1i*sqrt(-real(z));
end
end
