%% Groove-resolved radiation contributions along the physical eigenbranch
% This is the effective D(kx)a(kx) content used in the mechanism figure.
clear; clc;
rootDir = fileparts(mfilename('fullpath'));
outDir = fullfile(rootDir,'results','fig3_radiation_contributions_180k');
if ~exist(outDir,'dir'), mkdir(outDir); end
addpath(rootDir);

S = load(fullfile(rootDir,'results', ...
    'StrictRayleighBIC_180kHz_7p10deg_final.mat'));
cfg = S.cfg; cfg.solve_scattering = false;
kB = S.x(1); OB = 1-kB;
opB = S.op; zB = S.z(:);
CB = zB(opB.N+1:end); CB = CB/max(norm(CB),eps);
phi = dressed_basis(opB,zB);

% Resolve both the global evolution and the asymptotic common-zero region.
% Logarithmic offsets prevent the singular endpoint from being represented by
% only one neighboring sample, while the outer linear grid retains the full
% mechanism window used in the main figure.
dk = unique([logspace(-7,-3,49),linspace(1.05e-3,0.022,55)]).';
kLeft = kB-flipud(dk);
kRight = kB+dk;
left = one_side(flipud(kLeft),OB,cfg,CB,phi); left = flipud(left);
right = one_side(kRight,OB,cfg,CB,phi);

[b0B,bmB,aB,d0B,dmB,fidB] = contributions(opB,zB,cfg,phi);
endpoint = row_from_values(kB,OB,0,b0B,bmB,aB,d0B,dmB,fidB);
data = [left; endpoint; right];
writetable(data,fullfile(outDir,'radiation_contributions.csv'));
save(fullfile(outDir,'radiation_contributions.mat'),'data','S','cfg');

fprintf('Exported %d physical-branch points around kappa %.12f\n', ...
    height(data),kB);
fprintf('BIC n=0 sum %.3e, n=-1 sum %.3e\n', ...
    abs(sum(b0B)),abs(sum(bmB)));

function T = one_side(kappa,OmegaSeed,cfg,Cref,phi)
rows = cell(numel(kappa),1); seed = complex(OmegaSeed,1e-11);
for j=1:numel(kappa)
    p = ni2019_refine_outgoing_pole_kappa(cfg,kappa(j),seed, ...
        'OuterIterations',12,'Display','off', ...
        'FunctionTolerance',5e-12,'StepTolerance',5e-12);
    seed = p.Omega;
    z = p.mode_vector(:); C = z(p.operator.N+1:end);
    C = C/max(norm(C),eps);
    % Parallel-transport gauge, anchored to the strict BIC vector.
    phase = angle(Cref'*C); z = z*exp(-1i*phase);
    C = C*exp(-1i*phase); Cref = C;
    [b0,bm,a,d0,dm,fidelity] = contributions(p.operator,z,cfg,phi);
    rows{j} = row_from_values(kappa(j),p.Omega,p.Q,b0,bm,a, ...
        d0,dm,fidelity);
end
T = vertcat(rows{:});
end

function T = row_from_values(kappa,Omega,Q,b0,bm,a,d0,dm,fidelity)
Dm = [d0(:).';dm(:).'];
[~,Sm,Vm] = svd(Dm,'econ');
sigmaRatio = Sm(end,end)/max(Sm(1,1),eps);
darkVector = Vm(:,end);
alignment = abs(darkVector'*a(:))^2/ ...
    max(norm(darkVector)^2*norm(a(:))^2,eps);
darkMisalignment = max(0,1-min(1,alignment));
residual0 = abs(sum(b0))/max(sum(abs(b0)),eps);
residualm1 = abs(sum(bm))/max(sum(abs(bm)),eps);
T = table(kappa,real(Omega),imag(Omega),Q, ...
    real(a(1)),imag(a(1)),real(a(2)),imag(a(2)), ...
    real(d0(1)),imag(d0(1)),real(d0(2)),imag(d0(2)), ...
    real(dm(1)),imag(dm(1)),real(dm(2)),imag(dm(2)), ...
    fidelity(1),fidelity(2), ...
    real(b0(1)),imag(b0(1)),real(b0(2)),imag(b0(2)), ...
    real(sum(b0)),imag(sum(b0)),abs(sum(b0)), ...
    real(bm(1)),imag(bm(1)),real(bm(2)),imag(bm(2)), ...
    real(sum(bm)),imag(sum(bm)),abs(sum(bm)), ...
    sigmaRatio,darkMisalignment,residual0,residualm1, ...
    'VariableNames',{'kappa','Omega_real','Omega_imag','Q', ...
    'a_large_real','a_large_imag','a_small_real','a_small_imag', ...
    'd0_large_real','d0_large_imag','d0_small_real','d0_small_imag', ...
    'dm1_large_real','dm1_large_imag','dm1_small_real','dm1_small_imag', ...
    'basis_fidelity_large','basis_fidelity_small', ...
    'b0_large_real','b0_large_imag','b0_small_real','b0_small_imag', ...
    'b0_sum_real','b0_sum_imag','b0_sum_abs', ...
    'bm1_large_real','bm1_large_imag','bm1_small_real','bm1_small_imag', ...
    'bm1_sum_real','bm1_sum_imag','bm1_sum_abs', ...
    'radiation_sigma_ratio','dark_state_misalignment', ...
    'b0_normalized_residual','bm1_normalized_residual'});
end

function [b0,bm,a,d0,dm,fidelity] = contributions(op,z,cfg,phi)
C = z(op.N+1:end); scale=max(norm(C),eps); C=C/scale;
betaSin = op.beta_sin_normalized(:);
cosDepth = op.cos_depth_normalized(:);
orders = op.orders(:); kx = op.kx(:);
widths=cfg.widths(:).'; gaps=cfg.gaps(:).';
occupied=sum(widths)+sum(gaps); xleft=zeros(1,2);
xleft(1)=.5*(cfg.a-occupied); xleft(2)=xleft(1)+widths(1)+gaps(1);
b0 = channel_pair(0,C); bm = channel_pair(-1,C);

% Fixed linear projection onto the two dressed aperture modes extracted at
% the strict BIC.  Each dressed coordinate retains every transverse order.
surface=reshape(C.*op.cos_depth_normalized(:),op.K,op.L);
a=complex(zeros(1,2));
a(1)=phi(:,1)'*surface(:,1);
a(2)=phi(:,2)'*surface(:,2);
fidelity=zeros(1,2);
for ell=1:2
    fidelity(ell)=abs(a(ell))^2/max(norm(surface(:,ell))^2,eps);
end
d0 = channel_pair(0,reshape(phi./reshape(cosDepth,op.K,op.L),[],1));
dm = channel_pair(-1,reshape(phi./reshape(cosDepth,op.K,op.L),[],1));

    function b=channel_pair(order,coefficient)
        n=find(orders==order,1); b=complex(zeros(1,2));
        for ell=1:2
            for q=0:op.K-1
                j=(ell-1)*op.K+q+1;
                cmap=(widths(ell)/cfg.a)*projection( ...
                    -kx(n),q,xleft(ell),widths(ell));
                b(ell)=b(ell)+1i*cmap*betaSin(j)*coefficient(j);
            end
        end
    end
end

function phi=dressed_basis(op,z)
C=z(op.N+1:end); C=C/max(norm(C),eps);
surface=reshape(C.*op.cos_depth_normalized(:),op.K,op.L);
phi=complex(zeros(op.K,op.L));
for ell=1:op.L
    phi(:,ell)=surface(:,ell)/max(norm(surface(:,ell)),eps);
end
end

function value=projection(kx,q,xl,t)
alpha=q*pi/t;
value=exp(-1i*kx*xl)*.5*(expint_local(kx-alpha,t)+ ...
    expint_local(kx+alpha,t));
end

function value=expint_local(kappa,t)
z=kappa*t;
if abs(z)<1e-8
    value=1-1i*z/2-z^2/6;
else
    value=exp(-1i*z/2)*sin(z/2)/(z/2);
end
end
