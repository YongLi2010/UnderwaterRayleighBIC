function op = ni2019_full_eigen_operator(cfg)
%NI2019_FULL_EIGEN_OPERATOR Pole-free homogeneous modal-matching operator.
%   Unknowns are [A; C], where A are outgoing Floquet amplitudes and C are
%   groove cosine-mode amplitudes referenced at the rigid groove bottom:
%       p_g = C*cos(alpha*x)*cos(beta*(y+d)).
%   Keeping C avoids the tan(beta*d) poles introduced when groove variables
%   are eliminated in a scattering calculation.

if ~isfield(cfg,'N'), cfg.N=101; end
if ~isfield(cfg,'K'), cfg.K=10; end
lambda=cfg.lambda; a=cfg.a; k0=2*pi/lambda;
theta=deg2rad(cfg.theta_i_deg); kBloch=k0*sin(theta);
widths=cfg.widths(:).'; depths=cfg.depths(:).'; gaps=cfg.gaps(:).';
L=numel(widths); K=cfg.K; P=(cfg.N-1)/2;
orders=(-P:P).'; G=2*pi*orders/a; kx=kBloch+G;
ky=outgoing_sqrt(k0^2-kx.^2);

occupied=sum(widths)+sum(gaps);
if isfield(cfg,'x0'), x0=cfg.x0; else, x0=.5*(a-occupied); end
xleft=zeros(1,L); xleft(1)=x0;
for ell=2:L, xleft(ell)=xleft(ell-1)+widths(ell-1)+gaps(ell-1); end

nG=L*K;
B=complex(zeros(nG,cfg.N)); Cmap=complex(zeros(cfg.N,nG));
cosDepth=complex(zeros(nG,1)); betaSin=cosDepth;
verticalScale=zeros(nG,1); logVerticalScale=zeros(nG,1); betaList=cosDepth;
normCos=[1,.5*ones(1,K-1)];
row=0;
for ell=1:L
    t=widths(ell); d=depths(ell); xl=xleft(ell);
    for q=0:K-1
        row=row+1; alpha=q*pi/t; beta=groove_sqrt(k0^2-alpha^2);
        for n=1:cfg.N
            B(row,n)=projection(kx(n),q,xl,t)/normCos(q+1);
            Cmap(n,row)=(t/a)*projection(-kx(n),q,xl,t);
        end
        [cosDepth(row),betaSin(row),verticalScale(row), ...
            logVerticalScale(row)]=normalized_vertical(beta,k0,d);
        betaList(row)=beta;
    end
end

% Pressure continuity in every aperture mode, followed by global normal
% velocity continuity. With exp(+j*omega*t), rho*omega*v_y in a groove is
% -j*beta*C*sin(beta*d).
F=[B,-diag(cosDepth); diag(ky/k0),1i*Cmap*diag(betaSin)];
rowNorm=max(vecnorm(F,2,2),sqrt(eps));
Frow=F./rowNorm;
colNorm=max(vecnorm(Frow,2,1),sqrt(eps));
Fscaled=Frow./colNorm;

op=struct('F',F,'Fscaled',Fscaled,'orders',orders,'G',G,'kx',kx,'ky',ky, ...
    'k0',k0,'a',a,'lambda',lambda,'N',cfg.N,'K',K,'L',L, ...
    'xleft',xleft,'widths',widths,'depths',depths,'gaps',gaps, ...
    'cos_depth_normalized',cosDepth,'beta_sin_normalized',betaSin, ...
    'vertical_scale',verticalScale,'log_vertical_scale',logVerticalScale, ...
    'beta',betaList, ...
    'row_scale',rowNorm,'column_scale',colNorm, ...
    'n_floquet',cfg.N,'n_groove',nG);
end

function value=projection(kx,q,xl,t)
alpha=q*pi/t;
value=exp(-1i*kx*xl)*.5*(expint(kx-alpha,t)+expint(kx+alpha,t));
end
function value=expint(kappa,t)
z=kappa*t;
if abs(z)<1e-8
    value=1-1i*z/2-z^2/6;
else
    value=exp(-1i*z/2)*sin(z/2)/(z/2);
end
end
function value=outgoing_sqrt(z)
value=complex(zeros(size(z))); prop=real(z)>=0;
value(prop)=sqrt(real(z(prop))); value(~prop)=-1i*sqrt(-real(z(~prop)));
end
function value=groove_sqrt(z)
if real(z)>=0, value=sqrt(real(z)); else, value=1i*sqrt(-real(z)); end
end

function [cNormalized,sNormalized,scale,logScale]=normalized_vertical(beta,k0,d)
% Normalize cos(beta*d) and (beta/k0)sin(beta*d) without cosh overflow.
z=beta*d; growth=abs(imag(z));
if growth<300
    c=cos(z); s=(beta/k0)*sin(z);
    scale=max([abs(c),abs(s),sqrt(eps)]);
    cNormalized=c/scale; sNormalized=s/scale;
    logScale=log(scale);
    return;
end
decay=exp(-2*growth);
scaledCosh=.5*(1+decay);
scaledSinh=.5*sign(imag(z))*(1-decay);
x=real(z);
cScaled=cos(x)*scaledCosh-1i*sin(x)*scaledSinh;
sScaled=(beta/k0)*(sin(x)*scaledCosh+1i*cos(x)*scaledSinh);
scaledFloor=sqrt(eps)*exp(-growth);
scaledNorm=max([abs(cScaled),abs(sScaled),scaledFloor]);
cNormalized=cScaled/scaledNorm;
sNormalized=sScaled/scaledNorm;
logScale=growth+log(scaledNorm);
scale=exp(min(logScale,log(realmax)));
end
