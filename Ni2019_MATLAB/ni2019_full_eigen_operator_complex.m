function op=ni2019_full_eigen_operator_complex(cfg,k0,kBloch,targetOrder,targetKy)
%NI2019_FULL_EIGEN_OPERATOR_COMPLEX Full operator on a specified Riemann sheet.
%   k0 and kBloch have inverse-length units. targetKy explicitly selects the
%   sheet of the Rayleigh channel targetOrder. Open non-target channels are
%   analytically continued from positive real ky; closed channels are
%   continued from negative imaginary ky (decaying on the real axis).

if ~isfield(cfg,'N'), cfg.N=61; end
if ~isfield(cfg,'K'), cfg.K=10; end
a=cfg.a; widths=cfg.widths(:).'; depths=cfg.depths(:).'; gaps=cfg.gaps(:).';
L=numel(widths); K=cfg.K; P=(cfg.N-1)/2; orders=(-P:P).';
G=2*pi*orders/a; kx=kBloch+G; ky=complex(zeros(size(kx)));
for n=1:numel(orders)
    if orders(n)==targetOrder
        ky(n)=targetKy;
    else
        z=k0^2-kx(n)^2;
        zReference=real(k0)^2-kx(n)^2;
        if zReference>0
            ky(n)=sqrt(z);                 % outgoing open-channel sheet
            if real(ky(n))<0, ky(n)=-ky(n); end
        else
            ky(n)=-sqrt(z);                % decaying closed-channel sheet
            if imag(ky(n))>0, ky(n)=-ky(n); end
        end
    end
end

occupied=sum(widths)+sum(gaps);
if isfield(cfg,'x0'), x0=cfg.x0; else, x0=.5*(a-occupied); end
xleft=zeros(1,L); xleft(1)=x0;
for ell=2:L, xleft(ell)=xleft(ell-1)+widths(ell-1)+gaps(ell-1); end

nG=L*K; B=complex(zeros(nG,cfg.N)); Cmap=complex(zeros(cfg.N,nG));
cosDepth=complex(zeros(nG,1)); betaSin=cosDepth;
verticalScale=zeros(nG,1); logVerticalScale=zeros(nG,1); betaList=cosDepth;
normCos=[1,.5*ones(1,K-1)]; row=0;
for ell=1:L
    t=widths(ell); d=depths(ell); xl=xleft(ell);
    for q=0:K-1
        row=row+1; alpha=q*pi/t; beta=sqrt(k0^2-alpha^2);
        for n=1:cfg.N
            B(row,n)=projection(kx(n),q,xl,t)/normCos(q+1);
            Cmap(n,row)=(t/a)*projection(-kx(n),q,xl,t);
        end
        [cosDepth(row),betaSin(row),verticalScale(row), ...
            logVerticalScale(row)]=normalized_vertical(beta,k0,d);
        betaList(row)=beta;
    end
end

F=[B,-diag(cosDepth);diag(ky/k0),1i*Cmap*diag(betaSin)];
rowScale=max(vecnorm(F,2,2),sqrt(eps)); Frow=F./rowScale;
columnScale=max(vecnorm(Frow,2,1),sqrt(eps)); Fscaled=Frow./columnScale;
op=struct('F',F,'Fscaled',Fscaled,'orders',orders,'kx',kx,'ky',ky, ...
    'k0',k0,'kBloch',kBloch,'target_order',targetOrder,'target_ky',targetKy, ...
    'cos_depth_normalized',cosDepth,'beta_sin_normalized',betaSin, ...
    'vertical_scale',verticalScale,'log_vertical_scale',logVerticalScale, ...
    'beta',betaList, ...
    'row_scale',rowScale,'column_scale',columnScale,'N',cfg.N,'K',K,'L',L);
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


function [cNormalized,sNormalized,scale,logScale]=normalized_vertical(beta,k0,d)
% Scale out exp(abs(Im(beta*d))) before evaluating deep evanescent modes.
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
