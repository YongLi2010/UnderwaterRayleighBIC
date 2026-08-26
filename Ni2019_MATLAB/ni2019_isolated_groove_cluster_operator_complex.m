function op=ni2019_isolated_groove_cluster_operator_complex( ...
    widths,depths,gaps,K,k0,varargin)
%NI2019_ISOLATED_GROOVE_CLUSTER_OPERATOR_COMPLEX Open-groove cluster QNM.
%   High-order groove modes are matched to the continuous outgoing angular
%   spectrum above one infinite rigid baffle.  There is no periodic copy.

p=inputParser;
addParameter(p,'SMax',max(300,(K+12)*pi));
addParameter(p,'PanelsPerPi',4);
addParameter(p,'ReferenceLength',min(widths));
parse(p,varargin{:}); opt=p.Results;
widths=widths(:).'; depths=depths(:).'; gaps=gaps(:).';
L=numel(widths);
if numel(depths)~=L || numel(gaps)~=max(L-1,0)
    error('widths/depths/gaps dimensions are inconsistent.');
end
aRef=opt.ReferenceLength; occupied=sum(widths)+sum(gaps);
xleft=zeros(1,L); xleft(1)=-occupied/2;
for ell=2:L
    xleft(ell)=xleft(ell-1)+widths(ell-1)+gaps(ell-1);
end
z=k0*aRef;
if abs(imag(z))<1e-14
    zContour=z+1i*1e-12*max(abs(z),1);
else
    zContour=z;
end
[phi,wPhi]=composite_gauss(-pi/2,pi/2,16);
[r,wTail]=composite_gauss(0,sqrt(opt.SMax), ...
    ceil(opt.SMax/pi*opt.PanelsPerPi));
sC=zContour*sin(phi); sP=zContour+r.^2; sN=-zContour-r.^2;
[PC,MC]=cluster_projection(sC,aRef,widths,xleft,K);
[PP,MP]=cluster_projection(sP,aRef,widths,xleft,K);
[PN,MN]=cluster_projection(sN,aRef,widths,xleft,K);
kyP=sqrt(zContour^2-sP.^2); kyN=sqrt(zContour^2-sN.^2);
IC=transpose(PC)*((zContour*wPhi).*MC);
tailJacobian=2*r;
IP=transpose(PP)*((wTail.*tailJacobian.*zContour./kyP).*MP);
IN=transpose(PN)*((wTail.*tailJacobian.*zContour./kyN).*MN);

nG=L*K; normByMode=zeros(nG,1); widthByMode=zeros(nG,1);
beta=complex(zeros(nG,1)); cosDepth=beta; betaSin=beta;
verticalScale=zeros(nG,1); logVerticalScale=zeros(nG,1);
row=0;
for ell=1:L
    for q=0:K-1
        row=row+1; normByMode(row)=1/(1+(q>0));
        widthByMode(row)=widths(ell);
        beta(row)=sqrt(k0^2-(q*pi/widths(ell))^2);
        [cosDepth(row),betaSin(row),verticalScale(row), ...
            logVerticalScale(row)]=normalized_vertical( ...
            beta(row),k0,depths(ell));
    end
end
radiation=diag(1./normByMode)*(IC+IP+IN)* ...
    diag(widthByMode/aRef)/(2*pi);
F=-1i*radiation*diag(betaSin)-diag(cosDepth);
rowScale=max(vecnorm(F,2,2),sqrt(eps)); Frow=F./rowScale;
columnScale=max(vecnorm(Frow,2,1),sqrt(eps)); Fscaled=Frow./columnScale;
op=struct('F',F,'Fscaled',Fscaled,'radiation_matrix',radiation, ...
    'k0',k0,'widths',widths,'depths',depths,'gaps',gaps, ...
    'xleft',xleft,'K',K,'L',L,'beta',beta, ...
    'cos_depth_normalized',cosDepth,'beta_sin_normalized',betaSin, ...
    'vertical_scale',verticalScale,'log_vertical_scale',logVerticalScale, ...
    'row_scale',rowScale,'column_scale',columnScale, ...
    'reference_length',aRef,'quadrature_smax',opt.SMax, ...
    'quadrature_panels_per_pi',opt.PanelsPerPi);
end

function [Pplus,Pminus]=cluster_projection(s,aRef,widths,xleft,K)
L=numel(widths); Pplus=complex(zeros(numel(s),L*K));
Pminus=Pplus; column=0; kx=s/aRef;
for ell=1:L
    t=widths(ell); xl=xleft(ell);
    for q=0:K-1
        column=column+1;
        localPlus=.5*(stable_expint(kx*t-q*pi)+ ...
            stable_expint(kx*t+q*pi));
        localMinus=.5*(stable_expint(-kx*t-q*pi)+ ...
            stable_expint(-kx*t+q*pi));
        Pplus(:,column)=exp(-1i*kx*xl).*localPlus;
        Pminus(:,column)=exp(1i*kx*xl).*localMinus;
    end
end
end

function [s,w]=composite_gauss(left,right,nPanels)
xg=[-.989400934991650;-.944575023073233;-.865631202387832; ...
    -.755404408355003;-.617876244402644;-.458016777657227; ...
    -.281603550779259;-.095012509837637;.095012509837637; ...
    .281603550779259;.458016777657227;.617876244402644; ...
    .755404408355003;.865631202387832;.944575023073233; ...
    .989400934991650];
wg=[.027152459411754;.062253523938648;.095158511682493; ...
    .124628971255534;.149595988816577;.169156519395003; ...
    .182603415044924;.189450610455069;.189450610455069; ...
    .182603415044924;.169156519395003;.149595988816577; ...
    .124628971255534;.095158511682493;.062253523938648; ...
    .027152459411754];
edges=linspace(left,right,nPanels+1);
centers=.5*(edges(1:end-1)+edges(2:end)); half=.5*diff(edges);
s=reshape(ones(numel(xg),1)*centers+xg*half,[],1);
w=reshape(wg*half,[],1);
end

function value=stable_expint(z)
value=ones(size(z)); regular=abs(z)>=1e-8;
value(regular)=exp(-1i*z(regular)/2).* ...
    (sin(z(regular)/2)./(z(regular)/2));
zs=z(~regular); value(~regular)=1-1i*zs/2-zs.^2/6;
end

function [cNormalized,sNormalized,scale,logScale]=normalized_vertical(beta,k0,d)
z=beta*d; growth=abs(imag(z));
if growth<300
    c=cos(z); s=(beta/k0)*sin(z);
    scale=max([abs(c),abs(s),sqrt(eps)]);
    cNormalized=c/scale; sNormalized=s/scale; logScale=log(scale); return;
end
decay=exp(-2*growth); scaledCosh=.5*(1+decay);
scaledSinh=.5*sign(imag(z))*(1-decay); x=real(z);
cScaled=cos(x)*scaledCosh-1i*sin(x)*scaledSinh;
sScaled=(beta/k0)*(sin(x)*scaledCosh+1i*cos(x)*scaledSinh);
scaledFloor=sqrt(eps)*exp(-growth);
scaledNorm=max([abs(cScaled),abs(sScaled),scaledFloor]);
cNormalized=cScaled/scaledNorm; sNormalized=sScaled/scaledNorm;
logScale=growth+log(scaledNorm);
scale=exp(min(logScale,log(realmax)));
end
