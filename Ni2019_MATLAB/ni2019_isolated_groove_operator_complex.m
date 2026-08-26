function op=ni2019_isolated_groove_operator_complex(width,depth,K,k0,varargin)
%NI2019_ISOLATED_GROOVE_OPERATOR_COMPLEX Isolated groove in a rigid baffle.
%   The discrete Floquet radiation sum is replaced by its a->infinity
%   angular-spectrum integral, while all K transverse groove modes are
%   retained.  The exp(+i*omega*t) convention is used.

p=inputParser;
addParameter(p,'SMax',max(120,(K+12)*pi));
addParameter(p,'PanelsPerPi',8);
parse(p,varargin{:}); opt=p.Results;

z=k0*width;
% Continue the physical Sommerfeld contour together with the two branch
% points.  The central segment joins -z to +z; the two horizontal tails
% return to infinity.  This avoids the discontinuous channel-by-channel
% sheet switches produced by sampling the real kx axis at complex k0.
if abs(imag(z))<1e-14
    zContour=z+1i*1e-12*max(abs(z),1);
else
    zContour=z;
end
[phi,wPhi]=composite_gauss(-pi/2,pi/2,16);
[r,wTail]=composite_gauss(0,sqrt(opt.SMax), ...
    ceil(opt.SMax/pi*opt.PanelsPerPi));
sCentral=zContour*sin(phi);
sPositive=zContour+r.^2; sNegative=-zContour-r.^2;
PCentral=projection_matrix(sCentral,K);
PMinusCentral=projection_matrix(-sCentral,K);
PPositive=projection_matrix(sPositive,K);
PMinusPositive=projection_matrix(-sPositive,K);
PNegative=projection_matrix(sNegative,K);
PMinusNegative=projection_matrix(-sNegative,K);
kyPositive=sqrt(zContour^2-sPositive.^2);
kyNegative=sqrt(zContour^2-sNegative.^2);
integralCentral=transpose(PCentral)*((zContour*wPhi).*PMinusCentral);
tailJacobian=2*r;
integralPositive=transpose(PPositive)* ...
    ((wTail.*tailJacobian.*zContour./kyPositive).*PMinusPositive);
integralNegative=transpose(PNegative)* ...
    ((wTail.*tailJacobian.*zContour./kyNegative).*PMinusNegative);
normCos=[1,.5*ones(1,K-1)];
radiation=diag(1./normCos)* ...
    (integralCentral+integralPositive+integralNegative)/(2*pi);

q=(0:K-1).'; alpha=q*pi/width; beta=sqrt(k0^2-alpha.^2);
cosDepth=complex(zeros(K,1)); betaSin=cosDepth;
verticalScale=zeros(K,1); logVerticalScale=zeros(K,1);
for j=1:K
    [cosDepth(j),betaSin(j),verticalScale(j),logVerticalScale(j)]= ...
        normalized_vertical(beta(j),k0,depth);
end
F=-1i*radiation*diag(betaSin)-diag(cosDepth);
rowScale=max(vecnorm(F,2,2),sqrt(eps)); Frow=F./rowScale;
columnScale=max(vecnorm(Frow,2,1),sqrt(eps)); Fscaled=Frow./columnScale;
op=struct('F',F,'Fscaled',Fscaled,'radiation_matrix',radiation, ...
    'k0',k0,'width',width,'depth',depth,'K',K,'q',q, ...
    'beta',beta,'cos_depth_normalized',cosDepth, ...
    'beta_sin_normalized',betaSin,'vertical_scale',verticalScale, ...
    'log_vertical_scale',logVerticalScale,'row_scale',rowScale, ...
    'column_scale',columnScale,'quadrature_smax',opt.SMax, ...
    'quadrature_panels_per_pi',opt.PanelsPerPi);
end

function [s,w]=composite_gauss(left,right,nPanels)
% Sixteen-point Gauss-Legendre rule, composited uniformly in s.
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

function P=projection_matrix(s,K)
P=complex(zeros(numel(s),K));
for q=0:K-1
    P(:,q+1)=.5*(stable_expint(s-q*pi)+stable_expint(s+q*pi));
end
end

function value=stable_expint(z)
value=ones(size(z)); regular=abs(z)>=1e-8;
value(regular)=exp(-1i*z(regular)/2).* ...
    (sin(z(regular)/2)./(z(regular)/2));
zs=z(~regular);
value(~regular)=1-1i*zs/2-zs.^2/6;
end

function [cNormalized,sNormalized,scale,logScale]=normalized_vertical(beta,k0,d)
z=beta*d; growth=abs(imag(z));
if growth<300
    c=cos(z); s=(beta/k0)*sin(z);
    scale=max([abs(c),abs(s),sqrt(eps)]);
    cNormalized=c/scale; sNormalized=s/scale; logScale=log(scale);
    return;
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
