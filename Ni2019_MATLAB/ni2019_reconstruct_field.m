function field = ni2019_reconstruct_field(result,x,y,component)
%NI2019_RECONSTRUCT_FIELD Reconstruct pressure and normalized velocity.
%   FIELD = ... (RESULT,X,Y) evaluates the total field on meshgrid vectors
%   X and Y. Above the surface (y>=0), incident plus all retained Floquet
%   modes are used. Inside each groove, the cosine waveguide expansion is
%   used. Points in the rigid material below y=0 are returned as NaN.
%
%   FIELD = ... (RESULT,X,Y,'scattered') omits the incident wave above the
%   surface, matching the reflected-field panels in Figs. 2 and 5.
%
%   Velocity is normalized as rho*omega*v, so intensity can be plotted up
%   to the common positive factor 1/(rho*omega):
%       Ix ~ real(p*conj(vx_normalized))/2.

if nargin<4, component='total'; end
component = validatestring(component,{'total','scattered'});
[X,Y] = meshgrid(x,y);
p = nan(size(X)); vx = p; vy = p;
above = Y>=0;

if strcmp(component,'total')
    pa = exp(-1i*result.kx(result.orders==0).*X(above) + ...
             1i*result.ky_incident.*Y(above));
    vxa = result.kx(result.orders==0).*pa;
    vya = -result.ky_incident.*pa;
else
    pa = complex(zeros(nnz(above),1)); vxa=pa; vya=pa;
end
for n = 1:numel(result.orders)
    pr = result.A(n).*exp(-1i*result.kx(n).*X(above) - ...
                         1i*result.ky(n).*Y(above));
    pa = pa+pr;
    vxa = vxa+result.kx(n).*pr;
    vya = vya+result.ky(n).*pr;
end
p(above)=pa; vx(above)=vxa; vy(above)=vya;

for ell = 1:numel(result.widths)
    inside = Y<0 & Y>=-result.depths(ell) & ...
        X>=result.xleft(ell) & X<=result.xleft(ell)+result.widths(ell);
    if ~any(inside,'all'), continue; end
    u = X(inside)-result.xleft(ell);
    yy = Y(inside);
    pg = complex(zeros(size(u))); vxg=pg; vyg=pg;
    for q = 0:result.K-1
        alpha = q*pi/result.widths(ell);
        beta = groove_beta(result.k0^2-alpha^2);
        Pq = result.groove_surface_coefficients(q+1,ell);
        vertical = cos(beta*(yy+result.depths(ell)))/cos(beta*result.depths(ell));
        pg = pg + Pq*cos(alpha*u).*vertical;
        vxg = vxg - 1i*alpha*Pq*sin(alpha*u).*vertical;
        vyg = vyg - 1i*beta*Pq*cos(alpha*u).* ...
            sin(beta*(yy+result.depths(ell)))/cos(beta*result.depths(ell));
    end
    p(inside)=pg; vx(inside)=vxg; vy(inside)=vyg;
end

field = struct('x',x,'y',y,'X',X,'Y',Y,'p',p,'vx',vx,'vy',vy, ...
    'Ix',real(p.*conj(vx))/2,'Iy',real(p.*conj(vy))/2, ...
    'component',component);
end

function b = groove_beta(z)
if real(z)>=0, b=sqrt(real(z)); else, b=1i*sqrt(-real(z)); end
end
