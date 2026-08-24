function value=ni2019_route_objective(u,x,aPhysical,cWater,N,K,center,scale)
%NI2019_ROUTE_OBJECTIVE Negative anomalous-order efficiency for local search.
%   U is dimensionless. CENTER=[theta_deg,f_kHz] and SCALE gives the local
%   coordinate scale. A smooth quadratic penalty keeps the search in the
%   single-anomalous-order operating window.

theta=center(1)+scale(1)*u(1);
fKHz=center(2)+scale(2)*u(2);
penalty=0;
if theta<3.5, penalty=penalty+100*(3.5-theta)^2; end
if theta>8.0, penalty=penalty+100*(theta-8.0)^2; end
if fKHz<198.5, penalty=penalty+100*(198.5-fKHz)^2; end
if fKHz>202.0, penalty=penalty+100*(fKHz-202.0)^2; end

f=1e3*fKHz;
Omega=f*aPhysical/cWater;
fRayleigh=cWater/(aPhysical*(1+sind(theta)));
fOpposite=cWater/(aPhysical*(1-sind(theta)));
if f<=fRayleigh
    penalty=penalty+1e-4*(fRayleigh-f+1)^2;
end
if f>=fOpposite
    penalty=penalty+1e-4*(f-fOpposite+1)^2;
end

cfg=struct('a',1,'lambda',1/Omega,'theta_i_deg',theta, ...
    'depths',x(2:3),'widths',x(4:5),'gaps',x(6), ...
    'N',N,'K',K,'solve_scattering',true);
try
    R=ni2019_modal_solver(cfg);
    eta=R.eta(R.orders==-1);
    if ~isfinite(eta) || R.energy_error>1e-7
        value=10+penalty;
    else
        value=-eta+penalty;
    end
catch
    value=10+penalty;
end
end
