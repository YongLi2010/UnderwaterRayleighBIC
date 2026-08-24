function result=ni2019_refine_two_symmetric_groove_double_rayleigh_bic(cfg,y0,varargin)
%NI2019_REFINE_TWO_SYMMETRIC_GROOVE_DOUBLE_RAYLEIGH_BIC Refine an odd root.
%   R=NI2019_REFINE_TWO_SYMMETRIC_GROOVE_DOUBLE_RAYLEIGH_BIC(CFG,Y0)
%   refines the normalized geometry
%
%       Y=[d/a,w/a,g/a]
%
%   of two identical grooves centered in a period.  The frequency and Bloch
%   number are fixed at the simultaneous Rayleigh point
%
%       kappa=0,  Omega=a/lambda=1.
%
%   The full modal-matching operator is restricted to the mirror-odd space
%
%       A(-n)=-A(n),
%       C(q,2)=-(-1)^q C(q,1).
%
%   A(-1), A(0), and A(+1) are omitted identically.  Consequently the
%   optimization tests the strict double-Rayleigh BIC compatibility
%   condition without allowing an even branch or nonzero grazing pressure.
%
%   Name/value options:
%       Truncation              [N K], default [cfg.N cfg.K]
%       LowerBounds             lower bounds for [d/a,w/a,g/a]
%       UpperBounds             upper bounds for [d/a,w/a,g/a]
%       FillMax                 upper bound for 2*w/a+g/a, default .96
%       FillPenaltyWeight       least-squares penalty weight, default 1e5
%       MaxFunctionEvaluations  default 1800
%       MaxIterations           default 250
%       FunctionTolerance       default 1e-14
%       StepTolerance           default 1e-13
%       OptimalityTolerance     default 1e-13
%       FiniteDifferenceType    default 'central'
%       Display                 default 'off'

p=inputParser;
addParameter(p,'Truncation',[],@(x)isnumeric(x)&&isequal(size(x),[1 2]));
addParameter(p,'LowerBounds',[.03 .05 .01], ...
    @(x)isnumeric(x)&&numel(x)==3&&all(isfinite(x)));
addParameter(p,'UpperBounds',[1.20 .49 .40], ...
    @(x)isnumeric(x)&&numel(x)==3&&all(isfinite(x)));
addParameter(p,'FillMax',.96,@(x)isscalar(x)&&isfinite(x)&&x>0);
addParameter(p,'FillPenaltyWeight',1e5,@(x)isscalar(x)&&isfinite(x)&&x>0);
addParameter(p,'MaxFunctionEvaluations',1800,@(x)isscalar(x)&&x>=10);
addParameter(p,'MaxIterations',250,@(x)isscalar(x)&&x>=1);
addParameter(p,'FunctionTolerance',1e-14,@(x)isscalar(x)&&x>0);
addParameter(p,'StepTolerance',1e-13,@(x)isscalar(x)&&x>0);
addParameter(p,'OptimalityTolerance',1e-13,@(x)isscalar(x)&&x>0);
addParameter(p,'FiniteDifferenceType','central',@(x)ischar(x)||isstring(x));
addParameter(p,'Display','off',@(x)ischar(x)||isstring(x));
parse(p,varargin{:});
opt=p.Results;

if ~isfield(cfg,'a')||~isscalar(cfg.a)||~isfinite(cfg.a)||cfg.a<=0
    error('cfg.a must be a positive finite scalar.');
end
if isempty(opt.Truncation)
    if ~isfield(cfg,'N')||~isfield(cfg,'K')
        error('Supply Truncation or provide cfg.N and cfg.K.');
    end
    truncation=[cfg.N cfg.K];
else
    truncation=round(opt.Truncation);
end
if mod(truncation(1),2)~=1||truncation(1)<5||truncation(2)<2
    error('Truncation requires odd N>=5 and K>=2.');
end
if exist('lsqnonlin','file')~=2
    error('lsqnonlin from Optimization Toolbox is required.');
end

lb=opt.LowerBounds(:).';
ub=opt.UpperBounds(:).';
y0=y0(:).';
if any(lb>=ub)||any(y0<lb)||any(y0>ub)
    error('Bounds must be ordered and contain y0.');
end
if 2*y0(2)+y0(3)>opt.FillMax
    error('Initial geometry exceeds FillMax.');
end

% Fix the arbitrary singular-vector phase at one groove-coordinate entry.
[~,initialState]=projected_residual(y0,[]);
grooveStart=numel(initialState.positive_orders)+1;
[~,relativeAnchor]=max(abs(initialState.v(grooveStart:end)));
anchor=grooveStart+relativeAnchor-1;

historyY=zeros(0,3);
historyNorm=zeros(0,1);
lsqopt=optimoptions('lsqnonlin','Display',char(opt.Display), ...
    'MaxFunctionEvaluations',round(opt.MaxFunctionEvaluations), ...
    'MaxIterations',round(opt.MaxIterations), ...
    'FunctionTolerance',opt.FunctionTolerance, ...
    'StepTolerance',opt.StepTolerance, ...
    'OptimalityTolerance',opt.OptimalityTolerance, ...
    'FiniteDifferenceType',char(opt.FiniteDifferenceType));

[y,resnorm,residual,exitflag,output]=lsqnonlin(@objective,y0,lb,ub,lsqopt);
[compatibility,state]=projected_residual(y,anchor);
strict=ni2019_strict_rayleigh_operator(state.local_cfg, ...
    'TargetOrder',-1,'EnforceOtherThreshold',true);

result=struct();
result.y=y;
result.depth_over_a=y(1);
result.width_over_a=y(2);
result.gap_over_a=y(3);
result.depths_over_a=[y(1) y(1)];
result.widths_over_a=[y(2) y(2)];
result.gaps_over_a=y(3);
result.fill_fraction=2*y(2)+y(3);
result.kappa=0;
result.Omega=1;
result.theta_deg=0;
result.truncation=truncation;
result.compatibility_residual=compatibility;
result.compatibility_norm=norm(compatibility);
result.sigma_ratio=state.sigma_ratio;
result.sigma_min=state.sigma_min;
result.raw_residual=state.raw_residual;
result.odd_operator=state;
result.strict_operator=strict;
result.removed_orders=[-1 0 1];
result.mode=state.mode;
result.transverse=state.transverse;
result.resnorm=resnorm;
result.residual=residual;
result.exitflag=exitflag;
result.output=output;
result.anchor_index_odd=anchor;
result.history=struct('y',historyY,'compatibility_norm',historyNorm);
result.options=opt;

    function r=objective(yTrial)
        [rComplex,~]=projected_residual(yTrial,anchor);
        fillViolation=max(2*yTrial(2)+yTrial(3)-opt.FillMax,0);
        r=[real(rComplex);imag(rComplex); ...
            sqrt(opt.FillPenaltyWeight)*fillViolation];
        historyY(end+1,:)=yTrial(:).';
        historyNorm(end+1,1)=norm(r);
    end

    function [r,state]=projected_residual(yTrial,anchorIndex)
        local=cfg;
        local.lambda=cfg.a;
        local.theta_i_deg=0;
        local.depths=[yTrial(1) yTrial(1)]*cfg.a;
        local.widths=[yTrial(2) yTrial(2)]*cfg.a;
        local.gaps=yTrial(3)*cfg.a;
        local.N=truncation(1);
        local.K=truncation(2);
        local.solve_scattering=false;
        local=rmfield_if_present(local,'x0');
        op=ni2019_full_eigen_operator(local);

        if op.L~=2
            error('The local operator must contain exactly two grooves.');
        end
        kyTolerance=1e-8*max(abs(op.k0),1);
        finiteOpen=abs(imag(op.ky))<=kyTolerance&real(op.ky)>kyTolerance;
        threshold=abs(op.ky)<=kyTolerance;
        removedOrders=op.orders(finiteOpen|threshold);
        if ~isequal(removedOrders(:),[-1;0;1])
            error('Expected removed orders [-1 0 1], obtained [%s].', ...
                strtrim(sprintf('%d ',removedOrders)));
        end

        positiveOrders=(2:(op.N-1)/2).';
        nOddA=numel(positiveOrders);
        T=complex(zeros(op.N+op.n_groove,nOddA+op.K));
        for jj=1:nOddA
            idPlus=find(op.orders==positiveOrders(jj),1);
            idMinus=find(op.orders==-positiveOrders(jj),1);
            T(idPlus,jj)=1/sqrt(2);
            T(idMinus,jj)=-1/sqrt(2);
        end
        for q=0:op.K-1
            col=nOddA+q+1;
            T(op.N+q+1,col)=1/sqrt(2);
            T(op.N+op.K+q+1,col)=-(-1)^q/sqrt(2);
        end

        FoddFull=op.F*T;
        independentRows=[(1:op.K).';2*op.K+find(op.orders>=1)];
        Fodd=FoddFull(independentRows,:);
        rowScale=max(vecnorm(Fodd,2,2),sqrt(eps));
        Frow=Fodd./rowScale;
        columnScale=max(vecnorm(Frow,2,1),sqrt(eps));
        Fscaled=Frow./columnScale;
        if any(~isfinite(Fscaled(:)))
            error('Odd reduced operator contains NaN or Inf.');
        end
        [~,S,V]=svd(Fscaled,'econ');
        singularValues=diag(S);
        v=V(:,end);
        if ~isempty(anchorIndex)
            if anchorIndex>numel(v)||abs(v(anchorIndex))<100*eps
                error('The fixed phase anchor became singular during refinement.');
            end
            v=v*exp(-1i*angle(v(anchorIndex)));
        end
        sigmaMax=max(singularValues(1),eps);
        r=(Fscaled*v)/sigmaMax;

        uOdd=v./columnScale(:);
        z=T*uOdd;
        z=z/max(norm(z),eps);
        rawResidual=norm(op.F*z)/max(norm(z),eps);
        A=z(1:op.N);
        Cscaled=z(op.N+1:end);
        Cphysical=Cscaled./op.vertical_scale;
        CbyGroove=reshape(Cphysical,op.K,2);
        pressureSurface=Cscaled.*op.cos_depth_normalized;
        velocitySurface=Cscaled.*op.beta_sin_normalized;
        pressureByGroove=reshape(pressureSurface,op.K,2);
        velocityByGroove=reshape(velocitySurface,op.K,2);
        qNorm=sqrt(sum(abs(CbyGroove).^2,2));
        pressureQNorm=sqrt(sum(abs(pressureByGroove).^2,2));
        velocityQNorm=sqrt(sum(abs(velocityByGroove).^2,2));
        qFraction=qNorm.^2/max(sum(qNorm.^2),eps);
        pressureQFraction=pressureQNorm.^2/ ...
            max(sum(pressureQNorm.^2),eps);
        velocityQFraction=velocityQNorm.^2/ ...
            max(sum(velocityQNorm.^2),eps);
        mode=struct('z_full',z,'A',A,'C_scaled',Cscaled, ...
            'C_physical',Cphysical,'C_by_groove',CbyGroove, ...
            'surface_pressure_coefficients',pressureSurface, ...
            'surface_velocity_coefficients',velocitySurface);
        transverse=struct('q',(0:op.K-1).','norm_by_q',qNorm, ...
            'fraction_by_q',qFraction, ...
            'higher_order_fraction',sum(qFraction(2:end)), ...
            'q1_fraction',qFraction(2), ...
            'surface_pressure_norm_by_q',pressureQNorm, ...
            'surface_pressure_fraction_by_q',pressureQFraction, ...
            'surface_velocity_norm_by_q',velocityQNorm, ...
            'surface_velocity_fraction_by_q',velocityQFraction, ...
            'surface_pressure_higher_order_fraction', ...
            sum(pressureQFraction(2:end)), ...
            'surface_velocity_higher_order_fraction', ...
            sum(velocityQFraction(2:end)));
        state=struct('local_cfg',local,'operator',op,'T_odd',T, ...
            'F_odd_full_raw',FoddFull,'F_odd_raw',Fodd, ...
            'F_odd_scaled',Fscaled,'independent_rows',independentRows, ...
            'row_scale',rowScale,'column_scale',columnScale, ...
            'v',v,'u_odd',uOdd,'singular_values',singularValues, ...
            'sigma_min',singularValues(end), ...
            'sigma_ratio',singularValues(end)/sigmaMax, ...
            'raw_residual',rawResidual,'positive_orders',positiveOrders, ...
            'mode',mode,'transverse',transverse);
    end
end

function cfg=rmfield_if_present(cfg,name)
if isfield(cfg,name)
    cfg=rmfield(cfg,name);
end
end
