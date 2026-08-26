%% Full complex Bloch band continuation for the 180-kHz two-groove design
% All bands entering 0<Re(Omega)<1.20 at either Brillouin-zone edge are
% continued independently toward Gamma. Positive and negative kappa are
% calculated separately as a reciprocity audit.
clear; clc;
rootDir=fileparts(mfilename('fullpath'));
outDir=fullfile(rootDir,'results','full_complex_band_180k');
S=load(fullfile(rootDir,'results','StrictRayleighBIC_180kHz_7p10deg_final.mat'));
E=load(fullfile(outDir,'edge_discovery.mat'));
cfg=S.cfg; cfg.N=89; cfg.K=11; cfg.solve_scattering=false;

step=0.005;
kPos=(0.5:-step:0).';
kNeg=(-0.5:step:0).';
fprintf('Tracking %d positive-k and %d negative-k points, %d edge roots.\n', ...
    numel(kPos),numel(kNeg),numel(E.rp));

positive=track_family(cfg,kPos,E.rp,'positive');
save(fullfile(outDir,'full_complex_band_checkpoint.mat'),'positive','cfg','kPos');
negative=track_family(cfg,kNeg,E.rn,'negative');

% Reorder to increasing kappa and combine without duplicating Gamma.
[kp,ip]=sort(positive.kappa); [kn,in]=sort(negative.kappa);
positive=reorder_family(positive,ip,kp);
negative=reorder_family(negative,in,kn);
kappa=[negative.kappa;positive.kappa(2:end)];
Omega=[negative.Omega;positive.Omega(2:end,:)];
Q=[negative.Q;positive.Q(2:end,:)];
sigma_ratio=[negative.sigma_ratio;positive.sigma_ratio(2:end,:)];
raw_residual=[negative.raw_residual;positive.raw_residual(2:end,:)];
overlap=[negative.overlap;positive.overlap(2:end,:)];
exitflag=[negative.exitflag;positive.exitflag(2:end,:)];

% Reciprocity audit: compare independently calculated +/- kappa spectra
% after sorting each frequency row.  The assignment-independent error is
% useful even if two branches exchange labels at a crossing.
nHalf=numel(positive.kappa); reciprocityError=nan(nHalf,1);
for m=1:nHalf
    kpValue=positive.kappa(m);
    [~,j]=min(abs(negative.kappa+kpValue));
    op=sort(positive.Omega(m,:)); on=sort(negative.Omega(j,:));
    reciprocityError(m)=max(abs(op-on));
end

% Rayleigh-BIC endpoint and reciprocal partner are exact strict results.
kappaBIC=S.x(1); OmegaBIC=1-kappaBIC;

save(fullfile(outDir,'full_complex_band.mat'),'cfg','kappa','Omega','Q', ...
    'sigma_ratio','raw_residual','overlap','exitflag','positive','negative', ...
    'reciprocityError','kappaBIC','OmegaBIC','step','-v7.3');

fprintf('\nFull-band continuation complete.\n');
fprintf('  max sigma ratio = %.3e\n',max(sigma_ratio,[],'all','omitnan'));
fprintf('  max raw residual = %.3e\n',max(raw_residual,[],'all','omitnan'));
fprintf('  max reciprocity error = %.3e\n',max(reciprocityError,[],'omitnan'));
for j=1:size(Omega,2)
    fprintf('  band %d: ReOmega range [%.5f, %.5f], min ImOmega %.3e\n',j, ...
        min(real(Omega(:,j)),[],'omitnan'),max(real(Omega(:,j)),[],'omitnan'), ...
        min(abs(imag(Omega(:,j))),[],'omitnan'));
end

function family=track_family(cfg,kappa,edgeRoots,label)
nK=numel(kappa); nB=numel(edgeRoots);
Omega=complex(nan(nK,nB)); Q=nan(nK,nB); sigma=Q; raw=Q; ov=Q; flags=Q;
previous=edgeRoots(:).'; previous2=previous;
for j=1:nB
    Omega(1,j)=previous(j).Omega; Q(1,j)=previous(j).Q;
    sigma(1,j)=previous(j).sigma_ratio; raw(1,j)=previous(j).raw_residual;
    ov(1,j)=1; flags(1,j)=previous(j).exitflag;
end
for m=2:nK
    current=repmat(previous(1),1,nB);
    for j=1:nB
        if m==2
            predictor=previous(j).Omega;
        else
            predictor=previous(j).Omega+(previous(j).Omega-previous2(j).Omega);
        end
        seeds=[predictor,previous(j).Omega, ...
            predictor+0.004i,predictor-0.004i, ...
            predictor+0.012,predictor-0.012];
        [candidate,success]=best_candidate(cfg,kappa(m),seeds,previous(j),predictor);
        if ~success
            warning('%s band %d failed at kappa=%+.5f',label,j,kappa(m));
            candidate=previous(j); candidate.kappa=kappa(m);
            candidate.Omega=complex(NaN,NaN); candidate.Q=NaN;
            candidate.sigma_ratio=NaN; candidate.raw_residual=NaN;
        end
        current(j)=candidate;
        Omega(m,j)=candidate.Omega; Q(m,j)=candidate.Q;
        sigma(m,j)=candidate.sigma_ratio; raw(m,j)=candidate.raw_residual;
        flags(m,j)=candidate.exitflag;
        if success
            ov(m,j)=mode_overlap(previous(j),candidate);
        end
    end
    % Report accidental branch collapse; retain the data for later audit
    % rather than silently perturbing or swapping labels.
    for j=1:nB
        for ell=j+1:nB
            if isfinite(Omega(m,j)) && abs(Omega(m,j)-Omega(m,ell))<2e-5
                warning('%s branches %d/%d coalesced at kappa=%+.5f', ...
                    label,j,ell,kappa(m));
            end
        end
    end
    previous2=previous; previous=current;
    if mod(m-1,10)==0 || m==nK
        fprintf('  %-8s %3d/%3d, kappa=%+.3f\n',label,m,nK,kappa(m));
    end
end
family=struct('kappa',kappa,'Omega',Omega,'Q',Q,'sigma_ratio',sigma, ...
    'raw_residual',raw,'overlap',ov,'exitflag',flags);
end

function [best,success]=best_candidate(cfg,kappa,seeds,previous,predictor)
best=previous; success=false; bestScore=Inf;
for s=1:numel(seeds)
    try
        candidate=ni2019_refine_outgoing_pole_kappa(cfg,kappa,seeds(s), ...
            'OuterIterations',7,'Display','off');
    catch
        continue;
    end
    om=candidate.Omega;
    valid=isfinite(real(om)) && isfinite(imag(om)) && real(om)>0.005 && ...
        real(om)<1.35 && imag(om)>-2e-6 && imag(om)<1.0 && ...
        candidate.sigma_ratio<2e-6 && candidate.raw_residual<2e-3;
    if ~valid, continue; end
    overlap=mode_overlap(previous,candidate);
    score=abs(om-predictor)/0.05+0.65*(1-overlap);
    if score<bestScore
        best=candidate; bestScore=score; success=true;
    end
end
end

function value=mode_overlap(a,b)
ca=a.C_physical(:); cb=b.C_physical(:);
value=abs(ca'*cb)/max(norm(ca)*norm(cb),eps);
value=min(max(real(value),0),1);
end

function out=reorder_family(in,id,kappa)
out=in; out.kappa=kappa;
fields={'Omega','Q','sigma_ratio','raw_residual','overlap','exitflag'};
for j=1:numel(fields), out.(fields{j})=in.(fields{j})(id,:); end
end
