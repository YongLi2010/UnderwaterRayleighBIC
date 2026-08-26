function roots = ni2019_discover_outgoing_poles_kappa(cfg,kappa,varargin)
%NI2019_DISCOVER_OUTGOING_POLES_KAPPA Multistart pole discovery in a window.
p=inputParser;
addParameter(p,'OmegaRange',[0.02,1.20]);
addParameter(p,'RealSeeds',0.05:0.04:1.20);
addParameter(p,'ImagSeeds',[0.001,0.015,0.06,0.20]);
addParameter(p,'ComplexSeeds',complex([]));
addParameter(p,'MergeTolerance',4e-4);
addParameter(p,'SigmaTolerance',2e-7);
addParameter(p,'RawTolerance',2e-4);
addParameter(p,'Verbose',true);
parse(p,varargin{:}); opt=p.Results;

cartesianSeeds=complex(zeros(1,numel(opt.RealSeeds)*numel(opt.ImagSeeds)));
count=0;
for ir=1:numel(opt.RealSeeds)
    for ii=1:numel(opt.ImagSeeds)
        count=count+1;
        cartesianSeeds(count)=complex(opt.RealSeeds(ir),opt.ImagSeeds(ii));
    end
end
allSeeds=[opt.ComplexSeeds(:).',cartesianSeeds];
roots=struct([]); nTotal=numel(allSeeds); count=0;
if opt.Verbose
    fprintf('Pole discovery at kappa=%+.4f: %d seeds\n',kappa,nTotal);
end
for iseed=1:numel(allSeeds)
        count=count+1;
        seed=allSeeds(iseed);
        try
            candidate=ni2019_refine_outgoing_pole_kappa(cfg,kappa,seed, ...
                'OuterIterations',7,'Display','off');
        catch
            continue;
        end
        om=candidate.Omega;
        valid=isfinite(real(om)) && isfinite(imag(om)) && ...
            real(om)>opt.OmegaRange(1) && real(om)<opt.OmegaRange(2) && ...
            imag(om)>-1e-7 && imag(om)<0.9 && ...
            candidate.sigma_ratio<opt.SigmaTolerance && ...
            candidate.raw_residual<opt.RawTolerance;
        if ~valid, continue; end
        if isempty(roots)
            roots=candidate;
        else
            distance=arrayfun(@(r)abs(r.Omega-om),roots);
            [dmin,id]=min(distance);
            if dmin<opt.MergeTolerance
                if candidate.sigma_ratio<roots(id).sigma_ratio
                    roots(id)=candidate;
                end
            else
                roots(end+1)=candidate; %#ok<AGROW>
            end
        end
    if opt.Verbose && mod(iseed,24)==0
        fprintf('  %d/%d seeds, %d candidates\n',count,nTotal,numel(roots));
    end
end
if ~isempty(roots)
    [~,id]=sort(real([roots.Omega])); roots=roots(id);
end
if opt.Verbose
    fprintf('  accepted %d poles at kappa=%+.4f\n',numel(roots),kappa);
    for j=1:numel(roots)
        fprintf('    %2d Omega=% .8f%+.8fi Q=%.3g\n',j, ...
            real(roots(j).Omega),imag(roots(j).Omega),roots(j).Q);
    end
end
end
