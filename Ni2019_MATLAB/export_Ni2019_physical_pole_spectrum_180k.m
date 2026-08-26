%% Export a physically segmented outgoing-pole spectrum for Fig. 2b
% The nonlinear roots are reused from the independently seeded positive-k
% census.  This script does not force them into closed-system-like bands.
% Adjacent k points are joined only by a global one-to-one assignment that
% passes three hard gates: the same outgoing/decaying channel signature,
% strong groove-mode overlap, and a small complex-frequency displacement.
% The strict Rayleigh-BIC branch is exported separately elsewhere because
% its vanishing threshold amplitude allows it to cross a Rayleigh cut.
clear; clc;

rootDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(rootDir,'results','full_complex_band_180k');
D = load(fullfile(dataDir,'positive_pole_census.mat'), ...
    'rootsByK','kappa','cfg');

overlapGate = 0.95;
frequencyGate = 0.035;
unmatchedCost = 0.30;

pointRows = zeros(0,13);
pointIdByK = cell(numel(D.kappa),1);
pointCounter = 0;
for m = 1:numel(D.kappa)
    roots = D.rootsByK{m};
    ids = zeros(1,numel(roots));
    for r = 1:numel(roots)
        om = roots(r).Omega;
        if real(om)<=0 || real(om)>1.0005, continue; end
        pointCounter = pointCounter+1;
        ids(r) = pointCounter;
        [signature,openCount] = channel_signature(roots(r),D.kappa(m));
        cutDistance = min(abs(real(om)-abs(D.kappa(m)+(-1:1))));
        guided = double(openCount==0);
        pointRows(end+1,:) = [pointCounter,m,D.kappa(m),r, ... %#ok<SAGROW>
            real(om),imag(om),roots(r).Q,roots(r).sigma_ratio, ...
            roots(r).raw_residual,signature,openCount,guided,cutDistance];
    end
    pointIdByK{m} = ids;
end

segmentRows = zeros(0,14);
for m = 1:numel(D.kappa)-1
    left = D.rootsByK{m}; right = D.rootsByK{m+1};
    leftVisible = find(arrayfun(@(p) real(p.Omega)>0 && ...
        real(p.Omega)<=1.0005,left));
    rightVisible = find(arrayfun(@(p) real(p.Omega)>0 && ...
        real(p.Omega)<=1.0005,right));
    if isempty(leftVisible) || isempty(rightVisible), continue; end

    cost = inf(numel(leftVisible),numel(rightVisible));
    overlap = nan(size(cost)); deltaOmega = overlap;
    sameSignature = false(size(cost));
    for il = 1:numel(leftVisible)
        rl = left(leftVisible(il));
        [sigL,~] = channel_signature(rl,D.kappa(m));
        cL = rl.C_physical(:);
        for ir = 1:numel(rightVisible)
            rr = right(rightVisible(ir));
            [sigR,~] = channel_signature(rr,D.kappa(m+1));
            cR = rr.C_physical(:);
            overlap(il,ir) = abs(cL'*cR)/max(norm(cL)*norm(cR),eps);
            deltaOmega(il,ir) = abs(rr.Omega-rl.Omega);
            sameSignature(il,ir) = sigL==sigR;
            if sameSignature(il,ir) && overlap(il,ir)>=overlapGate && ...
                    deltaOmega(il,ir)<=frequencyGate
                cost(il,ir) = (1-overlap(il,ir)) + ...
                    0.12*(deltaOmega(il,ir)/frequencyGate)^2;
            end
        end
    end

    pairs = matchpairs(cost,unmatchedCost,'min');
    for p = 1:size(pairs,1)
        il = pairs(p,1); ir = pairs(p,2);
        if ~isfinite(cost(il,ir)), continue; end
        iRoot = leftVisible(il); jRoot = rightVisible(ir);
        id0 = pointIdByK{m}(iRoot); id1 = pointIdByK{m+1}(jRoot);
        if id0==0 || id1==0, continue; end
        om0 = left(iRoot).Omega; om1 = right(jRoot).Omega;
        [signature,~] = channel_signature(left(iRoot),D.kappa(m));
        segmentRows(end+1,:) = [id0,id1,D.kappa(m),D.kappa(m+1), ... %#ok<SAGROW>
            real(om0),imag(om0),real(om1),imag(om1), ...
            left(iRoot).Q,right(jRoot).Q,overlap(il,ir), ...
            deltaOmega(il,ir),signature,cost(il,ir)];
    end
end

points = array2table(pointRows,'VariableNames', ...
    {'point_id','k_index','kappa','root_index','Omega_real','Omega_imag', ...
    'Q','sigma_ratio','raw_residual','sheet_signature','open_count', ...
    'guided','cut_distance'});
segments = array2table(segmentRows,'VariableNames', ...
    {'point_id_0','point_id_1','kappa_0','kappa_1','Omega_real_0', ...
    'Omega_imag_0','Omega_real_1','Omega_imag_1','Q_0','Q_1', ...
    'mode_overlap','delta_Omega','sheet_signature','assignment_cost'});

writetable(points,fullfile(dataDir,'physical_outgoing_pole_points.csv'));
writetable(segments,fullfile(dataDir,'physical_outgoing_pole_segments.csv'));
save(fullfile(dataDir,'physical_outgoing_pole_spectrum.mat'), ...
    'points','segments','overlapGate','frequencyGate','unmatchedCost');

fprintf('Physical outgoing-pole spectrum exported.\n');
fprintf('  points: %d, accepted local segments: %d\n',height(points),height(segments));
fprintf('  minimum accepted overlap: %.8f\n',min(segments.mode_overlap));
fprintf('  maximum accepted |Delta Omega|: %.6g\n',max(segments.delta_Omega));
fprintf('  hard gates: overlap >= %.2f, |Delta Omega| <= %.3f\n', ...
    overlapGate,frequencyGate);

function [signature,openCount] = channel_signature(root,kappa)
% Encode the physical-boundary status of n=-1,0,+1 as a three-bit integer.
orders = -1:1;
isOpen = real(root.Omega).^2 > (kappa+orders).^2;
signature = sum(double(isOpen).*[1,2,4]);
openCount = sum(isOpen);
end
