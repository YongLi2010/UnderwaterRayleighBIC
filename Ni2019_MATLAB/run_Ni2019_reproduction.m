%% Reproduce diffraction efficiencies reported by Ni et al. (2019)
% H. Ni et al., Phys. Rev. B 100, 104104 (2019), Fig. 2.
% This script uses the published geometries and an independent
% modal-matching implementation. No optimization or COMSOL is required.

clear; close all; clc;

base.lambda = 1;           % normalized wavelength
base.theta_i_deg = 0;
base.N = 101;              % paper: N=101 Floquet orders
base.K = 10;               % paper: 10 groove modes

cases = published_cases(base);
results = cell(size(cases));

fprintf('Ni2019 modal-matching reproduction (A_inc = 1)\n');
fprintf('%-15s %8s %10s %10s %12s\n', ...
    'case','order','angle(deg)','eta','sum(eta)');
fprintf('%s\n', repmat('-',1,62));

for c = 1:numel(cases)
    results{c} = ni2019_modal_solver(cases(c).cfg);
    r = results{c};
    ids = find(r.is_propagating);
    for m = 1:numel(ids)
        id = ids(m);
        if m == 1
            label = cases(c).name;
        else
            label = '';
        end
        fprintf('%-15s %+8d %10.3f %10.6f %12s\n', label, ...
            r.orders(id), r.theta_deg(id), r.eta(id), ...
            ternary(m==1, sprintf('%.8f',r.total_efficiency), ''));
    end
    fprintf('  energy error = %.3e, cond = %.3e\n\n', ...
        r.energy_error, r.condition_number);
end

fprintf('Symmetric splitter ratios (negative order / positive order)\n');
for c = 3:4
    r = results{c};
    im = find(r.orders==-1,1); ip = find(r.orders==1,1);
    fprintf('  %-15s  |A_-1|/|A_+1| = %.6f,  eta_-1/eta_+1 = %.6f\n', ...
        cases(c).name, abs(r.A(im))/abs(r.A(ip)), r.eta(im)/r.eta(ip));
end
fprintf(['NOTE: the published "5:3" geometry reproduces an amplitude ratio near 5:3;\n' ...
         '      its normal-power ratio is therefore near 25:9, not 5:3.\n\n']);

% Compact efficiency comparison for all published structures.
figure('Color','w','Name','Ni2019 efficiency reproduction');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for c = 1:numel(cases)
    nexttile;
    r = results{c};
    id = r.is_propagating;
    bar(r.orders(id),r.eta(id),0.62,'FaceColor',[0.10 0.45 0.75]);
    xlabel('Diffraction order n'); ylabel('\eta_n');
    title(sprintf('%s,  \\Sigma\\eta=%.6f',cases(c).name,r.total_efficiency));
    ylim([0 1.05]); grid on; box on;
end

% Truncation check for the most extreme (-88 deg) example.
Nlist = [21 41 61 81 101 121];
Klist = [4 6 8 10 12];
etaTarget = nan(numel(Klist),numel(Nlist));
energyError = etaTarget;
for ik = 1:numel(Klist)
    for in = 1:numel(Nlist)
        cfg = cases(2).cfg;
        cfg.N = Nlist(in); cfg.K = Klist(ik);
        rr = ni2019_modal_solver(cfg);
        etaTarget(ik,in) = rr.eta(rr.orders==-1);
        energyError(ik,in) = rr.energy_error;
    end
end
figure('Color','w','Name','Ni2019 truncation convergence');
subplot(1,2,1);
plot(Nlist,etaTarget,'o-','LineWidth',1.2);
xlabel('Floquet truncation N'); ylabel('\eta_{-1}'); grid on;
legend(compose('K=%d',Klist),'Location','best'); title('-88 deg target efficiency');
subplot(1,2,2);
semilogy(Nlist,max(energyError,eps),'o-','LineWidth',1.2);
xlabel('Floquet truncation N'); ylabel('|1-\Sigma\eta_n|'); grid on;
legend(compose('K=%d',Klist),'Location','best'); title('Energy-conservation error');

assignin('base','Ni2019_cases',cases);
assignin('base','Ni2019_results',results);

function cases = published_cases(base)
% dx in the paper is the edge-to-edge distance between neighboring grooves.
theta = [81 88 72 72];
names = {'-81 deg','-88 deg','+/-72 deg 5:3','+/-72 deg 1:1'};
widthFrac = {[0.500 0.148],[0.085 0.527],[0.159 0.230],[0.094 0.229]};
depthFrac = {[0.399 0.121],[0.419 0.070],[0.229 0.170],[0.207 0.205]};
gapFrac = [0.111 0.101 0.176 0.206];
cases = struct('name',names,'cfg',cell(1,4));
for c = 1:4
    cfg = base;
    cfg.a = base.lambda/sind(theta(c));
    cfg.widths = widthFrac{c}*cfg.a;
    cfg.depths = depthFrac{c}*base.lambda;
    cfg.gaps = gapFrac(c)*cfg.a;
    cases(c).cfg = cfg;
end
end

function y = ternary(condition,a,b)
if condition, y = a; else, y = b; end
end
