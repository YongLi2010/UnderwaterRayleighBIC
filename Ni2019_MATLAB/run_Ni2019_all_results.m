%% Universal reproduction and paper-comparison program for Ni et al. (2019)
% Covers every distinct printed geometry in Figs. 2-6.

clear; close all; clc;
cases = ni2019_case_library('lambda',1,'N',101,'K',10);
results = cell(size(cases));
outputDir = fullfile(pwd,'results');
if ~exist(outputDir,'dir'), mkdir(outputDir); end

rowTemplate = struct('case_id','','figure','','order',NaN,'paper_angle_deg',NaN, ...
    'calculated_angle_deg',NaN,'angle_error_deg',NaN,'paper_abs_A',NaN, ...
    'calculated_abs_A',NaN,'amplitude_relative_error',NaN, ...
    'paper_efficiency',NaN,'calculated_efficiency',NaN,'total_efficiency',NaN, ...
    'energy_error',NaN);
nComparisonRows=sum(arrayfun(@(s)numel(s.target_orders),cases));
rows=repmat(rowTemplate,1,nComparisonRows);
rowIndex=0;

fprintf('\nNI2019: all printed geometries, direct modal-matching calculation\n');
fprintf('%-5s %-28s %8s %10s %10s %10s\n','ID','name','order','angle','|A|','eta');
fprintf('%s\n',repmat('-',1,80));

for c = 1:numel(cases)
    r = ni2019_modal_solver(cases(c).cfg);
    results{c} = r;
    prop = find(r.is_propagating);
    for k = 1:numel(prop)
        ii = prop(k);
        fprintf('%-5s %-28s %+8d %10.3f %10.5f %10.6f\n', ...
            cases(c).id,truncate(cases(c).name,28),r.orders(ii), ...
            r.theta_deg(ii),abs(r.A(ii)),r.eta(ii));
    end
    fprintf('      sum eta = %.12f, energy error = %.3e\n',r.total_efficiency,r.energy_error);

    for k = 1:numel(cases(c).target_orders)
        order = cases(c).target_orders(k);
        ii = find(r.orders==order,1);
        paperAngle = pick(cases(c).paper_angles_deg,k);
        paperA = pick(cases(c).paper_amplitudes,k);
        paperEta = pick(cases(c).paper_efficiencies,k);
        calcA = abs(r.A(ii));
        rowIndex=rowIndex+1;
        rows(rowIndex) = struct( ...
            'case_id',cases(c).id,'figure',cases(c).figure,'order',order, ...
            'paper_angle_deg',paperAngle,'calculated_angle_deg',r.theta_deg(ii), ...
            'angle_error_deg',r.theta_deg(ii)-paperAngle, ...
            'paper_abs_A',paperA,'calculated_abs_A',calcA, ...
            'amplitude_relative_error',relerr(calcA,paperA), ...
            'paper_efficiency',paperEta,'calculated_efficiency',r.eta(ii), ...
            'total_efficiency',r.total_efficiency,'energy_error',r.energy_error);
    end
end

comparison = struct2table(rows);
writetable(comparison,fullfile(outputDir,'Ni2019_paper_comparison.csv'));
save(fullfile(outputDir,'Ni2019_all_results.mat'),'cases','results','comparison');

% Summary plot: every propagating order for all geometries.
fig1 = figure('Color','w','Position',[100 100 1250 720]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for c = 1:numel(cases)
    nexttile; r=results{c}; id=r.is_propagating;
    bar(r.orders(id),r.eta(id),.65,'FaceColor',[.10 .45 .75]);
    ylim([0 1.05]); grid on; box on; xlabel('Order n'); ylabel('\eta_n');
    title(sprintf('%s: %s',cases(c).id,cases(c).name),'Interpreter','none');
end
exportgraphics(fig1,fullfile(outputDir,'Ni2019_all_efficiencies.png'),'Resolution',180);

% Reconstructed one-period fields corresponding to Figs. 2 and 5.
fig2 = figure('Color','w','Position',[100 100 1250 720]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for c = 1:numel(cases)
    r=results{c};
    ymax = 3*r.a; ymin = -max(r.depths);
    f = ni2019_reconstruct_field(r,linspace(0,r.a,280), ...
        linspace(ymin,ymax,480),'scattered');
    nexttile;
    imagesc(f.x/r.a,f.y/r.a,real(f.p)); axis xy tight;
    clim(max(abs(clim))*[-1 1]); colormap(turbo); box on;
    xlabel('x/a'); ylabel('y/a'); title(cases(c).id,'Interpreter','none');
end
exportgraphics(fig2,fullfile(outputDir,'Ni2019_reconstructed_fields.png'),'Resolution',180);

% Accuracy report for the values that are explicitly numerical in the paper.
hasA = ~isnan(comparison.paper_abs_A);
hasEta = ~isnan(comparison.paper_efficiency);
fprintf('\nExplicit paper-number comparison\n');
disp(comparison(hasA | hasEta,:));
fprintf('Mean amplitude relative error (where printed): %.4f %%\n', ...
    100*mean(comparison.amplitude_relative_error(hasA),'omitnan'));
fprintf('Maximum amplitude relative error (where printed): %.4f %%\n', ...
    100*max(comparison.amplitude_relative_error(hasA),[],'omitnan'));

% Automated diagnostics for the two inconsistencies exposed by direct use
% of the printed parameters.
r0 = results{5};
fprintf('\nPrinted-parameter diagnostics\n');
fprintf('F4a printed 0:1 case: eta(-1,0,+1) = %.6f, %.6f, %.6f.\n', ...
    r0.eta(r0.orders==-1),r0.eta(r0.orders==0),r0.eta(r0.orders==1));
r5 = results{6};
fprintf('F5 printed period gives order -1 angle %.6f deg (paper: -28.4 deg).\n', ...
    r5.theta_deg(r5.orders==-1));

assignin('base','Ni2019_cases',cases);
assignin('base','Ni2019_results',results);
assignin('base','Ni2019_comparison',comparison);

function value = pick(v,k)
if isempty(v), value=NaN; else, value=v(k); end
end
function e = relerr(a,b)
if isnan(b), e=NaN; else, e=abs(a-b)/max(abs(b),eps); end
end
function s = truncate(s,n)
if strlength(s)>n, s=extractBefore(string(s),n); end
end
