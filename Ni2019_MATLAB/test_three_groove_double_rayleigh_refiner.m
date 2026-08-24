function T=test_three_groove_double_rayleigh_refiner()
%TEST_THREE_GROOVE_DOUBLE_RAYLEIGH_REFINER Regression test for exact root.

x0=[.25 .65 .25 .10 .52 .10 .07 .07];
cfg=struct('a',1,'lambda',1,'theta_i_deg',0, ...
    'depths',x0(1:3),'widths',x0(4:6),'gaps',x0(7:8), ...
    'N',41,'K',3,'solve_scattering',false);
lb=[.08 .40 .08 .045 .38 .045 .025 .025];
ub=[.35 .90 .35 .18 .65 .18 .15 .15];
r=ni2019_refine_three_groove_double_rayleigh_bic(cfg,x0, ...
    'Truncation',[41 3],'LowerBounds',lb,'UpperBounds',ub, ...
    'FillMax',.94,'MaxFunctionEvaluations',1600, ...
    'MaxIterations',180,'Display','off');
d=ni2019_three_groove_parity_diagnostic(r.strict_operator);

assert(isequal(r.removed_orders(:),[-1;0;1]));
assert(r.sigma_ratio<1e-10);
assert(r.strict_operator.strict_residual<1e-9);
assert(all(r.strict_operator.mode.A( ...
    ismember(r.strict_operator.full_operator.orders,[-1 0 1]))==0));
assert(r.strict_operator.groove.pressure_proxy_fraction>.5);
assert(r.fill_fraction<.94 && all(r.x>lb+1e-6) && all(r.x<ub-1e-6));
assert(d.is_symmetric && strcmp(d.parity_label,'odd'));
assert(d.coefficient_residuals.odd<1e-6);

fprintf('three-groove double-Rayleigh refiner test passed.\n');
fprintf('  sigma ratio=%.3e, raw residual=%.3e, groove fraction=%.6f\n', ...
    r.sigma_ratio,r.strict_operator.strict_residual, ...
    r.strict_operator.groove.pressure_proxy_fraction);
fprintf('  parity=%s, odd residual=%.3e, removed orders:', ...
    d.parity_label,d.coefficient_residuals.odd);
fprintf(' %d',r.removed_orders); fprintf('\n  x=');
fprintf(' %.12g',r.x); fprintf('\n');

T=struct('result',r,'parity',d);
end
