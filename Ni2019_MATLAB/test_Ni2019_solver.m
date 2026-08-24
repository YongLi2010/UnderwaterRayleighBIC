%% Regression tests for all published Ni2019 geometries
clear; close all; clc;
cases=ni2019_case_library;
Ni2019_results=cell(1,numel(cases));
for c=1:numel(cases), Ni2019_results{c}=ni2019_modal_solver(cases(c).cfg); end

r81 = Ni2019_results{1};
r88 = Ni2019_results{2};
r53 = Ni2019_results{3};
r11 = Ni2019_results{4};

assert(r81.eta(r81.orders==-1) > 0.999, 'The -81 deg efficiency is too low.');
assert(r88.eta(r88.orders==-1) > 0.998, 'The -88 deg efficiency is too low.');
assert(abs(r11.eta(r11.orders==-1)-0.5) < 5e-4, 'The 1:1 splitter is not reproduced.');
assert(abs(r11.eta(r11.orders==+1)-0.5) < 5e-4, 'The 1:1 splitter is not reproduced.');

ampRatio53 = abs(r53.A(r53.orders==-1))/abs(r53.A(r53.orders==+1));
assert(abs(ampRatio53-5/3) < 0.02, 'Published 5:3 amplitude ratio is not reproduced.');
for c = 1:4
    assert(Ni2019_results{c}.energy_error < 1e-10, 'Energy conservation failed.');
end

r5=Ni2019_results{6};
assert(abs(abs(r5.A(r5.orders==-1))-.7562)/.7562 < .003, ...
    'Four-groove -1 amplitude does not match the paper.');
assert(abs(abs(r5.A(r5.orders==2))-1.7883)/1.7883 < .003, ...
    'Four-groove +2 amplitude does not match the paper.');
for c=1:numel(cases)
    assert(Ni2019_results{c}.energy_error<1e-10,'Energy conservation failed.');
end

fprintf('All Ni2019 regression tests passed.\n');
