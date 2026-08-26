%% Sub-microscale kappa refinement of the free-plate Rayleigh crossing
clear; clc;
rootDir = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(rootDir,'COMSOL_MATLAB','rounded_180k_results');
model = mphload(fullfile(outDir, ...
    'RayleighBIC_rounded_180k_coupled_free.mph'),'CoupledCrossingRefine');
model.study('std1').feature('eig').set('neigs',16);
model.study('std1').feature('eig').set('shift','f0');
a = model.param.evaluate('a');
c0 = model.param.evaluate('c0');

kappa = (0.1208550:0.0000005:0.1208620).';
n = numel(kappa);
frequency = complex(nan(n,1));
for i=1:n
    model.param.set('kxB',sprintf('2*pi*(%.16g)/a',kappa(i)));
    model.study('std1').run;
    f = mphglobal(model,'freq');
    candidate = find(real(f)>175e3 & real(f)<180e3);
    [~,local] = min(abs(real(f(candidate))-177.724e3));
    frequency(i)=f(candidate(local));
    fprintf('%2d/%2d kappa %.10f f %.9f %+.9fi kHz\n',i,n, ...
        kappa(i),real(frequency(i))/1e3,imag(frequency(i))/1e3);
end
Omega = real(frequency)*a/c0;
deltaOmega = Omega-(1-kappa);
Q = abs(real(frequency)./(2*imag(frequency)));
T = table(kappa,real(frequency),imag(frequency),Omega,deltaOmega,Q, ...
    'VariableNames',{'kappa','frequency_real_hz','frequency_imag_hz', ...
    'Omega','delta_Omega','Q'});
writetable(T,fullfile(outDir,'coupled_free_rayleigh_crossing_refined.csv'));

kRA = zeroCross(kappa,deltaOmega);
kIm = zeroCross(kappa,imag(frequency));
summary = table(kRA,kIm,abs(kRA-kIm), ...
    'VariableNames',{'kappa_Rayleigh','kappa_zero_imag','kappa_separation'});
writetable(summary,fullfile(outDir, ...
    'coupled_free_rayleigh_crossing_summary.csv'));
fprintf('kappa_Rayleigh = %.12f\n',kRA);
fprintf('kappa_zero_imag = %.12f\n',kIm);
fprintf('separation = %.3e\n',abs(kRA-kIm));


function root = zeroCross(x,y)
idx = find(y(1:end-1).*y(2:end)<=0,1,'first');
if isempty(idx)
    root=nan;
else
    root=x(idx)-y(idx)*(x(idx+1)-x(idx))/(y(idx+1)-y(idx));
end
end
