%% Re-optimize the Fig. 4(a) 0:1 target when using the printed parameters fails
clear; clc;
cases=ni2019_case_library('lambda',1,'N',61,'K',10);
printed=ni2019_modal_solver(cases(5).cfg);
fprintf('Printed geometry eta(-1,0,+1): %.6f %.6f %.6f\n', ...
    printed.eta(printed.orders==-1),printed.eta(printed.orders==0), ...
    printed.eta(printed.orders==1));

best=ni2019_optimize_geometry(cases(5).cfg,1,1,'NumStarts',8,'Display','off');
r=best.result;
fprintf('Re-optimized geometry (not claimed to be a correction of the paper):\n');
fprintf('  widths/a      = %s\n',mat2str(best.cfg.widths/best.cfg.a,7));
fprintf('  depths/lambda = %s\n',mat2str(best.cfg.depths/best.cfg.lambda,7));
fprintf('  gaps/a        = %s\n',mat2str(best.cfg.gaps/best.cfg.a,7));
fprintf('  eta(-1,0,+1)  = %.8f %.8f %.8f\n', ...
    r.eta(r.orders==-1),r.eta(r.orders==0),r.eta(r.orders==1));

validationCfg=best.cfg; validationCfg.N=101; validationCfg.K=10;
validation=ni2019_modal_solver(validationCfg);
fprintf('  N=101 validation = %.8f %.8f %.8f\n', ...
    validation.eta(validation.orders==-1),validation.eta(validation.orders==0), ...
    validation.eta(validation.orders==1));
best.validation_N101_K10=validation;
outputDir=fullfile(pwd,'results');
if ~exist(outputDir,'dir'), mkdir(outputDir); end
save(fullfile(outputDir,'Ni2019_F4a_reoptimized.mat'),'best');
assignin('base','Ni2019_F4a_reoptimized',best);
