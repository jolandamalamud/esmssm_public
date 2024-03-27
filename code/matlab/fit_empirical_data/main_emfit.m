option.filepath = '~/phd/projects/esmssm/data/data_prep/';
option.goalpath = '~/phd/projects/esmssm/results/modelfits/';
option.dataset = 'twindata';
option = set_universal_options(option);
fit = load_data(option);
%% 
Nsj = length(fit);
parfor k = 1:Nsj
    
    rng(1);
    suboption = set_data_options(option, fit{k});
    fprintf('%s: subject or subsystem %i \n', suboption.dataset, k);
    fprintf('**********************************\n');
    
    nettrue = struct();
    nettrue.x = fit{k}.data;
    nettrue.input = fit{k}.input;
    suboption.empirical = false;
    suboption.r = 2;
    netest = ssm_run_algorithm(suboption, nettrue, []);
    fit{k}.est = netest;
    
end

save([option.goalpath, option.dataset, '/ssmfit_r2.mat'], 'fit');