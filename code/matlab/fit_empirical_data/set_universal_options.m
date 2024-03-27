function option = set_universal_options(option)

    % which data set:
    % 'twindataset'
    % 'esm_merge2.0','esm_merge4.0','esm_merge7.0','esm_merge10.0'
    % 'leemput_patient'
    % 'leemput_control'
    % 'wichers'
    % 'delespaul'
%     option.dataset = ;
%     option.filepath = [];
    option.empirical = true;
 
    % model characteristics
    option.dense = false; % dense or sparse data
    option.zero_centered = false; % center data around 0?
    option.split = false;
    
    option.bias = true; % modelling with bias?
    option.input = true; % modelling with input?
    option.pca = false; % first do pca on data?
    option.smoothing = false; % first smooth data?
    option.S_fixed = false;
    
    option.opt = 2; % em 
    option.maxEMiter = 200;
    option.maxms = 1; % multistarts
    option.save_matfile = false; % save fit
    option.figure = false;
    
    option.k = 60;
    
end