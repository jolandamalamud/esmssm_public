function option = ssm_sim_set_options()

    option.opt = 2; 
    % 1=1 EM step with same initial conditions (checks algorithm), 
    % 2=run full EM with random initial conditions

    option.dimZ = 4; % number of hidden states
    option.dimX = 4; % number of observations
    option.dimInp = 16; % number of inputs
    option.T = 60; % number of timepoints
    option.zero_centered = false;
    option.bias = false; % bias or mean adjustment
    option.input = false; % input
    option.runs = 1000; % number of different simulations
    
    if option.dimZ ~= option.dimX % latent space
        option.B_full = false;
        option.B_fixed = false;
    else
        option.B_fixed = true; % direct mapping
    end
    option.S_fixed = false; % keep Sigma fixed
    option.sparse = true; % sparse data
    
    if option.sparse
        option.data_timepoints = true; % take timing from empirical data
        if option.data_timepoints == true
            option.days = 6;
        else
            option.gap = 10; % gap between timepoints
            option.T = option.T * 100; % * option.gap; % new number of timepoints
        end
    end
    
    option.maxms = 1;
    option.maxEMiter = 200; % maximum of iterations
    option.multistart = false; % multistarts of EM
    option.figure = false; % save figures
    option.trackEM = false; % save estimated parameter & state after every iteration
    option.save_matfile = false;
end