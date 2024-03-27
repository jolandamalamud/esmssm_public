%% ssm_master
% EM on simulated time series from a linear state space model:
% latent model: z_t = (A+W) * z_(t-1) + h + C * s_t + err1, err1 = N(0, S)
% observation model: x_t = z_t + err2, err2 = N(0, G)

% clear all; close all; addpath(genpath(pwd));

% set options
option = ssm_sim_set_options();

stats = struct;

s = 1; r = 1;
while s <= option.runs
    % for different settings same randomizations
    rng(r);
    
    % create ground truth for the ssm
    option.empirical = false;
    [nettrue, option] = ssm_ground_truth(option, s);
    
%     % use empirical parameter settings
%     [nettrue, option] = ssm_using_empirical_estimates(fit{s}, option);

    % simulate time series
    if option.sparse
        nettrue = ssm_simulate_sparse(nettrue, option, r);
    else
        nettrue= ssm_simluate(nettrue, option);
    end
    
    r = r + 1;
    
%     % assure stationarity of data
%     for ts = 1:size(nettrue.z, 1)
%         sta(ts) = adftest(nettrue.z(ts, :));
%     end
%     
%     if sum(sta) < size(nettrue.z, 1)
%         continue;
%     end
    option.nonNaNidx = ~isnan(nettrue.x(1,:));
    option.empirical = true;
    netest = ssm_run_algorithm(option, nettrue, []);
    
    if isempty(netest)
        continue
    else
        stats(s).netest = netest;
        stats(s).nettrue = nettrue;
        stats(s).option = option;
        s = s + 1;
    end
    
    if option.dimZ ~= option.dimX
        z_swapped = ssm_swap_latents(netest.z, nettrue.z);
        netest.z = z_swapped;
    end
    
end

% compute EM performance
option.state_stats = true;
statistics = ssm_statistic(option, stats);

if option.save_matfile
    if option.sparse; net = 'sparse'; else net = 'dense'; t = option.T; end
    if option.bias; b = '_bias'; else b = ''; end
    if option.input; i = '_input'; else i = ''; end
    save(['/matfiles/ssm_simulation_' net b i '_T' num2str(option.T)], ...
        'stats');
end

