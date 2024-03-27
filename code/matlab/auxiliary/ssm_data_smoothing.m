function smoothed_data = ssm_data_smoothing(data, scaling)
    %% smooth data to other time scale to reduce sparsity
    T = size(data, 2);
    
    for t = 1:floor(T / scaling)
        smoothed_data(:, t) = ...
            nanmean(data(:, (t-1) * scaling + 1 : t * scaling), 2); 
    end
    
    smoothed_data(:, t) = nanmean(data(:, (t-1) * scaling + 1 : T), 2); 
    
end