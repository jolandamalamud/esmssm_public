function [nettrue, option] = ssm_using_empirical_estimates(fit, option)

% using as parameter setting the estimated parameters from empirical data

    nettrue = fit.parameter;
    option.dimZ = size(fit.state, 1);
    option.dimX = size(fit.data, 1);
    option.dimInp = size(fit.input, 1);
    
    if option.sparse
        option.T = size(fit.data, 2);
        nettrue.input = fit.input;
       	option.timing = find(~isnan(fit.data(1,:)));
        option.gap = mean(option.timing(2:end) - option.timing(1:end-1));
    else
        option.T = size(fit.data(:, ~isnan(fit.data(1,:))),2);
        option.gap = 1;
        if option.input 
            nettrue.input = fit.input(:, ~isnan(fit.data(1,:))); 
        else
            nettrue.input = zeros(2, option.T);
        end
    end
    
%     % extend time
%     option.T = 200;
%     option.input = zeros(2, option.T);
%     % reduce noise
%     nettrue.G = 0.01 * eye(option.dimZ);
    
end

