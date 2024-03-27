function [nettrue, option] = ssm_ground_truth(option)

%% parameter initialization
%--------------------------------------------------------------------------

    nl = option.dimZ; % number hidden states
    ns = option.dimX; % number observations
    ni = option.dimInp; % number inputs
    
    if option.sparse % s network
        
            if option.data_timepoints
                schedule = esm_schedule(option.days, 9*60, 22*60, 10);
                data_points = schedule(:,2) - schedule(1,2) + 1;
                T = data_points(end); 
                option.tau = length(data_points);
                % gap is the mean difference between consecutive timings
                gap = mean(data_points(2:end) - data_points(1:end-1));
                if option.consider_gap == true
                    gap = 1;
                end
                option.timing = data_points;
            else
                T = option.T; % time steps
                gap = option.gap; % gap between time steps
                data_points = 1:gap:T;
                option.tau = length(data_points);
                option.timing = data_points;
            end
            
            if option.input
                for i = 1:ni
                    freq = 0.8;%abs(0.25 + 0.25*randn(1));
                    if freq > 1; freq = 1; end
                    inp(i,:) = real(rand(1,option.tau) > (1-freq));
                end
                
                % ensure full rank matrix
                while (rank(inp) ~= ni) 
                    for i = 1:ni
                        freq = 0.8;%abs(0.25 + 0.25*randn(1));
                        if freq > 1; freq = 1; end
                        inp(i,:) = real(rand(1,option.tau) > (1-freq));
                    end
                end
                input = nan(ni, T);
                input(:, data_points) = inp;
            end
           
    else % dense network
            
            T = option.T; % time steps
            gap = 1; % no gap
            option.tau = T;
            option.timing = [1:T];
            if option.input
                for i = 1:ni
                    freq = abs(0.25 + 0.25*randn(1));
                    if freq > 1; freq = 1; end
                    input(i,:) = real(rand(T) > (1-freq));
                end
                % ensure full rank matrix
                while (rank(input) ~= ni) 
                    for i = 1:ni
                        freq = abs(0.25 + 0.25*randn(1));
                        if freq > 1; freq = 1; end
                        input(i,:) = real(rand(T) > (1-freq));
                    end
                end
            end
    end
    
    option.T = T; option.gap = gap;
    
    if ~option.input
        input = zeros(ni, T);
    end
    
    nettrue = ssm_initialization(option);
  
    nettrue.input=input;
    nettrue.mu0=randi([1,7],ns,1);
    
    end