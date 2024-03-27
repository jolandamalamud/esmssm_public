function netest = ssm_run_algorithm_annealing(option, nettrue, s)

    X = nettrue.x;
    input = nettrue.input;
   
    init = ssm_initialization(option);

    if option.S_fixed 
        init.S = nettrue.S;
    end
    
    if option.zero_centered
        init.mu0 = zeros(option.dimZ, 1);   
    else
        init.mu0 = X(:, 1);
    end
    
    netest = init;



    %% annealing EM
    % 1) Sigma = diag(10^-i)
%     for i = 0:2
%         for t = 1:1
%             netest.S = 10^-i * eye(option.dimZ);
%             [Ez, Ezz, Ezz_1] = ssm_inference(netest, X, input);
%             netest = ssm_estimation(option, Ez, Ezz, Ezz_1, X, input);
%         end
%     end
    
    i = 1; 
    MLL = []; LLR = 1e8;
    tol = 1e-4;  % tolerance (EM)
    MaxIter = 1000;   % max iterations
    fixed_S = eye(option.dimZ);
    option.S_fixed = true;
    option.G_fixed = false;
    
    while (i<2 || LLR > tol * abs(MLL(1))) && i < MaxIter    

        % iterate over E- and M-step
        [Ez, Ezz, Ezz_1, MLL(i)] = ssm_inference(netest, X, input);
        netest = ssm_estimation(option, Ez, Ezz, Ezz_1, X, input);
        
        if option.S_fixed
            netest.S = fixed_S;
        end
        
        if i > 1; LLR = MLL(i) - MLL(i-1); else LLR = 1e8; end

        fprintf('it=%i LL=%g \r', i, MLL(i))
        i = i + 1;
    end
    
%     if any(diag(netest.G) > 10)
%         keyboard;
%     end
    
    fixed_G = netest.G;
    option.G_fixed = true;
    option.S_fixed = false;
    
    i = 1; 
    MLL = []; LLR = 1e8;
    tol = 1e-4;  % tolerance (EM)
    MaxIter = 1000;   % max iterations

    while (i<2 || LLR > tol * abs(MLL(1))) && i < MaxIter    

        % iterate over E- and M-step
        [Ez, Ezz, Ezz_1, MLL(i)] = ssm_inference(netest, X, input);
        netest = ssm_estimation(option, Ez, Ezz, Ezz_1, X, input);
        
        if option.G_fixed
            netest.G = fixed_G;
        end

        if i > 1; LLR = MLL(i) - MLL(i-1); else LLR = 1e8; end

        fprintf('it=%i LL=%g \r', i, MLL(i))
        i = i + 1;
    end
    
    netest.z = Ez;
    
end