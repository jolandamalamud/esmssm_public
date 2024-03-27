% running algortihm for smartphone-based sampling data as followed:
% z_t = A z_t-1 + W z_t-1 + h + e_t , e_t ~ N(0,S)
% x_t = B z_t + nu_t , nu_t ~ N(0,G)

function netest = ssm_run_algorithm(option, nettrue, init)

    opt = option.opt;
    X = nettrue.x;
    input = nettrue.input;

    %% Run Optimization
    %--------------------------------------------------------------------------
    if isequal(opt,1) % do E and M in one step with exact init conditions

        % E-step
        [Ez, Ezz, Ezz_1, MLL] = ssm_inference_KS(option, nettrue, X, input);

        % M-step
        netest = ssm_estimation(option, Ez, Ezz, Ezz_1, X, input, [], []);
        netest.z = Ez;
        keyboard;
        
    %-----------------------------------------------------------------------------

    elseif isequal(opt,2) %  do E and M in one step with random init conditions
        
        for ms = 1:option.maxms
            if isempty(init)
                % initialize with random conditions
                init = ssm_initialization(option);
                init.G = diag(nanvar(X, 0, 2));
                init.mu0 = X(:,1);

                if option.S_fixed 
                    init.S = nettrue.S;
                end
            end

            netest = init;

            % EM
            i = 1; 
            MLL = []; LLR = 1e8;
            tol = 1e-4;  % tolerance (EM)
            MaxIter = option.maxEMiter;   % max iterations

            while (i<3 || LLR > tol * abs(MLL(1)) || LLR < 0) && i < MaxIter    

                % iterate over E- and M-step
                [Ez, Ezz, Ezz_1, MLL(i)] = ssm_inference_KS(option, netest, X, input);

                if isnan(MLL(i))
                    netest = [];
                    return; 
                end

                netest = ssm_estimation(option, Ez, Ezz, Ezz_1, X, input, [], []);

                if option.S_fixed 
                    netest.S = nettrue.S;
                end


                if i > 2; LLR = MLL(i) - MLL(i-1); else LLR = 1e8; end

                fprintf('it=%i LL=%g ms=%g \r', i, MLL(i), ms)
                i = i + 1;

            end

            ll_new = MLL(end);

            if ms == 1 || ll_new > ll_old
                netest.z = Ez;
                netest.LL = MLL;
                netest.initialization = init;
                ll_old = ll_new;
            end
    end

end
