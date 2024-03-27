function ssm_emfit_fast_annealing(sj, option)
    rng(1)
    datafit = importdata(sprintf('%s/data_prep_cluster/%s/%s', ...
                option.filepath, option.foldername, option.moodfiles(sj).name));
    subject_fit = datafit;
    mood = datafit.data;
    inp = datafit.input;
    option.timing = find(~isnan(mood(1,:)));
    
    option.nonNaNidx = ~isnan(mood(1, :)); % non NaN entries
    option.dimZ = size(mood, 1); % number hidden states
    option.dimX = size(mood, 1); % number observations
    option.dimInp = size(inp, 1); % number observations
    option.timing = find(~isnan(mood(1,:))); % observation timings
    option.gap = mean(option.timing(2:end) - option.timing(1:end-1));
    option.T = size(mood,2);
    option.tau = size(mood(:,option.nonNaNidx), 2);
    
    for ms = 1:option.maxms
        % initialize parameters randomly
        init = ssm_initialization(option);
        if option.zero_centered
            init.mu0 = zeros(option.dimZ, 1);
        else
            init.mu0 = mood(:, 1);
        end
        init.G = diag(nanvar(mood, 0, 2));
        
        % ANNEALING EM
        tol = 1e-3; % tolerance (EM)
        MaxIter = 100;  % max iterations
        netest = init;
        
        for j = 0:3
            % 1) slowely reducing dynamics noise Sigma = diag(10^-i) and
            % Gamma estimated in first interation only, then fixed
            netest.S = 10^-j * eye(option.dimZ);
            
            i = 1; MLL = []; LLR = 1e8;
            while (i==1 || LLR > tol * abs(MLL(1)) || LLR < 0) && i < MaxIter  

                % iterate over E- and M-step
                [Ez, Ezz, Ezz_1, MLL(i)] = ssm_inference_KS(option, netest, mood, inp);
                if j == 0
                    netest = ssm_estimation(option, Ez, Ezz, Ezz_1, mood, inp, netest.S, []);
                else
                    netest = ssm_estimation(option, Ez, Ezz, Ezz_1, mood, inp, netest.S, netest.G);
                end

                if i > 1; LLR = MLL(i) - MLL(i-1); else LLR = 1e8; end

                fprintf('it=%i LL=%g \r', i, MLL(i))
                i = i + 1;
            end
        end

        ll_new = MLL(end);

        if ms == 1 || ll_new > ll_old
           state = Ez;
           parameter = netest;
           loglike = MLL;
           ll_old = ll_new;
        end
    end
    
    subject_fit.state = state;
    subject_fit.parameter = parameter;
    subject_fit.LL = loglike;
    subject_fit.option = option;

    save([option.goalpath '/fit_annealing_s' num2str(sj)], 'subject_fit');
end
