function ssm_data_singlesubject_cluster(sj, option)
    
    fit = struct;

    subject_mood = importdata(sprintf('%s/data_prep/%s/%s', ...
                option.filepath, option.foldername, option.moodfiles(sj).name));

    subject_timing = find(~isnan(subject_mood(1,:)));
    
    mood = subject_mood(:, 1:subject_timing(end));
    
    if option.zero_centered
        mood = mood - nanmean(mood);
    end
    
    option.gap = mean(subject_timing(2:end) - subject_timing(1:end-1));

    if option.input
        % load inputdata per subject
        inputdata = importdata(sprintf('%s/data_prep/%s/%s', ...
                option.filepath, option.foldername, option.inputfiles(sj).name));
        subject_input = inputdata(:, ~isnan(inputdata(1,:)));
        input = zeros(size(subject_input,1), subject_timing(end));
        input(:, subject_timing) = subject_input(:, 1:length(subject_timing));
    else
        input = zeros(2, subject_timing(end));
    end
    
    option.dimZ = size(mood, 1); % number hidden states
    option.dimX = size(mood, 1); % number observations
    option.dimInp = size(input, 1); % number observations
    
    % initialize parameters randomly
    init = ssm_initialization(option);
    
    if option.zero_centered
        init.mu0 = zeros(option.dimZ, 1);
    else
        init.mu0 = nanmean(mood, 2);
    end
    
    netest = init;

    %% annealing EM
    
    i = 1; 
    MLL = []; LLR = 1e8;
    tol = 1e-4;  % tolerance (EM)
    MaxIter = 1000;   % max iterations
    fixed_S = 10^-3 * eye(option.dimZ);
    option.S_fixed = true;
    option.G_fixed = false;
    
    while (i<3 || LLR > tol * abs(MLL(2)) || LLR < 0) && i < MaxIter    

        % iterate over E- and M-step
        [Ez, Ezz, Ezz_1, MLL(i)] = ssm_inference_cluster(netest, mood, input);
        netest = ssm_estimation_cluster(option, Ez, Ezz, Ezz_1, mood, input);
        
        if option.S_fixed
            netest.S = fixed_S;
        end
        
        if i > 2; LLR = MLL(i) - MLL(i-1); else LLR = 1e8; end

        fprintf('it=%i LL=%g \r', i, MLL(i))
        i = i + 1;
    end
    
    fixed_G = netest.G;
    option.G_fixed = true;
    option.S_fixed = false;
    
    i = 1; 
    MLL = []; LLR = 1e8;
    tol = 1e-4;  % tolerance (EM)
    MaxIter = 1000;   % max iterations

    while (i<3 || LLR > tol * abs(MLL(2)) || LLR < 0) && i < MaxIter    

        % iterate over E- and M-step
        [Ez, Ezz, Ezz_1, MLL(i)] = ssm_inference_cluster(netest, mood, input);
        netest = ssm_estimation_cluster(option, Ez, Ezz, Ezz_1, mood, input);
        
        if option.G_fixed
            netest.G = fixed_G;
        end

        if i > 2; LLR = MLL(i) - MLL(i-1); else LLR = 1e8; end

        fprintf('it=%i LL=%g \r', i, MLL(i))
        i = i + 1;
    end
    
    fit.parameter = netest;
    fit.initialization = init;
    fit.state = Ez;
    fit.data = mood;
    fit.input = input;
    fit.LL = MLL;
    
    save(sprintf('%s/data_fit/%s/ssm_fit_%s%s', option.filepath, ...
    option.foldername, option.foldername, sj), 'fit');

end