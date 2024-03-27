close all, clear all;
maxruns = 10;
parfor s = 1:maxruns
    option.dim=4; 
    option.days=5; 
    schedule = esm_schedule(option.days, 9*60, 22*60, 10);
    schedule = schedule(:,2) - schedule(1,2) + 1;
    option.T = schedule(end);
    A_diag = diag(rand(option.dim, 1));
    A_off = randn(option.dim);
    A = A_diag + A_off - diag(A_off);

    % check for stability, absolut eigenvalues < 1
    while find(abs(eig(A)) >= 1)
        A_diag = diag(rand(option.dim, 1));
        A_off = randn(option.dim);
        A = A_diag + A_off - diag(A_off);
    end
    
    nettrue = struct;
    nettrue.A = A;
    nettrue.S = eye(option.dim); 
    nettrue.mu0 = randi([1,7], option.dim, 1); 
    nettrue = ssm_simulate_var(nettrue, option);

    Mdl = varm(option.dim,1);
    sparse_data = nettrue.x(:, schedule);
    EstMdl = estimate(Mdl, sparse_data');

    % plot(sparse_data', 'o')

    %%
    sim_data = struct;
    sim_data.x = nan(option.dim, option.T);
    sim_data.x(:, schedule) = sparse_data;
    sim_data.data = sim_data.x;
    sim_data.input = zeros(2, option.T);

    option = set_universal_options(option);
    option = set_data_options(option, sim_data);
    option.input=false;
    option.bias = true;

    netest = ssm_run_algorithm(option, sim_data, []);

    fit{s}.true = nettrue;
    fit{s}.true.sparse = sim_data.x;
    fit{s}.est = netest;
    fit{s}.ar = EstMdl;
end

%%
for s =1:maxruns
    vvkf = get_eigvec(fit{s}.est.A+fit{s}.est.W);
    vvt = get_eigvec(fit{s}.true.A);
    vvvar = get_eigvec(fit{s}.ar.AR{1});
    vvs(:,:,s) = [vvt(:,1),vvvar(:,1),vvkf(:,1)]; 
end

%%
for s =1:maxruns
    timing = find(~isnan(fit{s}.true.sparse));
    gap = mean(timing(2:end)-timing(1:end-1));
    Akf = (fit{s}.est.A+fit{s}.est.W)^(gap);
    At = fit{s}.true.A;
    Avar = fit{s}.ar.AR{1};
    AA(:,:,s) = [reshape(At, 16, 1),reshape(Avar, 16, 1),reshape(Akf, 16, 1)]; 
end