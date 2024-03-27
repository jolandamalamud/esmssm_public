function nettrue = ssm_simulate_var(nettrue, option)

%% create ground truth data
%--------------------------------------------------------------------------
    
    m = option.dim; % number hidden states
    T = option.T;

    A=nettrue.A; S=nettrue.S; mu0=nettrue.mu0;
    
%% Run this network (create simulated states) and get observations
%--------------------------------------------------------------------------
    X = zeros(m, T); %observations
    X(:, 1) = mvnrnd(mu0, S);

    for t = 2 : T
        err = mvnrnd(zeros(m, 1), S)';
        X(:, t) = A * X(:, t-1) + err;
    end
    
    nettrue.x = X;
end