function nettrue = ssm_simulate_sparse(nettrue, option)

%% create ground truth data
%--------------------------------------------------------------------------
    
    nl = option.dimZ; % number hidden states
    ns = option.dimX; % number observations
    ni = option.dimInp; % number observations
    T = option.T;

    B=nettrue.B; G=nettrue.G; A=nettrue.A; W=nettrue.W; h=nettrue.h;
    C=nettrue.C; input=nettrue.input; S=nettrue.S; mu0=nettrue.mu0;
    
%% Run this network (create simulated states) and get observations
%--------------------------------------------------------------------------
    Z = zeros(nl, T); %latent states
    X = zeros(ns, T); %observations
    
    Z(:, 1) = mvnrnd(mu0 + C * input(:, 1), S)';
    X(:, 1) = mvnrnd(B * Z(:,1), G);

    for t = 2 : T
        errz = mvnrnd(zeros(nl, 1), S)';
        if ismember(t, option.timing)
            Z(:, t) = (A + W) * Z(:, t-1) + C * input(:, t) + h + errz;
            errx = mvnrnd(zeros(ns, 1), G)';
            X(:, t) = B * Z(:, t) + errx;
        else
            Z(:, t) = (A + W) * Z(:, t-1) + h + errz;
            X(:, t) = nan;
        end
    end
    
    nettrue.x = X; nettrue.z = Z;
    
%% true fixed point
nettrue.FP = ssm_fixed_point(nettrue);

%% plot observations and latent states
%--------------------------------------------------------------------------
    if option.figure
        
        figure(); hold on;
        subplot(1,2,1); plot(nettrue.x', '.'); title('observations');
        subplot(1,2,2); plot(nettrue.z', '.'); title('latent state'); 
        hold off;
        
    end
    
    
%% compute true controllablility
%--------------------------------------------------------------------------
%     nettrue.controllability = ctrb((nettrue.A + nettrue.W), nettrue.C);

end