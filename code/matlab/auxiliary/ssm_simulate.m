function nettrue = ssm_simulate(nettrue, option)

    nl = option.dimZ; % number hidden states
    ns = option.dimX; % number observations
    ni = option.dimInp; % number observations
    T = option.T;     % time steps
    
    B=nettrue.B; G=nettrue.G; A=nettrue.A; W=nettrue.W; h=nettrue.h;
    C=nettrue.C; input=nettrue.input; S=nettrue.S; mu0=nettrue.mu0;
    
%% Run this network (create simulated states) and get observations
%--------------------------------------------------------------------------

    Z = zeros(nl, T); %latent states
    X = zeros(ns, T); %observations

    Z(:, 1) = mvnrnd(mu0 + C * input(:, 1), S)';
    X(:, 1) = mvnrnd(B * Z(:,1), G);
    
%     for t = 1 : T - 1
%         Z(:, t+1) = mvnrnd((A + W) * Z(:, t) + h + C * input(:, t+1), S)';
%         X(:, t+1) = mvnrnd( B * Z(:, t+1), G)';
%     end
    
    for t = 1 : T - 1
        errz = mvnrnd(zeros(nl, 1), S)';
        errx = mvnrnd(zeros(ns, 1), G)';
        Z(:, t+1) = (A + W) * Z(:, t) + h + C * input(:, t+1) + errz;
        X(:, t+1) = B * Z(:, t+1) + errx;
    end
    
    nettrue.x = X; nettrue.z = Z; nettrue.T = T;

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
    
%% compute controllablility
%--------------------------------------------------------------------------
    
%     nettrue.controllability = ctrb((nettrue.A + nettrue.W), nettrue.C);

end