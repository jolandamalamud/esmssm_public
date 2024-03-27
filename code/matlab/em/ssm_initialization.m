function [nettrue] = ssm_initialization(option)
    % parameter initialization
    nl = option.dimZ; % number hidden states
    ns = option.dimX; % number observations
    ni = option.dimInp; % number observations
    
    if option.empirical == false
        
        % initialize dynamics matrix 
        % (for sparse data with eigenvalues close to 1)
        W = randn(nl) ./ nl; W = W - diag(diag(W));
        A = diag(rand(nl, 1));
        AW = real((A + W)^(1/option.gap));

        % check for stability, absolut eigenvalues < 1
        while find(abs(eig(AW)) >= 1)
            W = randn(nl)./ nl; W = W - diag(diag(W));
            A = diag(rand(nl, 1));
            AW = real((A + W)^(1/option.gap));
        end

        A = diag(diag(AW)); % autoregressive weights
        W = AW - A; % interaction weights

        if option.bias % bias
            h = randn(nl, 1) / option.gap; 
        else
            h = zeros(nl, 1);
        end

        C = randn(nl, ni)*0.1; % input weights
    
    elseif option.empirical == true
        
        A = eye(nl);
        W = zeros(nl);
        h = zeros(nl,1);
        C = zeros(nl, ni); 
        
    end
    
    S = eye(nl)./option.gap; % process noise variance
    B = eye(nl); % mapping observation to latents (usually identity)
    G = eye(ns); % observation noise variance

    nettrue.B=B; nettrue.G=G; nettrue.A=A; nettrue.W=W; nettrue.h=h;
    nettrue.C=C; nettrue.S=S;
    
end