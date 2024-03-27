function [LL, DL] = ssm_joint_likelihood_Classo_vec(C, X, inp, Z, A, h, mu0, S, G, r)
    
    [nl, T] = size(Z);
    tau = size(X(:,~isnan(X(1,:))),2);
    NoNaNidx = ~isnan(X(1,:));
    
    for t = 1:tau
        ii = inp(:,NoNaNidx); 
        inp_trans(:,:,t) = ssm_vector_transformation(ii(:,t),nl);
        inptimesC(:,t) = inp_trans(:,:,t) * C; 
    end

    inp(:,~NoNaNidx) = 0;
    inpC = zeros(nl, T);
    inpC(:,NoNaNidx) = inptimesC;
    ni = size(inp, 1);
    
    P1 = Z(:,1) + inpC(:,1); 
    LL1 = P1' * inv(diag(S)) * P1;
    
    P2 = Z(:,2:end) - (A * Z(:,1:end-1) + h + inpC(:,2:end));
    LL2 = sum(sum(P2 .* (inv(diag(S)) * P2)));
    
    P3 = X(:,~isnan(X(1,:))) - Z(:,~isnan(X(1,:))); 
    LL3 = sum(sum(P3 .* (inv(G) * P3)));
    LL = -1/2 * (T * log(2 * pi * det(diag(S))) + ...
        tau * log(2 * pi * det(G)) + LL1 + LL2 + LL3) + ...
        r * sum(abs(C));
    
    ss = diag(inv(diag(S))); ss = repmat(ss, 1, ni)'; ss = ss(:);
    
    DL = 2 * ss .* C * inp(:,1) * inp(:,1)' ...
        - 2 * inv(diag(S)) * (Z(:,1) - mu0) * inp(:,1)' ...
        + 2 * inv(diag(S)) * C * inp(:,2:end) * inp(:,2:end)' ...
        - 2 * inv(diag(S)) * (Z(:,2:end) - A * Z(:,1:end-1) - h) * inp(:,2:end)' ...
        + r * C ./ abs(C);
    %         + r * det(eye(nl) - C * ones(ni,nl)) * inv(eye(nl) - C * ones(ni,nl));
end