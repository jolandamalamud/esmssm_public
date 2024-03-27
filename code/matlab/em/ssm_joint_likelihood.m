function LL = ssm_joint_likelihood(X, inp, Z, par)

    [nl, T] = size(Z);
    tau = size(X,2);
    inp = inp(:,~isnan(inp(1,:)));

    P1 = Z(:,1) + par.C * inp(:,1); 
    LL1 = P1' * inv(par.S) * P1;
    
    P2 = Z(:,2:T) - ((par.A + par.W) * Z(:,1:T-1) + par.h + par.C * inp(:,2:T));
    LL2 = sum(sum(P2 .* (inv(par.S) * P2)));
    
    P3 = X(:,~isnan(X(1,:))) - Z(:,~isnan(X(1,:))); 
    LL3 = sum(sum(P3 .* (inv(par.G) * P3)));
    LL = -1/2 * (T * log(2 * pi * det(par.S)) + ...
        tau * log(2 * pi * det(par.G)) + LL1 + LL2 + LL3);
end