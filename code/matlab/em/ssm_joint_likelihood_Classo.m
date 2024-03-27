function [LL, DL] = ssm_joint_likelihood_Classo(C, option, X, inp, Z, A, h, mu0, S, G)
    
    iS = inv(S);
    iG = inv(G);
    pi_power = (2 * pi)^option.dimZ;
    detS = det(S);
    detG = det(G);
    C_inp_one = C * inp(:,1);
    A_times_Z = A * Z(:,1:end-1);
    C_times_inp = C * inp(:,2:end);
    
    P1 = Z(:,1) - (mu0 + C_inp_one); 
    LL1 = - 0.5 * P1' * iS * P1 - 0.5 * log(pi_power * detS);
    
    P2 = Z(:,2:end) - (A_times_Z + h + C_times_inp);
    LL2 = - 0.5 * sum(sum(P2 .* (iS * P2))) ...
          - (option.T-1)/2 * log(pi_power * detS);
    
    P3 = X(:,option.nonNaNidx) - Z(:,option.nonNaNidx); 
    LL3 = - 0.5 * sum(sum(P3 .* (iG * P3))) ...
          - option.tau/2 * log(pi_power * detG);
    
    
    LL =  LL1 + LL2 + LL3 - option.r(1) * sum(abs(C), 'all');
    
    DLC1 = - iS * ((C_inp_one - Z(:,1) + mu0) * inp(:,1)');
    
    DLC2 = - iS * ((C_times_inp + A_times_Z - Z(:,2:end)) ...
            * inp(:,2:end)' + h * sum(inp(:,2:end),2)');
                
    DLC3 = - option.r(1) * sign(C); %r * (det(eye(nl) - C * ones(ni,nl)) * inv(eye(nl) - C * ones(ni,nl))') * ones(ni,nl)';
    
    DLC = DLC1 + DLC2 + DLC3;
    
    LL = - LL;
    DL = - DLC;

end