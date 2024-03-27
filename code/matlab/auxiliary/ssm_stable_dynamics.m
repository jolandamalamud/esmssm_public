function stab_par = ssm_stable_dynamics(EL, ER)
    
    nl = size(ER, 1);
    np = length(ER(:)) - nl^2;
    HH = [];
    for i=1:nl, HH = blkdiag(HH,EL'); end
    Y = ER';
    hh = Y(1:end)';
    ii = -Inf * ones(nl, nl); ii_lb = ii - diag(diag(ii));
    ii_lb(isnan(ii_lb)) = -1;
    lb = [ii_lb(:); -Inf * ones(np, 1)];
    ii = Inf * ones(nl, nl); ii_ub = ii - diag(diag(ii));
    ii_ub(isnan(ii_ub)) = 1;
    ub = [ii_ub(:); Inf * ones(np, 1)];
    stab_par = quadprog(HH,-hh,[],[],[],[],lb,ub);

%         nl = size(A,1);
%         LB = -Inf * ones(nl, nl);
%         LB = LB - diag(diag(LB));
%         LB(isnan(LB)) = 0;
%         UB = Inf * ones(nl, nl);
%         UB = UB - diag(diag(UB));
%         UB(isnan(UB)) = 1;
% 
% 
%         CON = zeros(nl^2, nl^2);
%         [u, s, v] = svd(A);
%         i = 1;
%         while s(1) > 1 && i < 16
%             con = u(:,1) * v(:,1)';
%             CON(i,:) = con(:);
%             P = kron(eye(nl), Ezz);
%             a = quadprog(P, -Ezz_1(:), CON, ones(nl^2, 1), [], [], LB, UB);
%             A_stable = reshape(a, nl, nl);
%             [u, s, v] = svd(A_stable);
%             i = i + 1;
%         end
%        
end