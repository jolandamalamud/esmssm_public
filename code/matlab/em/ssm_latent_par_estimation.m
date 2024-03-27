function [mu0, A, h, C] = ssm_latent_par_estimation(option, T, X, z, inp, ...
    Ezz_1_2toT, Ez_1z_1_2toT, Ez_2toT, Ez_1_2toT, inpEz, inpEz_1, iinp_2toT, inpinp)
    
    fminopt=optimoptions('fminunc','display','off','GradObj','on',...
        'TolFun',1e-8,'algorithm','trust-region');

    A = Ezz_1_2toT / Ez_1z_1_2toT;
    
    if option.bias
        h = (Ez_2toT - A * Ez_1_2toT) / ((T-1));
    else 
        h = zeros(option.dimZ, 1);
    end
    
    if option.input
        C0 = (inpEz' - h * iinp_2toT' - A * inpEz_1') / inpinp;
    else
        C = zeros(option.dimZ, option.dimInp);
    end
    
    mu0 = (z(:,1)- C * inp(:,1));
    
    tol = 1e-4;
    i = 1;
    
    par(:, :, i) = [A, h, C, mu0];
    matrixdistance = 10;
    while matrixdistance > tol
        i = i + 1;
        
        A = (Ezz_1_2toT -  h * Ez_1_2toT' - C * inpEz_1) / Ez_1z_1_2toT;
   
       if option.bias
           h = (Ez_2toT - A * Ez_1_2toT - C * iinp_2toT) / (T-1);
       else 
           h = zeros(nl, 1);
       end
       
       if option.input
           C = (inpEz' - h * iinp_2toT' - A * inpEz_1' - mu0 * inp(:,1)') / inpinp;
           par.A = diag(A); par.W = A - diag(A); par.h = h; par.S = S; par.G = G;
           C = fminunc(@(x)ssm_joint_likelihood_Classo(x, X, inp, z, par, r), C0, fminopt);
       else
           C = zeros(option.dimZ, option.dimInp);
       end
        
        mu0 = (z(:,1)- C * inp(:,1));
        
        par(:, :, i) = [A, h, C, mu0];
        matrixdistance = sum(sum((abs(par(:, :, i) - par(:, :, i-1)))));
    end
    
end