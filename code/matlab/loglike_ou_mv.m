function ll = loglike_ou_mv(x, mu, sigma, lambda, dt)
    T = size(x,2);
    for t = 2:T
        sigma0 = sigma^2 * (1-exp(-2*lambda*dt(t-1))) / (2*lambda);
        sigma0 = sqrt(sigma0);

        prefactor = 1 / sqrt(2*pi*sigma0^2);

        f(t-1) = prefactor * exp(- (x(:,t) - x(:,t-1) * exp(-lambda*dt(t-1)) ...
            - mu * (1-exp(-lambda*dt(t-1))))^2 / (2*sigma0^2));
    end
    
    idx = f == 0;
    f(idx) = 10^-8;
    
    ll = -sum(log(f));
    
end