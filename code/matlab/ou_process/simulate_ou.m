function x = simulate_ou(x0, mu, sigma, lambda, t)
    % Simulates a sample path
    T = length(t);
    ns = length(x0);
    x = zeros(T,1);
    x(1) = x0;
    dt = t(2:end) - t(1:end-1);
    
    x = x0;
    
    for i = [2:T]
        exp_minus_lambda_deltat = exp(-lambda * dt(i-1));
        dWt = sqrt((1-exp(-2*lambda*dt(i-1)))/(2*lambda)) * randn;
        x(i) = exp_minus_lambda_deltat * x(i-1) + ...
            (1 - exp_minus_lambda_deltat) * mu + sigma * dWt;
    end
end