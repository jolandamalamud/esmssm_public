function [z_s, P, P_a, MLL] = ssm_inference_KS(option, parameter, X, inp)

    B=parameter.B; G=parameter.G; A=parameter.A; W=parameter.W; 
    h=parameter.h; C=parameter.C; S=parameter.S; mu0=parameter.mu0;

    LL = 0;
    
    z_p = zeros(option.dimZ, option.T);
    V_p = zeros(option.dimZ, option.dimZ, option.T);
    K = zeros(option.dimZ, option.dimZ, option.T);
    z_f = zeros(option.dimZ, option.T);
    V_f = zeros(option.dimZ, option.dimZ, option.T);
    J = zeros(option.dimZ, option.dimZ, option.T);
    z_s = zeros(option.dimZ, option.T);
    V_s = zeros(option.dimZ, option.dimZ, option.T);
    
%   kalman filter
    for t = 1:option.T
        if t == 1
            z_p(:, t) =  mu0 + C * inp(:, 1);
            V_p(:, :, t) = S;
        else
            z_p(:, t) = (A + W) * z_f(:, t-1) + h + C * inp(:, t);
            V_p(:, :, t) = (A + W) * V_f(:, :, t-1) * (A + W)' + S;
        end
        if ~ismember(option.timing, t)
            z_f(:, t) = z_p(:, t);
            V_f(:, :, t) = V_p(:, :, t);
        else
            K(:, :, t) = V_p(:, :, t) * inv(V_p(:, :, t) + G);
            z_f(:, t) = z_p(:, t) + K(:, :, t) * (X(:, t) - z_p(:, t));
            V_f(:, :, t) = V_p(:, :, t) - K(:, :, t) * V_p(:, :, t);
            % marginal loglikelihood p(X)
            H = log(det(G + V_p(:, :, t)));
            LL = LL + H + (X(:, t) - z_p(:, t))' * inv(V_p(:, :, t) + G)...
                * (X(:, t) - z_p(:, t));
        end
    end
    % marginal loglikelihood p(X)
    MLL = - 0.5 * option.tau * option.dimZ * log(2 * pi) - 0.5 * LL;
    
%   kalman smoother
    z_s(:, option.T) = z_f(:, option.T);
    V_s(:, :, option.T) = V_f(:, :, option.T);
    for t = option.T:-1:2
        J(:, :, t-1) = V_f(:, :, t-1) * (A + W)' * inv(V_p(:, :, t));
        if ismember(option.timing, t)
            z_s(:, t-1) = z_f(:, t-1) + J(:, :, t-1) * (z_s(:, t) - ...
                (A + W) * z_f(:, t-1)) + (J(:, :, t-1) * (A + W) - eye(option.dimZ)) ...
                * (V_f(:, :, t-1) * (A + W)' * inv(S) * h);
        else
            z_s(:, t-1) = z_f(:, t-1) + J(:, :, t-1) * (z_s(:, t) - ...
                            (A + W) * z_f(:, t-1)) + (J(:, :, t-1) * (A + W) - eye(option.dimZ)) ...
                            * (V_f(:, :, t-1) * (A + W)' * inv(S) * (h + C * inp(:, t)));
        end
        V_s(:, :, t-1) = V_f(:, :, t-1) + J(:, :, t-1) * (V_s(:, :, t) ...
                            - V_p(:, :, t)) * J(:, :, t-1)';
    end

    % one lag covariance
    for t = option.T:-1:2
        V_a(:, :, t) = V_s(:, :, t) * J(:, :, t-1)';
    end

    % set response
    for t = 1:option.T
        z(:, t) = z_s(:, t); % E[z(t)]
        P(:, :, t) = V_s(:, :, t) + z_s(:, t) * z_s(:, t)'; % E[z(t)z(t)']
        if ~(t == 1)
            P_a(:, :, t) = V_a(:, :, t) + z_s(:, t) * z_s(:, t-1)'; % E[z(t)z(t-1)']
        end
    end

end