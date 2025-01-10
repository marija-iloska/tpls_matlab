function [J] = pred_error_baselines(y, H, t, t0, var_y, theta)

% Get k
idx1 = find(theta ~= 0);


% theta is at t0

% Initialize
Dk = inv(H(1:t0, idx1)'*H(1:t0, idx1));
theta_k = Dk*H(1:t0, idx1)'*y(1:t0);
e = [];

for i = t0+1:t

    e(end+1) = y(i) - H(i,idx1)*theta_k;

    if (i == t)
        break
    end

    % Compute theta_(k+1, t-1), check Dk indices
    [theta_k, Dk] = time_update(y(i), H(i, idx1), theta_k, var_y, Dk);

end

% Predictive Residual error
J = sum(e.^2);

end