function [J, e_final] = pred_error_lasso(y, H, t, t0, var_y, theta)

% Get k
idx1 = find(theta ~= 0);


% theta is at t0

% Initialize
Dk = inv(H(1:t0, idx1)'*H(1:t0, idx1));
theta_k = Dk*H(1:t0, idx1)'*y(1:t0);
e = 0;
J = 0;

for i = t0+1:t

    e(end+1) = y(i) - H(i,idx1)*theta_k;


    if (i == t)
        break
    end

    % Compute theta_(k+1, t-1), check Dk indices
    [theta_k, Dk, ~] = time_update(y, H(1:i, idx1), i, theta_k, var_y, Dk, J);

end

% Predictive Residual error
J = sum(e.^2);

e_final = e(end);

end