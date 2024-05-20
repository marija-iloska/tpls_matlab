function [J, e, mse] = true_PE(y, H, t0, T, idx1, var_y)


% Initialize
Dp = inv(H(1:t0, idx1)'*H(1:t0, idx1));
theta_p = Dp*H(1:t0, idx1)'*y(1:t0);

e = [];
J = [];
mse = [];

for i = t0+1:T

    e(end+1) = y(i) - H(i,idx1)*theta_p;
    J(end+1) = sum(e.^2);

    mse(end+1) = var_y + var_y*H(i,idx1)*Dp*H(i,idx1)';

    if (i == T)
        break
    end

    % Compute theta_(k+1, t-1), check Dk indices
    [theta_p, Dp] = time_update(y(i), H(i, idx1), theta_p, var_y, Dp);


end

% Predictive Residual error
e_final = e(end);

end