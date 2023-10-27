function [J, e] = true_PE(y, H, t0, T, idx1, var_y)


% Initialize
Dk = inv(H(1:t0, idx1)'*H(1:t0, idx1));
theta_k = Dk*H(1:t0, idx1)'*y(1:t0);

e = [];
J = [];

for i = t0+1:T

    e(end+1) = y(i) - H(i,idx1)*theta_k;
    J(end+1) = sum(e.^2);

    if (i == T)
        break
    end

    % Compute theta_(k+1, t-1), check Dk indices
    [theta_k, Dk, ~] = time_update(y, H(1:i, idx1), i, theta_k, var_y, Dk, J);


end

% Predictive Residual error
e_final = e(end);

end