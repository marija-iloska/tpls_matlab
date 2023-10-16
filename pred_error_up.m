function [G, E, start] = pred_error_up(y, Hk, t, t0, var_y, J_old, start)

% Get k
k = length(Hk(1,1:end-1));


% Extract starting materials
[theta_test, Dk_test] = start{:};
Dk = inv(Hk(1:t0, 1:k)'*Hk(1:t0, 1:k));
theta_k = Dk*Hk(1:t0, 1:k)'*y(1:t0);

% Update and get new starts
% k+1  x  k+1 at t0
[theta_kk, Dkk, ~, ~] = ols_updates(y, Hk, k, 1, t0+1, Dk, theta_k);
start = {theta_kk, Dkk};

% Initialize
G = [];
THETA = [];
V = Hk(1:t0, 1:k)'*Hk(1:t0, k+1);

for i = t0+1:t


    % Compute Ai
    A = Hk(i, 1:k)*Dk*V - Hk(i,k+1);

    % Compute Gi
    G(end+1) = A*theta_kk(end);

    % Stack theta estimates into matrix
    THETA = [THETA; theta_k'];

    if (i == t)
        break
    end

    % Compute theta_(k+1, t-1), check Dk indices
    [theta_kk, Dkk, ~] = time_update(y, Hk(1:i, :), i, theta_kk, var_y, Dkk, J_old);
    [theta_k, Dk, ~] = time_update(y, Hk(1:i, 1:k), i, theta_k, var_y, Dk, J_old);

    % Update V
    V = V + Hk(i, 1:k)'*Hk(i, k+1);

end

% Predictive Residual error
E = y(t0+1:t) - sum( Hk(t0+1:t, 1:k).*THETA , 2 );



end