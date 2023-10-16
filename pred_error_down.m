function [G, E, start] = pred_error_down(y, H, t, t0, var_y, J_old, j, start)


% k+1 x k+1 at t0
% Extract starting materials
[theta_kk, Dkk] = start{:};

% Get k
k = length(theta_kk(1:end-1));



% Downdate and get new starts
% k x k  at t0
[theta_k_test, Dk_test, Hk, ~] = ols_downdates(theta_kk, Dkk, j, H, t);


Dk = inv(Hk(1:t0, 1:k)'*Hk(1:t0, 1:k));
theta_k = Dk*Hk(1:t0, 1:k)'*y(1:t0);

start = {theta_k, Dk};


% Initialize
G = [];
THETA = [];
V = Hk(t0, 1:k)'*Hk(t0, k+1);

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

% Compute predictive residual error
E = y(t0+1:t) - sum( Hk(t0+1:t, 1:k).*THETA , 2 );


end