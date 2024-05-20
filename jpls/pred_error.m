function [G, E] = pred_error(y, Hk, t, t0, var_y, J_old)

% Get k
k = length(Hk(1,1:end-1));


% k+1 x k+1 at t0-1
% Dkk_old = inv(Hk(1:t0-1, :)'*Hk(1:t0-1, :));
% theta_kk_old = Dkk_old*Hk(1:t0-1, :)'*y(1:t0-1);

% k x k at t0
Dk = inv(Hk(1:t0, 1:k)'*Hk(1:t0, 1:k));
theta_k = Dk*Hk(1:t0, 1:k)'*y(1:t0);

% k+1  x  k+1 at t0
[theta_kk, Dkk, ~, ~] = ols_updates(y, Hk, k, 1, t0+1, Dk, theta_k);

G = [];
THETA = [];

for i = t0+1:t


    % Compute Ai
    A = Hk(i, 1:k)*Dk*Hk(1:i-1, 1:k)'*Hk(1:i-1, k+1) - Hk(i,k+1);

    % Compute Gi
    G(end+1) = A*theta_kk(end);

    % Stack theta estimates into matrix
    THETA = [THETA; theta_k'];


    % Compute theta_(k+1, t-1), check Dk indices
    [theta_kk, Dkk] = time_update(y(i), Hk(i, :), theta_kk, var_y, Dkk);
    [theta_k, Dk] = time_update(y(i), Hk(i, 1:k), theta_k, var_y, Dk);

end

E = y(t0+1:t) - sum( Hk(t0+1:t, 1:k).*THETA , 2 );



end