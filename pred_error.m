function [G, E] = pred_error(y, Hk, t, t0, var_y, J_old, theta, D)

% Get k
k = length(Hk(1,1:end-1));


% k+1 x k+1 at t0-1
% Dkk_old = inv(Hk(1:t0-1, :)'*Hk(1:t0-1, :));
% theta_kk_old = Dkk_old*Hk(1:t0-1, :)'*y(1:t0-1);

% k x k at t0-1
Dk = inv(Hk(1:t0-1, 1:k)'*Hk(1:t0-1, 1:k));
theta_k = Dk*Hk(1:t0-1, 1:k)'*y(1:t0-1);

% k+1  x  k+1 at t0-1
[theta_kk, Dkk, ~, ~] = ols_updates(y, Hk, k, 1, t0, Dk, theta_k);

G = [];
THETA = [];

for i = t0:t


    % Compute Ai
    A = Hk(i, 1:k)*Dk*Hk(1:i-1, 1:k)'*Hk(1:i-1, k+1) - Hk(i,k+1);

    % Compute Gi
    G(end+1) = A*theta_kk(end);

    % Stack theta estimates into matrix
    THETA = [THETA; theta_k'];

    if (i == t)
         disp('stop')
    end

    % Compute theta_(k+1, t-1), check Dk indices
    [theta_kk, Dkk, ~] = time_update(y, Hk(1:i, :), i, theta_kk, var_y, Dkk, J_old);
    [theta_k, Dk, ~] = time_update(y, Hk(1:i, 1:k), i, theta_k, var_y, Dk, J_old);

end

E = y(t0:t) - sum( Hk(t0:t, 1:k).*THETA , 2 );

%theta'


end