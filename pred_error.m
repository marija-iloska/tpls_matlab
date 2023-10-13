function [G, E, Dkk_old, Dk_old] = pred_error(y, Hk, t, t0, var_y, J_old, theta, Dk)

% Get k
k = length(Hk(1,1:end-1));

% SCRIPT TO COMPUTE TRUE PREDICTIVE ERROR
% k+1 x k+1 at t0-1
% [theta_kk_old, Dkk_old] = Dk_jump(y, Hk, 2);
% [theta_k_old, Dk_old] = Dk_jump(y, Hk(:,1:k), 2);

Dkk_old = inv(Hk(1:t0-1, :)'*Hk(1:t0-1, :));
theta_kk_old = Dkk_old*Hk(1:t0-1, :)'*y(1:t0-1);

% k x k
Dk_old = inv(Hk(1:t0-1, 1:k)'*Hk(1:t0-1, 1:k));
theta_k_old = Dk_old*Hk(1:t0-1, 1:k)'*y(1:t0-1);

G = [];
V = [];
THETA = [];

for i = t0:t

    % Compute theta_(k+1, t-1), check Dk indices
    [theta_kk_new, Dkk_new, ~] = time_update(y, Hk(1:i, :), i, theta_kk_old, var_y, Dkk_old, J_old);
    [theta_k_new, Dk_new, ~] = time_update(y, Hk(1:i, 1:k), i, theta_k_old, var_y, Dk_old, J_old);


    % Compute Ai
    A = Hk(i, 1:k)*Dk_old*Hk(1:i-1, 1:k)'*Hk(1:i-1, k+1) - Hk(i,k+1);

    % Compute Gi
    G(end+1) = A*theta_kk_old(end);

    % Compute Vi
    %V(end+1) = (y(i) - Hk(i, 1:k)*theta_k_old)*G(end);

    % Stack theta estimates into matrix
    THETA = [THETA; theta_k_old'];

    % Reset old/new
    theta_k_old = theta_k_new;
    theta_kk_old = theta_kk_new;
    if (i == t)
        Dk_old
        Dk
        Dkk_old
    end
    Dk_old = Dk_new;
    Dkk_old = Dkk_new;


end

E = y(t0:t) - sum( Hk(t0:t, 1:k).*THETA , 2 );

%theta'


end