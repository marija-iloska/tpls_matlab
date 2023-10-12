function [G, V] = pred_error_down(y, Hk, t, t0, var_y, Dk_old, Dkk_old, theta_old)

% Get k
k = length(Hk(1,1:end-1));

G = [];
THETA = [];


for i = t:-1:t0+1

    % Compute new Dk(k, i-1) from old available Dk(k, i)
    Dk_new = Dk_old/inv(eye(k) - Hk(i,1:k)'*Hk(i,1:k)*Dk_old);

    % Compute new Dk(k+1, i-1) from old available Dk(k+1, i)
    Dkk_new = Dkk_old/inv(eye(k+1) - Hk(i,1:k+1)'*Hk(i,1:k+1)*Dkk_old);

    % Compute K gain
    xd = Dkk_old*Hk(i,1:k+1)';
    K = xd/(var_y  + Hk(i,1:k+1)*xd);

    % Compute new est theta(k+1, i-1) from theta(k+1, i)
    theta_new = (eye(k+1) - K*Hk(i, 1:k+1))/(theta_old' - K'*y(i));

    % Compute Ai
    A = Hk(i, 1:k)*Dk_old*Hk(1:i-1, 1:k)'*Hk(1:i-1, k+1) - Hk(i,k+1);

    % Compute Gi
    G(end+1) = A*theta_old(end);

    % Stack theta estimates into matrix
    THETA = [THETA; theta_old'];

    % Reset old/new
    theta_old = theta_new;
    Dk_old = Dk_new;
    Dkk_old = Dkk_new;


end

V = y(t0+1:t) - sum( Hk(t0+1:t, :).*THETA , 2 );


end