function [J, theta_k, Dk, Hk, Sigma] = initialize(y, H, t, k, var_y)


% Initialize first Hk
Hk = H(1:t,1:k);


% Initialize first Dk
Dk = inv(Hk'*Hk);

% Compute iniital estimate of theta_k
theta_k = Dk*Hk'*y(1:t);

% Initial covariance of data
Sigma = Dk/var_y;

% Initial error
J = sum((y(1:t+1) - H(1:t+1, 1:k)*theta_k).^2);

% Initialize first Hk
Hk = H(:,1:k);

end