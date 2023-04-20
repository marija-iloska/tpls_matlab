function [J, theta_k, Dk, Hk, Sigma] = initialize(y, H, t, var_y)


% Initialize first Hk
Hk = H(1:t,1);


% Initialize first Dk
Dk = 1/(Hk'*Hk);

% Compute iniital estimate of theta_k
theta_k = Dk*Hk'*y(1:2);

% Initial covariance of data
Sigma = Dk/var_y;


% Initial error
J = sum((y(1:2) - Hk*theta_k).^2);


end