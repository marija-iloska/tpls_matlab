function [J, theta_k, Dk, Hk, Sigma] = initialize_D(y, H, t, var_y)


% Initialize first Hk
Hk = H(1:t,:);


% Initialize first Dk
Dk = inv(Hk'*Hk);

% Compute iniital estimate of theta_k
theta_k = Dk*Hk'*y(1:t);

% Initial covariance of data
Sigma = Dk/var_y;

% Initial error
J = sum((y(1:t+1) - H(1:t+1, :)*theta_k).^2);


end