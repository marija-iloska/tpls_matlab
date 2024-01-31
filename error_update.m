function [mse] = error_update(H, t0, T, var_y, idx1)


% Initialize
Dp = inv(H(1:t0, idx1)'*H(1:t0, idx1));
k = length(idx1);

mse = 0;

for i = t0+1:T

    mse(end+1) = var_y + var_y*H(i,idx1)*Dp*H(i,idx1)';


    % Update covariance for new data
    Sigma = var_y*Dp;

    % Update gain
    K = Sigma*H(i,idx1)'/(var_y + H(i,idx1)*Sigma*H(i,idx1)');

    % Update covariance
    Sigma = (eye(k) - K*H(i,idx1))*Sigma;
    Dp = Sigma/var_y;


end






end