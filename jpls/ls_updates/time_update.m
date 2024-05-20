function [theta_k, Dk] = time_update(y, hk, theta_k, var_y, Dk)


    % Get current dimension
    k = length(theta_k);

    % Update covariance for new data
    Sigma = var_y*Dk;

    % Current error
    et = (y - hk*theta_k);

    % Update gain
    K = Sigma*hk'/(var_y + hk*Sigma*hk');

    % Update estimate
    theta_k = theta_k + K*et;

    % Update covariance
    Sigma = (eye(k) - K*hk)*Sigma;
    Dk = Sigma/var_y;

    


end
