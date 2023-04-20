function [theta_k, Sigma, J] = time_update(y, Hk, t, theta_k, var_y, Dk, J)


    k = length(theta_k);

    % Update covariance for new data
    Sigma = var_y*Dk;

    % Current error
    et = (y(t) - Hk(t,:)*theta_k);

    % Update gain
    K = Sigma*Hk(t,:)'/(var_y + Hk(t,:)*Sigma*Hk(t,:)');

    % Update estimate
    theta_k = theta_k + K*et;

    % Update covariance
    Sigma = (eye(k) - K*Hk(t,:))*Sigma;

    % Update error
    J = J + et^2;


end
