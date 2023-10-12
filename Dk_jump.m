function [theta_k, Dk] = Dk_jump(y, Hk, t)


    D1 = 1/(Hk(1:t, 1)'*Hk(1:t, 1));
    k = length(Hk(1,:));
    theta_k = D1*Hk(1:t,1)'*y(1:t);

    P = eye(t) - Hk(1:t,1)*D1*Hk(1:t,1)';

    % Reuseable terms
    temp = Hk(1:t,2:k)'*P;
    xd = temp*y(1:t);
    d =  D1*Hk(1:t, 1)'*Hk(1:t,2:k);


    % Jump m steps
    D22 = inv(Hk(1:t, 2:k)'*P*Hk(1:t,2:k));
    D12 = - d*D22;
    D21 = D12';
    D11 = D1 - D12*Hk(1:t, 2:k)'*Hk(1:t, 1)*D1;
    %D11 = D1 + D22*(d*d');

    % Update covariance for new data
    Dk = [D11 D12; D21 D22];

   % Update theta_k
    theta_k = [ theta_k - d*D22*xd;    D22*xd ];



end
