function [theta_k, Dk, Hk] = ols_updates(y, H, k, j, t, Dk, theta_k)

    % Current input data
    Hk = H(1:t, 1:k);

    Pk_norm = eye(t) - Hk*Dk*Hk';

    % Take the new observation column h(k+1)
    hk = H(1:t,k+j);

    % Reuseable terms
    temp = hk'*Pk_norm;
    xd = temp*y;
    d = Dk*Hk'*hk;

    % Compute terms of D(k+1)
    DK22 = 1 / ( temp*hk );
    DK12 = - d*DK22;
    DK21 = DK12';
    DK11 = Dk + DK22*(d*d');

    % Update D(k) to D(k+1)
    Dk = [ DK11, DK12 ;  DK21, DK22 ];

    % Update theta_k
    theta_k = [ theta_k - d*xd*DK22;    xd*DK22 ];

    % Update Hk
    Hk = [Hk, hk];



end
