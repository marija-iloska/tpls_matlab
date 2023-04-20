function [theta_k, Dk, Hk, J, H] = ols_updates(y, H, k, j, t, Dk, theta_k, J)

    % Current input data
    K = length(H(1,:));
    Hk = H(1:(t-1), 1:k);
    y = y(1:(t-1));

    Pk_norm = eye(t-1) - Hk*Dk*Hk';

    % Take the new observation column h(k+1)
    hk = H(1:(t-1),k+j);

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

    % Update Hk in time and k
    Hk = H(1:t, [1:k, k+j]);

    % Update original available data H (swap column order )
    H = H(:, [1:k, k+j, setdiff( (k+1):K, (k+j) ) ]);

    % Compute Jk ---> Jk+
    J = J - theta_k(end)*xd;



end
