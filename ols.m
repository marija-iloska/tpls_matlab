function [theta_k, Dk, Jk, error_store] = ols(y, H, epsilon, K)

% Initialize first Hk
Hk = H(:,1);

T = length(H(:,1));

% Initialize first Dk
Dk = 1/(Hk'*Hk);


% Compute iniital estimate of theta_k
theta_k = Dk*Hk'*y;


% Compute first Jmink 
Jk = sum( (y - Hk*theta_k).^2 );
error_store = Jk;


k = 1;
while (Jk > epsilon)

    % Order of model update
    k = k + 1;

    % Compute projection matrix Pk
    Pk_norm = eye(T) - Hk*Dk*Hk';

    % Take the new observation column h(k+1)
    hk = H(:,k);

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
    Jk = Jk - xd^2*DK22;
    error_store(k) = Jk;

    % Max observations
    if (k == K)
        break
    end


end



end