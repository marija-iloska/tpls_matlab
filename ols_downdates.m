function [theta_k, Dk, Hk, J, H] = ols_downdates(y, theta_k, Dk, min_k, H, t, t0, var_y, J_old)


% Code for D change (before update)
% Let be the column removed (or min_k)
dx = length(H(1,:));
K = length(Dk(1,:));
idx = 1:K;
idx = setdiff(idx, min_k);


% Update H
H = H(:, [idx, K+1:dx, min_k]);

% Hk update
Hk = H(1:t, 1:K);

Dknew = Dk(idx, idx);
Dknew(K, 1:K-1) = Dk(min_k, idx);
Dknew(:, K) = Dk([idx, min_k], min_k);

% Now Dk downdate code
DK11 = Dknew(1:K-1, 1:K-1);
DK12 = Dknew(1:K-1, K);
%DK21 = Dknew(K, 1:K-1);
DK22 = Dknew(K,K);

Dk = DK11 - DK12*DK12' /( DK22^2 );


idx = 1:length(theta_k);
idx = setdiff(idx, min_k);

% Ratio
ratio = DK12/DK22;


% Update rest
theta_k(idx) = theta_k(idx) - ratio*theta_k(min_k);

% Compute Jk ---> Jk+
%J_old = J_old + theta_k(min_k)^2/DK22;


% Rest
theta_k(min_k) = [];
k = length(theta_k);

% Compute predictive Jk ---> Jk-
%J =  (y(t) - Hk(t, :)*theta_k)^2;
%J = sum( (y(1:t) - Hk*theta_k).^2);

% SCRIPT TO COMPUTE TRUE PREDICTIVE ERROR
[G, V] = pred_error(y, Hk, t, t0, var_y, J_old);
J = J_old - ( G*G' - 2*G*V );




end