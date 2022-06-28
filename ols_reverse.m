function [theta_k, theta_clean, Dk_update] = ols_reverse(theta_clean, Dk, min_k, store, dx)


% Code for D change (before update)
% Let be the column removed (or min_k)
K = length(Dk(1,:));
idx = 1:K;
idx = setdiff(idx, min_k);

Dknew = Dk(idx, idx);
Dknew(K, 1:K-1) = Dk(min_k, idx);
Dknew(:, K) = Dk([idx, min_k], min_k);

% Now Dk downdate code
DK11 = Dknew(1:K-1, 1:K-1);
DK12 = Dknew(1:K-1, K);
%DK21 = Dknew(K, 1:K-1);
DK22 = Dknew(K,K);

Dk_update = DK11 - DK12*DK12' /( DK22^2 );


idx = 1:length(theta_clean);
idx = setdiff(idx, min_k);

% Ratio
ratio = DK12/DK22;

% Update rest
theta_clean(idx) = theta_clean(idx) - ratio*theta_clean(min_k);

% East
theta_clean(min_k) = [];


% Theta original
theta_k = zeros(dx,1);
theta_k(setdiff(1:dx,store)) = theta_clean;



end