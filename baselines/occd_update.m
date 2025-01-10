function [theta_est, rn, Rn] = occd_update(yn, Xn, rn, Rn, n, K, theta_est, all_but_j, var_y)

rn = rn + yn*Xn;
Rn = Rn + Xn'*Xn;


lambda = sqrt(2*var_y*n*log(K));

for j = 1:K
    rnp = rn(j) - sum(Rn(j,all_but_j{j})*theta_est(all_but_j{j}));

    theta_est(j) = sign(rnp)*max(abs(rnp) - lambda, 0)/Rn(j,j);

end

end


