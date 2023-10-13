function [theta_k, H, J,  Dk, k] = jump_up(y, dx, k, Dk, theta_k, J, H, t, t0, var_y)


for j = 1:(dx - k)

    % Update current theta by jth basis function
    [theta_temp, D_temp, Hk_temp,  Hnew] = ols_updates(y, H, k, j, t, Dk, theta_k);

    % Compute J(k,t) ---> J(k+1, t)
    [G, E] = pred_error(y, Hk_temp, t, t0, var_y, J, theta_temp, D_temp);
    %J_temp = J - (G*G' - 2*G*E);
    J_temp = J + (G*G' + 2*G*E);

                           
    % Corresponding variables store
    H_store{j} = Hnew;
    D_store{j} = D_temp;
    theta_store{j} = theta_temp;
    J_store(j) = J_temp;


end

% Choose min J to update
min_idx = find(J_store == min(J_store));


% Update all parameters
theta_k = theta_store{min_idx};
H  = H_store{min_idx};  % For optimization I could save indices here
J = J_store(min_idx);
Dk = D_store{min_idx};
k = length(theta_k);

end