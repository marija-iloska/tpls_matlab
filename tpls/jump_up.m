function [theta_k, idx_H, J,  Dk, k] = jump_up(y, dx, k, Dk, theta_k, J, Ht, t, t0, var_y)



parfor j = 1:(dx - k)

    % Update current theta by jth basis function
    [theta_store{j}, D_store{j}, Hk_temp,  idx_store{j}] = ols_updates(y, Ht, k, j, t, Dk, theta_k);

    % Compute PE J(k,t) ---> J(k+1, t)
    [G, E] = pred_error(y, Hk_temp, t, t0, var_y, J);
    J_store(j) = J + (G*G' + 2*G*E);

end

% Choose min J to update
min_idx = find(J_store == min(J_store));


% Update all parameters
theta_k = theta_store{min_idx};
idx_H  = idx_store{min_idx}; 
J = J_store(min_idx);
Dk = D_store{min_idx};

% Final dimension
k = k + 1;

end