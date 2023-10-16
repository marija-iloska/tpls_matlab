function [theta_k, idx_H, J,  Dk, k, start] = jump_up(y, dx, k, Dk, theta_k, J, H, t, t0, var_y, start)



for j = 1:(dx - k)

    % Update current theta by jth basis function
    [theta_store{j}, D_store{j}, Hk_temp,  idx_store{j}] = ols_updates(y, H, k, j, t, Dk, theta_k);

    % Compute PE J(k,t) ---> J(k+1, t)
    [~, ~, start_store{j}] = pred_error_up(y, Hk_temp, t, t0, var_y, J, start);
    [G, E] = pred_error(y, Hk_temp, t, t0, var_y, J);
    J_store(j) = J + (G*G' + 2*G*E);

end

% Choose min J to update
min_idx = find(J_store == min(J_store));


% Update all parameters
theta_k = theta_store{min_idx};
idx_H  = idx_store{min_idx};  % For optimization I could save indices here
J = J_store(min_idx);
Dk = D_store{min_idx};
start = start_store{min_idx};

% Final dimension
k = k + 1;

end