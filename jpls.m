function [theta_store, Hk, k_store, models_sorted, count_sorted, idx_orls, J_pred, J_incr, J_pred_tot, J_dec, J_curr] = jpls(y, H, dy, var_y, n, Nb, D)

% Store
H_true = H;
T = length(H(:,1));

% Initialize model order
k = dy;

%Initialize using t data points
[J, J_curr, theta_k, Dk, Hk, ~] = initialize(y, H, n, k, var_y);


[~, idx_sort] = sort(theta_k, 'descend');
H = H(:, idx_sort);
k = floor(dy/2);
[J, J_curr, theta_k, Dk, Hk,~] = initialize(y, H, n, k, var_y);

J_incr = J;
J_pred = 0;
J_pred_tot = 0;

J_up = 0;
J_do = 0;
J_st = 0;

M ={};
theta_store = {};
Dk_jump = {Dk, Dk, Dk};

% Start time loop
for t = n+1:T-1

    % Reset
    J_jump = {J, Inf, Inf};

    % STAY SAME
    J_jump{1} =  sum( (y(1:t) - H(1:t, 1:k)*theta_k).^2);
    Dk_jump{1} = Dk;
    k_jump{1} = k;
    H_jump{1} = H;
    theta_jump{1} = theta_k;

    % JUMP UP +
    if (dy > k)
        [theta_jump{2}, H_jump{2}, J_jump{2}, Dk_jump{2}, k_jump{2}] = jump_up(y, dy, k, Dk, theta_k, J, H, t);
    end

    % JUMP DOWN -
    if (k > 1)
        [theta_jump{3}, H_jump{3}, J_jump{3}, Dk_jump{3}, k_jump{3}] = jump_down(y, k, Dk, theta_k, J, H, t);
    end



    % Find Model with lowest PredError
    Js = [J_jump{1}, J_jump{2}, J_jump{3}];
    minJ = find(Js == min(Js));
    J_up(end+1) = J_jump{2};
    J_st(end+1) = J_jump{1};
    J_do(end+1) = J_jump{3};

    % Assign quantities to chosen model
    H = H_jump{minJ};
    k = k_jump{minJ};
    Dk = Dk_jump{minJ};
    theta_k = theta_jump{minJ};
    J = J_jump{minJ};

    % Update and store terms
    Hk = H(1:t, 1:k);
    k_store(t) = k;
    theta_store{end+1} = theta_k;
    J_pred(end+1) = J;
    J_incr(end+1) = J_incr(end) + (y(t) - H(t, 1:k)*theta_k)^2;
    J_pred_tot(end+1) = J_pred_tot(end) + J;
    J_curr(end+1) = (y(t) - H(t, 1:k)*theta_k)^2;


    % Check which model was selected at time t
    [~, idx_orls] = ismember(Hk(1,:), H_true(1,:));
    M{end+1} = [sort(idx_orls, 'ascend'), zeros(1, dy - length(idx_orls)) ];

    % TIME UPDATE
    [theta_k, Dk, ~] = time_update(y, Hk, t, theta_k, var_y, Dk, J);




end
J_dec = {J_st, J_up, J_do};

% Apply models
[models_sorted, count_sorted, idx_orls] = model_sorting(M, Nb, dy);



end