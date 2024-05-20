function [theta_k, Dk, k, hk, J] = model_update(y, Ht, theta_k, Dk, J, K, var_y, t0, t, idx_H)


    %% SETUP
    k = length(theta_k);

    % Update to J(k,t) from J(k,t-1)
    J = J + (y(t) - Ht(t, 1:k)*theta_k)^2; 

    % Collect current states theta_(k, t-1) J(k,t), Dk(k, t-1)
    stay = {theta_k, idx_H, J, Dk, k};
    
    
    % Reset PE storage
    J_jump = {J, Inf, Inf};


    %% MOVES 

    % STAY SAME
    [theta_jump{1}, idx_jump{1}, J_jump{1}, Dk_jump{1}, k_jump{1}] = stay{:};

    % JUMP UP +
    if (K > k)
        [theta_jump{2}, idx_jump{2}, J_jump{2}, Dk_jump{2}, k_jump{2}] = jump_up(y, K, k, Dk, theta_k, J, Ht, t, t0, var_y) ;
    end

    % JUMP DOWN -
    if (k > 1)
        [theta_jump{3}, idx_jump{3}, J_jump{3}, Dk_jump{3}, k_jump{3}] = jump_down(y, k, Dk, theta_k, J, Ht, t, t0, var_y, K);
    end


    %% CRITERION CHOICE
    % Find Model with lowest PredError
    Js = [J_jump{1}, J_jump{2}, J_jump{3}];
    minJ = find(Js == min(Js));


    % Assign quantities to chosen model: all(t-1)
    hk = Ht(t, idx_jump{minJ});
    k = k_jump{minJ};
    Dk = Dk_jump{minJ};
    theta_k = theta_jump{minJ};
    J = J_jump{minJ};





end