function [E_add, E_rmv] = expectations(y, H, t0, T, idx1, var_y, theta)


% Initialize
K = length(H(1,:));
idx0 = setdiff(1:K, idx1);
p = length(idx1);
pj = K - p;

% Sort H
H = H(:, [idx1, idx0]);
theta = theta(idx1);

% Get true estimates
Dp = inv(H(1:t0, 1:p)'*H(1:t0, 1:p));
theta_p = Dp*H(1:t0, 1:p)'*y(1:t0);



for t = t0+1:T

    
    % ADDITION  ===================================================
    for j = 1:pj
        % It updates DIM at t-1
        [~, Dpp, ~,~] = ols_updates(y, H, p, j, t, Dp, theta_p);

        % D(p+1, t-1)
        b_add = - Dpp(1:end-1,end)/Dpp(end, end);
        Q_add = H(t, p+j) - H(t, 1:p)*b_add;

        % Expectation E(p+1) - E(p)   single and batch
        E_add(t,j) = var_y*Q_add^2*Dpp(end,end);


    end

    % REMOVAL  ===================================================
    for j = 1:p

        idx = setdiff(1:p,j);
        
        % Get Dk tilde by swapping 
        Dp_swap = Dp(idx, idx);
        Dp_swap(p, 1:p-1) = Dp(j, idx);
        Dp_swap(:, p) = Dp([idx, j], j);
        
        % D(p, t-1)
        b_rmv = - Dp_swap(1:end-1,end)/Dp_swap(end, end);
        Q_rmv = H(t,j) - H(t, idx)*b_rmv;


        % Expectation E(p-1) - E(p)   single and batch
        E_rmv(t,j) = Q_rmv^2*(theta(j)^2 - var_y*Dp_swap(end,end));
     
    end


    % Compute theta_(k+1, t-1), check Dk indices
    [theta_p, Dp, ~] = time_update(y, H(1:t, 1:p), t, theta_p, var_y, Dp, 1);



end


end