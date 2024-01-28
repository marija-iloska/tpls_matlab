function [Es_add, Es_rmv, Eb_add, Eb_rmv] = expectations(y, H, t0, T, idx1, var_y, theta)


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
        Qb_add = H(1:t, p+j) - H(1:t, 1:p)*b_add;
        Q_add = Qb_add(end);


        % Expectation E(p+1) - E(p)   single and batch
        Eb_add(t,j) = var_y*Qb_add'*Qb_add*Dpp(end,end) - var_y; % - 2*var_y;
        Es_add(t,j) = var_y*Q_add^2*Dpp(end,end);


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
        Qb_rmv = H(1:t, j) - H(1:t, idx)*b_rmv;
        Q_rmv = Qb_rmv(end);
        QQ = Qb_rmv'*Qb_rmv;

        % Expectation E(p-1) - E(p)   single and batch
        Es_rmv(t,j) = Q_rmv^2*(theta(j)^2 - var_y*Dp_swap(end,end));
        
        % New Q
        q = H(t,j) - H(t, idx)*b_rmv;
        Qtilde = H(1:t-1, j) - H(1:t-1, idx)*b_rmv;
        T1 = q^2*(theta(j)^2 - var_y*Dp_swap(end,end));
        T2 =  Qtilde'*Qtilde*(theta(j)^2 + var_y*Dp_swap(end,end));
        Eb_rmv(t,j) =  2*var_y + T1 - T2 + 2*theta(j)^2/Dp_swap(end,end);
        %Eb_rmv(t,j) = QQ*(theta(j)^2  - var_y*Dp_swap(end,end)) + 2*var_y;

    end





    % Compute theta_(k+1, t-1), check Dk indices
    [theta_p, Dp, ~] = time_update(y, H(1:t, 1:p), t, theta_p, var_y, Dp, 1);



end


end