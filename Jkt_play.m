clear all
close all
clc

% Settings
var_y = 0.01; % Variance
p_s = 0.3;   % Sparsity percent
dx = 10;      % System dimension
T = 100;     % Time series length
r = 0.5;     % Range of input data H
rt = 4;      % Range of theta


% Current dimension dtheta
k = 1;

%Create data
[y, H, theta] = generate_data(T, dx, r,rt,  p_s, var_y);

H_true = H;


% Initialize first Hk
Hk = H(1:2,1);


% Initialize first Dk
Dk = 1/(Hk'*Hk);

% Compute iniital estimate of theta_k
theta_k = Dk*Hk'*y(1:2);

% Initial covariance of data
Sigma = Dk/var_y;


% Initial error
J = sum((y(1:2) - Hk*theta_k).^2);


for t = 3:T

    J_track(t) = J;
    % COMPUTE Mk estimate     k ---> k+
    clear H_store D1_store theta1_store J1_store
    if (dx > k)
        for j = 1:(dx - k)

            % Update current theta by jth basis function
            [theta_temp, D_temp, H_temp, J_temp, Hnew] = ols_updates(y, H, k, j, t, Dk, theta_k, J);

            % COMPUTE time estimate    t ---> t+
            [theta_temp, Sigma, J_temp] = time_update(y, H_temp, t, theta_temp, var_y, D_temp, J_temp);

            % Corresponding variables store
            H_store{j} = Hnew;
            D1_store{j} = D_temp;
            theta1_store{j} = theta_temp;
            J1_store(j) = J_temp;
            Sigma_store{j} = Sigma;
        end

        % Choose min J to update
        min_idx = find(J1_store == min(J1_store));

        % Update all parameters
        theta_k = theta1_store{min_idx};
        H  = H_store{min_idx};  % For optimization I could save indices here
        J = J1_store(min_idx);
        Sigma = Sigma_store{min_idx};
        Dk = Sigma/var_y;
        k = length(theta_k);
    else

        % COMPUTE time estimate    t ---> t+
        Hk = H(1:t, 1:k);
        [theta_k, Sigma, J] = time_update(y, Hk, t, theta_k, var_y, Dk, J);
        Dk = Sigma/var_y;
    end

end





[~, idx] = ismember(H_true(1,:), H(1,:))

theta'
theta_k(idx)'


