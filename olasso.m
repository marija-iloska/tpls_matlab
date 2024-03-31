function [theta_olin, idx_olin, J, plot_stats] = olasso(y, H, t0, epsilon, var_y, idx_h)

% Dimensions
T = length(y);
dy = length(H(1,:));

% Define initial batch
y0 = y(1:t0);
H0 = H(1:t0, :);

% Define initial batch terms
xy0 = H0'*y0;
xx0 = H0'*H0;

% EIG
a = eig(xx0);
step = 0.001*t0/max(real(a));

% Initial estimate
[B, STATS] = lasso(H0, y0, 'CV', 5);
theta_olin = B(:, STATS.IndexMinMSE);

% Initialize terms
xy = zeros(dy,1);
xx = zeros(dy,dy);

% theta at t0
e = [];
e_init = e;
J = [];

% For plotting
correct = 0;
incorrect = 0;
theta_store = [];

for t = t0+1:T

    % Pred Error
    [J(end+1), ~] = pred_error_lasso(y, H, t, t0, var_y, theta_olin, e_init);

    % Updates
    xx = xx + H(t,:)'*H(t,:);
    xy = xy + H(t,:)'*y(t);    
    [theta_olin, ~] = olin_lasso(xy0, xx0, xy, xx, theta_olin, epsilon, step, t0, t, dy);

    % Evaluate model
    idx_olin = find(theta_olin ~= 0)';
    correct(end+1) = sum(ismember(idx_olin, idx_h));
    incorrect(end+1) = length(idx_olin) - correct(end);

    %theta_store = [theta_store; theta_olasso'];
end

% Concatenate results
plot_stats = {correct, incorrect};


end