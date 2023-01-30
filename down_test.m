% OLS reverse test script 
clear all
close all
clc

% Order - recursive least squares example

% Settings
var_y = 1;
g = @(x) x;
p_s = 0.3;
dx = 10;
T = 100;
r = 2; % Range of input data H
rt = 5;  % Range of theta

%SSM
tr = @(coeff, states) coeff*g(states);

%Create data
[y, H, theta, a] = generate_data(T, dx, r, rt, p_s, var_y, tr, g);



% Call ORLS
epsilon = 0.01;
t = T-10;
[theta_k, Dk, Jk, error_store] = ols(y(1:t), H(1:t, :), epsilon, dx);

k = length(theta_k);

t = T-9;

for j = 1:k
    % Update current theta by remoivng jth basis function
    [theta_temp, D_temp, H_temp] = ols_downdates(theta_k, Dk, j, H, t);

    % Compute that probability
    log_pk0(j) = log_mvnpdf(y(1:t), H_temp, theta_temp, var_y);

    % Corresponding variables store
    H0_store{j} = H_temp;
    D0_store{j} = D_temp;
    theta0_store{j} = theta_temp;
end

% Compute probability of staying Mk
log_pk = log_mvnpdf(y(1:t), H(1:t, 1:k), theta_k, var_y);

% Normalize Probabilities
log_probs = [log_pk, log_pk0];
probs = exp(log_probs - max(log_probs));
probs = probs./sum(probs);


% Sample model
ki = datasample(1:length(probs), 1, 'Weights', probs);

m = ki - 1;

% If REMOVE COLUMN
if (ki > 1)
    % Update original H (swap column order)
    H = H(:, [setdiff(1:dx, m), m]);

    % Update rest
    Hk = H0_store{m};
    Dk = D0_store{m};
    theta_k = theta0_store{m};

end
