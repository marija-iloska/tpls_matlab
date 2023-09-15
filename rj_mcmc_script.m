clear all
close all
clc

% Create data
var_y = 1; % Variance
p_s = 0.6;    % Sparsity percent
K = 5;      % System dimension
N = 60;     % Time series length
r = 1;     % Range of input data H
rt = 1;      % Range of theta

[y, H, theta_true] = generate_data(N, K, r,rt,  p_s, var_y);
idx_h = find(theta_true ~= 0)';



% Settings for Algorithm

% Data partition
n = 20;

% Number of sweeps
Ns = 1020;
[M_max, theta_RJ, models_sorted, count_sorted, Nm] = rj_mcmc(y, H, n, Ns);

% LS estimator for reference
theta_LS = inv(H'*H)*H'*y;

% Legend Label
for m = 1:Nm
    str{m} = join(['M_{', num2str(m), '}']);
end




% Bar plot
figure(1)
b = bar(count_sorted/Ns, 'FaceColor', 'flat');
ylabel('Number of Visits')
title('Models visited','FontSize',20)
set(gca, 'FontSize', 20); %, 'XTicklabel', str)
grid on



