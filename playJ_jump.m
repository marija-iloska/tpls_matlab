clc 
close all
clear ll

% Write a code that computes jumping up, staying same, jumping down Jkt

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

% Store 
H_true = H;

% Initialize using t data points
t = 2;
[J, theta_k, Dk, Hk, Sigma] = initialize(y, H, t, var_y);

% Start time loop
for t = 2:T

    % JUMP UP +
    % Model compute

    % Time compute


    % JUMP DOWN -
    % Model compute

    % Time compute

    % ELSE 
    % Time compute

end









