clear all
close all
clc


% example_code.m
% This script is an example code on how to run TPLS.
% The code is easy to use and follows the format:

% TPLS_________________________________________
% Initialize --> model(t0)
% for t = t0+1, ... T
%   receive new data D(t) = { y(t), H(t,:) }
%   results(t) = model_update( model(t-1), D(t) )
%   model(t) = time_update( results(t), D(t) )
% end


% Some Important Variables -----------------------------------------

% y, H  - Output data, Feature matrix
% theta - True model parameter
% theta_est - the last parameter estimate
% idx_est - set of the indices of the selected features
% J - the predictive error of the last estimate
% correct - how many features selected were in the true model
% incorrect - how many features selected were not in the true model



% Paths to access functions from other folders
function_paths = [genpath('tpls/'), genpath('util/'), genpath('baselines/')];

% Add the paths
addpath(function_paths)


%% GENERATE SYNTHETIC DATA
% Settings
var_y = 0.1;            % Observation noise Variance
ps = 6;                 % Number of 0s in theta
K = 12;                 % Number of available features
var_features =  1;      % Range of input data H
var_theta = 0.5;        % Variance of theta
T = 300;                 % Number of data points
p = K - ps;             % True model dimension

% Initial batch of data
t0 = K+1;

%Create data
[y, H, theta] = generate_data(T, K, var_features, var_theta,  ps, var_y);
idx_h = find(theta ~= 0)';

%% INITIALIZE

% Get estimate using half the total features
k = floor(K/2);
[~, ~, theta_k, Dk, ~,~] = initialize(y, H, t0, k, var_y);

% Initialize variables
J_pred = [];
J = 0;

% Model storage
correct = zeros(1,T-t0);
incorrect = zeros(1,T-t0);

% Set of all feature indices
idx_H = 1:K;
idx_tpls_all = idx_H;


%% JPLS LOOP
% Start time loop
for t = t0+1:T

    % Data
    Ht = H(1:t, idx_tpls_all);
    yt = y(1:t);

    % TPLS one instant
    [theta_k, Dk, k, hk, J] = model_update(yt, Ht, theta_k, Dk, J, K, var_y, t0, t, idx_H);


    % STORE PREDICTIVE ERROR 
    J_pred(end+1) = J;

     % Check which model was selected to update feature order
    [~, idx_tpls_all] = ismember(hk, H(2,:));

    % Features used
    idx_est = idx_tpls_all(1:k);


    % TIME UPDATE    
    % theta(k,t) <-- theta(k,t-1) and Dk(t) <-- Dk(t-1)
    [theta_k, Dk] = time_update(y(t), H(t, idx_est), theta_k, var_y, Dk);



    % EVALUATION
    correct(t-t0) = sum(ismember(idx_est, idx_h));
    incorrect(t-t0) = length(idx_est) - correct(t-t0);


end

% Concatenate results
plot_stats = {correct, incorrect};


% Results for barplots
[correct, incorrect] = plot_stats{:};
tpls_features = [correct;  incorrect];

% Estimate of theta (final model (or it can be most visited model) )
%theta_est = theta_k;



%% PLOTTING

% Colors, FontSize, Linewidths
load plot_settings.mat

% Time range to plot
time_plot = t0+1:T;

% Create figure 
figure('Renderer', 'painters', 'Position', [200 300 1000 400])

% Features Bar plot
subplot(1,2,1)
formats = {fsz, fszl, lwdt, c_tpls, c_inc, c_true, ''};
bar_plots(tpls_features, t0+1, T, p, K, formats)

% Predictive Error plot
subplot(1,2,2)
plot(time_plot, J_pred, 'Color', c_tpls, 'LineWidth', lwd)
xlim([t0+1, T])
set(gca, 'FontSize', 15)
legend('J_{TPLS}', 'FontSize', fszl)
xlabel('Time', 'FontSize', fsz)
ylabel('Predictive Error', 'FontSize', fsz)
grid on
sgtitle('\bf{TPLS}', 'FontSize',fsz)


