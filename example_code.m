clear all
close all
clc

% example_code.m
% This script is an example code on how to run JPLS.
% The code is easy to use and follows the format:
% results = jpls(arguments)


% Some Important Variables -----------------------------------------

% y, H  - Output data, Feature matrix
% theta - True model parameter
% theta_est - the last parameter estimate
% idx_est - set of the indices of the selected features
% J - the predictive error of the last estimate
% correct - how many features selected were in the true model
% incorrect - how many features selected were not in the true model



% Paths to access functions from other folders
function_paths = [genpath('jpls/'), genpath('util/'), ...
    genpath('predictive_error/'), genpath('baselines/')];

% Add the paths
addpath(function_paths)

%% Main Script

% Settings
var_y = 1;            % Observation noise Variance
ps = 4;                 % Number of 0s in theta
K = 10;                 % Number of available features
var_features =  1;      % Range of input data H
var_theta = 0.5;        % Variance of theta
T = 150;                 % Number of data points
p = K - ps;             % True model dimension

% Initial batch of data
t0 = K+1;

%Create data
[y, H, theta] = generate_data(T, K, var_features, var_theta,  ps, var_y);
idx_h = find(theta ~= 0)';


% JPLS =================================================================
[theta_est, idx_est, J, plot_stats] = jpls(y, H, K, var_y, t0, idx_h);


% Results for barplots
[correct, incorrect] = plot_stats{:};
jpls_features = [correct;  incorrect];



%% PLOTTING

% Colors, FontSize, Linewidths
load plot_settings.mat

% Time range to plot
time_plot = t0+1:T;

% Create figure 
figure('Renderer', 'painters', 'Position', [200 300 1000 400])

% Features Bar plot
subplot(1,2,1)
formats = {fsz, fszl, lwdt, c_jpls, c_inc, c_true, ''};
bar_plots(jpls_features, t0+1, T, p, K, formats)

% Predictive Error plot
subplot(1,2,2)
plot(time_plot, J, 'Color', c_jpls, 'LineWidth', lwd)
xlim([t0+1, T])
set(gca, 'FontSize', 15)
legend('J_{JPLS}', 'FontSize', fszl)
xlabel('Time', 'FontSize', fsz)
ylabel('Predictive Error', 'FontSize', fsz)
grid on
sgtitle('\bf{JPLS}', 'FontSize',fsz)



