clear all
close all
clc

% Paths to access functions from other folders
function_paths = [genpath('tpls/'), genpath('util/'), ...
    genpath('predictive_error/'), genpath('baselines/')];

% Add the paths
addpath(function_paths)

% Main Script ==========================================================

% Settings
var_y = 1;            % Observation noise Variance
ps = 5;                 % Number of 0s in theta
K = 12;                 % Number of available features
var_features =  1;      % Range of input data H
var_theta = 1;        % Variance of theta
T = 200;                 % Number of data points
p = K - ps;             % True model dimension

% SYNTHETIC DATA =================================================================
%Create data
[y, H, theta] = generate_data(T, K, var_features, var_theta,  ps, var_y);

% Indices of true features
idx_h = find(theta ~= 0)';

% Pad original true indices for comparison later
idx_h_padded = [idx_h zeros(1, K - length(idx_h))];


% OLASSO params
epsilon = 1e-7;

% Initial batch of datad
t0 = K+1;

% rjMCMC params
n = round(0.2*T);
Ns = 700;
Nb = 50;


%% CALL METHODS


% TPLS =================================================================
[theta_tpls, idx_tpls, J, plot_stats, idx_tpls_store] = tpls(y, H, K, var_y, t0, idx_h);

% Results for plotting
[tpls_correct, tpls_incorrect] = plot_stats{:};
J_tpls = J;


% OlinLASSO =================================================================
[theta_olin, idx_olin, J, plot_stats, idx_olin_store] = olasso(y, H, t0, epsilon, var_y, idx_h);

% Results for plotting
[olin_correct, olin_incorrect] = plot_stats{:};
J_olin = J;

% OCCD =================================================================
[theta_occd, idx_occd, J, plot_stats, idx_occd_store] = occd(y, H, t0, var_y, idx_h);

% Results for plotting
[occd_correct, occd_incorrect] = plot_stats{:};
J_occd = J;


% RJ MCMC =================================================================
% Data partition and Number of sweeps
[idx_mcmc, theta_RJ, plot_stats, J] = rj_mcmc(y, H, n, Ns, Nb, idx_h, var_y, t0);

% Results for plotting
[mcmc_correct, mcmc_incorrect] = plot_stats{:};
J_mcmc = J;




%% GROUND TRUTHS ===========================================================
% GENIE
[J_true, ~] = true_PE(y, H, t0, T, idx_h, var_y);

% SUPER GENIE
e_super = y(t0+1:end) - H(t0+1:end,:)*theta;
J_super = cumsum(e_super.^2);



% SINGLE EXPECTATIONS =====================================================
[E_add, E_rmv] = expectations(y, H, t0, T, idx_h, var_y, theta);


% BARS (for statistical performance)
tpls_features = [tpls_correct;  tpls_incorrect];
olin_features = [olin_correct;  olin_incorrect];
occd_features = [occd_correct;  occd_incorrect];
mcmc_features = [mcmc_correct;  mcmc_incorrect];




%% FIGURE 1: TPLS  vs  OLinLASSO  vs  RJMCMC  barplots
% Specific run with feature bar plots

% Colors, FontSizes, Linewidths
load plot_settings.mat
fszl= 10;

% Time range to plot
time_plot = t0+1:T;


% BAR PLOTS SPECIFIC RUN =========================================
figure('Renderer', 'painters', 'Position', [900 100 1000 900])

% TPLS
subplot(3,2,1)
formats = {fsz, fszl, fsz_title, lwdt, c_tpls, c_inc, c_true, 'TPLS', 'Time'};
bar_plots(tpls_features, t0+1, T, p, K, formats)

% OLinLASSO
subplot(3,2,2)
formats = {fsz, fszl, fsz_title, lwdt, c_olin, c_inc, c_true, 'OLinLASSO', 'Time'};
bar_plots(olin_features, t0+1, T, p, K, formats)

% OCCD
subplot(3,2,3)
formats = {fsz, fszl, fsz_title, lwdt, c_occd, c_inc, c_true, 'OCCD-TWL', 'Time'};
bar_plots(occd_features, t0+1, T, p, K, formats)

% RJMCMC
subplot(3,2,4)
formats = {fsz, fszl, fsz_title, lwdt, c_mcmc, c_inc, c_true, 'RJMCMC', 'Iteration'};
bar_plots(mcmc_features, 1, Ns, p, K, formats)


% PREDICTIVE ERRORs ===============================================

% Difference
subplot(3,2,6)
hold on
plot(time_plot, J_occd - J_true, 'Color', c_occd, 'LineWidth', lwd)
plot(time_plot, J_olin - J_true, 'Color', c_olin, 'LineWidth', lwd)
plot(time_plot, J_mcmc - J_true, 'Color', c_mcmc, 'LineWidth', lwd)
plot(time_plot, J_tpls - J_true, 'Color', c_tpls, 'LineWidth', lwd)
yline(0, 'Color',c_true, 'linewidth', 2)
hold off
xlim([t0+1, T])
set(gca, 'FontSize', 15)
title('Relative', 'FontSize', fsz_title)
legend('\Delta J_{OCCD}','\Delta J_{OLin}', '\Delta J_{RJMCMC}', '\Delta J_{TPLS}', 'FontSize', fszl-2)
xlabel('Time', 'FontSize', fsz)
ylabel('Predictive Error Difference', 'FontSize', fsz)
grid on
box on

% Raw
subplot(3,2,5)
hold on
plot(time_plot, J_occd, 'Color', c_occd, 'LineWidth', lwd)
plot(time_plot, J_olin,  'Color', c_olin, 'LineWidth', lwd)
plot(time_plot, J_mcmc,  'Color', c_mcmc, 'LineWidth', lwd)
plot(time_plot, J_tpls,  'Color', c_tpls, 'LineWidth', lwd)
plot(time_plot, J_true,  'Color', c_true, 'LineWidth', lwd)
plot(time_plot, J_super, 'Color', c_true, 'LineWidth', lwd, 'LineStyle','--')
hold off
xlim([t0+1, T])
set(gca, 'FontSize', 15)
legend('J_{OCCD}','J_{OLin}', 'J_{RJMCMC}', 'J_{TPLS}',  'J_{GENIE}', 'J_{TRUTH}',  'FontSize', fszl-2)
title('Predictive Error', 'FontSize', fsz_title)
ylabel('Predictive Error ', 'FontSize', fsz)
xlabel('Time', 'FontSize', fsz)
grid on
box on
 
lwd = lwd+0.5;

%% FIGURES 4 and 5:  PLOT EXPECTATIONS

% Import colors
load colors.mat
col{7} = [0,0.5,0];
title_str = 'INSTANT';
y_str = '\Delta_n';

time_plot = time_plot(1:100);
% ADD A FEATURE =======================================================
figure('Renderer', 'painters', 'Position', [0 100 450 900])
subplot(4,1,1)
formats = {fsz, lwd, col, 'Extra Feature', '\Delta^{+m}_n', ''};
expectation_plots(E_add(time_plot,:), time_plot, K-p,  formats)

subplot(4,1,2)
formats = {fsz, lwd, col, '', '\Sigma \Delta^{+m}_n', ''};
expectation_plots(cumsum(E_add(time_plot,:)), time_plot, K-p,  formats)

% REMOVING A FEATURE ==================================================
%col{7} = [0,0,0];
subplot(4,1,3)
formats = {fsz, lwd, col, 'Removed Feature', '\Delta^{-m}_n', ''};
expectation_plots(E_rmv(time_plot,:), time_plot, p,  formats)

subplot(4,1,4)
formats = {fsz, lwd, col, '', '\Sigma \Delta^{-m}_n', 'n'};
expectation_plots(cumsum(E_rmv(time_plot,:)), time_plot, p,  formats)


