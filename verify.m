
%% Main Code
clear all
close all
clc

% Settings
var_y = 1;   % Variance
ps = 5;     % Sparsity percent
K = 10;      % System dimension
r =  3;       % Range of input data H
rt = 0.5;      % Range of theta
T = 5000;
p = K - ps;
t0 = K+1;

% JPLS params
Tb = 5;


% Parallel runs
R = 100;

% Initialize arrays
idx0 = datasample(1:K, ps, 'replace', false);

%Create data
[y_clean, H, theta] = fix_data(T, K, r, rt, idx0);
idx_h = find(theta ~= 0)';

fsz = 20;
lwd = 3;

tic
parfor run = 1:R

    % Add noise to data
    noise = normrnd(0, var_y^0.5, T,1);
    y = y_clean + noise;


    % Pad original true indices for comparison
    idx_h_padded = [idx_h zeros(1, K - length(idx_h))];


%     % PJ ORLS___________________________________________________
%     tic
%     [theta_jpls, H_jpls, ~,  error_stats, plot_stats] = jpls(y, H, K, var_y, t0, Tb, idx_h);
%     [jpls_missing, jpls_correct, jpls_wrong] = plot_stats{:};
%     [J_pred, e] = error_stats{:};
%     toc


    % MSE SUPER GENIE
    mse_super(run, :) = cumsum(noise(t0+1:T).^2);

    % BARS
    %jpls_f(run, :, :) = [jpls_correct;  jpls_wrong; jpls_missing]; 

%     figure(1)
%     subplot(2,4, run)
%     jb = bar(t0:T, squeeze(jpls_f(run, :,:)), 'stacked', 'FaceColor', 'flat', 'FaceAlpha', 1);
%     jb(1).CData = [0.7, 0, 0];
%     jb(2).CData = [0,0,0];
%     jb(3).CData = [0.6, 0.6, 0.6];
%     hold on
%     yline(K-ps, 'Color', 'b', 'LineWidth', 5)
%     ylim([0, K])
%     set(gca, 'FontSize', 15)
%     legend('Correct', 'Incorrect', 'Missing', 'True Order', 'FontSize', 10)
%     title('JPLS', 'FontSize', 15)
%     ylabel('Number of Features ', 'FontSize', fsz)
%     xlabel('Time', 'FontSize', 15)


end
toc 


% GENIE 
[genie] = error_update(H, t0, T, var_y, idx_h);

mse_genie = cumsum(genie);


figure;
plot(mse_genie,  'k','Linewidth',3)
hold on
for r = 1:R
    plot(mse_super(r,:), 'Linewidth',1)
    hold on
end
plot(mse_genie,  'y','Linewidth',3)
hold on
set(gca, 'FontSize', 20)
ylabel('Cumulative Noise Squared ', 'FontSize', 15)
xlabel('Time', 'FontSize', 15)
legend('MSE Genie', '\Sigma E(\eta_t)^2', 'FontSize', 20, 'location', 'northwest')
title('Same System - Noise Realizations', 'FontSize', 15)
grid on






