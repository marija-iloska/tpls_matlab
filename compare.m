
%% Main Code
clear all
close all
clc

% Settings
var_y = 1;   % Variance
ps = 9;     % Sparsity percent
dy = 15;      % System dimension
r =  1;       % Range of input data H
rt = 0.5;      % Range of theta
T = 500;
p = dy - ps;

% OLASSO params 
epsilon = 1e-7;
t0 = dy+1;

% JPLS params
Tb = t0;
init = t0;

% rjMCMC params
n = round(0.2*T);
Ns = 500;
Nb = 1;

% Parallel runs
R = 1;

% Initialize arrays
% time_mcmc = zeros(R);
% time_jpls = zeros(R);
% time_olin = zeros(R);
mcmc_run = zeros(R);
jpls_run = zeros(R);
olin_run = zeros(R);



tic
for run = 1:R

    %Create data
    [y, H, theta] = generate_data(T, dy, r, rt,  ps, var_y);
    idx_h = find(theta ~= 0)';


    % Pad original true indices for comparison
    idx_h_padded = [idx_h zeros(1, dy - length(idx_h))];


    % PJ ORLS___________________________________________________
    tic
    [theta_jpls, H_jpls, model_stats,  error_stats, plot_stats] = jpls(y, H, dy, var_y, init, Tb, idx_h);
    [jpls_missing, jpls_correct, jpls_wrong] = plot_stats{:};
    [models_jpls, count_jpls, idx_jpls, idx_store] = model_stats{:};
    [J_pred, e] = error_stats{:};
    toc
    Jpred_jpls(run,:) = J_pred;
    e_jpls(run,:) = e;


    % Olin LASSO___________________________________________________
    tic
    [theta_olin, idx_olin, models_olin, count_olin, e, J_pred, olin_correct, olin_wrong, olin_missing] = olasso(y, H, t0, epsilon, var_y, idx_h);
    toc
    Jpred_olin(run,:) = J_pred;
    e_olin(run,:) = e;



    % RJ MCMC ___________________________________________________
    % Data partition and Number of sweeps
%     [idx_mcmc, theta_RJ, models_mcmc, count_mcmc, Nm, mcmc_stats, ~] = rj_mcmc(y, H, n, Ns, Nb, idx_h, var_y);
%     [mcmc_missing, mcmc_correct, mcmc_wrong] =mcmc_stats{:};
%     [J_mcmc, ~] = true_PE(y, H, t0, T, idx_mcmc, var_y);

    % GENIE 
    [J_true(run,:), e_true(run,:)] = true_PE(y, H, t0, T, idx_h, var_y);

    % SUPER GENIE
    e_super(run,:) = y(t0+1:end) - H(t0+1:end,:)*theta;
    J_super(run,:) = cumsum(e_super(run,:).^2);


    % SINGLE EXPECTATIONS
%     [Es_add, Es_rmv, Eb_add, Eb_rmv] = expectations(y, H, t0, T, idx_h, var_y, theta);


    % BARS
    jpls_f(run, :, :) = [jpls_correct;  jpls_wrong]; % jpls_missing]; 
    olin_f(run, :, :) = [olin_correct;  olin_wrong]; % olin_missing]; 
%     mcmc_f(run, :, :) = [mcmc_correct;  mcmc_wrong]; % mcmc_missing]; 



end
toc 

jpls_features = squeeze(mean(jpls_f,1));
olin_features = squeeze(mean(olin_f,1));
% mcmc_features = squeeze(mean(mcmc_f,1));


e_olin = mean(e_olin, 1);
e_jpls = mean(e_jpls, 1);
e_true = mean(e_true, 1);
J_jpls = mean(Jpred_jpls, 1);
J_olin = mean(Jpred_olin, 1);
J_true = mean(J_true,1);
J_super = mean(J_super,1);


% Anything below 5
jpls_run(jpls_run > 4) = 5;
olin_run(olin_run > 4) = 5;
%mcmc_run(mcmc_run > 4) = 5;

% For Labels
str_dy = num2str(dy);
str_k = num2str(dy - ps);
str_T = num2str(T);
str_v = num2str(var_y);
str_R = num2str(R);

% filename = join(['Results/T', str_T, '_K', str_dy, '_k', str_k, '_v', str_v, ...
%     '_R', str_R, '.mat']);
% 
% save(filename)

%% EXPECTATIONS

fsz = 20;
lwd = 3;

time_plot = t0+1:T;

% PLOT EXPECTATIONS
figure;
subplot(2,1,1)
for j = 1:dy-p
    plot(time_plot,Es_add(t0+1:end,j), 'Linewidth',2)
    hold on
end
yline(0, 'k', 'linewidth',1)
set(gca, 'FontSize', 15)
title('Signle Instant', 'FontSize', 15)
ylabel('MSE diff ', 'FontSize', fsz)
xlabel('Time', 'FontSize', 15)
xlim([t0+1, T])


subplot(2,1,2)
for j = 1:dy-p
    plot(time_plot, cumsum(Es_add(t0+1:end,j)), 'Linewidth',2)
    hold on
end
yline(0, 'k', 'linewidth',1)
set(gca, 'FontSize', 15)
title('In Batch', 'FontSize', 15)
ylabel('MSE diff ', 'FontSize', fsz)
xlabel('Time', 'FontSize', 15)
xlim([t0+1, T])

sgtitle('Extra feature:  E_{+j,n} - E_{p,n}', 'fontsize', 15)


figure;
subplot(2,1,1)
for j = 1:p
    plot(time_plot,Es_rmv(t0+1:end,j), 'Linewidth',2)
    hold on
end
yline(0, 'k', 'linewidth',1)
set(gca, 'FontSize', 15)
title('Single Instant', 'FontSize', 15)
ylabel('MSE diff', 'FontSize', fsz)
xlabel('Time', 'FontSize', 15)
xlim([t0+1, T])


subplot(2,1,2)
for j = 1:p
    plot(time_plot, cumsum(Es_rmv(t0+1:end,j)), 'Linewidth',2)
    hold on
end
yline(0, 'k', 'linewidth',1)
set(gca, 'FontSize', 15)
title('In Batch ','FontSize', 15)
ylabel('MSE diff', 'FontSize', fsz)
xlabel('Time', 'FontSize', 15)
xlim([t0+1, T])

sgtitle('Removed feature:  E_{-j,n} - E_{p,n}', 'fontsize', 15)




%% GENIE SUPER GENIE TEST


c_true = [248, 257, 133]/256;
c_olin = [2, 209, 123]/256;
c_olin = [0, 0.8, 0];
c_jpls = [252, 10, 67]/256;
c_jpls = [0.8, 0, 0];
c_mcmc = [43, 115, 224]/256;
c_true = [0,0,0];
c_inc = [0.4, 0.4, 0.4];
% 
% figure;
% 
% subplot(1,2,1)
% plot(J_jpls, 'Color', c_jpls, 'LineWidth', lwd)
% hold on
% plot(J_olin, 'Color', c_olin, 'LineWidth', lwd)
% hold on
% plot(J_mcmc(end)*ones(1,length(J_jpls)), 'Color', c_mcmc, 'LineWidth', lwd)
% hold on
% plot(J_true, 'Color', c_true, 'LineWidth', lwd)
% hold on
% plot(J_super, 'Color', [0, 0, 0], 'LineWidth', lwd)
% hold on
% set(gca, 'FontSize', 15)
% legend('J_{JPLS}', 'J_{OLinLASSO}', 'J_{RJMCMC}', 'J_{TRUE FEATURES}', 'J_{GROUND TRUTH}', 'FontSize', 15)
% ylabel('Predictive Error ', 'FontSize', fsz)
% xlabel('Time', 'FontSize', 15)


% % Create figure name
% filename = join(['figs_spec/genie_T', str_T, '_tr',num2str(rt), '_K', str_dy, '_k', str_k, '_v', str_v, ...
%     '_R', str_R, '.eps']);
% 
% % Save figure
% print(gcf, filename, '-depsc2', '-r300');





% FEATURES
fsz = 12;
fszl = 10;
lwd = 2;
lwdt = 4;
figure;

% JPLS features BAR plots
subplot(3,2,1)
jb = bar(t0:T, jpls_features', 'stacked', 'FaceColor', 'flat', 'FaceAlpha', 1);
jb(1).CData = c_jpls;
jb(2).CData = c_inc;
% jb(3).CData = [0.6, 0.6, 0.6];
hold on
yline(dy-ps, 'Color', c_true, 'LineWidth', lwdt)
ylim([0, dy])
xlim([t0+1, T])
set(gca, 'FontSize', 15)
legend('Correct', 'Incorrect', 'True Dim', 'FontSize', fszl)
title('JPLS', 'FontSize', 15)
ylabel('Number of Features ', 'FontSize', fsz)
xlabel('Time', 'FontSize', 15)


% OLinLASSO features BAR plots
%subplot(1,3,2)
subplot(3,2,3)
ob = bar(t0:T, olin_features', 'stacked', 'FaceColor', 'flat', 'FaceAlpha', 1);
hold on
ob(1).CData = c_olin;
ob(2).CData =  c_inc;
% ob(3).CData = [0.6, 0.6, 0.6];
yline(dy-ps, 'Color', c_true, 'LineWidth', lwdt)
ylim([0, dy])
xlim([t0+1, T])
set(gca, 'FontSize', 15)
legend('Correct', 'Incorrect', 'True Dim', 'FontSize', fszl)
title('OLinLASSO', 'FontSize', 15)
ylabel('Number of Features ', 'FontSize', fsz)
xlabel('Time', 'FontSize', 15)


% PREDICTIVE ERROR Plots
%subplot(1,3,3);
subplot(3,2,4)
time_plot = t0+1:T;
plot(time_plot, J_olin - J_true, 'Color', c_olin, 'LineWidth', lwd)
hold on
% plot(time_plot, J_mcmc - J_true, 'Color', c_mcmc, 'LineWidth', lwd)
hold on
plot(time_plot, J_jpls - J_true,  'Color', c_jpls, 'LineWidth', lwd)
hold on
yline(0, 'Color',c_true, 'linewidth', lwdt)
xlim([t0+1, T])
set(gca, 'FontSize', 15)
title('Relative', 'FontSize', 15)
legend('J_{OLin} - J_{FEATURES}', 'J_{RJMCMC} - J_{FEATURES}', 'J_{JPLS} - J_{FEATURES}', 'FontSize', fszl)
xlabel('Time', 'FontSize', fsz)
ylabel('Predictive Error Difference', 'FontSize', fsz)
grid on


subplot(3,2,2)
plot(time_plot,J_olin, 'Color', c_olin, 'LineWidth', lwd)
hold on
% plot(time_plot, J_mcmc, 'Color', c_mcmc, 'LineWidth', lwd)
hold on
plot(time_plot,J_jpls, 'Color', c_jpls, 'LineWidth', lwd)
hold on
plot(time_plot, J_true, 'Color', c_true, 'LineWidth', lwd)
hold on
plot(time_plot,J_super, 'Color', [0, 0, 0], 'LineWidth', lwd, 'LineStyle','--')
hold on
xlim([t0+1, T])
set(gca, 'FontSize', 15)
legend('J_{OLinLASSO}', 'J_{RJMCMC}', 'J_{JPLS}',  'J_{FEATURES}', 'J_{TRUTH}',  'FontSize', fszl)
title('Predictive Error', 'FontSize', 15)
ylabel('Predictive Error ', 'FontSize', fsz)
xlabel('Time', 'FontSize', fsz)
grid on


% subplot(3,1,3)
% rj = bar(mcmc_features', 'stacked', 'FaceColor', 'flat', 'FaceAlpha', 1);
% hold on
% rj(1).CData = c_mcmc;
% rj(2).CData =  c_inc;
% % rj(3).CData = [0.6, 0.6, 0.6];
% yline(dy-ps, 'Color', c_true, 'LineWidth', lwdt)
% ylim([0, dy])
% set(gca, 'FontSize', 15)
% legend('Correct', 'Incorrect', 'True Dim', 'FontSize', fszl)
% title('RJMCMC', 'FontSize', 15)
% ylabel('Number of Features ', 'FontSize', fsz)
% xlabel('Sweeps', 'FontSize', 15)

% Create title string
title_str = join(['\sigma^2 = ', str_v, ...
    '     {\bf h_k } ~ N( {\bf 0}, ', num2str(r), '{\bf I} )  ' , '  {\bf \theta } ~ N( {\bf 0} ', num2str(rt), '{\bf I} ) ']) ; 

% Make super title
%sgtitle(title_str, 'FontSize', 20)



%% SAVE FIGURE
% % Create figure name
% filename = join(['figs_spec/stack_T', str_T, '_tr',num2str(rt), '_K', str_dy, '_k', str_k, '_v', str_v, ...
%     '_R', str_R, '.eps']);
% 
% % Save figure
% print(gcf, filename, '-depsc2', '-r300');


