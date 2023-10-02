clear all
close all
clc

R = 1000;
var_y = 0.001;   % Variance
ps = 1;     % Sparsity percent
dy = 7;      % System dimension


str_dy = num2str(dy);
str_k = num2str(dy - ps);
str_v = num2str(var_y);
str_R = num2str(R);

all_T = [30, 60, 120, 240, 480];
len_T = length(all_T);

count = 0;

for t = 1:len_T

    % Get correct Data
    T = all_T(t);
    str_T = num2str(T);

    % COnstruct full filename to load
    filename = join(['Results/T', str_T, '_K', str_dy, '_k', str_k, '_v', str_v, ...
        '_R', str_R, '.mat']);
     str_title = join(['T', str_T, ' K', str_dy, ' k', str_k]);
     str_T_title = join(['T = ', str_T]);

    % Load data
    load(filename)

    % Variables of interest
    time_mcmc = time_mcmc(:,1);
    time_orls = time_orls(:,1);
    run_time_ratio(t) = mean(time_mcmc./time_orls);
    avg_time_mcmc(t) = mean(time_mcmc);
    avg_time_orls(t) = mean(time_orls);

    mcmc{t} = mcmc_run(:,1);
    orls{t} = orls_run(:,1);

    % Number of times chosen correctly 
    mcmc_first(t) = sum(mcmc{t} == 1);
    orls_first(t) = sum(orls{t} == 1);

    mcmc_second(t) = sum(mcmc{t} == 2);
    orls_second(t) = sum(orls{t} == 2);

    mcmc_third(t) = sum(mcmc{t} == 3);
    orls_third(t) = sum(orls{t} == 3);

%     mcmc_(t) = sum(mcmc{t} == 3);
%     orls_third(t) = sum(orls{t} == 3);


     count = count + 1;
%     figure(1);
%     subplot(len_T, 2, count)
%     histogram(orls{t}, 'FaceColor', [176, 123, 173]/256, 'FaceAlpha', 0.8, ...
%         'EdgeColor', [80, 0, 110]/256, 'EdgeAlpha', 1, 'LineWidth', 1.5)
%     orls_x = get(gca, 'xTick');
%     xticks(unique(round(orls_x)));
%     ylim([0,R])
%     xlim([0,5.5])
%     hold on
%     set(gca, 'FontSize', 15)
%     title(str_T_title, 'FontSize', 15)
%     ylabel('Runs', 'FontSize', 15)
%     xlabel('Rank of Correct Model', 'FontSize', 15)
%     legend('pjORLS', 'FontSize', 15)
%     grid on
% 
%     count  = count + 1;
% 
% 
%     subplot(len_T, 2, count)
%     histogram(mcmc{t}, 'FaceColor', [9, 173, 168]/256, 'FaceAlpha', 0.8, ...
%         'EdgeColor', [31, 61, 60]/256, 'EdgeAlpha', 1, 'LineWidth', 1.5)
%     mcmc_x = get(gca, 'xTick');
%     xticks(unique(round(mcmc_x)));
%     ylim([0,R])
%     xlim([0,5.5])
%     set(gca, 'FontSize', 15)
%     title(str_T_title, 'FontSize', 15)
%     ylabel('Runs', 'FontSize', 15)
%     xlabel('Rank of Correct Model', 'FontSize', 15)
%     legend('rjMCMC', 'FontSize', 15)
%     grid on





    figure(2)
    subplot(len_T, 1, t)
    histogram(orls{t}, 'FaceColor', [176, 123, 173]/256, 'FaceAlpha', 0.5, ...
        'EdgeColor', [80, 0, 110]/256, 'EdgeAlpha', 1)
    orls_x = unique(round(get(gca, 'xTick')));
    ylim([0,R])
    hold on
    histogram(mcmc{t}, 'FaceColor', [9, 173, 168]/256, 'FaceAlpha', 0.3, ...
        'EdgeColor', [31, 61, 60]/256, 'EdgeAlpha', 1)
    mcmc_x = unique(round(get(gca, 'xTick')));
    x = {orls_x, mcmc_x};
    len = [length(orls_x), length(mcmc_x)];
    idx = find(max(len) == len);
    set(gca, 'FontSize', 15)  
    xticks(x{idx(1)});
    title(str_title ,  'FontSize', 15)
    ylabel('Runs', 'FontSize', 15)
    xlabel('Rank of Correct Model', 'FontSize', 15)
    legend('pjORLS', 'rjMCMC', 'FontSize',15)
    grid on




end

% str_dim = join(['K =', str_dy, ', k = ', str_k]);
% sgtitle(str_dim, 'FontSize', 15)

% fig_name = join(['figs/', 'K', str_dy', '_k', str_k, '.eps' ]);
% print(gcf, fig_name, '-depsc2', '-r300');
% 

% load run_time.mat
% lwd = 3;
% figure;
% plot(all_T, run_time_ratioK5_k3, 'Color', 'k', 'Linewidth', lwd, 'LineStyle', '--')
% hold on
% plot(all_T, run_time_ratioK7_k3, 'Color', 'r', 'Linewidth', lwd, 'LineStyle', '--')
% hold on
% plot(all_T, run_time_ratioK7_k6, 'Color', 'b', 'Linewidth', lwd, 'LineStyle', '--')
% set(gca, 'FontSize', 15)
% xlabel('Number of Observations T', 'FontSize', 20)
% ylabel('t_{rjMCMC} / t_{pjORLS}', 'FontSize', 25)
% title('Average Run Time Ratio', 'FontSize', 15)
% legend('K5 k3', 'K7 k3', 'K7 k6', 'FontSize', 15)
% grid on
% 
% 
% 
% figure;
% plot(all_T, orls_first, 'Color', [176, 123, 173]/256, 'linewidth', lwd, 'LineStyle',':')
% hold on
% plot(all_T, orls_second, 'Color', [176, 123, 173]/256, 'linewidth', lwd, 'LineStyle','--')
% hold on
% plot(all_T, orls_third, 'Color', [176, 123, 173]/256, 'linewidth', lwd, 'LineStyle','-')
% hold on
% plot(all_T, mcmc_first, 'Color', [9, 173, 168]/256, 'linewidth', lwd, 'LineStyle',':')
% hold on
% plot(all_T, mcmc_second,'Color', [9, 173, 168]/256, 'linewidth', lwd, 'LineStyle','--')
% hold on
% plot(all_T, mcmc_third,'Color', [9, 173, 168]/256, 'linewidth', lwd, 'LineStyle','-')
% grid on
% ylim([0, R])
% set(gca, 'FontSize', 15)
% xlabel('Number of Observations T', 'FontSize', 20)
% ylabel('Number of Runs', 'FontSize', 20)
% legend('pjORLS first rank', 'pjORLS second rank', 'pjORLS third rank', 'rjMCMC first rank',  'rjMCMC second rank', 'rjMCMC third rank', 'FontSize',15)

% run_time_ratioK5_k3 = run_time_ratio;
% 
% 
% save('run_time.mat', 'run_time_ratioK7_k3', 'run_time_ratioK7_k6', 'run_time_ratioK5_k3')
% 
