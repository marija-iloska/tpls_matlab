clear all
close all
clc

% Settings
var_y = 0.001;   % Variance
ps = 2;     % Sparsity percent
dy = 5;      % System dimension
r = 1;       % Range of input data H
rt = 2;      % Range of theta
T = 240;
n = round(0.3*T);
Ns = 3000;
Nb = 1300;
Tb = 40;

R = 1000;

parpool(25)

time_mcmc = zeros(R);
time_orls = zeros(R);
orls_run = zeros(R);
mcmc_run = zeros(R);


tic
parfor run = 1:R

    %Create data
    [y, H, theta] = generate_data(T, dy, r, rt,  ps, var_y);
    idx_h = find(theta ~= 0)';

    % Pad original true indices for comparison
    idx_h_padded = [idx_h zeros(1, dy - length(idx_h))];


    % RJ MCMC ___________________________________________________
    % Data partition and Number of sweeps
    tic
    [idx_mcmc, theta_RJ, models_mcmc, count_mcmc, Nm] = rj_mcmc(y, H, n, Ns, Nb);
    time_mcmc(run) = toc;


    % Check through all models
    idx_corr_mcmc = 0;
    for m = 1:length(models_mcmc(:,1))
        if (sum(models_mcmc(m,:) == idx_h_padded ) == dy)
            idx_corr_mcmc = m;
        end
    end


    % PJ ORLS___________________________________________________
    tic
    init = dy + 1;
    [theta_k, Hk, k_store, k_mode, models_orls, count_orls, idx_orls] = pj_orls(y, H, dy, var_y, init, Tb);
    time_orls(run) = toc;



    % Check through all models
    idx_corr_orls = 0;
    for m = 1:length(models_orls(:,1))
        if (sum(models_orls(m,:) == idx_h_padded ) == dy)
            idx_corr_orls = m;
        end
    end

    orls_run(run) = idx_corr_orls;
    mcmc_run(run) = idx_corr_mcmc;
end
toc

% Anything below 5
orls_run(orls_run > 4) = 5;
mcmc_run(mcmc_run > 4) = 5;


% Average run time ratio
avg_time = mean(time_mcmc./time_orls);


% 
% figure;
% histogram(orls_run, 'FaceColor', [176, 123, 173]/256, 'FaceAlpha', 0.5, ...
%     'EdgeColor', [80, 0, 110]/256, 'EdgeAlpha', 1)
% orls_x = get(gca, 'xTick');
% xticks(unique(round(orls_x)));
% ylim([0,R])
% hold on
% histogram(mcmc_run, 'FaceColor', [9, 173, 168]/256, 'FaceAlpha', 0.3, ...
%     'EdgeColor', [31, 61, 60]/256, 'EdgeAlpha', 1)
% set(gca, 'FontSize', 15)
% %title('pjORLS', 'FontSize', 15)
% ylabel('Percentage', 'FontSize', 15)
% xlabel('Rank of Correct Model', 'FontSize', 15)
% legend('pjORLS', 'rjMCMC', 'FontSize',15)
% grid on
% 
% 
% figure;
% subplot(2,1,1)
% histogram(orls_run, 'FaceColor', [176, 123, 173]/256, 'FaceAlpha', 0.8, ...
%     'EdgeColor', [80, 0, 110]/256, 'EdgeAlpha', 1, 'LineWidth', 1.5)
% orls_x = get(gca, 'xTick');
% xticks(unique(round(orls_x)));
% ylim([0,R])
% xlim([0,5.5])
% hold on
% set(gca, 'FontSize', 15)
% title('pjORLS', 'FontSize', 15)
% ylabel('Percentage', 'FontSize', 15)
% xlabel('Rank of Correct Model', 'FontSize', 15)
% grid on
% 
% 
% subplot(2,1,2)
% histogram(mcmc_run, 'FaceColor', [9, 173, 168]/256, 'FaceAlpha', 0.8, ...
%     'EdgeColor', [31, 61, 60]/256, 'EdgeAlpha', 1, 'LineWidth', 1.5)
% mcmc_x = get(gca, 'xTick');
% xticks(unique(round(mcmc_x)));
% ylim([0,R])
% xlim([0,5.5])
% set(gca, 'FontSize', 15)
% title('rjMCMC', 'FontSize', 15)
% ylabel('Percentage', 'FontSize', 15)
% xlabel('Rank of Correct Model', 'FontSize', 15)
% grid on
% 
% 
% 
% 
% figure;
% histogram(mcmc_run, 'FaceColor', [9, 173, 168]/256, 'FaceAlpha', 0.5, ...
%     'EdgeColor', [31, 61, 60]/256, 'EdgeAlpha', 1, 'LineWidth', 3)
% mcmc_x = get(gca, 'xTick');
% xticks(unique(round(mcmc_x)));
% ylim([0,R])
% set(gca, 'FontSize', 15)
% title('rjMCMC', 'FontSize', 15)
% ylabel('Percentage', 'FontSize', 15)
% xlabel('Rank of Correct Model', 'FontSize', 15)
% grid on
% 
% figure;
% histogram(orls_run, 'FaceColor', [176, 123, 173]/256, 'FaceAlpha', 0.5, ...
%     'EdgeColor', [80, 0, 110]/256, 'EdgeAlpha', 1, 'LineWidth', 1)
% orls_x = get(gca, 'xTick');
% xticks(unique(round(orls_x)));
% ylim([0,R])
% set(gca, 'FontSize', 15)
% title('pjORLS', 'FontSize', 15)
% ylabel('Percentage', 'FontSize', 15)
% xlabel('Rank of Correct Model', 'FontSize', 15)
% grid on
R = length(orls_run);

str_dy = num2str(dy);
str_k = num2str(dy - ps);
str_T = num2str(T);
str_v = num2str(var_y);
str_R = num2str(R);

filename = join(['Results/T', str_T, '_K', str_dy, '_k', str_k, '_v', str_v, ...
    '_R', str_R, '.mat']);

save(filename)



% Bar plot
% figure;
% b_orls = bar(count_orls/(T-Tb), 'FaceColor', 'flat');
% ylabel('Number of Visits')
% title('ORLS Models visited ','FontSize',20)
% set(gca, 'FontSize', 20);
% grid on
% if (idx_corr_orls==0)
%     text(1,0.1, 'True Model NOT visited', 'FontSize', 15)
% else
%     b_orls.CData(idx_corr_orls,:) = [0.5, 0, 0];
% end


% Bar plot
% figure;
% b_mcmc = bar(count_mcmc/Ns, 'FaceColor', 'flat');
% ylabel('Number of Visits')
% title('RJMCMC Models visited','FontSize',20)
% set(gca, 'FontSize', 20);
% grid on
% b_mcmc.CData(idx_corr_mcmc,:) = [0, 0, 0];




% filename = 'TestORLS.eps'; % join(['figs23/pjorls', num2str(run), '.eps']);
% print(gcf, filename, '-depsc2', '-r300');


% filename = 'TestMCMC.eps';  % join(['figs23/rjmcmc', num2str(run), '.eps']);
% print(gcf, filename, '-depsc2', '-r300');

