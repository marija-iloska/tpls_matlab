clear all
close all
clc

% Settings
var_y = 0.01;   % Variance
ps = 3;     % Sparsity percent
dy = 6;      % System dimension
r = 2;       % Range of input data H
rt = 0.5;      % Range of theta
T = 300;

% OLASSO params
epsilon = 1e-7;
t0 = 50;

% JPLS params
Tb = 50;

% rjMCMC params
n = round(0.2*T);
Ns = 2000;
Nb = 1000;

% Parallel runs
R = 1;

% Initialize arrays
time_mcmc = zeros(R);
time_orls = zeros(R);
time_olasso = zeros(R);
orls_run = zeros(R);
mcmc_run = zeros(R);
olin_run = zeros(R);


for run = 1:R

    %Create data
    [y, H, theta] = generate_data(T, dy, r, rt,  ps, var_y);
    idx_h = find(theta ~= 0)';

    % Pad original true indices for comparison
    idx_h_padded = [idx_h zeros(1, dy - length(idx_h))];


    % RJ MCMC ___________________________________________________
    % Data partition and Number of sweeps
%     tic
%     [idx_mcmc, theta_RJ, models_mcmc, count_mcmc, Nm] = rj_mcmc(y, H, n, Ns, Nb);
%     time_mcmc(run) = toc;
% 
% 
%     % Check through all models
%     idx_corr_mcmc = 0;
%     for m = 1:length(models_mcmc(:,1))
%         if (sum(models_mcmc(m,:) == idx_h_padded ) == dy)
%             idx_corr_mcmc = m;
%         end
%     end


    % PJ ORLS___________________________________________________
    tic
    init = dy + 1;
    [theta_k, Hk, k_store, k_mode, models_orls, count_orls, idx_orls, J] = pj_orls(y, H, dy, var_y, init, Tb);
    toc
    time_orls(run) = toc;



    % Check through all models
    idx_corr_orls = 0;
    for m = 1:length(models_orls(:,1))
        if (sum(models_orls(m,:) == idx_h_padded ) == dy)
            idx_corr_orls = m;
        end
    end



    % Olin LASSO
    tic
    [theta_olasso, idx_olasso, models_olasso, count_lasso, J_lasso] = olasso(y, H, t0, epsilon);
    toc
    time_olasso(run) = toc;
    
    % Check through all models
    idx_corr_olasso = 0;
    for m = 1:length(models_olasso(:,1))
        if (sum(models_olasso(m,:) == idx_h_padded ) == dy)
            idx_corr_olasso = m;
        end
    end



    orls_run(run) = idx_corr_orls;
    %mcmc_run(run) = idx_corr_mcmc;
    olin_run(run) = idx_corr_olasso;
end


% Anything below 5
orls_run(orls_run > 4) = 5;
%mcmc_run(mcmc_run > 4) = 5;


% Average run time ratio
%avg_time = mean(time_mcmc./time_orls);


R = length(orls_run);

str_dy = num2str(dy);
str_k = num2str(dy - ps);
str_T = num2str(T);
str_v = num2str(var_y);
str_R = num2str(R);

% filename = join(['Results/T', str_T, '_K', str_dy, '_k', str_k, '_v', str_v, ...
%     '_R', str_R, '.mat']);
% 
% save(filename)


% Bar plot
figure;
subplot(2,2,1)
per_lasso = count_lasso/sum(count_lasso);
b_lasso = bar(per_lasso, 'FaceColor', 'flat');
ylabel('Number of Visits')
title('OLinLASSO Models visited ','FontSize',20)
set(gca, 'FontSize', 20); 
grid on
if (idx_corr_olasso == 0)
    text(1, 0.5*max(per_lasso), 'True Model NOT visited', 'FontSize', 15)
else
    b_lasso.CData(idx_corr_olasso,:) = [0, 0.5, 0];
end


% Bar plot
subplot(2,2, 2)
per_orls = count_orls/sum(count_orls);
b_orls = bar(per_orls, 'FaceColor', 'flat');
%ylim([0, 0.4])
ylabel('Number of Visits')
title('ORLS Models visited ','FontSize',20)
set(gca, 'FontSize', 20);
grid on
if (idx_corr_orls==0)
    text(1,0.5*max(per_orls), 'True Model NOT visited', 'FontSize', 15)
else
    b_orls.CData(idx_corr_orls,:) = [0.5, 0, 0];
end

subplot(2,2,3)
plot(J_lasso)
hold on
plot(J)


% Bar plot
% subplot(1,3, 3)
% b_mcmc = bar(count_mcmc/sum(count_mcmc), 'FaceColor', 'flat');
% %ylim([0, 0.4])
% ylabel('Number of Visits')
% title('RJMCMC Models visited','FontSize',20)
% set(gca, 'FontSize', 20);
% grid on
% b_mcmc.CData(idx_corr_mcmc,:) = [0, 0, 0];
