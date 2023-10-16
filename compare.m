
%% Main Code
clear all
close all
clc

% Settings
var_y = 0.1;   % Variance
ps = 5;     % Sparsity percent
dy = 10;      % System dimension
r =  2;       % Range of input data H
rt = 1;      % Range of theta
T = 400;

% OLASSO params
epsilon = 1e-7;
t0 = 20;

% JPLS params
Tb = 20;
init = t0;

% rjMCMC params
% n = round(0.2*T);
% Ns = 2000;
% Nb = 1000;

% Parallel runs
R = 1;

% Initialize arrays
% time_mcmc = zeros(R);
% mcmc_run = zeros(R);

time_jpls = zeros(R);
jpls_run = zeros(R);
time_olin = zeros(R);
olin_run = zeros(R);


for run = 1:R

    %Create data
    [y, H, theta] = generate_data(T, dy, r, rt,  ps, var_y);
    idx_h = find(theta ~= 0)';


    % Pad original true indices for comparison
    idx_h_padded = [idx_h zeros(1, dy - length(idx_h))];


    % PJ ORLS___________________________________________________
    tic
    [theta_jpls, H_jpls,  models_jpls, count_jpls, idx_jpls, J, e] = jpls(y, H, dy, var_y, init, Tb);
    toc
    time_jpls(run) = toc;
    J_jpls(run,:) = J;
    e_jpls(run,:) = e;

    % Check through all models
    idx_corr_jpls = 0;
    for m = 1:length(models_jpls(:,1))
        if (sum(models_jpls(m,:) == idx_h_padded ) == dy)
            idx_corr_jpls = m;
        end
    end
    best_jpls = models_jpls(1,:);



    % Olin LASSO___________________________________________________
    tic
    [theta_olasso, idx_olin, models_olin, count_olin, e, J] = olasso(y, H, t0, epsilon);
    toc
    time_olin(run) = toc;
    J_olin(run,:) = J;
    e_olin(run,:) = e;

    % Check through all models
    idx_corr_olin = 0;
    for m = 1:length(models_olin(:,1))
        if (sum(models_olin(m,:) == idx_h_padded ) == dy)
            idx_corr_olin = m;
        end
    end
    best_olin = models_olin(1,:);


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


    % Store model ranks
    jpls_run(run) = idx_corr_jpls;
    olin_run(run) = idx_corr_olin;
    %mcmc_run(run) = idx_corr_mcmc;
end



% Anything below 5
jpls_run(jpls_run > 4) = 5;
olin_run(olin_run > 4) = 5;
%mcmc_run(mcmc_run > 4) = 5;


str_dy = num2str(dy);
str_k = num2str(dy - ps);
str_T = num2str(T);
str_v = num2str(var_y);
str_R = num2str(R);

% filename = join(['Results/T', str_T, '_K', str_dy, '_k', str_k, '_v', str_v, ...
%     '_R', str_R, '.mat']);
% 
% save(filename)


%% PLOTS 

% RESIDUAL PREDICTIVE ERROR PLOT
fsz = 20;
figure;
plot(t0+1:T, mean(e_olin,1), 'Color', [0, 0.5, 0], 'LineWidth', 0.5)
hold on
plot(t0+1:T, mean(e_jpls,1), 'Color', [0.5, 0, 0], 'LineWidth', 1)
hold on
text(t0+2, 0.5*max(e_jpls),  't_0',   'Color' , [0, 0, 0],'FontSize', 15)
hold on
xline(t0, 'Color', [0, 0, 0])
xlabel('Time', 'FontSize', fsz)
ylabel('Residual Predictive Error', 'FontSize', fsz)
legend('OLinLASSO','JPLS',  'FontSize',15); %, 'Location','northwest')



title_str = join(['\sigma^2_y = ', str_v, ...
    ',  h ~ N( 0, ', num2str(r), 'I ), ' , '  theta ~ N( 0, ', num2str(rt), 'I ) ']) ; %, ' K = ', str_dy, ',  p = ', str_k ]);

sgtitle(title_str, 'FontSize', 15)



%% 
% filename = join(['figsPE/K', str_dy, '_k', str_k, '_v', str_v, '_h', num2str(r), '.eps']);
% print(gcf, filename, '-depsc2', '-r300');



%%
% 
% figure;
% range = 50 : 60;
% range = range + 230;
% plot(mean(J_dec{1}(range), 1), 'LineWidth', 2)
% hold on
% plot(mean(J_dec{2}(range), 1), 'LineWidth', 2)
% hold on
% plot(mean(J_dec{3}(range), 1), 'LineWidth', 2)
% title('DEC', 'FontSize', 15)
% legend('STAY', 'UP', 'DOWN','FontSize', 15)


% Bar plot
% subplot(1,3, 3)
% b_mcmc = bar(count_mcmc/sum(count_mcmc), 'FaceColor', 'flat');
% %ylim([0, 0.4])
% ylabel('Number of Visits')
% title('RJMCMC Models visited','FontSize',20)
% set(gca, 'FontSize', 20);
% grid on
% b_mcmc.CData(idx_corr_mcmc,:) = [0, 0, 0];


% % Bar plot
% figure;
% subplot(1,3,1)
% per_lasso = count_lasso/sum(count_lasso);
% b_lasso = bar(per_lasso, 'FaceColor', 'flat');
% ylim([0, 0.5])
% ylabel('Number of Visits')
% title('OLinLASSO Models visited ','FontSize',20)
% set(gca, 'FontSize', 20); 
% grid on
% if (idx_corr_olasso == 0)
%     text(1, 0.5*max(per_lasso), 'True Model NOT visited', 'FontSize', 15)
% else
%     b_lasso.CData(idx_corr_olasso,:) = [0, 0.5, 0];
% end
% 
% 
% % Bar plot
% subplot(1,3, 2)
% per_orls = count_orls/sum(count_orls);
% b_orls = bar(per_orls, 'FaceColor', 'flat');
% ylim([0, 0.5])
% ylabel('Number of Visits')
% title('JPLS Models visited ','FontSize',20)
% set(gca, 'FontSize', 20);
% grid on
% if (idx_corr_orls==0)
%     text(1,0.5*max(per_orls), 'True Model NOT visited', 'FontSize', 15)
% else
%     b_orls.CData(idx_corr_orls,:) = [0.5, 0, 0];
% end



% filename = join(['figs/OLinLASSO/T', str_T, '_K', str_dy, '_k', str_k, '_v', str_v, ...
%     '_R', str_R, '.eps']);
% 
% print(gcf, filename, '-depsc2', '-r300');

