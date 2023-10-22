
%% Main Code
clear all
close all
clc

% Settings
var_y = 1;   % Variance
ps = 10;     % Sparsity percent
dy = 20;      % System dimension
r =  1;       % Range of input data H
rt = 0.5;      % Range of theta
T = 300;

% OLASSO params
epsilon = 1e-7;
t0 = 30;

% JPLS params
Tb = 5;
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
    [theta_jpls, H_jpls,  models_jpls, count_jpls, idx_jpls, e, J_pred, jpls_correct] = jpls(y, H, dy, var_y, init, Tb, idx_h);
    toc
    time_jpls(run) = toc;
    %J_jpls(run,:) = J_total;
    Jpred_jpls(run,:)=J_pred;
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
    [theta_olasso, idx_olin, models_olin, count_olin, e, J_pred, olin_correct] = olasso(y, H, t0, epsilon, var_y, idx_h);
    toc
    time_olin(run) = toc;
    %J_olin(run,:) = J_total;
    Jpred_olin(run,:) = J_pred;
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


    [J_true, e_true] = true_PE(y, H, t0, T, idx_h, var_y);

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
lwd = 1.5;
e_lim = max(max([e_olin, e_jpls]));
figure;
subplot(4,1,1)
plot(t0+1:T, mean(e_olin,1), 'Color', [0, 0.5, 0], 'LineWidth', 0.5)
hold on
plot(t0+1:T, mean(e_jpls,1), 'Color', [0.5, 0, 0], 'LineWidth', 1)
hold on
text(t0+2, 0.8*e_lim,  't_0',   'Color' , [0, 0, 0],'FontSize', 15)
hold on
xline(t0, 'Color', [0, 0, 0], 'linewidth',1)
xlabel('Time', 'FontSize', fsz)
ylabel('e_{k,t}', 'FontSize', fsz)
legend('OLinLASSO','JPLS',  'FontSize',15); %, 'Location','northwest')


subplot(4,1,2)
plot(t0+1:T, mean(e_olin.^2,1), 'Color', [0, 0.5, 0], 'LineWidth', 0.5)
hold on
plot(t0+1:T, mean(e_jpls.^2,1), 'Color', [0.5, 0, 0], 'LineWidth', 1)
hold on
text(t0+2, 1.5*e_lim,  't_0',   'Color' , [0, 0, 0],'FontSize', 15)
hold on
xline(t0, 'Color', [0, 0, 0], 'linewidth',1)
xlabel('Time', 'FontSize', fsz)
ylabel('e_{k,t}^2', 'FontSize', fsz)
legend('OLinLASSO', 'JPLS',  'FontSize',15); %, 'Location','northwest')


sgtitle('Residual Predictive Error')


subplot(4,1,3)
plot(t0+1:T, mean(Jpred_olin,1), 'Color', [0, 0.5, 0], 'LineWidth', lwd)
hold on
plot(t0+1:T, mean(Jpred_jpls,1), 'Color', [0.5, 0, 0], 'LineWidth', lwd)
hold on
plot(t0+1:T, J_true(2:end), 'k', 'LineStyle','--', 'LineWidth', lwd)
hold on
text(t0+2, 0.8*e_lim,  't_0',   'Color' , [0, 0, 0],'FontSize', 15)
hold on
xline(t0, 'Color', [0, 0, 0], 'linewidth',3)
xlabel('Time', 'FontSize', fsz)
ylabel('J_{k,t}', 'FontSize', fsz)
legend('OLinLASSO','JPLS', 'True Model PE',  'FontSize',15, 'Location','northwest')

subplot(4,1,4)
plot(t0+1:T, cumsum(Jpred_olin,2), 'Color', [0, 0.5, 0], 'LineWidth', lwd)
hold on
plot(t0+1:T, cumsum(Jpred_jpls,2), 'Color', [0.5, 0, 0], 'LineWidth', lwd)
hold on
plot(t0+1:T, cumsum(J_true(2:end)), 'k', 'LineStyle','--', 'LineWidth', lwd)
hold on
text(t0+2, 0.8*e_lim,  't_0',   'Color' , [0, 0, 0],'FontSize', 15)
hold on
xline(t0, 'Color', [0, 0, 0], 'linewidth',3)
xlabel('Time', 'FontSize', fsz)
ylabel('Cumulative J_{k,t}', 'FontSize', fsz)
legend('OLinLASSO','JPLS',  'FontSize',15, 'Location','northwest')

%sgtitle('Predictive Error', 'FontSize', fsz)

title_str = join(['\sigma^2_y = ', str_v, ...
    ',  h ~ N( 0, ', num2str(r), 'I ), ' , '  theta ~ N( 0, ', num2str(rt), 'I ) ']) ; %, ' K = ', str_dy, ',  p = ', str_k ]);

sgtitle(title_str, 'FontSize', 15)



%% 
filename = join(['figsPE/K', str_dy, '_k', str_k, '_v', str_v, '_h', num2str(r), '.eps']);
print(gcf, filename, '-depsc2', '-r300');



%%


% figure;
% plot(t0+1:T, mean(J_olin,1), 'Color', [0, 0.5, 0], 'LineWidth', 0.5)
% hold on
% plot(t0+1:T, mean(J_jpls,1), 'Color', [0.5, 0, 0], 'LineWidth', 1)
% hold on
% text(t0+2, 0.8*e_lim,  't_0',   'Color' , [0, 0, 0],'FontSize', 15)
% hold on
% xline(t0, 'Color', [0, 0, 0], 'linewidth',1)
% xlabel('Time', 'FontSize', fsz)
% ylabel('Cumulative J_{k,t}', 'FontSize', fsz)
% legend('OLinLASSO','JPLS',  'FontSize',15); %, 'Location','northwest')




% Bar plot
% subplot(1,3, 3)
% b_mcmc = bar(count_mcmc/sum(count_mcmc), 'FaceColor', 'flat');
% %ylim([0, 0.4])
% ylabel('Number of Visits')
% title('RJMCMC Models visited','FontSize',20)
% set(gca, 'FontSize', 20);
% grid on
% b_mcmc.CData(idx_corr_mcmc,:) = [0, 0, 0];


% Percent Visit
per_orls = count_jpls/sum(count_jpls);
per_lasso = count_olin/sum(count_olin);
y_lim = max(max([per_orls, per_lasso])) + 0.1;


% Bar plot
figure;
subplot(1,2,1)
b_lasso = bar(per_lasso, 'FaceColor', 'flat');
ylim([0, y_lim])
ylabel('Percent Visits')
title('OLinLASSO','FontSize',20)
set(gca, 'FontSize', 20); 
grid on
if (idx_corr_olin == 0)
    text(1, 0.5*max(per_lasso), 'True Model NOT visited', 'FontSize', 15)
else
    b_lasso.CData(idx_corr_olin,:) = [0, 0.5, 0];
end


% Bar plot
subplot(1,2, 2)

b_orls = bar(per_orls, 'FaceColor', 'flat');
ylim([0, y_lim])
ylabel('Percent Visits')
title('JPLS ','FontSize',20)
set(gca, 'FontSize', 20);
grid on
if (idx_corr_jpls==0)
    text(1,0.5*max(per_orls), 'True Model NOT visited', 'FontSize', 15)
else
    b_orls.CData(idx_corr_jpls,:) = [0.5, 0, 0];
end

sgtitle('Models Visited', 'FontSize', 15)

%% 
filename = join(['figs/OLinLASSO/T', str_T, '_K', str_dy, '_k', str_k, '_v', str_v, ...
    '_R', str_R, '.eps']);

print(gcf, filename, '-depsc2', '-r300');


figure;
yline(dy, 'Color', 'k', 'LineWidth', 1.5)
hold on
yline(dy-ps, 'Color', 'b', 'LineWidth', 1.5)
hold on
plot(t0+1:T-1, jpls_correct(2:end), '.', 'Color', [0.5, 0, 0], 'MarkerSize', 15)
hold on
plot(t0+1:T-1, olin_correct(2:end), '.', 'Color', [0, 0.5, 0], 'MarkerSize', 15)
ylim([0, dy+1])
legend('Available Features', 'True Number of Features', 'JPLS', 'OLinLASSO', 'FontSize', 15)
ylabel('Correct Number of Features', 'FontSize', 15)
xlabel('Time', 'FontSize', 15)


