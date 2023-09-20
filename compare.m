clear all
close all
clc

% Settings
var_y = 0.001;   % Variance
ps = 3;     % Sparsity percent
dy = 5;      % System dimension
T = 300;      % Time series length
r = 1;       % Range of input data H
rt = 2;      % Range of theta
n = round(0.3*T);
Ns = 2000;
Nb = 1000;
Tb = 30;


%Create data
[y, H, theta] = generate_data(T, dy, r, rt,  ps, var_y);      
idx_h = find(theta ~= 0)';
init = 2;


% RJ MCMC ___________________________________________________
% Data partition and Number of sweeps
tic
[idx_mcmc, theta_RJ, models_mcmc, count_mcmc, Nm] = rj_mcmc(y, H, n, Ns, Nb);
toc


% Pad original true indices for comparison
idx_h_padded = [idx_h zeros(1, dy - length(idx_h))];

% Check through all models
for m = 1:length(models_mcmc(:,1))
    if (sum(models_mcmc(m,:) == idx_h_padded ) == dy)
        idx_corr_mcmc = m;
    end
end



% filename = 'TestORLS.eps'; % join(['figs23/pjorls', num2str(run), '.eps']);
% print(gcf, filename, '-depsc2', '-r300');
% 


% Bar plot
figure;
b_mcmc = bar(count_mcmc/Ns, 'FaceColor', 'flat');
ylabel('Number of Visits')
title('RJMCMC Models visited','FontSize',20)
set(gca, 'FontSize', 20);
grid on
b_mcmc.CData(idx_corr_mcmc,:) = [0, 0, 0];

% filename = 'TestMCMC.eps';  % join(['figs23/rjmcmc', num2str(run), '.eps']);
% print(gcf, filename, '-depsc2', '-r300');


% PJ ORLS___________________________________________________
tic
[theta_k, Hk, k_store, k_mode, models_orls, count_orls, idx_orls] = pj_orls(y, H, dy, var_y, init, Tb);
toc 

[~, idx_orls_last] = ismember(Hk(1,:), H(1,:));
idx_orls_last = sort(idx_orls_last, 'ascend');




% Check through all models
for m = 1:length(models_orls(:,1))
    if (sum(models_orls(m,:) == idx_h_padded ) == dy)
        idx_corr_orls = m;
    end
end


% Bar plot
figure;
b_orls = bar(count_orls/(T-3), 'FaceColor', 'flat');
ylabel('Number of Visits')
title('ORLS Models visited ','FontSize',20)
set(gca, 'FontSize', 20); 
grid on
b_orls.CData(idx_corr_orls,:) = [0.5, 0, 0];



