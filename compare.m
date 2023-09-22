clear all
close all
clc

% Settings
var_y = 0.001;   % Variance
ps = 1;     % Sparsity percent
dy = 5;      % System dimension
T = 500;      % Time series length
r = 1;       % Range of input data H
rt = 2;      % Range of theta
n = round(0.3*T);
Ns = 2000;
Nb = 1000;
Tb = 350;


%Create data
[y, H, theta] = generate_data(T, dy, r, rt,  ps, var_y);      
idx_h = find(theta ~= 0)';
init = 10;


% RJ MCMC ___________________________________________________
% Data partition and Number of sweeps
% tic
% [idx_mcmc, theta_RJ, models_mcmc, count_mcmc, Nm] = rj_mcmc(y, H, n, Ns, Nb);
% toc
% 
% 
% Pad original true indices for comparison
idx_h_padded = [idx_h zeros(1, dy - length(idx_h))];

% % Check through all models
% for m = 1:length(models_mcmc(:,1))
%     if (sum(models_mcmc(m,:) == idx_h_padded ) == dy)
%         idx_corr_mcmc = m;
%     end
% end



% Bar plot
% figure;
% b_mcmc = bar(count_mcmc/Ns, 'FaceColor', 'flat');
% ylabel('Number of Visits')
% title('RJMCMC Models visited','FontSize',20)
% set(gca, 'FontSize', 20);
% grid on
% b_mcmc.CData(idx_corr_mcmc,:) = [0, 0, 0];
% 


% PJ ORLS___________________________________________________
tic
[theta_k, Hk, k_store, k_mode, models_orls, count_orls, idx_orls] = pj_orls(y, H, dy, var_y, init, Tb);
toc 



% Check through all models
idx_corr_orls = [];
for m = 1:length(models_orls(:,1))
    if (sum(models_orls(m,:) == idx_h_padded ) == dy)
        idx_corr_orls = m;
    end
end


% Bar plot
figure;
b_orls = bar(count_orls/(T-Tb), 'FaceColor', 'flat');
ylabel('Number of Visits')
title('ORLS Models visited ','FontSize',20)
set(gca, 'FontSize', 20); 
grid on
if (isempty(idx_corr_orls))
    text(1,0.1, 'True Model NOT visited', 'FontSize', 15)
else
    b_orls.CData(idx_corr_orls,:) = [0.5, 0, 0];
end



% % PJ ORLS___________________________________________________
% swap = datasample(1:dy, dy, 'replace', false);
% tic
% [theta_k, Hk, k_store, k_mode, models_swap, count_swap, idx_swap] = pj_orls(y, H(:, swap), dy, var_y, init, Tb);
% toc 
% 
% 
% 
% % Check through all models
% idx_corr_orls_swap = [];
% for m = 1:length(models_swap(:,1))
%     if (sum(models_swap(m,:) == idx_h_padded ) == dy)
%         idx_corr_orls_swap = m;
%     end
% end
% 
% 
% % Bar plot
% figure;
% b_orls = bar(count_swap/(T-Tb), 'FaceColor', 'flat');
% ylabel('Number of Visits')
% title('ORLS SWAP ','FontSize',20)
% set(gca, 'FontSize', 20); 
% grid on
% if (isempty(idx_corr_orls_swap))
%     text(1,0.1, 'True Model NOT visited', 'FontSize', 15)
% else
%     b_orls.CData(idx_corr_orls_swap,:) = [0.5, 0, 0];
% end

% filename = 'TestORLS.eps'; % join(['figs23/pjorls', num2str(run), '.eps']);
% print(gcf, filename, '-depsc2', '-r300');
% 

% filename = 'TestMCMC.eps';  % join(['figs23/rjmcmc', num2str(run), '.eps']);
% print(gcf, filename, '-depsc2', '-r300');

