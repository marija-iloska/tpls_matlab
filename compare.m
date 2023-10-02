clear all
close all
clc

% Settings
var_y = 0.01;   % Variance
ps = 4;     % Sparsity percent
dy = 7;      % System dimension
T = 200;      % Time series length
r = 3;       % Range of input data H
rt = 2;      % Range of theta

% ORLS params
Tb = 50;

% OLASSO params
epsilon = 1e-7;
t0 = dy + 3 ;

% RJMCMC params
n = round(0.2*T);
Ns = 2000;
Nb = 1000;




%Create data
[y, H, theta] = generate_data(T, dy, r, rt,  ps, var_y);      
idx_h = find(theta ~= 0)';
init = dy + 1;


% RJ MCMC ___________________________________________________
%Data partition and Number of sweeps
tic
[idx_mcmc, theta_RJ, models_mcmc, count_mcmc, Nm] = rj_mcmc(y, H, n, Ns, Nb);
toc


%Pad original true indices for comparison
idx_h_padded = [idx_h zeros(1, dy - length(idx_h))];

% Check through all models
for m = 1:length(models_mcmc(:,1))
    if (sum(models_mcmc(m,:) == idx_h_padded ) == dy)
        idx_corr_mcmc = m;
    end
end



% Bar plot
figure;
b_mcmc = bar(count_mcmc/sum(count_mcmc), 'FaceColor', 'flat');
ylabel('Number of Visits')
title('RJMCMC Models visited','FontSize',20)
set(gca, 'FontSize', 20);
grid on
b_mcmc.CData(idx_corr_mcmc,:) = [0, 0, 0];



% Olin LASSO
tic
[theta_olasso, idx_olasso, models_olasso, count_lasso] = olasso(y, H, t0, epsilon);
toc

% Check through all models
idx_corr_olasso = [];
for m = 1:length(models_olasso(:,1))
    if (sum(models_olasso(m,:) == idx_h_padded ) == dy)
        idx_corr_olasso = m;
    end
end



% Bar plot
figure;
per_lasso = count_lasso/sum(count_lasso);
b_orls = bar(per_lasso, 'FaceColor', 'flat');
ylabel('Number of Visits')
title('OLinLASSO Models visited ','FontSize',20)
set(gca, 'FontSize', 20); 
grid on
if (isempty(idx_corr_olasso))
    text(1, 0.5*max(per_lasso), 'True Model NOT visited', 'FontSize', 15)
else
    b_orls.CData(idx_corr_olasso,:) = [0, 0.5, 0];
end


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
per_orls = count_orls/sum(count_orls);
b_orls = bar(per_orls, 'FaceColor', 'flat');
ylabel('Number of Visits')
title('ORLS Models visited ','FontSize',20)
set(gca, 'FontSize', 20); 
grid on
if (isempty(idx_corr_orls))
    text(1, 0.5*max(per_orls), 'True Model NOT visited', 'FontSize', 15)
else
    b_orls.CData(idx_corr_orls,:) = [0.5, 0, 0];
end





