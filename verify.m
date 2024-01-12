
%% Main Code
clear all
close all
clc

% Settings
var_y = 1;   % Variance
ps = 3;     % Sparsity percent
K = 8;      % System dimension
r =  1;       % Range of input data H
rt = 0.5;      % Range of theta
T = 1000;
p = K - ps;
t0 = K+1;

% JPLS params
Tb = 5;
init = t0;


% Parallel runs
R = 48;

% Initialize arrays
jpls_run = zeros(R);
idx0 = datasample(1:K, ps, 'replace', false);

%Create data
[y_clean, H, theta] = fix_data(T, K, r, rt, idx0);
idx_h = find(theta ~= 0)';

tic
parfor run = 1:R

    % Add noise to data
    noise = mvnrnd(zeros(T,1), var_y*eye(T))';
    y = y_clean + noise;
    noise_super(run, :) = noise;


    % Pad original true indices for comparison
    idx_h_padded = [idx_h zeros(1, K - length(idx_h))];


    % PJ ORLS___________________________________________________
    tic
    [theta_jpls, H_jpls, model_stats,  error_stats, plot_stats] = jpls(y, H, K, var_y, init, Tb, idx_h);
    [jpls_missing, jpls_correct, jpls_wrong] = plot_stats{:};
    [models_jpls, count_jpls, idx_jpls] = model_stats{:};
    [J_pred, e] = error_stats{:};
    toc
    Jpred_jpls(run,:) = J_pred;
    e_jpls(run,:) = e;

    % GENIE 
    [J_true(run,:), e_true(run,:), noise_genie(run,:)] = true_PE(y, H, t0, T, idx_h, var_y);

    % SUPER GENIE
    e_super(run,:) = y(t0+1:end) - H(t0+1:end,:)*theta;
    J_super(run,:) = cumsum(e_super(run,:).^2);

    % MSE SUPER GENIE
    mse_super(run, :) = cumsum(noise.^2);
    mse_genie(run, :) = cumsum(noise_genie(run,:));


    % MSE GENIE 



    % BARS
    jpls_f(run, :, :) = [jpls_correct;  jpls_wrong; jpls_missing]; 



end
toc 

jpls_features = squeeze(mean(jpls_f,1));

% e_jpls = mean(e_jpls, 1);
% e_true = mean(e_true, 1);
% J_jpls = mean(Jpred_jpls, 1);
% J_true = mean(J_true,1);
% J_super = mean(J_super,1);


% Anything below 5
jpls_run(jpls_run > 4) = 5;

% For Labels
str_dy = num2str(K);
str_k = num2str(K - ps);
str_T = num2str(T);
str_v = num2str(var_y);
str_R = num2str(R);

close all
plot(cumsum(noise_super(1,:)))

figure;
subplot(1,2,1)
for r = 1:R
    plot(mse_super(r,:))
    hold on
end

subplot(1,2,2)
for r = 1:R
    plot(mse_genie(r,:))
    hold on
end


fsz = 20;
lwd = 3;

figure;
subplot(1,2,1)
for r = 1:R
    plot(J_super(r,:))
    hold on
end
plot(mean(J_super, 1), 'k', 'Linewidth', 1)
ylabel('J_{SUPER} ', 'FontSize', fsz)
xlabel('Time', 'FontSize', 15)


subplot(1,2,2)
for r = 1:R
    plot(J_true(r,:))
    hold on
end
plot(mean(J_true, 1), 'k', 'Linewidth', 1)
ylabel('J_{GENIE} ', 'FontSize', fsz)
xlabel('Time', 'FontSize', 15)

save('error_plots.mat')

