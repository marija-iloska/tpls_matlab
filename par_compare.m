clear all
close all
clc

% Settings
var_y = 1;   % Variance
ps = 12;     % Sparsity percent
dy = 15;      % System dimension
r =  0.5;       % Range of input data H
rt = 0.5;      % Range of theta
T = 6000;
D = 20;

% OLASSO params
epsilon = 1e-7;
t0 = 100;

% JPLS params
Tb = 300;
init = t0;


% Parallel runs
R = 32;

% Initialize arrays
time_orls = zeros(R);
time_olasso = zeros(R);
orls_run = zeros(R);
olin_run = zeros(R);

J_oi = zeros(R, T-t0);
J_ol = zeros(R, T-t0);

tic
parfor run = 1:R

    %Create data
    [y, H, theta] = generate_data(T, dy, r, rt,  ps, var_y);
    idx_h = find(theta ~= 0)';


    % Pad original true indices for comparison
    idx_h_padded = [idx_h zeros(1, dy - length(idx_h))];


    % JPLS ___________________________________________________
    tic
    [theta_k, Hk, k_store, models_orls, count_orls, idx_orls, J_incr] = pj_orls(y, H, dy, var_y, init, Tb, D);
    toc
    time_orls(run) = toc;
    %J_orls(run, :) = J;
    J_oi(run,:) = J_incr;



    % Check through all models
    idx_corr_orls = 0;
    for m = 1:length(models_orls(:,1))
        if (sum(models_orls(m,:) == idx_h_padded ) == dy)
            idx_corr_orls = m;
        end
    end
    %best_orls = models_orls(1,:);



    % Olin LASSO
    tic
    [theta_olasso, idx_olasso, models_olasso, count_lasso, J, J_incr] = olasso(y, H, t0, epsilon);
    toc
    time_olasso(run) = toc;
%     J_lasso(run,:) = J;
    J_ol(run, :) = J_incr;
    
    % Check through all models
    idx_corr_olasso = 0;
    for m = 1:length(models_olasso(:,1))
        if (sum(models_olasso(m,:) == idx_h_padded ) == dy)
            idx_corr_olasso = m;
        end
    end
    %best_lasso = models_olasso(1,:);


    orls_run(run) = idx_corr_orls;
    olin_run(run) = idx_corr_olasso;
end
toc




% Anything below 5
orls_run(orls_run > 4) = 5;




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

%% PLOTS 



% fsz = 20;
% figure;
% plot(init+1:T-1, mean(J_orls,1), 'Color', [0.5, 0, 0], 'LineWidth', 2)
% hold on
% plot(t0+1:T-1, mean(J_lasso,1), 'Color', [0, 0.5, 0], 'LineWidth', 2)
% set(gca, 'FontSize',15)
% xlabel('Time', 'FontSize', fsz)
% ylabel('Predictive Error', 'FontSize', fsz)
% legend('JPLS', 'OLinLASSO', 'FontSize',fsz)

fsz = 15;
figure;
plot(init+1:T, mean(J_oi,1), 'Color', [0.5, 0, 0], 'LineWidth', 2)
hold on
plot(t0+1:T, mean(J_ol,1), 'Color', [0, 0.5, 0], 'LineWidth', 2)
hold on
xline(t0, 'Color', [0, 0.5, 0])
hold on
text(t0+2, 0.5*max(mean(J_oi,1)),  't_0',   'Color' , [0, 0.5, 0],'FontSize', 15)
hold on
xline(init, 'Color', [0.5, 0, 0])
hold on
text(init+2, 0.5*max(mean(J_ol,1)), 't_0', 'Color' , [0.5, 0, 0],  'FontSize', 15)
set(gca, 'FontSize',15)
xlabel('Time', 'FontSize', fsz)
ylabel('Predictive Error', 'FontSize', fsz)
legend('JPLS', 'OLinLASSO', 'FontSize',fsz)


% filename = join(['figs/OLinLASSO/T', str_T, '_K', str_dy, '_k', str_k, '_v', str_v, ...
%     '_R', str_R, '.eps']);
% 
% print(gcf, filename, '-depsc2', '-r300');
% 

