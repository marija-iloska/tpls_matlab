clear all
close all
clc

% Settings
var_y = 0.1; % Variance
ps = 0.3;    % Sparsity percent
dy = 7;      % System dimension
T = 50;     % Time series length
r = 0.5;     % Range of input data H
rt = 5;      % Range of theta


R = 100;
tic
parfor run = 1 : R

    %Create data
    [y, H, theta] = generate_data(T, dy, r, rt,  ps, var_y);

    % Predictive jump ORLS
    [theta_k, Hk, k_store, k_mode] = pj_orls(y, H, dy, var_y);

    % Evaluate
    [dk, dk_mode, dk_est, check_mode, check, over, under, up, down] = eval_orls(theta, k_store, T-1);

    % Collect statistics
    count(run) = check;
    count_mode(run) = check_mode;
    count_over(run) = over;
    count_under(run) = under;
    count_up(run) = up;
    count_down(run) = down;

end
toc

c = sum(count);
cm = sum(count_mode);
un = sum(count_under);
ov = sum(count_over);
up = sum(count_up);
down = sum(count_down);



% Decorations
sz = 40;

% Reds
rp = [212, 19, 19]/256;
rs = [94, 4, 4]/256;

% Greens
gp = [87, 194, 105]/ 256;
gs = [5, 102, 37]/256;

% Blues
bp = [145, 165, 235]/256;
bs = [18, 22, 148]/256;

% Tir
tp = [30, 214, 187]/256;
ts = [14, 117, 102]/256;

% Magentas
ms = [120, 16, 71]/256;
mp = [214, 118, 169]/256;

% yp = Hk*theta_k;


% Bar plot
figure(1)
b = bar([cm, ov, un, up, down]*100/R, 'FaceColor', 'flat');
ylabel('Percent')
title(join(['R = ', num2str(R), ' runs']))
set(gca, 'FontSize', 20, 'xticklabel', {'Correct', 'Over by 1', 'Under by 1', 'Over', 'Under'})
b.CData(1,:) = gs;
b.CData(2,:) = rs;
b.CData(3,:) = bs;
b.CData(4,:) = ms;
b.CData(5,:) = [0, 0, 0];
ylim([0,100])
grid on


% figure(2)
% yline(dk, 'Color', gp, 'LineWidth', 3, 'LineStyle', '-')
% hold on
% scatter(1:T, k_store(1:T), sz,'filled', 'MarkerFaceColor', rs, 'Linewidth', 1);
% set(gca,'FontSize', 15, 'Linewidth',1)
% xlabel('Time', 'FontSize', 20)
% ylabel('Model Order','FontSize', 20)
% title('Convergence', 'FontSize',20)
% ylim([0, dy])
% legend('True', 'Estimate', 'FontSize', 15)








