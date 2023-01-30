clear all
close all
clc

% Settings
var_y = 0.1; % Variance
g = @(x) x;  % Transition function
p_s = 0.8;   % Sparsity percent
dx = 10;      % System dimension
T = 1000;     % Time series length
r = 0.5;     % Range of input data H
rt = 5;      % Range of theta

%SSM
tr = @(coeff, states) coeff*g(states);


R = 2;
tic
for run = 1 : R
    %Create data
    [y, H, theta, a] = generate_data(T, dx, r, rt, p_s, var_y, tr, g);

    % Start Time
    [theta_k, Hk, check, check_mode, over, under, up, down] = orls_jump(y, H, dx, var_y, theta);
    count(run) = check;
    count_mode(run) = check_mode;
    count_over(run) = over;
    count_under(run) = under;
    %count_other(run) = other;
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


% Bar plot
b = bar([c, ov, un, up, down]/10,'FaceColor', 'flat');
ylabel('Percent')
title('R = 1000 runs')
set(gca, 'FontSize', 20, 'xticklabel', {'Correct', 'Overestimate by 1', 'Underestimate by 1', 'Over', 'Under'})
b.CData(1,:) = gs;
b.CData(2,:) = rs;
b.CData(3,:) = bs;
b.CData(4,:) = ms;
b.CData(5,:) = [0, 0, 0];






% yline(dk, 'Color', gp, 'LineWidth',2, 'LineStyle', '-')
% hold on
% scatter(1:T, k_store(1:T), sz,'filled', 'MarkerFaceColor', rs, 'Linewidth', 1);
% set(gca,'FontSize',15, 'Linewidth',1)
% xlabel('Time', 'FontSize', 30)
% ylabel('Model Choice','FontSize', 30)
% title('Convergence', 'FontSize',30)
% ylim([0, dx])








