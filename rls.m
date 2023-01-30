clear all
close all
clc

% Sparse Bayesian Learning Demo

% Settings
var_x = 0.1;
var_y = 0.00001;
g = @(x) 1./(1 + exp(-x));
p_s = 0.3;
dx = 10;
T = 100;
r = 1;  % Deciding factor

% SSM
tr = @(coeff, states) coeff*g(states);
obs =@(coeff, states) coeff*states;

% Create data
[A, C, x, y_obs] = generate_states(T, dx, r, p_s, var_x, var_y, tr, obs);


% Initialie
ct = zeros(dx,T);
lambda = 1e-8;
P = eye(dx)/lambda;
beta = 0.999;


% Choose j
j = datasample(1:dx, 1);
cj = C(j,:)';
G = g(x);
count = zeros(1,dx);

% Start RLS
for t = 2:T

    % error
    e = x(j,t) - G(:, t-1)' * ct(:,t-1);
    
    % Term
    zt = P*G(:, t-1);
    
    % Kalman gain
    kt = zt/(beta + G(:,t-1)'*zt);

    % Update params
    ct(:,t) = ct(:,t-1) + kt*e;


    % Update P
    P = P/beta - kt*zt'/beta;


    % Get MSE
    mse(t-1) = sum( (cj - ct(:,t)).^2 )/dx;

end

cj';
count;

figure(1)
idx = find(cj == 0);
figure(1)
bar(count, 'FaceColor', 'k')
xlabel('Indices', 'FontSize',20)
ylabel('Number of samples selected 0', 'FontSize',20)
hold on
scatter(idx,zeros(1,length(idx)), 120, 'filled', 'r', 'Linewidth', 2)


figure(2)
k = datasample(1:dx, 1);
plot(cj(k)*ones(1,T), 'LineWidth',2)
hold on
plot(ct(k,:), 'k', 'LineStyle','--', 'LineWidth',2)
legend('True Cjk', 'Estimate Cjk')
title('Convergence')



