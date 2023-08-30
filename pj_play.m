clc
close all
clear all

% Write a code that computes jumping up, staying same, jumping down Jkt

% Settings
var_y = 0.1; % Variance
p_s = 0.5;   % Sparsity percent
dx = 7;      % System dimension
T = 100;     % Time series length
r = 0.5;     % Range of input data H
rt = 5;      % Range of theta

% Current dimension dtheta
k = 1;

%Create data
[y, H, theta] = generate_data(T, dx, r,rt,  p_s, var_y);

% Store
H_true = H;

% Initialize using t data points
t = 2;
[J, theta_k, Dk, Hk, Sigma] = initialize(y, H, t, k, var_y);
count1= 0;
count0 = 0;
Jup_track(1:2) = 0;
Jdown_track(1:2) = 0;
% Start time loop
tic
for t = 3:T-1

    % JUMP UP +
    J_up = Inf;
    if (dx > k)
        [theta_up, H_up, J_up, Sigma_up, Dk_up, k_up] = jump_up(y, dx, k, Dk, theta_k, J, H, t, var_y) ;
        Jup_track(t) = J_up;
    else
        Jup_track(t) = Jup_track(t-1);
    end

    % JUMP DOWN -
   J_down = Inf;
   if (k > 1)
        [theta_down, H_down, J_down, Sigma_down, Dk_down, k_down] = jump_down(y, k, Dk, theta_k, J, H, t, var_y);
        Jdown_track(t) = J_down;
   else
       Jdown_track(t) = Jdown_track(t-1);
   end

   % STAY SAME
   [theta_stay, Sigma_stay, J_stay] = time_update(y, Hk, t, theta_k, var_y, Dk, J);
   Dk_stay = Sigma_stay/var_y;

   J_track(t) = J_stay;
   Js = [J_stay, J_up, J_down];
   %minJ = find(Js == min(Js));
   Ws = exp(-(Js-min(Js)));
   %Ws = 1./Js;
   Ws = Ws./sum(Ws);
   minJ = datasample(1:3, 1, 'Weights', Ws);

  if (minJ == 3)
      k = k_down;
      H = H_down;
      Hk = H(1:t+1, 1:k);
      theta_k = theta_down;
      J = J_down;
      Sigma = Sigma_down;
      Dk = Dk_down;
  elseif (minJ == 2)
      H = H_up;
      k = k_up;
      Hk = H(1:t+1, 1:k);
      theta_k = theta_up;
      J = J_up;
      Sigma = Sigma_up;
      Dk = Dk_up;
  else
      theta_k = theta_stay;
      Sigma = Sigma_stay;
      J = J_stay;
      Dk = Sigma/var_y;
      Hk = H(1:t+1, 1:k);
  end
  k_track(t) = k;

end
toc


[~, idx] = ismember(H_true(1,:), H(1,:));
 

mode(k_track)
% mode(k_track(end-200:end))

theta'
theta_k'

lwd = 1;
fsz = 20;

theta_length = theta;
theta_length(theta_length==0) = [];
k_truth = length(theta_length);

tN = T-1;
figure(1)
plot(k_track(1:tN), '.', 'Color', 'm', 'MarkerSize',10)
hold on
yline(k_truth, 'k', 'linewidth', 1)
set(gca, 'FontSize', 16)
ylabel('Model Order k', 'FontSize', fsz)
xlabel('Time', 'FontSize', fsz)
ylim([0,dx])
legend('Estimate', 'True Order', 'FontSize', fsz)

tN = 60;
tN = T-1;
figure(2)
sz = 10;
plot(J_track(1:tN), '.', 'Color', 'k', 'MarkerSize', sz)
hold on
plot(Jup_track(1:tN), '.', 'Color', 'r','MarkerSize', sz)
hold on
plot(Jdown_track(1:tN), '.', 'Color', 'b','MarkerSize', sz)
hold on


plot(J_track(1:tN),  'Color', 'k', 'LineWidth', lwd)
hold on
plot(Jup_track(1:tN),  'Color', 'r', 'LineWidth', lwd)
hold on
plot(Jdown_track(1:tN),  'Color', 'b', 'LineWidth', lwd)
set(gca, 'FontSize', 16)
xlabel('Time', 'FontSize', fsz)
ylabel('Error', 'FontSize', fsz)
legend('J_{STAY}', 'J_{UP}', 'J_{DOWN}', 'FontSize', fsz)

