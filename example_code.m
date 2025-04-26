clear all
close all
clc


% example_code.m
% This script is an example code on how to run TPLS.
% The code is easy to use and follows the format:

% TPLS_________________________________________
% Initialize --> model(n0)
% for n = n0+1, ... N
%   receive new data D(n) = { y(n), H(n,:) }
%   compute pred error of models above <-- recursive PLS
%   compute pred error of models below <-- recursive PLS
%   move to model with lowest pred error <-- aORLS or dORLS
%   model(n) = time_update( model(n-1), D(n) ) <-- RLS 
% end


% Some Important Variables -----------------------------------------

% y, H  - Output data, Feature matrix
% theta - True model parameter
% idx - set of the indices of the selected features
% J - the predictive error of the last estimate
% correct - how many features selected were in the true model
% incorrect - how many features selected were not in the true model
% S_features_used - array with indices of features in current model
% S_features_used - array with indices of features not in the current model



% Paths to access functions from other folders
function_paths = [genpath('tpls/'), genpath('util/'), genpath('baselines/')];

% Add the paths
addpath(function_paths)


%% GENERATE SYNTHETIC DATA
% Settings
var_y = 0.1;            % Observation noise Variance
ps = 10;                 % Number of 0s in theta
K = 20;                 % Number of available features
var_features = 1;      % Range of input data H
var_theta = 1;        % Variance of theta
N = 300;                 % Number of data points
p = K - ps;             % True model dimension


%Create data
[y, H, theta] = generate_data(N, K, var_features, var_theta,  ps, var_y);
idx = find(theta ~= 0)';


%% INITIALIZE

% Initial batch of data (choose any that initially satisfy k < n0 )
n0 = K+1;

% Get estimate using any random k < n0 number of features 
k = 5;

S_features_used = datasample(1:K, k, 'replace', false);
S_features_unused = setdiff(1:K, S_features_used);

% Initial batch of data
y0 = y(1:n0);
H0 = H(1:n0,:);

[theta_k, Dk] = initialize(y0, H0(:,S_features_used));
model_at_n0 = {theta_k, Dk};

% Model storage
correct = zeros(1,N-n0);
incorrect = zeros(1,N-n0);

% Predictive errors init
J = 0;
J_pred = J;
J_above = Inf;
J_below = Inf;


%% JPLS LOOP
% Start time loop
tic
for n = n0+1:N

    % Data
    Hn = H(1:n, S_features_used);
    yn = y(1:n);

    % CURRENT MODEL PE
    J = J + (y(n) - Hn(n,:)*theta_k)^2;

    % MODEL ABOVE PE
    if k < K
        rm_idx = nan;
        parfor m = 1:K-k
             difference = pred_error_model_above(yn, [Hn H(1:n, S_features_unused(m))], n, n0, model_at_n0);
             J_above(m) = J + difference;     
        end
    end

   % MODEL BELOW PE
   if k > 1
    parfor m = 1:k
         difference = pred_error_model_below(yn, Hn, n, n0, m, model_at_n0);
         J_below(m) = J + difference;
    end
   end

    % Find min PE
    m_above  = find(J_above == min(J_above));
    m_below = find(J_below == min(J_below));
    minJ = min([J_above(m_above), J_below(m_below), J]);

    % IF BELOW
    if minJ == J_below(m_below)


        % Update base model
        [theta0, D0] = descendingORLS(model_at_n0{1}, model_at_n0{2}, m_below);
        model_at_n0 = {theta0, D0};


        % Move to neighbor model below
        [theta_k, Dk] = descendingORLS(theta_k, Dk, m_below);

        % Update predictive error
        J = J_below(m_below);

        % Update feature sets
        S_features_unused(end+1) = S_features_used(m_below);
        S_features_used(m_below) = [];

        % Update model dimension
        k = k - 1;

    % IF ABOVE
    elseif minJ == J_above(m_above)     

        % Update base model
        [theta0, D0] = ascendingORLS(y0, H0(:,S_features_used), H0(:, S_features_unused(m_above)), n0, model_at_n0{2}, model_at_n0{1});
        model_at_n0 = {theta0, D0};

        % Move to neighbor model above
        [theta_k, Dk] = ascendingORLS(yn(1:n-1), Hn(1:n-1,:), H(1:n-1, S_features_unused(m_above)), n-1, Dk, theta_k);
        
        % Update predictive error
        J = J_above(m_above);

        % Update feature sets
        S_features_used(end+1) = S_features_unused(m_above);
        S_features_unused(m_above) = [];

        % Update model dimension
        k = k + 1;

    end

    % Place holders
    J_above = Inf;
    J_below = Inf;
    
    % Store predictive error
    J_pred(end+1) = J;


    % TIME UPDATE    
    % theta(k,n) <-- theta(k,n-1) and Dk(n) <-- Dk(n-1)
    [theta_k, Dk] = RLS(y(n), H(n,S_features_used), theta_k, Dk);

    % EVALUATION
    correct(n-n0) = sum(ismember(S_features_used, idx));
    incorrect(n-n0) = length(S_features_used) - correct(n-n0);


end
toc
% Concatenate results
plot_stats = {correct, incorrect};


% Results for barplots
[correct, incorrect] = plot_stats{:};
tpls_features = [correct;  incorrect];



%% PLOTTING

% Colors, FontSize, Linewidths
load plot_settings.mat

% Time range to plot
time_plot = n0:N;

% Create figure 
figure('Renderer', 'painters', 'Position', [200 300 1000 400])

% Features Bar plot
subplot(1,2,1)
formats = {fsz, fszl, fsz, lwdt, c_tpls, c_inc, c_true, '', 'n'};
bar_plots(tpls_features, n0+1, N, p, K, formats)

% Predictive Error plot
subplot(1,2,2)
plot(time_plot, J_pred, 'Color', c_tpls, 'LineWidth', lwd)
xlim([n0+1, N])
set(gca, 'FontSize', 15)
legend('J_{TPLS}', 'FontSize', fszl)
xlabel('Time', 'FontSize', fsz)
ylabel('Predictive Error', 'FontSize', fsz)
grid on
sgtitle('\bf{TPLS}', 'FontSize',fsz)
hold on

