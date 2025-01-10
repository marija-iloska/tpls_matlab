function [theta_occd, idx_occd, J, plot_stats, idx_store] = occd(y, H, t0, var_y, idx_h)

% Dimensions
T = length(y);
K = length(H(1,:));

% Initial batch start
theta_occd = zeros(K,1);

% Denominators for each feature
for j = 1:K
    % Indexes of all elements except jth
    all_but_j{j} = setdiff(1:K, j);
end

rn = zeros(1,K);
Rn = zeros(K,K);
J = [];
correct = [];
incorrect = [];
idx_store = {};

for t = 1:T
    % Receive new data point Xn, yn
    Xn = H(t,:);
    yn = y(t);

    [theta_occd, rn, Rn] = occd_update(yn, Xn, rn, Rn, t, K, theta_occd, all_but_j, var_y);

    if t>t0
        [J(end+1)] = pred_error_baselines(y, H, t, t0, var_y, theta_occd);
        idx_occd = find(theta_occd ~= 0)';
        idx_store{end+1} = idx_occd;
    
        % EVALUATION
        correct(t-t0) = sum(ismember(idx_occd, idx_h));
        incorrect(t-t0) = length(idx_occd) - correct(t-t0);
    end

end

%J = [];


% for t = t0+1:T
% 
%     % CALL occd
%     % Receive new data point Xn, yn
%     Xn = H(t,:);
%     yn = y(t);
% 
%     rn = rn + yn*Xn;
%     Rn = Rn + Xn'*Xn;
% 
%     [theta_occd, rn, Rn] = occd_update(yn, Xn, rn, Rn, t, K, theta_occd, all_but_j, var_y);
%     [J(end+1)] = pred_error_baselines(y, H, t, t0, var_y, theta_occd);
%     idx_occd = find(theta_occd ~= 0)';
% 
%     % EVALUATION
%     correct(t-t0) = sum(ismember(idx_occd, idx_h));
%     incorrect(t-t0) = length(idx_occd) - correct(t-t0);
% 
% end


% Concatenate results
plot_stats = {correct, incorrect};

end