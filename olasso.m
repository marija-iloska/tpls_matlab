function [theta_olasso, idx_sorted, models_sorted, count_sorted, J_pred] = olasso(y, H, t0, epsilon)

% Dimensions
T = length(y);
dy = length(H(1,:));

% Define initial batch
y0 = y(1:t0);
H0 = H(1:t0, :);


% Define initial batch terms
xy0 = H0'*y0;
xx0 = H0'*H0;

% EIG
a = eig(xx0);
step = 0.01*t0/max(real(a));

% Initial estimate
[B, STATS] = lasso(H0, y0, 'CV', 5);
theta_olasso = B(:, STATS.IndexMinMSE);

% Initialize terms
xy = zeros(dy,1);
xx = zeros(dy,dy);

J_pred = zeros(1,T-1);

for t = t0+1:T-1

    % Updates
    xx = xx + H(t,:)'*H(t,:);
    xy = xy + H(t,:)'*y(t);
    
    [theta_olasso, loss{t}] = olin_lasso(xy0, xx0, xy, xx, theta_olasso, epsilon, step, t0, t, dy);
    J_pred(t) = J_pred(t-1) + (y(t+1) - H(t+1,:)*theta_olasso)^2;
    idx = find(theta_olasso ~= 0)';
    M{t-t0} = [idx, zeros(1, dy - length(idx))];
end


% Find unique models
models = unique(cell2mat(M'), 'rows');
Nm = length(models(:,1));
count = zeros(1,Nm);

% For each unique model
for m = 1:Nm

    % Current model to check for
    Mk = models(m,:);

    % Check all sweeps
    for s = 1:length(M)
        if (sum(Mk == M{s}) == dy)
            count(m) = count(m) + 1;
        end     
    end
end


% Sort count and find max
[count_sorted, idx_sorted] = sort(count, 'descend');
models_sorted = models(idx_sorted,:);
idx_sorted = nonzeros(models_sorted(1,:))';


end