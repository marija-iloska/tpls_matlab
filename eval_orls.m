function [dk, dk_mode, dk_est, check_mode, check, over_est, under_est, up_est, down_est] = eval_orls(theta, k_store, T)

% Get dk true order of theta
theta(theta==0) = [];
dk = length(theta);


% Initialize counts of over/under estimates
over_est = 0;
under_est= 0;
down_est = 0;
up_est = 0;

% Take last 10 % of dk estimates
last_T = round(0.15*T);

% Final sample estimate
dk_est = k_store(T);

% Mode of last 10%
dk_mode = mode(k_store((T-last_T):T));

% Is the last sample correct
check = (dk == dk_est);

% Is the mode correct
check_mode = (dk == dk_mode);

dk_est = dk_mode;
% Counting the estimates
if (dk_est == dk + 1)
    over_est = 1;
    under_est = 0;
    down_est = 0;
    up_est = 0;
elseif (dk_est == dk - 1)
    under_est = 1;
    over_est = 0;
    %other_est = 0;
    down_est = 0;
    up_est = 0;
elseif (dk_est == dk)
    down_est = 0;
    up_est = 0;
    under_est = 0;
    over_est = 0;
elseif (dk_est > dk + 1)
    up_est = 1;
    down_est = 0;
    under_est = 0;
    over_est = 0;
elseif (dk_est < dk - 1)
    down_est = 1;
    up_est = 0;
    under_est = 0;
    over_est = 0;
end

% Final sample estimate
dk_est = k_store(T);

end