clear all
close all
clc

% Settings
var_x = 0.1;
g = @(x) x;
p_s = 0.3;
dx = 15;
T = 50;
r = 0.5; % Range of input data H
rt = 2;  % Range of theta

%SSM
tr = @(coeff, states) coeff*g(states);

%Create data
[y, H, theta, a] = generate_data(T, dx, r, rt, p_s, var_x, tr, g);


% Initialize
yt = y(1);

% Initialize first Hk
Hk = H(1,1);

T = length(H(:,1));

% Initialize first Dk
Dk = 1/(Hk'*Hk);

% Compute iniital estimate of theta_k
theta_k = Dk*Hk'*y;

% Initial covariance of data
Sigma = Dk/var_x;

% First gain
K = 0.5*Dk*Hk;

% Initial theta update
theta_k = Dk*Hk*yt;

k = 1;
Hk = H(1:2, k);
idx = 1:k;
% Time loop
for t = 2:T-1

    % Compute probabilities 
    % yt becomes available
    yt(t) = y(t);
    if (t == 2)
        yt = yt';
    end
    k = length(theta_k);

    % For Mk - to stay
    log_pk = -0.5/var_x * sum((yt - Hk*theta_k).^2);

    % For M k+1
    % First update theta_k
    [theta_k1, Dk1] = ols_updates(yt, H, k, t, Dk, theta_k);
    log_pk1 = -0.5/var_x*sum((yt - H(1:t, [idx, k+1])*theta_k1).^2);

    % Store Dks
    Dk_store{1} = Dk;
    Dk_store{2} = Dk1;
    % Store thetas
    theta_store = {theta_k, theta_k1};

    % 
    if (k > 2)
        % For loop for all ks 
        theta_clean = theta_k;
        for tk = 1:k
            % Downdate theta_k
            [theta_temp, Dk_update] = ols_downdates(theta_k, Dk, tk);
            log_ptk(tk) = -0.5/var_x*sum((yt - H(1:t, setdiff(idx, tk))*theta_temp).^2);
            theta_store{tk + 2} = theta_temp;
            Dk_store{tk + 2} = Dk_update;
        end
    
        % Normalize probs
        probs = exp( [log_pk log_pk1 log_ptk] - max( [log_pk log_pk1 log_ptk] ) );

        % Choose model index
        ki = datasample(1:(k+2), 1, 'Weights', probs);
        clear log_ptk

    else

        % Theta store
        theta_store = {theta_k, theta_k1};
    
        % Normalize probs
        probs = exp( [log_pk log_pk1] - max( [log_pk log_pk1] ) );

        % Choose model index
        ki = datasample([1,2], 1, 'Weights', probs);

    end

    % Choose theta
    theta_k = theta_store{ki};
    Dk = Dk_store{ki};
    Sigma = var_x*Dk;

    % Update ki and time in H
    if (ki > 2)
        idx = setdiff(idx, ki-2);
        Hk = H(1:t, idx);
    elseif (ki == 2)
        idx = [idx k+1];
        Hk = H(1:t, idx);
    end
    k = length(theta_k);

    % Update time parameters

    % Update gain
    K = Sigma*Hk(t,:)'/(var_x + Hk(t,:)*Sigma*Hk(t,:)');

    % Update estimate
    theta_k = theta_k + K*(yt(t) - Hk(t,:)*theta_k);

    % Update covariance
    Sigma = (eye(k) - K*Hk(t,:))*Sigma;

    % Update Hk
    Hk = H(1:t+1, idx);

    

end










