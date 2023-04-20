clear all
close all
clc

% Settings
var_y = 0.001; % Variance
p_s = 0;   % Sparsity percent
dx = 50;      % System dimension
T = 100;     % Time series length
r = 0.5;     % Range of input data H
rt = 3;      % Range of theta


% Current dimension dtheta
k = 1;

% n - number of jumps we want
n1 = 20;
n = dx - (k + n1);


%Create data
[y, H, theta] = generate_data(T, dx, r,rt,  p_s, var_y);

H_true = H;


% Initialize first Hk
Hk = H(:,1);
T = length(H(:,1));

% Initialize first Dk
Dk = 1/(Hk'*Hk);

% Compute iniital estimate of theta_k
theta_k = Dk*Hk'*y;

% Initial covariance of data
Sigma = Dk/var_y;


% Initial theta update
theta_k = Dk*Hk'*y;




% ORDER UPDATE by + 1
% j is an index that iterates through the next unused basis functions in H
Dk_temp = Dk;
theta_temp = theta_k;
H_temp = H_true;
% Update n steps
tic
for j = 1:n1
    % Update current theta by jth basis function
    [theta_temp, Dk_temp, H_temp] = ols_updates_n(y, H, k, j, T, Dk_temp, theta_temp);
                                                 
end
toc

k = n1+1;
Dk = Dk_temp;
theta_k = theta_temp;
for j = 1:n
    % Update current theta by jth basis function
    [theta_temp, Dk_temp, H_temp] = ols_updates_n(y, H, k, j, T, Dk_temp, theta_temp);
                                                 
end

% Store res
Dkn_seq = Dk_temp;
Hkn_seq = H_temp;
theta_kn_seq = theta_temp;


% % ORDER UPDATE by + n
tic
% Obtain D22kn
% - need Pk
% - need Hn
% - need inverse algorithm - first try with normal inverse low dim
Hk = H(:, 1:k);
Pk = eye(T) - Hk*Dk*Hk';
Hn = H(:, k+1:k+n);
D22kn = inv(Hn'*Pk*Hn);
D22test = inverse_D22(Pk, Hn, n);



% Obtain D12kn and D21kn
% - need Hn
% - need Hk
% - need Dk
% - need D22kn
xd = Dk*Hk'*Hn;
D12kn = -xd*D22kn;
D21kn = D12kn';


% Obtain D11kn
% - need Hn
% - need Hk
% - need Dk
% - need D22kn
D11kn = Dk - D12kn*xd';

% Construct Dkn
Dkn = [D11kn, D12kn; D21kn, D22kn];


% Obtain theta_kn
% - need Hn
% - need D12kn
% - need D22kn
% - need Pk
% - need y
theta_kn = [theta_k + D12kn*Hn'*Pk*y;  D22kn*Hn'*Pk*y ];
toc


Dsub = [D11kn, D12kn; D21kn, D22test];
theta_sub = [theta_k + D12kn*Hn'*Pk*y;  D22test*Hn'*Pk*y ];



sum(abs(y - H*theta_kn_seq))/T
sum(abs(y - H*theta_kn))/T
sum(abs(y - H*theta_sub))/T
sum(abs(y - H*theta))/T





theta_kn_seq'
theta_kn'
theta_sub'
theta'




