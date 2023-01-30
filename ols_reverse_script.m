clear all
close all
clc

% Order - recursive least squares example

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



% Call ORLS
epsilon = 0.01;
[theta_k, Dk, Jk, error_store] = ols(y, H, epsilon, dx);




% Minimum possible error;
Jmin  = sum( (y - H*theta).^2 )


% Choose one to remove ( Think about different ways of choosing)
theta_clean = theta_k;

% Fix this
store = 0;

Dk_update = Dk;

for i = 1:5 

[~,  min_k] = sort( abs(theta_clean) ) ;
min_k = min_k(1);
min0 = find(theta_k == theta_clean(min_k));
store(i) = min0;


[theta_k, theta_clean, Dk_update] = ols_reverse(theta_clean, Dk_update, min_k, store, dx);

end

J = sum( (y - H*theta_k).^2)

a_est = (theta_k ~=0);

wrong = sum(a_est ~= a);




