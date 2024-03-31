function [y, H, theta] = fix_data(T, dy, var_h, rt,  idx)


% Generate random theta in the range between -rt, rt
theta = normrnd(0, rt, dy, 1);

% Set chosen indices to 0s
theta(idx) = 0;

% Create basis functions and data
H = normrnd(0, var_h, T, dy);

y = H*theta;


end
