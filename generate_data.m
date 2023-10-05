function [y, H, theta] = generate_data(T, dy, var_h, rt,  p_s, var_y)


% Choose random indices to be 0s
j = datasample(1:dy, p_s, 'replace', false);

% Generate random theta in the range between -rt, rt
theta = normrnd(0, rt, dy, 1);
%unifrnd(-rt, rt, dy, 1);

% Set chosen indices to 0s
theta(j) = 0;

% Create basis functions and data
H = sin(normrnd(0, var_h, T, dy));
%H = normrnd(0, var_h, T, dy);

y = H*theta;
y = y + mvnrnd(zeros(T,1), var_y*eye(T))';



end
