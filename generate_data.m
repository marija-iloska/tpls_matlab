function [y, H, theta, a] = generate_data(T, dx, r,rt,  p_s, var_x, tr, g)


% Initialize coefficient and adjacency matrices
C = unifrnd(-r, r , dx, dx);
% A = ones(dx, dx);


% for j = 1 : dx
%     idx = datasample(1:dx, round(p_s*dx));
%     A(j,idx) = 0;
% end


% C = C.*A;


% Generate the data
%x(:,1) = 0.5*rand(dx, 1);

% for t = 2:T
%     x(:,t) = tr(C, x(:,t-1));       
% end
% 
% H = x(:, 1:T-1)';
% 
j = datasample(1:dx, round(p_s*dx));
% theta = C(j,:)';

theta = unifrnd(-rt,rt, dx, 1);

theta(j) = 0;

H = unifrnd(-r, r, T-1, dx);

y = H*theta + mvnrnd(zeros(T-1,1), var_x*eye(T-1))';

a = (theta ~=0);


end
