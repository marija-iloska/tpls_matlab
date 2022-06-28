function [A, C, x, y] = generate_states(T, dx, r, p_s, var_x, var_y, tr, obs)


% Initialize coefficient and aadjacency matrices
C = unifrnd(-r, r , dx, dx);
A = ones(dx, dx);


for j = 1 : dx
    idx = datasample(1:dx, round(p_s*dx));
    A(j,idx) = 0;
end


C = C.*A;


% Generate the data
x(:,1) = 0.5*rand(dx, 1);
y(:,1) = x(:,1) + mvnrnd(zeros(dx,1), var_y*eye(dx))';

for t = 2:T
    x(:,t) = tr(C, x(:,t-1)) + mvnrnd(zeros(dx,1), var_x*eye(dx))';
    y(:,t) = obs(1, x(:,t)) + mvnrnd(zeros(dx,1), var_y*eye(dx))';       
end

end
