clear all
close all
clc

% Order - recursive least squares example

% Settings
var_x = 0.01;
g = @(x) x;
p_s = 0.3;
dx = 20;
T = 50;
r = 0.5;  % Deciding factor
rt = 5;

%SSM
tr = @(coeff, states) coeff*g(states);

%Create data
[y, H, theta, a] = generate_data(T, dx, r, rt, p_s, var_x, tr, g);

% T = 99;
% t = 0:1:T-1;
% 
% theta = [1.2, 0.03, 2.8, -4]';
% 
% dx = length(theta);
% 
% H = ones(T, dx);
% 
% H(:,2) = t';
% H(:,3) = t'.^2;
% H(:,4) = sin(t);
% 
% y = H*theta +  mvnrnd(zeros(T,1), var_x*eye(T))';



% Initialize first Hk
Hk = H(:,1);

T = length(H(:,1));

% Initialize first Dk
Dk = 1/(Hk'*Hk);


% Compute iniital estimate of theta_k
theta_k = Dk*Hk'*y;


% Compute first Jmink 
Jk = sum( (y - Hk*theta_k).^2 );
error_store = Jk;

epsilon = 0.01;
k = 1;
count = zeros(1,dx);
while (Jk > epsilon)

    % Order of model update
    k = k + 1;

    % Compute projection matrix Pk
    Pk_norm = eye(T) - Hk*Dk*Hk';

    % Take the new observation column h(k+1)
    hk = H(:,k);

    % Reuseable terms
    temp = hk'*Pk_norm;
    xd = temp*y;
    d = Dk*Hk'*hk;

    % Compute terms of D(k+1)
    DK22 = 1 / ( temp*hk );
    DK12 = - d*DK22;
    DK21 = DK12';
    DK11 = Dk + DK22*(d*d');

    % Update D(k) to D(k+1)
    Dk = [ DK11, DK12 ;  DK21, DK22 ];

    % Update theta_k
    theta_k = [ theta_k - d*xd*DK22;    xd*DK22 ];

    % Update Hk
    Hk = [Hk, hk];
    Jk = Jk - xd^2*DK22;
    error_store(k) = Jk;

    if (abs(Jk - error_store(k-1)) <= 1)
        theta_k(end) = 0;
        count(k) = count(k) + 1;

    end

    % Max observations
    if (k == dx)
        break
    end


end


% Minimum possible error;
Jmin  = sum( (y - H*theta).^2 );

plot(error_store)
de = error_store(1:end-1) - error_store(2:end);
%plot(de)


theta'

theta_k'

%
