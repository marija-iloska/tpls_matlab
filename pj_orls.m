function [theta_store, Hk, k_store, k_mode, models_sorted, count_sorted, idx_orls, J_pred, J_incr] = pj_orls(y, H, dy, var_y, n, Nb, D)

% Store
H_true = H;
T = length(H(:,1));
K = length(H(1,:));

% len = K/D;
% 
% 
% for d = 1:D
%     range{d} = d*len - len + 1 : d*len; 
%     [J(d), theta_d, Dk, Hk, ~] = initialize_D(y, H(:, range{d}), n, var_y);
%     theta_D{d} = theta_d;
%     Hd{d} = Hk;
%     Dd{d} = Dk;
% end
% minD = find(J == min(J));
% theta_k = theta_D{d};
% Hk = Hd{d};
% Dk = Dd{d};
% H = H(: , [range{minD}, setdiff(1:K, range{minD})]);
% k = len;



% Initialize model order 
k = dy;

%Initialize using t data points
[J, theta_k, Dk, Hk, ~] = initialize(y, H, n, k, var_y);


[~, idx_sort] = sort(theta_k, 'descend');
H = H(:, idx_sort);
k = floor(dy/2);
[J, theta_k, Dk, Hk,~] = initialize(y, H, n, k, var_y);

J_pred = [];
J_incr = J;

M ={};
theta_store = {};

% Start time loop
for t = n+1:T-1

    % JUMP UP +
    J_up = Inf;
    if (dy > k)
        [theta_up, H_up, J_up, Dk_up, k_up] = jump_up(y, dy, k, Dk, theta_k, J, H, t) ;        
    end

    % JUMP DOWN -
   J_down = Inf;
   if (k > 1)
       [theta_down, H_down, J_down, Dk_down, k_down] = jump_down(y, k, Dk, theta_k, J, H, t);
   end

   % STAY SAME
   J_stay =  sum( (y(1:t) - H(1:t, 1:k)*theta_k).^2);
   %J_stay = (y(t) - H(t, 1:k)*theta_k)^2;
   Dk_stay = Dk;

   % Compute weights based on errors
   Js = [J_stay, J_up, J_down];
   minJ = find(Js == min(Js));
   

  if (minJ == 3)
      k = k_down;
      H = H_down;
      theta_k = theta_down;
      J = J_down;
      Dk = Dk_down;
  elseif (minJ == 2)
      H = H_up;
      k = k_up;
      theta_k = theta_up;
      J = J_up;
      Dk = Dk_up;
  else
      theta_k = theta_k;
      J = J_stay;
      Dk = Dk_stay;
  end

  Hk = H(1:t+1, 1:k);
  k_store(t) = k;
  theta_store{end+1} = theta_k;

  % Check which model was selected at time t
  [~, idx_orls] = ismember(Hk(1,:), H_true(1,:));
  M{end+1} = [sort(idx_orls, 'ascend'), zeros(1, dy - length(idx_orls)) ];

  % TIME UPDATE
  J_pred(end+1) = J;
  J_incr(end+1) = J_incr(end) + (y(t) - H(t, 1:k)*theta_k)^2;
  [theta_k, Sigma, ~] = time_update(y, Hk, t, theta_k, var_y, Dk, J); 
  Dk = Sigma/var_y;


end

% Apply models
M_burn = M;
M_burn(1:Nb) = [];

% Find unique models
models = unique(cell2mat(M_burn'), 'rows');
Nm = length(models(:,1));
count = zeros(1,Nm);

% For each unique model
for m = 1:Nm

    % Current model to check for
    Mk = models(m,:);

    % Check all sweeps
    for s = 1: length(M) - Nb
        if (sum(Mk == M_burn{s}) == dy)
            count(m) = count(m) + 1;
        end     
    end

end

% Sort count and find max
[count_sorted, idx_sorted] = sort(count, 'descend');
models_sorted = models(idx_sorted,:);
idx_orls = nonzeros(models_sorted(1,:))';


% Get mode
k_mode = [mode(k_store), mode(k_store(end - round(0.25*T) : end))];



end