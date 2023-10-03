function [theta_store, Hk, k_store, k_mode, models_sorted, count_sorted, idx_orls, J_pred] = pj_orls(y, H, dy, var_y, n, Nb)

% Store
H_true = H;
T = length(H(:,1));

% Initialize model order 
k = dy;

% Initialize using t data points
[J, theta_k, Dk, Hk, ~] = initialize(y, H, n, k, var_y);
Jup_track(1:n) = 0;
Jdown_track(1:n) = 0;

[~, idx_sort] = sort(theta_k, 'descend');
H = H(:, idx_sort);
k = floor(dy/2);
[J, theta_k, Dk, Hk, ~] = initialize(y, H, n, k, var_y);
J_pred = J;

% Start time loop
for t = n+1:T-1

    % JUMP UP +
    J_up = Inf;
    if (dy > k)
        [theta_up, H_up, J_up, Dk_up, k_up] = jump_up(y, dy, k, Dk, theta_k, J, H, t) ;
        if (isinf(J_up) == 1)
            Jup_track(t) = nan;
        else
            Jup_track(t) = J_up;
        end
    else
        Jup_track(t) = Jup_track(t-1);
    end

    % JUMP DOWN -
   J_down = Inf;
   if (k > 1)
       [theta_down, H_down, J_down, Dk_down, k_down] = jump_down(y, k, Dk, theta_k, J, H, t, var_y);
       if (isinf(J_down) == 1)
            Jdown_track(t)= nan;
       else
           Jdown_track(t) = J_down;
       end
   else
       Jdown_track(t) = Jdown_track(t-1);
   end

   % STAY SAME
   J_stay =  J +  (y(t) - H(t, 1:k)*theta_k)^2;
   Dk_stay = Dk;

   % Compute weights based on errors
   J_track(t) = J_stay;
   Js = [J_stay, J_up, J_down];
   Ws = exp(-(Js-min(Js)));
   Ws = Ws./sum(Ws);
   %minJ = datasample(1:3, 1, 'Weights', Ws);
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
  theta_store{t} = theta_k;

  % Check which model was selected at time t
  [~, idx_orls] = ismember(Hk(1,:), H_true(1,:));
  M{t-2} = [sort(idx_orls, 'ascend'), zeros(1, dy - length(idx_orls)) ];

  % TIME UPDATE
  [~, Sigma, ~] = time_update(y, Hk, t, theta_k, var_y, Dk, J);
  Dk = Sigma/var_y;
  J_pred(t) = J;


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