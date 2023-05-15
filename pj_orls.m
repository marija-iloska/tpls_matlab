function [theta_k, Hk, k_store, k_mode] = pj_orls(y, H, dy, var_y)

% Store
H_true = H;
T = length(H(:,1));

% Initialize model order 
k = 1;

% Initialize using t data points
t = 2;
[J, theta_k, Dk, Hk, Sigma] = initialize(y, H, t, k, var_y);
Jup_track(1:2) = 0;
Jdown_track(1:2) = 0;

% Start time loop
for t = 3:T-1

    % JUMP UP +
    J_up = Inf;
    if (dy > k)
        [theta_up, H_up, J_up, Sigma_up, Dk_up, k_up] = jump_up(y, dy, k, Dk, theta_k, J, H, t, var_y) ;
        Jup_track(t) = J_up;
    else
        Jup_track(t) = Jup_track(t-1);
    end

    % JUMP DOWN -
   J_down = Inf;
   if (k > 1)
       [theta_down, H_down, J_down, Sigma_down, Dk_down, k_down] = jump_down(y, k, Dk, theta_k, J, H, t, var_y);
       Jdown_track(t) = J_down;
   else
       Jdown_track(t) = Jdown_track(t-1);
   end

   % STAY SAME
   [theta_stay, Sigma_stay, J_stay] = time_update(y, Hk, t, theta_k, var_y, Dk, J);
   Dk_stay = Sigma_stay/var_y;

   % Compute weights based on errors
   J_track(t) = J_stay;
   Js = [J_stay, J_up, J_down];
   Ws = exp(-(Js-min(Js)));
   Ws = Ws./sum(Ws);
   minJ = datasample(1:3, 1, 'Weights', Ws);

  if (minJ == 3)
      k = k_down;
      H = H_down;
      Hk = H(1:t+1, 1:k);
      theta_k = theta_down;
      J = J_down;
      Sigma = Sigma_down;
      Dk = Dk_down;
  elseif (minJ == 2)
      H = H_up;
      k = k_up;
      Hk = H(1:t+1, 1:k);
      theta_k = theta_up;
      J = J_up;
      Sigma = Sigma_up;
      Dk = Dk_up;
  else
      theta_k = theta_stay;
      Sigma = Sigma_stay;
      J = J_stay;
      Dk = Sigma/var_y;
      Hk = H(1:t+1, 1:k);
  end
  k_store(t) = k;

end

% Get mode
k_mode = [mode(k_store), mode(k_store(end - 0.25*T : end))];

end