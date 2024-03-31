function [M_max, theta_RJ, plot_stats, J] = rj_mcmc(y, H, n, Ns, Nb, idx, var_y)

% Get length of data
N = length(y);
K = length(H(1,:));

% Partition data
ye = y(1:n);
He = H(1:n, :);

% Track of indices in use
S = 1:K;
Mj = datasample(S, 1);
S = setdiff(S, Mj);

% Current model params
pj = 1;
s = 1;

correct = [];
incorrect = [];

% Start sweep
while s <= Ns

    if (isempty(S))
        % Record obtained model
        M{s} = [Mj, zeros(1, K-length(Mj))];
        pj = length(Mj);

        % Reset S and start a new sweep
        S = 1:K;
        s = s + 1;
        correct(end+1) = sum(ismember(Mj, idx));
        incorrect(end+1) = length(Mj) - correct(end);

        [Jtemp, ~] = true_PE(y, H, n, length(y), Mj, var_y);
        J(end+1) = Jtemp(end);
    end
    
    % Choose a random predictor
    k = datasample(S, 1);

    % Check if in current model
    if (ismember(k, Mj))

        % Sort
        Mk = sort( setdiff(Mj, k), "ascend");

        % Compute predictive density
        pk = pj - 1;
        [dens_k] = pdf_compute(N, n, pk, Mk, He, H, y, ye);
        [dens_j] = pdf_compute(N, n, pj, Mj, He, H, y, ye);

        % Accept / Reject
        if log(rand) < dens_k - dens_j           
            % New model
            Mj = Mk;
        end

    else
        % Sort
        Mk = sort( [Mj, k], "ascend");

        % Compute predictive density
        pk = pj + 1;
        [dens_k] = pdf_compute(N, n, pk, Mk, He, H, y, ye);
        [dens_j] = pdf_compute(N, n, pj, Mj, He, H, y, ye);

        % Accept / Reject
        if log(rand) < dens_k - dens_j           
            % New model
            Mj = Mk;
        end

    end

    % Update S
    S = setdiff(S, k);

end

% Apply burn-in
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
    for s = 1:(Ns - Nb)
        if (sum(Mk == M_burn{s}) == K)
            count(m) = count(m) + 1;
        end     
    end

end

% Sort count
[~, idx_sorted] = sort(count, 'descend');
models_sorted = models(idx_sorted,:);
M_max = nonzeros(models_sorted(1,:))';
plot_stats = {correct, incorrect};

% Get theta estimate of final winning model
theta_RJ = inv(H(:,M_max)'*H(:,M_max))*H(:,M_max)'*y;

% Compute predictive error for winning model
[J, ~] = true_PE(y, H, t0, T, M_max, var_y);


end