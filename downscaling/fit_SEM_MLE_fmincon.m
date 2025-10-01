function [lambda_hat, sigma2_hat, eps_mat, W_best, threshold_best] = fit_SEM_MLE_fmincon(residual_mat, metaTable, proj)

[y_proj, x_proj] = projfwd(proj, metaTable.lat, metaTable.lon);
coords = [x_proj, y_proj];
D = squareform(pdist(coords));
nStations = size(residual_mat,1);
nTime = size(residual_mat,2);

% Keep only calibrated stations
residual_mat(~metaTable.flag_cal,:) = NaN;
validStations = find(metaTable.flag_cal);

% Try multiple thresholds to pick W (same as before)
lb_distance = 25000;
ub_distance = 100000;
thresholds = linspace(lb_distance, ub_distance, 4);

best_LL = -Inf;

penalty_weight = 0.1;

for th_idx = 1:length(thresholds)
    th = thresholds(th_idx);

    % Gaussian decay spatial weights
    h = th;  % You can tune this percentile as needed
    W = exp(-(D.^2) / (2 * h^2));
    W(eye(nStations) == 1) = 0;  % No self-weight
    
    % Normalize rows
    row_sums = sum(W, 2);
    W = W ./ row_sums;
    W(isnan(W)) = 0;

    eigvals = eig(W);
    rho = max(abs(eigvals));
    if rho >= 1
        W = W / (rho + 1e-2);  % small buffer below 1
    end

    % Objective function: negative log-likelihood profiled over sigma2
    negLogLik = @(lambda) profileNegLogLik(lambda, residual_mat, W, validStations, nTime, penalty_weight);

    % Bounds and options
    lb = 0;
    ub = 0.9;
    options = optimoptions('fmincon', 'Display','off', 'Algorithm','interior-point');

    % Optimize lambda
    [lambda_opt, fval] = fmincon(negLogLik, 0.5, [], [], [], [], lb, ub, [], options);

    if -fval > best_LL
        best_LL = -fval;
        lambda_hat = lambda_opt;
        W_best = W;
        threshold_best = th;
    end
end

% Compute final eps_mat and sigma2
A_best = eye(nStations) - lambda_hat * W_best;
eps_mat = NaN(nStations, nTime);
for t = 1:nTime
    res_t = residual_mat(:,t);
    valid = ~isnan(res_t);
    A_sub = A_best(valid, valid);
    eps_mat(valid,t) = A_sub * res_t(valid);
end

sigma2_hat = nanmean(eps_mat.^2, 2);

end

function nLL = profileNegLogLik(lambda, residual_mat, W, validStations, nTime, alpha)
    nStations = length(validStations);
    A = eye(size(W,1)) - lambda * W;

    eps_mat = NaN(size(residual_mat));
    for t = 1:nTime
        res_t = residual_mat(:,t);
        valid = ~isnan(res_t);
        A_sub = A(valid, valid);
        eps_mat(valid,t) = A_sub * res_t(valid);
    end

    sigma2 = nanmean(eps_mat.^2, 2);
    if any(sigma2 <= 0)
        nLL = Inf; % Invalid solution
        return;
    end

    % Compute logdet(A) on calibrated stations
    A_sub = A(validStations, validStations);
    % Use logdet via chol if possible:
    [~, p] = chol(A_sub);
    if p == 0
        U = chol(A_sub);
        logdetA = 2*sum(log(diag(U)));
    else
        % fallback:
        logdetA = log(abs(det(A_sub)) + eps);
    end

    term1 = - (nTime/2)*sum(log(sigma2(validStations)));
    term2 = nTime * logdetA;

    % Quadratic form
    quad = 0;
    for t = 1:nTime
        et = eps_mat(validStations,t);
        quad = quad + sum((et.^2) ./ sigma2(validStations));
    end
    term3 = -0.5 * quad;

    logLik = term1 + term2 + term3;

    % Ridge-like penalty for lambda near 1
    % penalty = alpha * lambda^2 / (1 - lambda + eps);
    penalty = -alpha * log(1 - lambda + eps);

    nLL = -(logLik - penalty); % Penalized negative log-likelihood
end
