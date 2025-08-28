function [u_mat] = inverse_SEM_season(eps_ARMA, W_mam, W_jja, W_son, W_djf, lambda)

% Validate dimensions
[nStations, nTime] = size(eps_ARMA);

% Identity matrix
I = eye(nStations);

% Precompute A matrices for each season
A_mam = I - lambda * W_mam;
A_jja = I - lambda * W_jja;
A_son = I - lambda * W_son;
A_djf = I - lambda * W_djf;

% Output matrix
u_mat = NaN(nStations, nTime);

% Loop over time
for t = 1:nTime
    month_t = mod(t-1,12)+1;

    % Determine season
    if ismember(month_t, [3, 4, 5])
        A = A_mam;
    elseif ismember(month_t, [6, 7, 8])
        A = A_jja;
    elseif ismember(month_t, [9, 10, 11])
        A = A_son;
    else % DJF
        A = A_djf;
    end

    % Subset and solve
    eps_t = eps_ARMA(:, t);
    valid = ~isnan(eps_t);
    A_sub = A(valid, valid);
    eps_sub = eps_t(valid);
    
    % Backsolve
    u_mat(valid, t) = A_sub \ eps_sub;
end

end
