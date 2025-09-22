function [u_mat] = inverse_SEM(eps_ARMA, W, lambda)

% Validate dimensions
[nStations, nTime] = size(eps_ARMA);

% Identity matrix
I = eye(nStations);

% Precompute A matrices for each season
A = I - lambda * W;

% Output matrix
u_mat=A\eps_ARMA;

end
