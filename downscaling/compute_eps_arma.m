function eps_arma = compute_eps_arma(metaTable, list_ARMA, time_ESM, start_time, p, q)
% compute_epsilon: Simulates ARMA residuals for stations using ESM grid
%
% Inputs:
% - metaTable: Table with metadata including flags and coordinates
% - list_ARMA: Cell arrays of ARMA models for each grid point
% - time_ESM: Time vector (general)
% - start_time: Time index to start forecasting from
% - p, q: ARMA(p,q) model orders
%
% Outputs:
% - epsilon: Simulated residuals using list_ARMA
% - epsilon_lr: Simulated residuals using list_ARMA_lr

nStations = height(metaTable);
nTime = numel(time_ESM);

% Initialize outputs
eps_arma = zeros(nStations, nTime);
eps_arma(:,1:(start_time-1)) = 0;

% Loop over stations
metaTable_flag_cal = metaTable.flag_cal;
parfor i_ID = 1:nStations
    if metaTable_flag_cal(i_ID)
        model_pcr = list_ARMA{i_ID};

        % Extract ARMA model variances
        sigma = sqrt(model_pcr.Variance);
        
        % Simulate white noise upfront
        eta = sigma * randn(1, nTime);

        % Assuming AR and MA coefficients:
        phi = [model_pcr.AR{:}];            % p-length vector
        theta = [model_pcr.MA{:}];          % q-length vector
        
        % Initialization
        eps_tmp = zeros(1, nTime);
        
        for t = start_time:nTime
            ar_term = 0;
            ma_term = 0;
            for j = 1:p
                if t-j > 0
                    ar_term = ar_term + phi(j) * eps_tmp(t-j);
                end
            end
            for j = 1:q
                if t-j > 0
                    ma_term = ma_term + theta(j) * eta(t-j);
                end
            end
            eps_tmp(t) = ar_term + ma_term + eta(t);  % Proper ARMA simulation
        end
        eps_arma(i_ID, :) = eps_tmp;    
    end
end

end
