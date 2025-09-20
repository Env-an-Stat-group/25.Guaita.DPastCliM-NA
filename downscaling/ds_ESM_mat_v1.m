function [EOds_hat_mat, PI_mat] = ...
    ds_ESM_mat_v1(tgt_ESM_mth, tgt_lon, tgt_lat, time_ESM_mth, lm_list, mu_gO_local, O_t_local, flag_cal, ...
    eof_all, n_pc_array, metaTable, name_var)
%DS_ESM_MAT_V1  Downscale ESM data to station level using EOFs + regression.
%
% This function performs statistical downscaling of Earth System Model (ESM)
% data to station observations using Empirical Orthogonal Functions (EOFs)
% as predictors and station-specific linear regression models.
%
% Inputs:
%   - tgt_ESM_mth : Monthly ESM data (gridPoints x timesteps).
%   - tgt_lon     : Vector of longitudes for ESM grid points.
%   - tgt_lat     : Vector of latitudes for ESM grid points.
%   - time_ESM_mth: Vector of time steps (months).
%   - lm_list     : Cell array of station-level linear regression models.
%   - mu_gO_local : Table with station-specific means of transformed obs (for inverse transform).
%   - O_t_local   : Table with station-specific translation terms (for inverse transform).
%   - flag_cal    : Calibration flag (used in trend computation).
%   - eof_all     : Matrix of EOF loadings (gridPoints x PCs).
%   - n_pc_array  : Indices of selected PCs to project onto.
%   - metaTable   : Metadata table with station info (IDs, lon, lat, calibration flag).
%   - name_var    : Variable name ('tas' for temperature, 'pr' for precipitation).
%
% Outputs:
%   - EOds_hat_mat : Expected downscaled values at stations (stations x timesteps).
%   - PI_mat       : Prediction intervals (stations x timesteps x 2).
%
% Notes:
%   * Each station has its own regression model on selected PCs.
%   * Prediction intervals are estimated via Monte Carlo sampling of LR residuals.
%   * Inverse transform depends on 'name_var' (identity for tas, log/exp for pr).

%% prepare metadata
% Combine local station data (mean values and scaling factors) with metaTable
metaTable = innerjoin(metaTable, mu_gO_local, 'Keys', 'ID');  % Join mu_gO_local to metaTable by station ID
metaTable = innerjoin(metaTable, O_t_local, 'Keys', 'ID');  % Join O_t_local to metaTable by station ID

%% calculate the moving average to tranfer the trend given by the ESM to the downscaled Y

[mu_gO_trend_local,mu_mov_ESM_trend] = compute_mu_gO_trend_local_v1(tgt_ESM_mth, tgt_lon, tgt_lat, time_ESM_mth, flag_cal, name_var, metaTable);

%% Prepare predictands

% Select the EOFs based on the indices in n_pc_array
eof_subset = eof_all(:, n_pc_array);  % Subset the EOFs based on the selected PCs

% Normalize the target ESM data by subtracting the mean for each grid point
SX = tgt_ESM_mth - mu_mov_ESM_trend;  % ESM anomalies (centered data)
X = (SX' * eof_subset);  % Project the anomalies onto the selected EOFs

% Initialize matrices for the downscaled values, expected values, and prediction intervals
EOds_hat_mat = nan(height(metaTable),size(X, 1));  % Expected (mean) downscaled values at stations
PI_mat = nan(height(metaTable),size(X, 1),2);  % Expected (mean) downscaled values at stations

%% Get the mean prediction and the Prediction intervals

% Loop over each station in the metadata table
for i_ID = 1:height(metaTable)
    if metaTable.flag_cal(i_ID)
        lm_tmp = lm_list{i_ID};
        if not(isempty(lm_tmp))
            % Predict expected value on regression scale
            [ESgO_hat, ~] = predict(lm_tmp, X);  

            % Inverse transform expected value
            sigma_hat = sqrt(lm_tmp.MSE);  % regression standard deviation

            % get prediction intervals with montecarlo
            nSim=1000;
            nTime = size(X,1);
            % Create simulated noise
            u_sim = sigma_hat * randn(nTime, nSim);
            switch name_var
                case 'tas'
                    Ohat_sim = ESgO_hat + mu_gO_trend_local(i_ID,:)' + u_sim - metaTable.O_t(i_ID);
                    % Compute PI quantiles along the simulation dimension
                    PI_inf_tmp = quantile(Ohat_sim', 0.025)';  
                    PI_sup_tmp = quantile(Ohat_sim', 0.975)';  
                    % invert ESgO_hat and get the predicted mean
                    EOds_hat_tmp = inverseTransformMean(ESgO_hat, mu_gO_trend_local(i_ID,:)', metaTable.O_t(i_ID), sigma_hat^2, name_var);
                case 'pr'
                    Ohat_sim = exp(ESgO_hat + mu_gO_trend_local(i_ID,:)' + u_sim) - metaTable.O_t(i_ID);
                    Ohat_sim(Ohat_sim<0)=0;
                    % Compute PI quantiles along the simulation dimension
                    PI_inf_tmp = quantile(Ohat_sim', 0.025)';  
                    PI_sup_tmp = quantile(Ohat_sim', 0.975)';  
                    % invert ESgO_hat and get the predicted mean
                    EOds_hat_tmp = inverseTransformMean(ESgO_hat, mu_gO_trend_local(i_ID,:)', metaTable.O_t(i_ID), sigma_hat^2, name_var);
                    % correct negative values
                    PI_inf_tmp(PI_inf_tmp<0)=0;
                    PI_sup_tmp(PI_sup_tmp<0)=0;
                    EOds_hat_tmp(EOds_hat_tmp<0)=0;
            end
            % store results
            PI_mat(i_ID, :, 1) = PI_inf_tmp;  
            PI_mat(i_ID, :, 2) = PI_sup_tmp;              
            EOds_hat_mat(i_ID, :) = EOds_hat_tmp;
        end
    end
end

end
