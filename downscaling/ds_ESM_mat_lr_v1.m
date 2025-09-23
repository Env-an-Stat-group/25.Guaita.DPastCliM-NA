function [EOds_hat_mat_lr, PI_mat_lr] = ...
    ds_ESM_mat_lr_v1(tgt_ESM_mth, tgt_lon, tgt_lat, time_ESM_mth, lm_list_lr, mu_gO_local, O_t_local, flag_cal,...
    metaTable, name_var)
%DS_ESM_MAT_LR_V1  Local downscaling of ESM data using linear regression.
%
% This function applies station-specific linear regression (LR) models to
% downscale Earth System Model (ESM) output to station-level data. It
% produces both the expected downscaled values and prediction intervals
% (PIs). 
%
% Inputs:
%   - tgt_ESM_mth : Matrix of monthly ESM data (gridPoints x timesteps).
%   - tgt_lon     : Vector of longitudes for ESM grid points.
%   - tgt_lat     : Vector of latitudes for ESM grid points.
%   - time_ESM_mth: Vector of time steps (months).
%   - lm_list_lr  : Cell array of station-level linear regression models.
%   - mu_gO_local : Table with station-specific means of transformed obs (for inverse transform).
%   - O_t_local   : Table with station-specific translation terms (for inverse transform).
%   - flag_cal    : Calibration flag (used in trend computation).
%   - metaTable   : Metadata table containing station info (IDs, lon, lat, calibration flag).
%   - name_var    : Variable name ('tas' for temperature, 'pr' for precipitation).
%
% Outputs:
%   - EOds_hat_mat_lr : Expected downscaled values at stations (stations x timesteps).
%   - PI_mat_lr       : Prediction intervals (stations x timesteps x 2).
%
% Notes:
%   * For each station, the closest ESM grid point is used as predictor.
%   * Prediction intervals are estimated via Monte Carlo sampling of LR residuals.
%   * The inverse transformation depends on 'name_var' (identity for tas, 
%     log/exp transform with shift for pr).

%% Data Preprocessing: Prepare metadata and projection
% Join mu_gO_local and O_t_local data to metaTable for inverse transformations
metaTable = innerjoin(metaTable, mu_gO_local, 'Keys', 'ID');  % Add mean values for each station
metaTable = innerjoin(metaTable, O_t_local, 'Keys', 'ID');  % Add scaling factors for each station

%% calculate the moving average to tranfer the trend given by the ESM to the downscaled Y
[mu_gO_trend_local,~] = compute_mu_gO_trend_local_v1(tgt_ESM_mth, tgt_lon, tgt_lat, time_ESM_mth, flag_cal, name_var, metaTable);

%% Data Preprocessing: Prepare predictand

% define predictand
X = tgt_ESM_mth;  

% Initialize output matrices for downscaled values, prediction intervals, and confidence intervals
EOds_hat_mat_lr = nan(height(metaTable), size(X, 2));  % Stores expected downscaled values at stations
PI_mat_lr = nan(height(metaTable), size(X, 2),2);  % Stores expected downscaled values at stations

%% Downscale ESM at each station using linear regression models

% Create a grid of longitudes and latitudes
[tgt_longrid, tgt_latgrid] = ndgrid(tgt_lon, tgt_lat);

% Loop through each station and apply linear regression to downscale the data
for i_ID = 1:height(metaTable)
    if metaTable.flag_cal(i_ID)  % Only process stations marked for calibration
        lm_tmp = lm_list_lr{i_ID};  % Retrieve the linear regression model for the current station
        
        if not(isempty(lm_tmp))
            % Compute distances to the fixed coordinate (target ESM grid point)
            % Using Euclidean distance here; an alternative could be Haversine distance for geodesic measurements
            distances = sqrt((tgt_longrid - metaTable.lon(i_ID)).^2 + (tgt_latgrid - metaTable.lat(i_ID)).^2);  
            
            % Find the index of the closest ESM grid point to the station
            [~, min_idx] = min(distances(:));  % Find the grid point closest to the station
            
            % Perform predictions using the linear regression model
            % The SX(min_idx,:) is the subset of ESM data corresponding to the closest grid point
            [ESgO_hat, ~] = predict(lm_tmp, X(min_idx,:)');  % Predict the expected value

            % Inverse transform expected value
            sigma_hat = sqrt(lm_tmp.MSE);  % regression standard deviation

            % get prediction intervals with montecarlo
            nSim=1000;
            nTime = size(X,2);

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
            PI_mat_lr(i_ID, :, 1) = PI_inf_tmp;  
            PI_mat_lr(i_ID, :, 2) = PI_sup_tmp;              
            EOds_hat_mat_lr(i_ID, :) = EOds_hat_tmp;            
        end
    end
end

end
