function [dsEValue_mat_lr, PI_mat_lr] = ...
    ds_ESM_mat_lr_v1(tgt_ESM_mth, tgt_lon, tgt_lat, time_ESM_mth, lm_list_lr, mu_Y_local, Y_t_local, flag_cal,...
    metaTable, name_var)
% This function performs local downscaling of ESM data to station-level data
% using linear regression models. It produces both the predicted values and
% prediction intervals for each station and outputs downscaled maps using
% Gaussian Process Regression (GPR).
%
% Inputs:
%   - tgt_ESM_mth: Matrix of ESM data (grid points x timesteps)
%   - tgt_lon: Vector of longitudes corresponding to ESM grid points
%   - tgt_lat: Vector of latitudes corresponding to ESM grid points
%   - time_mth: Vector of time steps (months)
%   - lm_list_lr: Cell array of linear regression models for each station
%   - mu_Y_local: Table of mean values for each station (for inverse transformation)
%   - Y_t_local: Table of scaling factors (for inverse transformation)
%   - metaTable: Metadata table containing station info (IDs, locations, calibration flags)
%   - sigmanoise: Noise level used in GPR to regularize the model
%   - opt_gpr: Options for the Gaussian Process Regression model (e.g., kernel function)
%   - name_var: The name of the target variable for inverse transformations and labeling
%
% Outputs:
%   - dsValue_map_lr: Downscaled values at each grid point over time (maps)
%   - dsEValue_map_lr: Expected downscaled values at each grid point over time (maps)
%   - dsPIsup_map_lr: Upper bounds of prediction intervals at each grid point
%   - dsPIinf_map_lr: Lower bounds of prediction intervals at each grid point
%   - dsValue_mat_lr: Downscaled values at stations (stations x timesteps)
%   - dsEValue_mat_lr: Expected downscaled values at stations (stations x timesteps)

%% Data Preprocessing: Prepare metadata and projection
% Join mu_Y_local and Y_t_local data to metaTable for inverse transformations
metaTable = innerjoin(metaTable, mu_Y_local, 'Keys', 'ID');  % Add mean values for each station
metaTable = innerjoin(metaTable, Y_t_local, 'Keys', 'ID');  % Add scaling factors for each station

%% calculate the moving average to tranfer the trend given by the ESM to the downscaled Y
[mu_Y_trend_local,mu_mov_ESM_trend] = compute_mu_Y_trend_local_v1(tgt_ESM_mth, tgt_lon, tgt_lat, time_ESM_mth, flag_cal, name_var, metaTable);

%% Data Preprocessing: Prepare predictand

% Normalize the target ESM data by subtracting the mean (to focus on anomalies)
SX = tgt_ESM_mth - mu_mov_ESM_trend;  % ESM anomalies (centered data)

% Initialize output matrices for downscaled values, prediction intervals, and confidence intervals
dsEValue_mat_lr = nan(height(metaTable), size(SX, 2));  % Stores expected downscaled values at stations
PI_mat_lr = nan(height(metaTable), size(SX, 2),2);  % Stores expected downscaled values at stations

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
            [ESY_hat, ~] = predict(lm_tmp, SX(min_idx,:)');  % Predict the expected value and confidence intervals
            [~, SPI] = predict(lm_tmp, SX(min_idx,:)', 'prediction', 'observation');  % Predict prediction intervals (SPI)
                
            % Inverse transformation of predicted values to match the original variable scale
            dsEValue_mat_lr(i_ID,:) = inversetransform(ESY_hat, mu_Y_trend_local(i_ID,:)', metaTable.Y_t(i_ID), name_var);  % Expected downscaled values
            PI_mat_lr(i_ID, :, 1) = inversetransform(SPI(:, 1), mu_Y_trend_local(i_ID,:)', metaTable.Y_t(i_ID), name_var);  % Lower bound of prediction interval
            PI_mat_lr(i_ID, :, 2) = inversetransform(SPI(:, 2), mu_Y_trend_local(i_ID,:)', metaTable.Y_t(i_ID), name_var);  % Upper bound of prediction interval
        end
    end
end

% fix precipitation impossible values
if strcmp(name_var,'pr')
    flag_neg = dsEValue_mat_lr<0;
    dsEValue_mat_lr(flag_neg)=prctile(dsEValue_mat_lr,25,'all')*rand(1,sum(flag_neg,'all'));
end

end
