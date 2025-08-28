function [dsEValue_mat, PI_mat] = ...
    ds_ESM_mat_v1(tgt_ESM_mth, tgt_lon, tgt_lat, time_ESM_mth, lm_list, mu_Y_local, Y_t_local, flag_cal, ...
    eof_all, n_pc_array, metaTable, name_var)
% This function performs statistical downscaling of Earth System Model (ESM) data to local stations.
% It combines EOFs (Empirical Orthogonal Functions), regression models, and Gaussian Process Regression (GPR)
% to produce downscaled maps of the target variable along with prediction intervals.
%
% Inputs:
%   - tgt_ESM_mth: Monthly ESM data (grid points x timesteps)
%   - tgt_lon: Vector of longitudes corresponding to ESM grid points
%   - tgt_lat: Vector of latitudes corresponding to ESM grid points
%   - time_mth: Vector of time steps (#months since 0 CE)
%   - lm_list: Cell array of pre-trained linear regression models for each station
%   - mu_Y_local: Table of mean values of the observed variable for each station (for inverse transformation)
%   - Y_t_local: Table of scaling factors (e.g., standard deviation) for each station (for inverse transformation)
%   - flag_cal: flag for calibration years
%   - eof_all: Matrix of all available EOFs (Empirical Orthogonal Functions)
%   - n_pc_array: Array specifying the selected PCs (Principal Components) to project the data onto
%   - metaTable: Metadata table containing information about each station (IDs, location, calibration flags, etc.)
%   - sigmanoise: Noise level to regularize the Gaussian Process Regression (GPR)
%   - opt_gpr: Options for the Gaussian Process Regression model (e.g., kernel function, options)
%   - name_var: Name of the target variable for inverse transformations and labeling
%
% Outputs:
%   - dsValue_map: Downscaled values at each grid point over time
%   - dsEValue_map: Expected downscaled values (mean prediction)
%   - dsPIsup_map: Upper bounds of prediction intervals
%   - dsPIinf_map: Lower bounds of prediction intervals
%   - dsValue_mat: Downscaled values at station locations (stations x timesteps)
%   - dsEValue_mat: Expected downscaled values at station locations (stations x timesteps)

%% prepare metadata
% Combine local station data (mean values and scaling factors) with metaTable
metaTable = innerjoin(metaTable, mu_Y_local, 'Keys', 'ID');  % Join mu_Y_local to metaTable by station ID
metaTable = innerjoin(metaTable, Y_t_local, 'Keys', 'ID');  % Join Y_t_local to metaTable by station ID

%% calculate the moving average to tranfer the trend given by the ESM to the downscaled Y

[mu_Y_trend_local,mu_mov_ESM_trend] = compute_mu_Y_trend_local_v1(tgt_ESM_mth, tgt_lon, tgt_lat, time_ESM_mth, flag_cal, name_var, metaTable);

%% Prepare predictands

% Select the EOFs based on the indices in n_pc_array
eof_subset = eof_all(:, n_pc_array);  % Subset the EOFs based on the selected PCs

% Normalize the target ESM data by subtracting the mean for each grid point
SX = tgt_ESM_mth - mu_mov_ESM_trend;  % ESM anomalies (centered data)
X = (SX' * eof_subset);  % Project the anomalies onto the selected EOFs

% Initialize matrices for the downscaled values, expected values, and prediction intervals
dsEValue_mat = nan(height(metaTable),size(X, 1));  % Expected (mean) downscaled values at stations
PI_mat = nan(height(metaTable),size(X, 1),2);  % Expected (mean) downscaled values at stations

%% Downscale data at each station
% Loop over each station in the metadata table
for i_ID = 1:height(metaTable)
    if metaTable.flag_cal(i_ID)  % Only process stations marked for calibration
        lm_tmp = lm_list{i_ID};  % Retrieve the linear regression model for the current station

        if not(isempty(lm_tmp))
            % Predict the downscaled values from the EOF-projected PCs
            [ESY_hat, ~] = predict(lm_tmp, X);  % Predict the expected value (ESY_hat) and confidence intervals (SCI)
            [~, SPI] = predict(lm_tmp, X, 'prediction', 'observation');  % Predict prediction intervals (SPI)
    
            % Apply inverse transformations to get downscaled values back to the original scale
            dsEValue_mat(i_ID, :) = inversetransform(ESY_hat, mu_Y_trend_local(i_ID,:)', metaTable.Y_t(i_ID), name_var);  % Expected downscaled values
            PI_mat(i_ID, :, 1) = inversetransform(SPI(:, 1), mu_Y_trend_local(i_ID,:)', metaTable.Y_t(i_ID), name_var);  % Lower bound of prediction interval
            PI_mat(i_ID, :, 2) = inversetransform(SPI(:, 2), mu_Y_trend_local(i_ID,:)', metaTable.Y_t(i_ID), name_var);  % Upper bound of prediction interval

            % Apply inverse transformations to get downscaled values back to the original scale
            dsEValue_mat(i_ID, :) = inversetransform(ESY_hat, mu_Y_trend_local(i_ID,:)', metaTable.Y_t(i_ID), name_var);  % Expected downscaled values
        end
    end
end


% fix precipitation impossible values
if strcmp(name_var,'pr')
    flag_neg = dsEValue_mat<0;
    dsEValue_mat(flag_neg)=prctile(dsEValue_mat,25,'all')*rand(1,sum(flag_neg,'all'));
end

end
