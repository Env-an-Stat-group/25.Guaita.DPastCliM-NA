function [obsTable_mth, rmse_array, r2_array, lm_list] =  ...
    recal_predval_v2(n_pc_array, pc_all, eof_all, ...
    tgt_ESM_mth, mu_ESM, ...
    obsTable_mth, metaTable, time_mthcal, time_mth, name_var, opt_spatial, opt_ECIPI, n_err_iter)
% This function calibrates the downscaling model based on the selected
% principal components (PCs)
% and applies it to the entire dataset. It calculates RMSE, R2 values for validation,
% and interpolates missing data using Gaussian Process Regression (GPR) for non-calibration stations.

% Input Arguments:
% - n_pc_array: Array of indices indicating which principal components to use.
% - pc_all: Matrix containing all the principal components.
% - eof_all: Matrix of Empirical Orthogonal Functions (EOFs) used for dimensionality reduction.
% - tgt_ESM_mth: The Earth System Model (ESM) variable field for the target variable (predictor).
% - mu_ESM: The mean of the ESM over the calibration period.
% - obsTable_mth: Table of observations for each station.
% - metaTable: Metadata for the stations.
% - time_mthcal: The months selected for the calibration period.
% - time_mth: All the selected months.
% - name_var: Name of the variable (either 'tas' or 'pr').
% - opt_spatial: Boolean flag indicating if spatial interpolation is to be performed.
% - opt_ECIPI: Boolean flag for calculating prediction intervals and confidence intervals using GPR.
% - n_err_iter: Number of error iterations for model validation.

% Output:
% - obsTable_mth: Updated observation table with downscaled values and uncertainties.
% - rmse_array: Array of RMSE values for model validation.
% - r2_array: Array of R2 values for model validation.
% - lm_list: List of linear regression models (one for each station).

%% Initializing and Preprocessing Data
tic

% precompute variables
nStations = height(metaTable);
nObs = height(obsTable_mth);

% Precompute indices for each station for efficient lookups
flag_indices = cell(nStations, 1);
for i_ID = 1:nStations
    flag_indices{i_ID} = ismember(obsTable_mth.ID, metaTable.ID(i_ID));  % Indices for each station
end

% Select the principal components (PCs) based on the provided indices (n_pc_array)
pc_subset = pc_all(:, n_pc_array);  % Subset of PCs to use in the model
eof_subset = eof_all(:, n_pc_array);  % Corresponding subset of EOFs for dimensionality reduction

% Initialize an empty list to store linear regression models for each station
lm_list = cell(nStations, 1);

%% Calibration: Fit Linear Models for Each Station

% Parallel loop to calibrate the model for each station in metaTable
obsTable_mth_const = parallel.pool.Constant(obsTable_mth);
flag_indices_const = parallel.pool.Constant(flag_indices);

parfor i_ID = 1:nStations
    if metaTable.flag_cal(i_ID)  % Only calibrate if the station is marked for calibration
        % Define the data for the current station (station-specific anomalies)
        flag_ID = flag_indices_const.Value{i_ID} & obsTable_mth_const.Value.flag_cal;
        obsTable_tmp = obsTable_mth_const.Value(flag_ID, :);  % Filter for calibration data
        obsTable_tmp = sortrows(obsTable_tmp, "month_since_0CE", "ascend");  % Sort by time
        
        % Select the relevant time steps for calibration
        flag_time_tmp = ismember(time_mthcal, obsTable_tmp.month_since_0CE);
    
        % Define the response variable as the anomaly in observations (SY)
        Y = obsTable_tmp.gO - obsTable_tmp.mu_gO;
        X = pc_subset(flag_time_tmp, :);

        % Fit a linear regression model without intercept (downscaling model)
        lm_list{i_ID} = fitlm(X, Y, 'Intercept', false);
    end
end

%% Predictions: Apply Models to All Time Steps
tic
% Initialize column to store downscaled anomalies (predictions)
obsTable_mth.ESgOds_hat(:) = nan;
obsTable_mth.SgOds_hat(:) = nan;

% Flag for calibration stations (those included in the calibration phase)
flag_calstat = ismember(obsTable_mth.ID, metaTable.ID(metaTable.flag_cal));

% Calculate RMSE and R2 for model validation
Value = obsTable_mth.Value;
Value_val = Value(obsTable_mth.flag_val & flag_calstat);  % Model selection data (validation + calibration)

% Initialize arrays for RMSE and R2 values across error iterations
rmse_array = nan(1, n_err_iter);
r2_array = nan(1, n_err_iter);

% Project the ESM data onto the EOFs to generate the predictor matrix (X)
SM = tgt_ESM_mth - mu_ESM;  % ESM anomalies (centered)
X = (SM' * eof_subset);  % Project the ESM anomalies onto the EOF space

%% Error Iterations: Generate Downscaled Predictions

% Parallel loop over error iterations (multiple realizations)
lm_list_const = parallel.pool.Constant(lm_list);

% Precompute time index mappings once per station
ind_mth_all = cell(nStations, 1);
for i_ID = 1:nStations
    if metaTable.flag_cal(i_ID)
        idx = flag_indices{i_ID};
        [~, ind_mth_all{i_ID}] = ismember(obsTable_mth.month_since_0CE(idx), time_mth);
    end
end
ind_mth_const = parallel.pool.Constant(ind_mth_all);

% Precompute reusable constants
flag_ms = obsTable_mth.flag_val & flag_calstat;
n_val = sum(flag_ms);
mu_gO = obsTable_mth.mu_gO;
O_t = obsTable_mth.O_t;

% Preallocate outputs
rmse_array = nan(n_err_iter, 1);
r2_array = nan(n_err_iter, 1);

parfor i_erriter = 1:n_err_iter
    SgOds_hat = nan(nObs, 1);  % Preallocate once per iteration

    % Loop over calibration stations
    for i_ID = 1:nStations
        if metaTable.flag_cal(i_ID)
            lm_tmp = lm_list_const.Value{i_ID};
            if ~isempty(lm_tmp)
                idx = flag_indices_const.Value{i_ID};
                ind_mth = ind_mth_const.Value{i_ID};

                % Predict anomalies
                SgOds_hat_tmp = random(lm_tmp, X);  % Realization
                SgOds_hat(idx) = SgOds_hat_tmp(ind_mth);
            end
        end
    end

    % Inverse transform to original scale
    Ods_hat = inversetransform(SgOds_hat, mu_gO, O_t, name_var);

    % Extract validation data only once
    Ods_hat_val = Ods_hat(flag_ms);

    % RMSE
    rmse_array(i_erriter) = rmse(Ods_hat_val, Value_val, 'omitnan');

    % R² and adjusted R²
    SSE = sum((Value_val - Ods_hat_val).^2, 'omitnan');
    SST = sum((Value_val - mean(Value_val, 'omitnan')).^2, 'omitnan');
    R2_val = 1 - SSE / SST;

    r2_array(i_erriter) = 1 - (n_val / (n_val - length(n_pc_array))) * (1 - R2_val);
end

%% Final Predictions: Obtain Realizations and Confidence Intervals

% Loop over each station and obtain final predictions
for i_ID = 1:nStations
    if metaTable.flag_cal(i_ID)  % Only process calibrated stations
        lm_tmp = lm_list{i_ID};

        if not(isempty(lm_tmp))
            % Generate mean predictions and confidence intervals
            SgOds_hat_tmp = random(lm_tmp, X);  % Realized anomalies for the station
    
            % Update obsTable_mth with the predictions and intervals
            flag_ID = ismember(obsTable_mth.ID, metaTable.ID(i_ID));
            [~, ind_mth] = ismember(obsTable_mth.month_since_0CE(flag_ID), time_mth);
            obsTable_mth.SgOds_hat(flag_ID) = SgOds_hat_tmp(ind_mth);  % Realized predictions

            if opt_ECIPI
                [ESgOds_hat, ~] = predict(lm_tmp, X);%[ESgOds_hat, SCI] = predict(lm_tmp, X);  % Mean prediction and confidence intervals (SCI)
                %[~, SPI] = predict(lm_tmp, X, 'prediction', 'observation');  % Prediction intervals (SPI)
                obsTable_mth.ESgOds_hat(flag_ID) = ESgOds_hat(ind_mth);  % Mean prediction for the anomalies
                %obsTable_mth.CI_infSY(flag_ID) = SCI(ind_mth, 1);  % Lower bound of the confidence interval
                %obsTable_mth.CI_supSY(flag_ID) = SCI(ind_mth, 2);  % Upper bound of the confidence interval
                %obsTable_mth.PI_infSY(flag_ID) = SPI(ind_mth, 1);  % Lower bound of the prediction interval
                %obsTable_mth.PI_supSY(flag_ID) = SPI(ind_mth, 2);  % Upper bound of the prediction interval
            end
        end
    end
end

%% Inverse Transform: Convert Anomalies Back to Original Values

% Convert both the mean predictions and realizations back to the original scale
obsTable_mth.Ods_hat = inversetransform(obsTable_mth.SgOds_hat, obsTable_mth.mu_gO, obsTable_mth.O_t, name_var);
if opt_ECIPI
    %obsTable_mth.EOds_hat = inversetransform(obsTable_mth.ESgOds_hat, obsTable_mth.mu_gO, obsTable_mth.O_t, name_var);
    %obsTable_mth.CI_inf = inversetransform(obsTable_mth.CI_infSY, obsTable_mth.mu_gO, obsTable_mth.O_t, name_var);
    %obsTable_mth.CI_sup = inversetransform(obsTable_mth.CI_supSY, obsTable_mth.mu_gO, obsTable_mth.O_t, name_var);
    %obsTable_mth.PI_inf = inversetransform(obsTable_mth.PI_infSY, obsTable_mth.mu_gO, obsTable_mth.O_t, name_var);
    %obsTable_mth.PI_sup = inversetransform(obsTable_mth.PI_supSY, obsTable_mth.mu_gO, obsTable_mth.O_t, name_var);
end

%% Interpolation: Use GPR to Fill Missing Values for Non-Calibration Stations

if opt_spatial
    % Initialize results for each time step (for interpolation purposes)
    time_groups = findgroups(obsTable_mth.month_since_0CE);
    obsTable_results = cell(max(time_groups), 1);
    %obsTable_results = cell(1, length(time_mth));
    obsTable_mth_const = parallel.pool.Constant(obsTable_mth);

    parfor i_time = 1:max(time_groups)
        % Extract data for the current time step
        flag_time = (time_groups == i_time);
        obsTable_tmp = obsTable_mth_const.Value(flag_time, :);

        % Filter out rows with missing values in `EOds_hat`
        flag_sample = not(isnan(obsTable_tmp.Ods_hat));
        if sum(not(flag_sample))
            % fill missing data
            x = obsTable_tmp.lat(flag_sample);  
            y = obsTable_tmp.lon(flag_sample);  
            z1 = obsTable_tmp.Ods_hat(flag_sample);  

            % get unique points
            [coords,ia,~] = unique([x y],'stable','rows');
            x = coords(:,1);
            y = coords(:,2);
            z1 = z1(ia);
    
            F1 = scatteredInterpolant(x, y, double(z1), 'natural', 'nearest'); % method: 'linear', 'nearest', 'natural'

            % Predict missing values using GPR
            missing_x = obsTable_tmp.lat(not(flag_sample));  % Latitude of stations with missing data
            missing_y = obsTable_tmp.lon(not(flag_sample));  % Longitude of stations with missing data
            
            obsTable_tmp.Ods_hat(not(flag_sample)) = single(F1([missing_x, missing_y]));
            
            % If ECIPI option is enabled, calculate confidence and prediction intervals using GPR
            if opt_ECIPI
                %z2 = obsTable_tmp.EOds_hat(flag_sample);  % Mean predicted values
                %z2_CIinf = obsTable_tmp.CI_inf(flag_sample);  % Confidence interval lower bound
                %z2_CIsup = obsTable_tmp.CI_sup(flag_sample);  % Confidence interval upper bound
                %z_PIinf = obsTable_tmp.PI_inf(flag_sample);  % Prediction interval lower bound
                %z_PIsup = obsTable_tmp.PI_sup(flag_sample);  % Prediction interval upper bound

                %z2 = z2(ia);
                %z2_CIinf = z2_CIinf(ia);
                %z2_CIsup = z2_CIsup(ia);
                %z_PIinf = z_PIinf(ia);
                %z_PIsup = z_PIsup(ia);

                %F2 = scatteredInterpolant(x, y, double(z2), 'natural', 'nearest'); 
                %F2_CIinf = scatteredInterpolant(x, y, double(z2_CIinf), 'natural', 'nearest');
                %F2_CIsup = scatteredInterpolant(x, y, double(z2_CIsup), 'natural', 'nearest'); 
                %F1_PIinf = scatteredInterpolant(x, y, double(z_PIinf), 'natural', 'nearest'); 
                %F1_PIsup = scatteredInterpolant(x, y, double(z_PIsup), 'natural', 'nearest'); 

                % Predict missing values for confidence and prediction intervals
                %obsTable_tmp.EOds_hat(not(flag_sample)) = single(F2([missing_x, missing_y]));
                %obsTable_tmp.CI_inf(not(flag_sample)) = single(F2_CIinf([missing_x, missing_y]));
                %obsTable_tmp.CI_sup(not(flag_sample)) = single(F2_CIsup([missing_x, missing_y]));
                %obsTable_tmp.PI_inf(not(flag_sample)) = single(F1_PIinf([missing_x, missing_y]));
                %obsTable_tmp.PI_sup(not(flag_sample)) = single(F1_PIsup([missing_x, missing_y]));
            end
        end
        % Store the updated data for this time step
        obsTable_results{i_time} = obsTable_tmp;
    end

    % Combine results
    obsTable_mth = vertcat(obsTable_results{:});
end

%% Residuals: Compute differences between observations and predictions
obsTable_mth.u(:) = nan;   % residuals from ESgOds_hat (deterministic fit)

for i_ID = 1:nStations
    if metaTable.flag_cal(i_ID)
        % Indices for this station
        flag_ID = ismember(obsTable_mth.ID, metaTable.ID(i_ID));
        
        % Observed anomalies (SY)
        SgO = obsTable_mth.gO(flag_ID) - obsTable_mth.mu_gO(flag_ID);
        
        obsTable_mth.u(flag_ID) = SgO - obsTable_mth.ESgOds_hat(flag_ID);
    end
end

toc
end
