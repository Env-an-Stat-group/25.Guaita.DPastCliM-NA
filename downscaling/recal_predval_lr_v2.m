function [obsTable_mth, rmse_array_lr, r2_array_lr, lm_list] =  ...
    recal_predval_lr_v2(obsTable_mth, metaTable, time_mth, name_var, n_err_iter)
% calibrate the linear regression model for downscaling, and apply it to the whole dataset
% obsTable_mth = table containing all observations
% metaTable = metadata for the stations
% time_mth = all the selected months (starting from January 0 CE)
% name_var = name of the variable (originally either 'tas' or 'pr')
% n_err_iter = number of error iterations for model validation
tic

% Initialize an empty cell array to store linear regression models for each station
lm_list = cell(height(metaTable),1);

% Recalibrate the model for each station in metaTable
obsTable_mth_const = parallel.pool.Constant(obsTable_mth);

parfor i_ID = 1:height(metaTable)
    if metaTable.flag_cal(i_ID)  % Check if the station is marked for calibration
        % Identify observations for the current station
        flag_ID = ismember(obsTable_mth_const.Value.ID, metaTable.ID(i_ID));  % Find rows corresponding to the current station
        obsTable_tmp = obsTable_mth_const.Value(flag_ID & obsTable_mth_const.Value.flag_cal,:);  % Extract the relevant data for calibration
        obsTable_tmp = sortrows(obsTable_tmp, "month_since_0CE", "ascend");  % Sort data by time

        % Calculate anomalies for the predictor (ESM data) and the response (observations)
        SX = obsTable_tmp.ValueESM - obsTable_tmp.mu_ESM;  % Anomalies of ESM data
        SY = obsTable_tmp.Y - obsTable_tmp.mu_Y;  % Anomalies of observed data

        % Fit a linear model to the data (no intercept)
        lm_list{i_ID} = fitlm(SX, SY, 'Intercept', false);  % Fit the linear regression model for each station
    end
end

%%
% After calibration, use the linear models to make predictions for all time steps

% Initialize the downscaled values column with NaN (for anomaly predictions)
obsTable_mth.dsESY_lr(:) = nan;   

% Precompute indices for each station (for efficient lookups)
flag_indices = cell(height(metaTable), 1);
for i_ID = 1:height(metaTable)
    flag_indices{i_ID} = ismember(obsTable_mth.ID, metaTable.ID(i_ID));  % Indices for the current station
end

% Flag for calibration stations
flag_calstat = ismember(obsTable_mth.ID, metaTable.ID(metaTable.flag_cal));  

% Calculate RMSE and R2 adjusted distributions for model validation
Value = obsTable_mth.Value;
Value_val = Value(obsTable_mth.flag_val & flag_calstat);  % Values for validation (excluding calibration stations)

rmse_array_lr = nan(1, n_err_iter);  % Initialize array for RMSE values
r2_array_lr = nan(1, n_err_iter);  % Initialize array for R2 values

% Parallel loop over error iterations
obsTable_mth_const = parallel.pool.Constant(obsTable_mth);
lm_list_const = parallel.pool.Constant(lm_list);
parfor i_erriter = 1:n_err_iter
    % Temporary storage for predictions
    dsSY_lr_tmp = nan(height(obsTable_mth_const.Value), 1);

    % Inner loop over each station in metaTable
    for i_ID = 1:height(metaTable)
        if metaTable.flag_cal(i_ID)  % Only process stations marked for calibration
            lm_tmp = lm_list_const.Value{i_ID};  % Get the linear model for the current station
            idx = flag_indices{i_ID};  % Precomputed indices for the station

            % Project ESM data anomalies (SX) onto the linear model
            obsTable_tmp = obsTable_mth_const.Value(idx, :);
            SX = obsTable_tmp.ValueESM - obsTable_tmp.mu_ESM;  % Anomalies of ESM data

            % Make predictions using the linear model
            dsSY_lr_tmp(idx) = random(lm_tmp, SX);  % Get predicted anomalies for the station
        end
    end

    % Perform inverse transform to convert anomalies back to original scale
    dsValue = inversetransform(dsSY_lr_tmp, obsTable_mth_const.Value.mu_Y, obsTable_mth_const.Value.Y_t, name_var);
    dsValue_val = dsValue(obsTable_mth_const.Value.flag_val & flag_calstat);  % Predicted values for validation

    % Calculate RMSE and R2 adjusted for this error iteration
    rmse_array_lr(i_erriter) = rmse(dsValue_val, Value_val,'omitnan');  % RMSE calculation
    SSE = sum((Value_val - dsValue_val).^2,'omitnan');  % Sum of squared errors
    SST = sum((Value_val - mean(Value_val,'omitnan')).^2,'omitnan');  % Total sum of squares
    r2_array_lr(i_erriter) = 1 - SSE / SST;  % R2 calculation
end

%%
% Final step: Generate the realizations (predictions) for each station
for i_ID = 1:height(metaTable)
    if metaTable.flag_cal(i_ID)
        lm_tmp = lm_list{i_ID};  % Get the linear model for the current station

        % Project ESM data anomalies (SX) onto the model to get predictions
        flag_ID = ismember(obsTable_mth.ID, metaTable.ID(i_ID));  % Find rows for the current station
        obsTable_tmp = obsTable_mth(flag_ID, :);
        SX = obsTable_tmp.ValueESM - obsTable_tmp.mu_ESM;  % Anomalies of ESM data

        % Make predictions
        ESY_hat = predict(lm_tmp, SX);  % Mean predicted anomalies
        SY_hat = random(lm_tmp, SX);  % Realization of the predicted anomalies

        % Update the observation table with the predicted values
        obsTable_mth.dsESY_lr(flag_ID) = ESY_hat;  % Store mean anomaly predictions
        obsTable_mth.dsSY_lr(flag_ID) = SY_hat;  % Store realization anomaly predictions
    end
end

%%
% Transform back the downscaled anomalies to the original scale (both mean prediction and realizations)
obsTable_mth.dsEValue_lr    = inversetransform(obsTable_mth.dsESY_lr, obsTable_mth.mu_Y, obsTable_mth.Y_t, name_var);
obsTable_mth.dsValue_lr     = inversetransform(obsTable_mth.dsSY_lr, obsTable_mth.mu_Y, obsTable_mth.Y_t, name_var);

%%
% Interpolate NaN values for non-calibration stations using Gaussian Process Regression (GPR)
% Some stations may have missing downscaled values; we use GPR to interpolate them

% Preallocate a cell array to hold results for each time step
obsTable_results = cell(1, length(time_mth));

obsTable_mth_const = parallel.pool.Constant(obsTable_mth);
parfor i_time = 1:length(time_mth)
    % Extract data for the current time step
    flag_time = ismember(obsTable_mth_const.Value.month_since_0CE, time_mth(i_time));  % Find rows for the current time step
    obsTable_tmp = obsTable_mth_const.Value(flag_time, :);

    % Filter out rows with missing values in `dsEValue_lr`
    flag_sample = not(isnan(obsTable_tmp.dsEValue_lr));  % Find rows with valid predictions
    x = obsTable_tmp.lat(flag_sample);  % Latitude of stations with valid predictions
    y = obsTable_tmp.lon(flag_sample);  % Longitude of stations with valid predictions
    z1 = obsTable_tmp.dsEValue_lr(flag_sample);  % Mean predicted values
    z2 = obsTable_tmp.dsValue_lr(flag_sample);  % Realization values

    % get unique points
    [coords,ia,~] = unique([x y],'stable','rows');
    x = coords(:,1);
    y = coords(:,2);
    z1 = z1(ia);
    z2 = z2(ia);

    F1 = scatteredInterpolant(x, y, double(z1), 'natural', 'nearest'); 
    F2 = scatteredInterpolant(x, y, double(z2), 'natural', 'nearest'); 

    % Predict missing values using the GPR models
    missing_x = obsTable_tmp.lat(not(flag_sample));  % Latitude for missing stations
    missing_y = obsTable_tmp.lon(not(flag_sample));  % Longitude for missing stations

    % Predict values for missing stations
    obsTable_tmp.dsEValue_lr(not(flag_sample)) = single(F1([missing_x missing_y]));
    obsTable_tmp.dsValue_lr(not(flag_sample)) = single(F2([missing_x missing_y]));

    % Store the result for this iteration
    obsTable_results{i_time} = obsTable_tmp;
end

% Combine the results back into `obsTable_mth` after interpolation
for i_time = 1:length(time_mth)
    flag_time = ismember(obsTable_mth.month_since_0CE, time_mth(i_time));  % Find rows for the current time step
    obsTable_mth(flag_time, :) = obsTable_results{i_time};  % Update the main table with the interpolated data
end

toc

end
