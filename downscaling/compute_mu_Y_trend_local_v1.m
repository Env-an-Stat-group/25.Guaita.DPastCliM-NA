function [mu_Y_trend_local,mu_mov_ESM_trend] = compute_mu_Y_trend_local_v1(tgt_ESM_mth, tgt_lon, tgt_lat, time_ESM_mth, flag_cal, name_var, metaTable)
    % Ensure flag_cal is the same length as time_ESM_mth
    if length(time_ESM_mth) > length(flag_cal)
        flag_cal = [false(1, length(time_ESM_mth) - length(flag_cal)), flag_cal];
    end
    
    % Extract ESM data for the calibration period
    tgt_ESM_mth_tmp = reshape(tgt_ESM_mth, length(tgt_lon), length(tgt_lat), [], length(flag_cal));
    tgt_ESM_mth_cal = reshape(tgt_ESM_mth_tmp(:,:,:,flag_cal), size(tgt_ESM_mth,1), []);
    clear tgt_ESM_mth_tmp;
    
    % Create a grid of longitudes and latitudes
    [tgt_longrid, tgt_latgrid] = ndgrid(tgt_lon, tgt_lat);
    
    % Number of years for the moving average
    n_year_mov = 30;
    
    % Compute trends based on variable type
    switch name_var
        case 'tas'
            mu_cal_ESM   = mean(tgt_ESM_mth_cal, 2, 'omitnan');
            mu_mov_ESM_trend = movmean(tgt_ESM_mth, n_year_mov, 2, 'omitnan');
            mu_adjtrend  = mu_mov_ESM_trend - mu_cal_ESM;
        case 'pr'
            M_t = 1.1 + min(tgt_ESM_mth_cal, [], 2);
            mu_cal_ESM   = mean(log(tgt_ESM_mth_cal + M_t), 2, 'omitnan');
            mu_mov_ESM_trend = movmean(log(tgt_ESM_mth + M_t), n_year_mov, 2, 'omitnan');
            mu_adjtrend  = mu_mov_ESM_trend - mu_cal_ESM;
    end
    
    % Reshape adjustment trend
    mu_adjtrend = reshape(mu_adjtrend, length(tgt_lon), length(tgt_lat), []);
    
    % Initialize output
    mu_Y_trend_local = nan(height(metaTable), length(time_ESM_mth));
    
    % Interpolate trend adjustments for each station location
    for i_time = 1:length(time_ESM_mth)
        f_int = griddedInterpolant(tgt_longrid, tgt_latgrid, mu_adjtrend(:,:,i_time), 'nearest', 'nearest');
        for i_ID = 1:height(metaTable)
            mu_adjtrend_tmp = f_int(metaTable.lon(i_ID), metaTable.lat(i_ID));
            mu_Y_trend_local(i_ID, i_time) = metaTable.mu_Y(i_ID) + mu_adjtrend_tmp;
        end
    end
end