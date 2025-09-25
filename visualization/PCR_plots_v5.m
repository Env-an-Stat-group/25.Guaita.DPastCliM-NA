% plots after downscaling


clear
close all

%% parameters

disp('setting initial parameters...')

rng(812)

% parameters
path_main = 'C:\Users\guait\Università Cattolica del Sacro Cuore\FIS-AMB-ECOFIS - Documenti\GUAITA\PALEON\downscaling\';%'/data/pguaita/downscaling/';
addpath(genpath(fullfile(path_main,'matlab_code_git')));
name_model = 'MPI-ESM1-2-LR'; % model name
name_var = 'tas'; % variable name
name_experiment = 'past2k';
% starting and ending year for plotting
year_start = 0;
year_start_hist = 1875;
year_end = 2014;
year_base_ESM_hist = 1850; % first year in the experiment
path_fig = fullfile(path_main,['downscaling_output_' name_model],'figures_PCR');
path_file = fullfile(path_main,['downscaling_output_' name_model]);
path_downmodel = fullfile(path_main,['downscaling_models_' name_model]);
path_shp_file = fullfile(path_main,'/matlab_code_git/visualization/world_borders/ne_10m_admin_0_countries.shp'); 
suffix = '_Hartfordtest';

disp(suffix)
disp(name_var)
%% parameters that you most likely should not change

switch name_var
    case 'tas'
        mth_interval = (1:12)';
    case 'pr'
        mth_interval = (1:12)';
end

ind_season = {1:12,3:5,6:8,9:11,[12 1 2]}; % define month indexes for seasons
i_mth_txt_loop = {'Yearly','MAM','JJA','SON','DJF'};
flag_loop = {'flag_year','flag_mam','flag_jja','flag_son','flag_djf'};

%% Figures parameters

% figure colorbar limits
switch name_var
        case 'pr'
            var_title = 'pr';
            % absolute value parameters
            abs_limit = [0 7];
            abs_step = 1;
            n_color_abs = 14;
            palette_abs = @(m) Tol_seq_iridescent(m);

            % bias parameters
            bias_limit = [-2.5 2.5];
            bias_step = 1;
            n_color_bias = 10;
            palette_bias = @(m) Tol_div_BuRd(m);

            % eof parameters
            eof_limit = [-0.05 0.05];
            eof_step = 0.02;
            n_color_eof = 10;
            palette_eof = @(m) Tol_div_BuRd(m);

            % range parameters
            range_limit = [0 10];
            range_step = 2;
            n_color_range = 10;
            palette_range = @(m) Tol_seq_smoothrainbow(m);

        case 'tas'
            var_title = 'tas';
            % absolute value parameters
            abs_limit = [-25 35];
            abs_step = 10;
            n_color_abs = 12;
            palette_abs = @(m) Tol_seq_YlOrBr(m);
            
            % bias parameters
            bias_limit = [-5 5];
            bias_step = 2;
            n_color_bias = 10;
            palette_bias = @(m) Tol_div_BuRd(m);

            % range parameters
            range_limit = [0 10];
            range_step = 2;
            n_color_range = 10;
            palette_range = @(m) Tol_seq_smoothrainbow(m);

            % eof parameters
            eof_limit = [-0.05 0.05];
            eof_step = 0.02;
            n_color_eof = 10;
            palette_eof = @(m) Tol_div_BuRd(m);

end
% define time windows for maps
time_bound = [0 100; 950 1050; 1900 2000];

res_fig = '-r500';
res_fig_low = '-r100';
res_plot= '-r300';

disp('     Done.')

%% setup variables

disp('Defining initial variables...')

if not(exist(path_fig,'dir'))
    mkdir(path_fig);
end

% load target grid resolution
load(fullfile(path_main,['static_maps/downscaling_grid' suffix '.mat']));

% filters for the domain 
lim_lat = [min(lat),max(lat)];
lim_lon = [min(lon),max(lon)];
filter_lat = [min(lat),max(lat)];
filter_lon = [min(lon),max(lon)];

% filter the target grid
flag_tgt_lat = lat>= filter_lat(1) & lat<=filter_lat(2);
flag_tgt_lon = lon>= filter_lon(1) & lon<=filter_lon(2);
lat = lat(flag_tgt_lat);
lon = lon(flag_tgt_lon);

% define lon and lat grid
[longrid,latgrid] = ndgrid(lon,lat);

% load landsea mask
landseamask = ncread(fullfile(path_main,'static_maps','IMERG_land_sea_mask.nc'),'landseamask');
flag_land = not(landseamask == 100);
lat_land = ncread(fullfile(path_main,'static_maps','IMERG_land_sea_mask.nc'),'lat');
lon_land = ncread(fullfile(path_main,'static_maps','IMERG_land_sea_mask.nc'),'lon');
[flag_land,lon_land,lat_land] = cmipcoord_2_map(flag_land,lon_land,lat_land);
% regrid land-sea mask
[lonland_grid,latland_grid] = ndgrid(lon_land,lat_land);
f_int_land = griddedInterpolant(lonland_grid,latland_grid,double(flag_land),'nearest');
flag_land = logical(f_int_land(longrid,latgrid));

clear f_int_land landseamask lonland_grid latland_grid lon_land lat_land

disp('     Done.')

%% load downscaled data and observed data (gridded)

disp('loading downscaled data and observed data (gridded)...')

% load raw data
load(fullfile(path_file,[name_var '_' name_model '_' name_experiment '_raw_PCRdownscaled' suffix]))
load(fullfile(path_file,[name_var '_' name_model '_' name_experiment '_raw_LRdownscaled' suffix]))

clear PI_mat PI_mat_lr

% get the predicted values by adding u_mat to the predicted average
switch name_var
    case 'tas'
        Ods_hat_mat = EOds_hat_mat + u_mat;
        Ods_hat_mat_lr = EOds_hat_mat_lr + u_mat_lr;
    case 'pr'
        Ods_hat_mat = pr_realizations(EOds_hat_mat,u_mat,...
            path_main,name_var,suffix,name_model);
        Ods_hat_mat_lr = pr_realizations(EOds_hat_mat_lr,u_mat_lr,...
            path_main,name_var,suffix,name_model);
end

clear EOds_hat_mat EOds_hat_mat_lr epsilon 

% load observations
load(fullfile(path_file,[name_var '_observations' suffix]))

clear obsValue_map

% load and aggregate gridded maps
list_dir_PCR    = dir(fullfile(fullfile(path_file,[name_var '_' name_model '_' name_experiment '_PCRdownscaled' suffix '_*.mat'])));
list_dir_LR    = dir(fullfile(fullfile(path_file,[name_var '_' name_model '_' name_experiment '_LRdownscaled' suffix '_*.mat'])));
list_dir_PI_PCR = dir(fullfile(fullfile(path_file,[name_var '_' name_model '_' name_experiment '_PI_PCRdownscaled' suffix '_*.mat'])));

for i_file = 1:length(list_dir_PCR)
    load(fullfile(list_dir_PCR(i_file).folder,list_dir_PCR(i_file).name));
    load(fullfile(list_dir_LR(i_file).folder,list_dir_LR(i_file).name));
    load(fullfile(list_dir_PI_PCR(i_file).folder,list_dir_PI_PCR(i_file).name));
    if i_file==1
        Ods_hat_map_tmp = Ods_hat_map;
        Ods_hat_map_lr_tmp = Ods_hat_map_lr;
        PI_map_tmp = PI_map;
        time_array_tmp = time_array;
    else
        Ods_hat_map_tmp = cat(3,Ods_hat_map_tmp,Ods_hat_map);
        Ods_hat_map_lr_tmp = cat(3,Ods_hat_map_lr_tmp,Ods_hat_map_lr);
        PI_map_tmp = cat(3,PI_map_tmp,PI_map);
        time_array_tmp = cat(1,time_array_tmp, time_array);
    end
end

% finalize variables
Ods_hat_map = Ods_hat_map_tmp;
Ods_hat_map_lr = Ods_hat_map_lr_tmp;
PI_map = PI_map_tmp;
time_array = time_array_tmp;

clear *_tmp

% sort them according to time
[time_array,idx_time] = sort(time_array);
Ods_hat_map = Ods_hat_map(:,:,idx_time);
PI_map = PI_map(:,:,idx_time,:);

disp('     Done.')

% define season flags:

disp('define season flags')

flag_year= true(length(time_ESM),1);
flag_mam = false(length(time_ESM),1);
flag_jja = false(length(time_ESM),1);
flag_son = false(length(time_ESM),1);
flag_djf = false(length(time_ESM),1);
flag_mam([3:12:end 4:12:end 5:12:end]) = true;
flag_jja([6:12:end 7:12:end 8:12:end]) = true;
flag_son([9:12:end 10:12:end 11:12:end]) = true;
flag_djf([12:12:end 1:12:end 2:12:end]) = true;

disp('     Done.')

%% Compute SSIM of mean fields (only over valid timesteps)

disp('Compute SSIM of mean fields...')

% Preallocate SSIM arrays
n_time = size(Ods_hat_map, 3);
ssim_pcr     = nan(n_time, 1);
ssim_lr     = nan(n_time, 1);
ssim_esm      = nan(n_time, 1);

% Identify valid time steps (e.g. those with at least 50 valid pixels)
valid_times = false(n_time,1);
for t = 1:n_time
    D = obsValue_test_map(:,:,t);
    if sum(~isnan(D),'all') > 50
        valid_times(t) = true;
    end
end

% Compute mean fields only over valid time steps
A_mean = mean(Ods_hat_map(:,:,valid_times), 3, 'omitnan');       % PCR
B_mean = mean(Ods_hat_map_lr(:,:,valid_times), 3, 'omitnan');    % LR
C_mean = mean(tgt_ESM(:,:,valid_times), 3, 'omitnan');           % ESM
D_mean = mean(obsValue_test_map(:,:,valid_times), 3, 'omitnan');      % OBS

% Mask of valid pixels across all mean fields
valid_mask = ~isnan(A_mean) & ~isnan(B_mean) & ~isnan(C_mean) & ~isnan(D_mean);

% SSIM dynamic range
dr_mean = max(D_mean(valid_mask)) - min(D_mean(valid_mask));

% Compute SSIMs on the mean fields
ssim_mean_pcr = ssim(A_mean(valid_mask), D_mean(valid_mask), 'DynamicRange', dr_mean);
ssim_mean_lr  = ssim(B_mean(valid_mask), D_mean(valid_mask), 'DynamicRange', dr_mean);
ssim_mean_esm = ssim(C_mean(valid_mask), D_mean(valid_mask), 'DynamicRange', dr_mean);

% Display
fprintf('Mean field SSIMs:\n');
fprintf('  PCR vs OBS: %.4f\n', ssim_mean_pcr);
fprintf('  LR  vs OBS: %.4f\n', ssim_mean_lr);
fprintf('  ESM vs OBS: %.4f\n', ssim_mean_esm);

close all

%% aggregate observation table (which contains also the downscaled product)

disp('Aggregating obs table')
tic

for i_mth_interval = 1:size(mth_interval,1)
    
    i_mth = mth_interval(i_mth_interval,:);    
    i_mth_txt = [int2str(mth_interval(i_mth_interval,1)) '-' int2str(mth_interval(i_mth_interval,end))];
    disp(['doing month ' i_mth_txt])

    tgt_ESM_mth = reshape(squeeze(tgt_ESM(:,:,i_mth,:)),[],size(tgt_ESM,4)*length(i_mth));

    % load downscaling temperature model
    path_loadfile_tmp = fullfile(path_downmodel,[name_var '_PCR_models_mth=' i_mth_txt suffix]);
    load(path_loadfile_tmp)
    path_loadfile_tmp = fullfile(path_downmodel,[name_var '_PCR_original_data_mth=' i_mth_txt suffix]);
    load(path_loadfile_tmp)
    clear path_loadfile_tmp  

    % aggregating the observation table
    if i_mth_interval==1
        obsTable = obsTable_mth;
    else
        obsTable = vertcat(obsTable,obsTable_mth);
    end

end

clear obsTable_mth lm_list eof_all ...
    mu_Y_local Y_t_local pc_all *_mth

%% create observation matrix (which is used for comparison with the
% final downscaled product) - this is taking only the non-spatialized
% testing values

disp(['create observation matrix (which is used for comparison with the final downscaled product) ' ...
    '- this is taking only the non-spatialized testing values'])

% Preallocate
obsTable_mat       = nan(size(Ods_hat_mat));
ESM_mat            = nan(size(Ods_hat_mat));
obsTable_mat_test  = nan(size(Ods_hat_mat));
ESM_mat_test       = nan(size(Ods_hat_mat));

% Map obsTable.ID to metaTable row indices
[~, row_idx] = ismember(obsTable.ID, metaTable.ID);

% Map obsTable.month_since_0CE to time_ESM column indices
[~, col_idx] = ismember(obsTable.month_since_0CE, time_ESM);

% Compute linear indices
linear_idx = sub2ind(size(Ods_hat_mat), row_idx, col_idx);

% Fill in matrices
obsTable_mat(linear_idx) = obsTable.Value;
ESM_mat(linear_idx)      = obsTable.M;

% Filter test entries only
test_mask     = obsTable.flag_test;
row_idx_test  = row_idx(test_mask);
col_idx_test  = col_idx(test_mask);
linear_idx_test = sub2ind(size(Ods_hat_mat), row_idx_test, col_idx_test);

obsTable_mat_test(linear_idx_test) = obsTable.Value(test_mask);
ESM_mat_test(linear_idx_test)      = obsTable.M(test_mask);

clear obsTable metaTable

disp('     Done.')

%% calculate performance metrics and compile table

disp('calculate performance metrics and compile table')

% Initialize output containers
metrics = {'MeanBias', 'StdBias', 'RMSE', 'R2', 'MAE', 'KGE', 'SSIM_avg_time', 'SSIM_mean_field'};
methods = {'PCR', 'LR', 'ESM'};
nMetrics = length(metrics);
nMethods = length(methods);
nSeasons = length(i_mth_txt_loop);

% Preallocate results array for seasonal metrics
all_results = zeros(nMetrics * nSeasons, nMethods);
row_labels = cell(nMetrics * nSeasons, 1);

for i_season = 1:nSeasons
    % Select current time mask
    flag_season = eval(flag_loop{i_season});

    % Apply mask to time dimension (assumes time is 2nd dim of input matrices)
    y_true = obsTable_mat_test(:, flag_season);
    y_pcr  = Ods_hat_mat(:, flag_season);
    y_lr   = Ods_hat_mat_lr(:, flag_season);
    y_esm  = ESM_mat_test(:, flag_season);

    % Flatten to vectors
    y_true = y_true(:);
    y_pcr  = y_pcr(:);
    y_lr   = y_lr(:);
    y_esm  = y_esm(:);

    % Remove NaNs
    valid_mask = ~(isnan(y_true) | isnan(y_pcr) | isnan(y_lr) | isnan(y_esm));
    y_true = y_true(valid_mask);
    y_pcr  = y_pcr(valid_mask);
    y_lr   = y_lr(valid_mask);
    y_esm  = y_esm(valid_mask);

    % Biases
    bias_pcr = y_pcr - y_true;
    bias_lr  = y_lr - y_true;
    bias_esm = y_esm - y_true;

    % Compute KGE for each method
    kge_pcr = kge(y_true, y_pcr);
    kge_lr  = kge(y_true, y_lr);
    kge_esm = kge(y_true, y_esm);

    % Compute SSIM for this season - full vectors (average over time / realizations)
    ssim_season_avgtime = [...
        mean(ssim_pcr(flag_season),'omitnan'), ...
        mean(ssim_lr(flag_season),'omitnan'), ...
        mean(ssim_esm(flag_season),'omitnan')];

    % Compute mean fields only over valid time steps
    A_mean = mean(Ods_hat_map(:,:,valid_times & flag_season), 3, 'omitnan');       % PCR
    B_mean = mean(Ods_hat_map_lr(:,:,valid_times & flag_season), 3, 'omitnan');    % LR
    C_mean = mean(tgt_ESM(:,:,valid_times & flag_season), 3, 'omitnan');           % ESM
    D_mean = mean(obsValue_test_map(:,:,valid_times & flag_season), 3, 'omitnan');      % OBS
    
    % Mask of valid pixels across all mean fields
    valid_mask = ~isnan(A_mean) & ~isnan(B_mean) & ~isnan(C_mean) & ~isnan(D_mean);
    
    % SSIM dynamic range
    dr_mean = max(D_mean(valid_mask)) - min(D_mean(valid_mask));
    
    % Compute SSIMs on the mean fields
    ssim_mean_pcr = ssim(A_mean(valid_mask), D_mean(valid_mask), 'DynamicRange', dr_mean);
    ssim_mean_lr  = ssim(B_mean(valid_mask), D_mean(valid_mask), 'DynamicRange', dr_mean);
    ssim_mean_esm = ssim(C_mean(valid_mask), D_mean(valid_mask), 'DynamicRange', dr_mean);

    ssim_season_meanfield = [ssim_mean_pcr,ssim_mean_lr,ssim_mean_esm];

    % Metrics for this season
    results = [
        mean(bias_pcr), mean(bias_lr), mean(bias_esm);
        std(bias_pcr),  std(bias_lr),  std(bias_esm);
        sqrt(mean(bias_pcr.^2)), sqrt(mean(bias_lr.^2)), sqrt(mean(bias_esm.^2));
        1 - sum((y_true - y_pcr).^2)/sum((y_true - mean(y_true)).^2), ...
        1 - sum((y_true - y_lr).^2)/sum((y_true - mean(y_true)).^2), ...
        1 - sum((y_true - y_esm).^2)/sum((y_true - mean(y_true)).^2);
        mean(abs(bias_pcr)), mean(abs(bias_lr)), mean(abs(bias_esm));
        kge_pcr, kge_lr, kge_esm;
        ssim_season_avgtime;
        ssim_season_meanfield
    ];

    % Store in master results array
    row_start = (i_season - 1) * nMetrics + 1;
    all_results(row_start:row_start + nMetrics - 1, :) = results;

    % Generate row labels
    for j = 1:nMetrics
        row_labels{row_start + j - 1} = sprintf('%s_%s', metrics{j}, i_mth_txt_loop{i_season});
    end
end

% ---- Add monthly climatology difference and trend (only for yearly) ----
% Extract yearly mask

% Compute monthly climatology: mean per month over all sites and times
clim_true = zeros(12,1);
clim_pcr = zeros(12,1);
clim_lr = zeros(12,1);
clim_esm = zeros(12,1);

for m = 1:12
    flag_month = false(size(flag_year));
    flag_month(m:12:end) = true; % select all times for month m

    y_t = obsTable_mat_test(:, flag_month);
    y_p = Ods_hat_mat(:, flag_month);
    y_l = Ods_hat_mat_lr(:, flag_month);
    y_e = ESM_mat_test(:, flag_month);

    % Remove NaNs
    valid_mask = ~(isnan(y_t) | isnan(y_p) | isnan(y_l) | isnan(y_e));
    y_t = y_t(valid_mask);
    y_p = y_p(valid_mask);
    y_l = y_l(valid_mask);
    y_e = y_e(valid_mask);

    clim_true(m) = mean(y_t);
    clim_pcr(m) = mean(y_p);
    clim_lr(m) = mean(y_l);
    clim_esm(m) = mean(y_e);
end

% Monthly climatology difference = mean absolute difference over months
clim_diff_pcr = mean(abs(clim_pcr - clim_true));
clim_diff_lr = mean(abs(clim_lr - clim_true));
clim_diff_esm = mean(abs(clim_esm - clim_true));

% Trend analysis: simple linear regression slope of yearly data
time_idx = find(flag_year);
time_vec = time_idx; % use time indices as x values

obsTable_mat_tmp = squeeze(mean(reshape(obsTable_mat,size(obsTable_mat,1),12,[]),2,'omitnan'));
Ods_hat_mat_tmp = squeeze(mean(reshape(Ods_hat_mat,size(Ods_hat_mat,1),12,[]),2,'omitnan'));
Ods_hat_mat_lr_tmp = squeeze(mean(reshape(Ods_hat_mat_lr,size(Ods_hat_mat_lr,1),12,[]),2,'omitnan'));
ESM_mat_tmp = squeeze(mean(reshape(ESM_mat,size(ESM_mat,1),12,[]),2,'omitnan'));

% homogenize Nans across datasets for trend calculation
flag_nan_tmp = (isnan(obsTable_mat_tmp) | isnan(Ods_hat_mat_tmp) | isnan(Ods_hat_mat_lr_tmp) | isnan(ESM_mat_tmp));
obsTable_mat_tmp(flag_nan_tmp) = nan;
Ods_hat_mat_tmp(flag_nan_tmp) = nan;
Ods_hat_mat_lr_tmp(flag_nan_tmp) = nan;
ESM_mat_tmp(flag_nan_tmp) = nan;

mean_true = mean(obsTable_mat_tmp, 1,'omitnan')';
mean_pcr = mean(Ods_hat_mat_tmp, 1,'omitnan')';
mean_lr = mean(Ods_hat_mat_lr_tmp, 1,'omitnan')';
mean_esm = mean(ESM_mat_tmp, 1,'omitnan')';

% Remove NaNs for trend calculation
valid_t = ~isnan(mean_true) & ~isnan(mean_pcr) & ~isnan(mean_lr) & ~isnan(mean_esm);

% Fit linear trends
trend_true = polyfit(time_vec(valid_t), mean_true(valid_t), 1);
trend_pcr  = polyfit(time_vec(valid_t), mean_pcr(valid_t), 1);
trend_lr   = polyfit(time_vec(valid_t), mean_lr(valid_t), 1);
trend_esm  = polyfit(time_vec(valid_t), mean_esm(valid_t), 1);

% Slopes as trend metric
trend_slope_true = trend_true(1);
trend_slope_pcr  = trend_pcr(1);
trend_slope_lr   = trend_lr(1);
trend_slope_esm  = trend_esm(1);

% ---- Add monthly climatology difference and trend (only for yearly) ----

extra_metrics = {... 
    'ClimDiff', ...
    'TrendSlopeRatio_(model/obs)', ...
    'TrendSlopeRatio_(ds/esm)'};

% 1. Monthly climatology difference
clim_diff_row = [clim_diff_pcr, clim_diff_lr, clim_diff_esm];

% 2. Trend slope ratios
trend_model_obs = [trend_slope_pcr/trend_slope_true, ...
                   trend_slope_lr/trend_slope_true, ...
                   trend_slope_esm/trend_slope_true];

trend_ds_esm = [trend_slope_pcr/trend_slope_esm, ...
                trend_slope_lr/trend_slope_esm, ...
                trend_slope_esm/trend_slope_esm];

% Combine all extra rows
extra_results = [
    clim_diff_row;
    trend_model_obs;
    trend_ds_esm];

% Append extra metrics rows
all_results = [all_results; extra_results];

for i = 1:length(extra_metrics)
    row_labels{end+1} = extra_metrics{i};
end

% Create table
result_table_final = array2table(all_results, ...
    'RowNames', row_labels, ...
    'VariableNames', methods);

disp(result_table_final)

% Define filename with path and .csv extension
filename = fullfile(path_fig, ['performance_table_' name_var '.csv']);

% Save table as CSV
writetable(result_table_final, filename, 'WriteRowNames', true);

clear y_* valid_mask *_tmp Ods_hat_map_lr

disp('    Done.')

%% plots

disp('preparing large variables for plotting...')

Ods_hat_4d     = single(reshape(Ods_hat_map,length(lon),length(lat),12,[]));
clear Ods_hat_map 
PI_map_5d     = single(reshape(PI_map,length(lon),length(lat),12,2015,2));
clear PI_map
obsValue_test_4d    = single(reshape(obsValue_test_map,length(lon),length(lat),12,[]));
clear obsValue_test_map 
ESM_4d    = single(reshape(tgt_ESM,length(lon),length(lat),12,[]));
clear tgt_ESM 

% filter on land
flag_land_4d = repmat(flag_land,1,1,12,size(Ods_hat_4d,4));

% mask out sea nodes
Ods_hat_4d(not(flag_land_4d)) = nan;
obsValue_test_4d(not(flag_land_4d)) = nan;
ESM_4d(not(flag_land_4d)) = nan;

clear flag_land_4d

flag_land_5d = repmat(flag_land,1,1,12,size(Ods_hat_4d,4),2);
PI_map_5d(not(flag_land_5d)) = nan;

clear flag_land_5d

close all

disp('     Done.')

%% mean biases for PCR and for ESM over testing dataset

for i_season = 1:length(ind_season)
    i_mth_txt = i_mth_txt_loop{i_season};
    ind_season_tmp = ind_season{i_season};
    
    % Compute seasonal averages
    dsVal_tmp       = squeeze(Ods_hat_4d(:,:,ind_season_tmp,:));
    obsVal_test_tmp = squeeze(obsValue_test_4d(:,:,ind_season_tmp,:));
    ESM_tmp         = squeeze(ESM_4d(:,:,ind_season_tmp,:));
    
    text_title = [i_mth_txt ' ' name_var ' MB'];
    
    % --- PCR ---
    PCR_plot = mean(dsVal_tmp - obsVal_test_tmp,[3 4],'omitnan');
    save_name = fullfile(path_fig, [text_title '_PCR' suffix]);
    
    geosurfm_meansdlabel(PCR_plot, flag_land, lat, lon, bias_limit, ...
        save_name, unit_var, [text_title ' (PCR)'], ...
        lim_lat, lim_lon, bias_step, palette_bias(n_color_bias), path_shp_file, res_fig);

    % --- ESM ---
    ESM_plot = mean(ESM_tmp - obsVal_test_tmp,[3 4],'omitnan');
    save_name_ESM = fullfile(path_fig, [text_title '_ESM' suffix]);
    
    geosurfm_meansdlabel(ESM_plot, flag_land, lat, lon, bias_limit, ...
        save_name_ESM, unit_var, [text_title ' (ESM)'], ...
        lim_lat, lim_lon, bias_step, palette_bias(n_color_bias), path_shp_file, res_fig);
end

close all

clear *_tmp

disp ('mb for PCR')

%% percentile differences

for i_season = 1:length(ind_season)
    i_mth_txt = i_mth_txt_loop{i_season};
    ind_season_tmp = ind_season{i_season};
    dsVal_tmp       = squeeze(Ods_hat_4d(:,:,ind_season_tmp,:));
    obsVal_test_tmp = squeeze(obsValue_test_4d(:,:,ind_season_tmp,:));
    ESM_tmp         = squeeze(ESM_4d(:,:,ind_season_tmp,:));
    
    %10th percentile
    text_title = [i_mth_txt ' ' name_var ' 10th p. diff.'];    
    % --- PCR ---
    PCR_plot = prctile(dsVal_tmp,10,[3 4]) - prctile(obsVal_test_tmp,10,[3 4]);
    save_name = fullfile(path_fig, [text_title '_PCR' suffix]);
    
    geosurfm_meansdlabel(PCR_plot, flag_land, lat, lon, bias_limit, ...
        save_name, unit_var, [text_title ' (PCR)'], ...
        lim_lat, lim_lon, bias_step, palette_bias(n_color_bias), path_shp_file, res_fig);

    % --- ESM ---
    ESM_plot = prctile(ESM_tmp,10,[3 4]) - prctile(obsVal_test_tmp,10,[3 4]);
    save_name_ESM = fullfile(path_fig, [text_title '_ESM' suffix]);
    
    geosurfm_meansdlabel(ESM_plot, flag_land, lat, lon, bias_limit, ...
        save_name_ESM, unit_var, [text_title ' (ESM)'], ...
        lim_lat, lim_lon, bias_step, palette_bias(n_color_bias), path_shp_file, res_fig);

    %90th percentile
    text_title = [i_mth_txt ' ' name_var ' 90th p. diff.'];    
    % --- PCR ---
    PCR_plot = prctile(dsVal_tmp,90,[3 4]) - prctile(obsVal_test_tmp,90,[3 4]);
    save_name = fullfile(path_fig, [text_title '_PCR' suffix]);
    
    geosurfm_meansdlabel(PCR_plot, flag_land, lat, lon, bias_limit, ...
        save_name, unit_var, [text_title ' (PCR)'], ...
        lim_lat, lim_lon, bias_step, palette_bias(n_color_bias), path_shp_file, res_fig);

    % --- ESM ---
    ESM_plot = prctile(ESM_tmp,90,[3 4]) - prctile(obsVal_test_tmp,90,[3 4]);
    save_name_ESM = fullfile(path_fig, [text_title '_ESM' suffix]);
    
    geosurfm_meansdlabel(ESM_plot, flag_land, lat, lon, bias_limit, ...
        save_name_ESM, unit_var, [text_title ' (ESM)'], ...
        lim_lat, lim_lon, bias_step, palette_bias(n_color_bias), path_shp_file, res_fig);
    
end    

clear *_tmp

close all

disp ('bias of extremes')

%% absolute field and anomalies

for i_season = 1:length(ind_season)
    i_mth_txt = i_mth_txt_loop{i_season};
    ind_season_tmp = ind_season{i_season};
    flag_season = eval(flag_loop{i_season});
    dsVal_tmp       = squeeze(Ods_hat_4d(:,:,ind_season_tmp,:));
    dsVal_longterm = mean(dsVal_tmp,[3 4]);
    PI_range_tmp     = squeeze(PI_map_5d(:,:,ind_season_tmp,:,2)-PI_map_5d(:,:,ind_season_tmp,:,1));
    PI_range_longerm     = mean(PI_map_5d(:,:,ind_season_tmp,:,2)-PI_map_5d(:,:,ind_season_tmp,:,1),[3 4]);

    % Define figure save path
    save_name = fullfile(path_fig, [name_var '_abs field_' i_mth_txt suffix]);
    geosurfm_meansdlabel(dsVal_longterm, flag_land, lat, lon, abs_limit, ...
        save_name, unit_var, [i_mth_txt ' mean PCR-downscaled ' name_var ' (2k years)'], ...
        lim_lat, lim_lon, abs_step, palette_abs(n_color_abs), path_shp_file, res_fig);
    plot_anomaly(dsVal_tmp-dsVal_longterm, time_bound, year_start, flag_land, ...
        lat, lon, bias_limit/2, bias_step/2, palette_bias, path_shp_file, res_fig, ...
        name_var, i_mth_txt, unit_var, lim_lat, lim_lon, n_color_bias, path_fig, suffix)
    save_name = fullfile(path_fig, [name_var '_PI range_' i_mth_txt suffix]);
    geosurfm_meansdlabel(PI_range_longerm, flag_land, lat, lon, abs_limit*2, ...
        save_name, unit_var, [i_mth_txt ' mean PCR-downscaled 95% PI for ' name_var ' (2k years)'], ...
        lim_lat, lim_lon, abs_step*2, palette_abs(n_color_abs), path_shp_file, res_fig);
    plot_range_anomaly_PI(PI_range_tmp-PI_range_longerm, time_bound, year_start, flag_land, ...
        lat, lon, bias_limit/2, bias_step/2, palette_bias, path_shp_file, res_fig, ...
        name_var, i_mth_txt, unit_var, lim_lat, lim_lon, n_color_bias, path_fig, suffix)
    % ESM plots
    save_name = fullfile(path_fig, [name_var '_ESM field_' i_mth_txt suffix]);
    geosurfm_meansdlabel(dsVal_longterm, flag_land, lat, lon, abs_limit, ...
        save_name, unit_var, [i_mth_txt ' mean ESM ' name_var ' (2k years)'], ...
        lim_lat, lim_lon, abs_step, palette_abs(n_color_abs), path_shp_file, res_fig);
    plot_anomaly_ESM(dsVal_tmp-dsVal_longterm, time_bound, year_start, flag_land, ...
        lat, lon, bias_limit/2, bias_step/2, palette_bias, path_shp_file, res_fig, ...
        name_var, i_mth_txt, unit_var, lim_lat, lim_lon, n_color_bias, path_fig, suffix)

end

close all

disp ('absolute fields')

disp('Finished.')
