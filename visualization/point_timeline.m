% plots after downscaling


clear
close all

%% parameters

disp('setting initial parameters...')

rng(812)

% parameters
path_main = '/data/pguaita/downscaling/';
addpath(genpath(fullfile(path_main,'matlab_code_git')));
name_model = 'MPI-ESM1-2-LR'; % model name
name_var = 'pr'; % variable name
name_experiment = 'past2k';
n_min_yr = 100; % minimum number of years for stations
path_fig = fullfile(path_main,['downscaling_output_' name_model],'figures_PCR');
path_file = fullfile(path_main,['downscaling_output_' name_model]);
path_obs = fullfile(path_main,'obs_data');
path_downmodel = fullfile(path_main,['downscaling_models_' name_model]);
path_shp_file = fullfile(path_main,'/matlab_code_git/visualization/world_borders/ne_10m_admin_0_countries.shp'); 
suffix = '_NA_020';

% define subdomain for station selection (each row is a different
% subdomain)
lon_subdom = [-73.1 -71.1];%; -94.8 -92.8; -96.2 -94.2; -120.8 -110.8];
lat_subdom = [41.5 43.5];%; 43.7 45.7; 46.2 48.2; 29.9 37.9];

%% load grid and define limits
load(fullfile(path_main, ['static_maps/downscaling_grid' suffix '.mat']));
lat = lat(lat >= min(lat) & lat <= max(lat));
lon = lon(lon >= min(lon) & lon <= max(lon));
lim_lat = [min(lat), max(lat)];
lim_lon = [min(lon), max(lon)];

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

%% load observation table and downscaled matrix (mean value and prediction intervals)

disp('load observation table and downscaled station data (including PI)')

load(fullfile(path_obs,['GHCNm_' name_var suffix '.mat']))
load(fullfile(path_file,[name_var '_' name_model '_' name_experiment '_raw_PCRdownscaled' suffix '.mat']))

%load metadata
metaTable = readtable(fullfile(path_file,[name_var '_metadata_table' suffix '.csv']));

disp('     Done.')

% get the predicted values by adding u_mat to the predicted average
switch name_var
    case 'tas'
        Ods_hat_mat = EOds_hat_mat + u_mat; %EOds_hat_mat + u_mat;
    case 'pr'
        Ods_hat_mat = pr_realizations(EOds_hat_mat,u_mat,...
            path_main,name_var,suffix,name_model);
end

clear EOds_hat_mat u_mat

%% filter stations for sufficient data coverage

disp('filter stations for sufficient data coverage')

% Preallocate arrays
nStations = height(metaTable);
length_year = zeros(nStations,1);

% Loop only once over stations to compute timeseries length
for i_ID = 1:nStations
    flag_ID = ismember(obsTable.ID, metaTable.ID(i_ID));
    obsTable_tmp = obsTable(flag_ID, :);
    % filter on months 1950-1955 and 1995-2000
    % flag_mth_2000 = obsTable_tmp.month_since_0CE>=23920 & obsTable_tmp.month_since_0CE<24000;  
    % obsTable_tmp = obsTable_tmp(flag_mth_2000,:);
    length_year(i_ID) = floor(height(obsTable_tmp) / 12);
end

% Identify short stations
idx_remove = (length_year < n_min_yr);

% Count removed
n_removed = sum(idx_remove);

disp(['removing ' int2str(n_removed) '/' int2str(nStations) ...
    ' timeseries shorter than ' int2str(n_min_yr) ' years.']);

% remove from obsTable
obsTable = obsTable(~ismember(obsTable.ID, metaTable.ID(idx_remove)), :);

% Remove from tables and matrices in one shot
metaTable(idx_remove,:) = [];
Ods_hat_mat(idx_remove,:) = [];
PI_mat(idx_remove,:,:) = [];

disp('     Done.')

%% load observations
load(fullfile(path_file,[name_var '_observations' suffix '.mat']),'tgt_ESM')

%% select stations to plot
close all

% filter metaTable to select ID in the selected subdomain
for i_sd = 1:size(lon_subdom,1)
    flag_subdom_tmp = lon_subdom(i_sd,1)<metaTable.lon & metaTable.lon<lon_subdom(i_sd,2) & ...
        lat_subdom(i_sd,1)<metaTable.lat & metaTable.lat<lat_subdom(i_sd,2);
    if i_sd == 1
        flag_subdom = flag_subdom_tmp;
    else
        flag_subdom = flag_subdom | flag_subdom_tmp;
    end
end

metaTable = metaTable(flag_subdom,:);
Ods_hat_mat = Ods_hat_mat(flag_subdom,:);
PI_mat = PI_mat(flag_subdom,:,:);

disp(['left with ' int2str(height(metaTable)) '/' int2str(nStations) ...
    ' after selecting subdomain(s).']);

switch name_var
    case 'pr'
        list_point = metaTable.ID; %{'USC00193052','USC00190666','USC00190120','USC00190190','USC00067432','USC194105','USC00049452'};%metaTable.ID;
    case 'tas'
        list_point = {'USW00014740','USW00093193'};%metaTable.ID;
end

%% plot map with stations

plot_station_map_long(metaTable, lim_lat, lim_lon, path_shp_file, ...
    fullfile(path_fig, ['longstations_' name_var '.png']), [name_var ' GHCN-m stations covering most of 1950-1955 or 1995-2000']); % >' int2str(n_min_yr) ' years' ]);

%% Loop over points
for i_station = 1:length(list_point)
    ID_stat = list_point{i_station};    
    flag_ID_meta = ismember(metaTable.ID, ID_stat);
    metaTable_tmp = metaTable(flag_ID_meta,:);
    flag_ID_obsTable = ismember(obsTable.ID, ID_stat);
    obsTable_tmp = obsTable(flag_ID_obsTable,:);

    % Find closest grid node
    [~, lat_id] = min(abs(lat - metaTable_tmp.lat));  % nearest latitude index
    [~, lon_id] = min(abs(lon - metaTable_tmp.lon));  % nearest longitude index

    % Convert to linear index if your data are 2D [lat x lon x time]
    node_id = sub2ind([length(lat), length(lon)], lat_id, lon_id);
        
    % Apply region mask: keep only rows (grid cells) inside the region
    dsVal_tmp       = Ods_hat_mat(flag_ID_meta, :);
    dsPIinf_tmp     = PI_mat(flag_ID_meta, :, 1);
    dsPIsup_tmp     = PI_mat(flag_ID_meta, :, 2);
    obs_value_tmp = obsTable_tmp.Value;
    time_obs_tmp = obsTable_tmp.month_since_0CE; 
    [time_obs_tmp, idx_sort] = sort(time_obs_tmp);
    obs_value_tmp = obs_value_tmp(idx_sort);
    ESM_tmp         = squeeze(tgt_ESM(lon_id,lat_id, :));
    
    % add station name to output filename
    suffix_region = [ID_stat ' (' num2str(metaTable_tmp.lon) '°,' num2str(metaTable_tmp.lat) '°)'];
        
    %% Call plotting function with region data

    % yearly mean
    plot_downscaled_timeline_point_year(name_var, unit_var, time_ESM, time_obs_tmp, ...
        obs_value_tmp, ESM_tmp, dsVal_tmp, ...
        path_fig, suffix_region, name_model, res_plot)

    % time slices
    plot_downscaled_timeline_point(name_var, unit_var, time_ESM, time_obs_tmp, ...
        obs_value_tmp, ESM_tmp, dsVal_tmp, dsPIinf_tmp, dsPIsup_tmp, ...
        path_fig, suffix_region, name_model, res_plot, [250 255], true)

    plot_downscaled_timeline_point(name_var, unit_var, time_ESM, time_obs_tmp, ...
        obs_value_tmp, ESM_tmp, dsVal_tmp, dsPIinf_tmp, dsPIsup_tmp, ...
        path_fig, suffix_region, name_model, res_plot, [950 955], false)

    plot_downscaled_timeline_point(name_var, unit_var, time_ESM, time_obs_tmp, ...
        obs_value_tmp, ESM_tmp, dsVal_tmp, dsPIinf_tmp, dsPIsup_tmp, ...
        path_fig, suffix_region, name_model, res_plot, [1250 1255], false)

    plot_downscaled_timeline_point(name_var, unit_var, time_ESM, time_obs_tmp, ...
        obs_value_tmp, ESM_tmp, dsVal_tmp, dsPIinf_tmp, dsPIsup_tmp, ...
        path_fig, suffix_region, name_model, res_plot, [1995 2000], false)
    
end
