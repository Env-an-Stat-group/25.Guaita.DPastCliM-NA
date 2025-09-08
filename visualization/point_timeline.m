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
name_var = 'tas'; % variable name
name_experiment = 'past2k';
n_min_yr = 120; % minimum number of years for stations
path_fig = fullfile(path_main,['downscaling_output_' name_model],'figures_PCR');
path_file = fullfile(path_main,['downscaling_output_' name_model]);
path_obs = fullfile(path_main,'obs_data');
path_downmodel = fullfile(path_main,['downscaling_models_' name_model]);
path_shp_file = fullfile(path_main,'/matlab_code_git/visualization/world_borders/ne_10m_admin_0_countries.shp'); 
suffix = '_NA_020';

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

dsValue = dsEValue_mat + u_mat;

clear dsEValue_mat u_mat

%% filter stations for sufficient data coverage

disp('filter stations for sufficient data coverage')

% Preallocate arrays
nStations = height(metaTable);
length_year = zeros(nStations,1);

% Loop only once over stations to compute timeseries length
for i_ID = 1:nStations
    flag_ID = ismember(obsTable.ID, metaTable.ID(i_ID));
    obsTable_tmp = obsTable(flag_ID, :);
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
dsValue(idx_remove,:) = [];
PI_mat(idx_remove,:,:) = [];

disp('     Done.')

%% load observations
load(fullfile(path_file,[name_var '_observations' suffix '.mat']),'tgt_ESM')

%% plot map with stations

plot_station_map_long(metaTable, lim_lat, lim_lon, path_shp_file, ...
    fullfile(path_fig, ['longstations_' name_var '.png']), [name_var ' GHCN-m stations >' int2str(n_min_yr) ' years' ]);

%% select stations to plot
close all

switch name_var
    case 'pr'
        list_point = {'USC00043747','USC00200032'};
    case 'tas'
        list_point = {'USC00351862','USW00014742'};
end

% Loop over points
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
    dsVal_tmp       = dsValue(flag_ID_meta, :);
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

    plot_downscaled_timeline_point(name_var, unit_var, time_ESM, time_obs_tmp, ...
        obs_value_tmp, ESM_tmp, dsVal_tmp, dsPIinf_tmp, dsPIsup_tmp, ...
        path_fig, suffix_region, name_model, res_plot, [250 255])

    plot_downscaled_timeline_point(name_var, unit_var, time_ESM, time_obs_tmp, ...
        obs_value_tmp, ESM_tmp, dsVal_tmp, dsPIinf_tmp, dsPIsup_tmp, ...
        path_fig, suffix_region, name_model, res_plot, [800 805])

    plot_downscaled_timeline_point(name_var, unit_var, time_ESM, time_obs_tmp, ...
        obs_value_tmp, ESM_tmp, dsVal_tmp, dsPIinf_tmp, dsPIsup_tmp, ...
        path_fig, suffix_region, name_model, res_plot, [1295 1300])

    plot_downscaled_timeline_point(name_var, unit_var, time_ESM, time_obs_tmp, ...
        obs_value_tmp, ESM_tmp, dsVal_tmp, dsPIinf_tmp, dsPIsup_tmp, ...
        path_fig, suffix_region, name_model, res_plot, [1750 1755])

    plot_downscaled_timeline_point(name_var, unit_var, time_ESM, time_obs_tmp, ...
        obs_value_tmp, ESM_tmp, dsVal_tmp, dsPIinf_tmp, dsPIsup_tmp, ...
        path_fig, suffix_region, name_model, res_plot, [1945 1950])

    plot_downscaled_timeline_point(name_var, unit_var, time_ESM, time_obs_tmp, ...
        obs_value_tmp, ESM_tmp, dsVal_tmp, dsPIinf_tmp, dsPIsup_tmp, ...
        path_fig, suffix_region, name_model, res_plot, [1995 2000])

    plot_downscaled_timeline_point(name_var, unit_var, time_ESM, time_obs_tmp, ...
        obs_value_tmp, ESM_tmp, dsVal_tmp, dsPIinf_tmp, dsPIsup_tmp, ...
        path_fig, suffix_region, name_model, res_plot, [2000 2005])
    
end
