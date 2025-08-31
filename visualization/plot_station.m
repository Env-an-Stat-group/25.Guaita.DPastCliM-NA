clear
close all

%% parameters
disp('setting initial parameters...')

rng(812)

% paths and settings
path_main = '/data/pguaita/downscaling/';
addpath(genpath(fullfile(path_main,'matlab_code_git')));
name_model = 'MPI-ESM1-2-LR';
suffix = '_NA_020';
path_fig = fullfile(path_main, ['downscaling_output_' name_model], 'figures_PCR');
path_file = fullfile(path_main, ['downscaling_output_' name_model]);
path_shp_file = fullfile(path_main,'/matlab_code_git/visualization/world_borders/ne_10m_admin_0_countries.shp'); 

%% load grid and define limits
load(fullfile(path_main, ['static_maps/downscaling_grid' suffix '.mat']));
lat = lat(lat >= min(lat) & lat <= max(lat));
lon = lon(lon >= min(lon) & lon <= max(lon));
lim_lat = [min(lat), max(lat)];
lim_lon = [min(lon), max(lon)];

%% load metadata
meta_tas = readtable(fullfile(path_file, ['tas_metadata_table' suffix '.csv']));
meta_pr  = readtable(fullfile(path_file, ['pr_metadata_table' suffix '.csv']));

%% create output folder
if ~exist(path_fig, 'dir')
    mkdir(path_fig)
end

%% plot temperature stations
plot_station_map_with_density(meta_tas, lim_lat, lim_lon, path_shp_file, ...
    fullfile(path_fig, 'temperature_station_locations.png'), 'tas GHCN-m stations');

%% plot precipitation stations
plot_station_map_with_density(meta_pr, lim_lat, lim_lon, path_shp_file, ...
    fullfile(path_fig, 'precipitation_station_locations.png'), 'pr GHCN-m Stations');
