clear
close all

%% parameters
disp('setting initial parameters...')

rng(812)

% paths and settings
path_main = 'C:\Users\guait\OneDrive - Università Cattolica del Sacro Cuore\PALEON\downscaling\'; %'/data/pguaita/downscaling/';
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
close all

        star_list = {'CA001100120', 'USC00026796', ...
            'USC00047738', 'USW00023188', 'USW00024216','USC00053005',...
            'USC00198367','USC00200146','USC00351862','USW00014742', 'USW00023234'};

plot_station_map_with_density(meta_tas, lim_lat, lim_lon, path_shp_file, ...
    fullfile(path_fig, 'temperature_station_locations.png'), 'tas GHCN-m stations', ...
    star_list);

%% plot precipitation stations

close all

        star_list = {'USC00043747','USC00353827','USC00456096','USC00456678',...
            'USC00252205','USC00096335','USC00200032','USC00200230','USC00080478',...
            'USC00176430'
        };


plot_station_map_with_density(meta_pr, lim_lat, lim_lon, path_shp_file, ...
    fullfile(path_fig, 'precipitation_station_locations.png'), 'pr GHCN-m Stations',...
    star_list);
