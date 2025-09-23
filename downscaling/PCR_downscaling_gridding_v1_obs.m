%PCR downscaling - code description
%
% PR GUAITA - Jan 2025

clear
close all

%% parameters (that you should most likely change)

disp('Starting downscaling...')
rng(812)

% parameters
path_main = '/data/pguaita/downscaling/';
addpath(genpath(fullfile(path_main,'matlab_code_git')));
name_model = 'MPI-ESM1-2-LR'; % model name
name_var = 'pr'; % variable name
name_experiment = 'past2k';
% starting and ending year for the PMIP expeirment considered.
% weird enough, but the past experiments start counting years from 7000 CE
year_start_experiment = 7000;
year_end_experiment = 8850;
year_base_ESM_experiment_model = 7000; %first year in the experiment
year_base_ESM_experiment_real = 0; %first year in the experiment
% starting and ending year for the historical period
year_start_hist = 1850;
year_end_hist = 2014;
year_base_ESM_hist = 1850; % first year in the experiment
path_ESM = fullfile(path_main,'CMIP_data');
path_output = fullfile(path_main, ['downscaling_output_' name_model]);
path_fig = fullfile(path_main,['downscaling_models_' name_model],'figures_PCR');
path_downmodel = fullfile(path_main,['downscaling_models_' name_model]);
path_shp_file = fullfile(path_main,'/matlab_code_git/visualization/world_borders/ne_10m_admin_0_countries.shp'); 
suffix = '_NA_020';

disp('obs gridding')
disp(suffix)
disp(name_var)
%% parameters that you most likely should not change

% define projection
proj = projcrs(5070);  % NAD83 / Conus Albers

switch name_var
    case 'tas'
        mth_interval = (1:12)';
        unit_var = '°C';
    case 'pr'
        mth_interval = (1:12)';
        unit_var = 'mm/day';
end

%% setup variables

disp('Setting up variables...')

if not(exist(path_fig,'dir'))
    mkdir(path_fig)
end

% load target grid resolution
load(fullfile(path_main,['static_maps/downscaling_grid' suffix '.mat']));

% filters for the domain 
lim_lat = [min(lat),max(lat)];
lim_lon = [min(lon),max(lon)];
filter_lat = [min(lat),max(lat)];
filter_lon = [min(lon),max(lon)];

tgt_lon = lon;
tgt_lat = lat;

% filter the target grid
flag_tgt_lat = tgt_lat>= filter_lat(1) & tgt_lat<=filter_lat(2);
flag_tgt_lon = tgt_lon>= filter_lon(1) & tgt_lon<=filter_lon(2);
tgt_lat = tgt_lat(flag_tgt_lat);
tgt_lon = tgt_lon(flag_tgt_lon);

% define lon and lat grid
[tgt_longrid,tgt_latgrid] = ndgrid(tgt_lon,tgt_lat);

% load landsea mask
landseamask = ncread(fullfile(path_main,'static_maps','IMERG_land_sea_mask.nc'),'landseamask');
flag_land = not(landseamask == 100);
lat_land = ncread(fullfile(path_main,'static_maps','IMERG_land_sea_mask.nc'),'lat');
lon_land = ncread(fullfile(path_main,'static_maps','IMERG_land_sea_mask.nc'),'lon');
[flag_land,lon_land,lat_land] = cmipcoord_2_map(flag_land,lon_land,lat_land);
% regrid land-sea mask
[lonland_grid,latland_grid] = ndgrid(lon_land,lat_land);
f_int_land = griddedInterpolant(lonland_grid,latland_grid,double(flag_land),'nearest');
flag_land = logical(f_int_land(tgt_longrid,tgt_latgrid));


%% load and stack ESM data for the specified model

disp('load and stack ESM data for the specified model')
tic

% load nc files for the historical experiment
file_list_ESM = dir(fullfile(path_ESM,[name_var '_Amon_' name_model '_historical_*.nc']));

[tgt_ESM_hist, time_ESM_hist] = ...
    process_ESM_Data_v1(file_list_ESM, path_ESM, name_model, name_var, ...
        year_start_hist, year_end_hist, year_base_ESM_hist, filter_lat, filter_lon, tgt_lat, tgt_lon);

% load nc files for the past-time experiment
file_list_ESM = dir(fullfile(path_ESM,[name_var '_Amon_' name_model '_' name_experiment '_*.nc']));

[tgt_ESM_past, time_ESM_past] = ...
    process_ESM_Data_v1(file_list_ESM, path_ESM, name_model, name_var, ...
        year_start_experiment, year_end_experiment, year_base_ESM_experiment_model, filter_lat, filter_lon, tgt_lat, tgt_lon);

% the time array for the past ESM is still referred to the experiment time
% (e.g. 7000 CE), which is not the real year.
% correct the time array for the past experiment so that it matches the
% actual real time array (starting from the real year for the experiment, 
% e.g. 850 CE for past1000, or 0 CE for past2k). This is till in #months
% since 0 CE.
time_ESM_past = time_ESM_past + (year_base_ESM_experiment_real-year_base_ESM_experiment_model)*12;

% concatenate the two ESM matrixes and the two time arrays
tgt_ESM = cat(3,tgt_ESM_hist,tgt_ESM_past);
time_ESM = cat(1,time_ESM_hist,time_ESM_past);
% sort them by time
[time_ESM,ind_timesort] = sort(time_ESM);
tgt_ESM=tgt_ESM(:,:,ind_timesort);

clear tgt_ESM_* time_ESM_*

tgt_ESM = reshape(tgt_ESM,length(tgt_lon),length(tgt_lat),12,[]);
time_ESM = reshape(time_ESM,12,[]);

disp(['Done in ' num2str(round(toc/60,1)) ' mins'])

%% aggregate observation table (which contains also the downscaled product)

disp('Aggregating obs table')
tic

for i_mth_interval = 1:size(mth_interval,1)
    
    i_mth = mth_interval(i_mth_interval,:);    
    i_mth_txt = [int2str(mth_interval(i_mth_interval,1)) '-' int2str(mth_interval(i_mth_interval,end))];
    disp(['doing month ' i_mth_txt])

    tgt_ESM_mth = reshape(squeeze(tgt_ESM(:,:,i_mth,:)),[],size(tgt_ESM,4)*length(i_mth));
    time_ESM_mth = reshape(squeeze(time_ESM(i_mth,:)),[],size(tgt_ESM,4)*length(i_mth));

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

disp(['Done in ' num2str(round(toc/60,1)) ' mins'])

%% create observation maps

disp('create observation maps')
tic

% initialize variable to store final maps
obsValue_map        = nan(length(tgt_lon)*length(tgt_lat),numel(time_ESM));
obsValue_test_map   = nan(length(tgt_lon)*length(tgt_lat),numel(time_ESM));
flag_nantest = false(1,numel(time_ESM));

% size of iterations for the parfor loop
n_iter = 48;
for i_iter = 1:ceil(numel(time_ESM)/n_iter)
    idx_ini = n_iter*(i_iter-1)+1;
    idx_fin = min(numel(time_ESM),n_iter*i_iter);
    idx_iter= idx_ini:idx_fin;
    obsTable_slice = obsTable(ismember(obsTable.month_since_0CE,time_ESM(idx_iter)),:);
    parfor i_time = idx_iter
        flag_time = ismember(obsTable_slice.month_since_0CE,time_ESM(i_time));
        obsTable_tmp = obsTable_slice(flag_time,:);
        obsTable_tmp_test = obsTable_tmp(obsTable_tmp.flag_test,:);
    
        % observations
        x_obs = obsTable_tmp.lon;
        y_obs = obsTable_tmp.lat;
        z_Value = obsTable_tmp.Value;    
    
        x_obs_test = obsTable_tmp_test.lon;
        y_obs_test = obsTable_tmp_test.lat;
        z_Value_test = obsTable_tmp_test.Value;
    
    
        if height(obsTable_tmp)>height(metaTable)*0.1
            [y_proj, x_proj] = projfwd(proj, y_obs, x_obs);  % first lat, then lon in degrees
            [ygrid_proj, xgrid_proj] = projfwd(proj, tgt_latgrid(:), tgt_longrid(:)); 
            F_inf = scatteredInterpolant(x_proj, y_proj, double(z_Value), 'natural', 'nearest');
            obsValue_map(:,i_time) = F_inf(xgrid_proj, ygrid_proj);
        end
        if height(obsTable_tmp_test)>height(metaTable)*0.05
            [y_proj, x_proj] = projfwd(proj, y_obs_test, x_obs_test);  % first lat, then lon in degrees
            [ygrid_proj, xgrid_proj] = projfwd(proj, tgt_latgrid(:), tgt_longrid(:)); 
            F_inf = scatteredInterpolant(x_proj, y_proj, double(z_Value_test), 'natural', 'nearest');
            obsValue_test_map(:,i_time) = F_inf(xgrid_proj, ygrid_proj);
        else
            flag_nantest(i_time) = true;
        end
    end
    disp(['iteration ' int2str(i_iter) ' #timestep(' int2str(idx_ini) '-' int2str(idx_fin) ') done'])
    disp(['---> Done. Been running for ' num2str(round(toc/60,1)) ' mins'])
end

disp(['Done in ' num2str(round(toc/60,1)) ' mins'])

% convert output variables to single float

disp('converting output to single float')
obsValue_map = single(round(obsValue_map,2));
obsValue_test_map = single(round(obsValue_test_map,2));

%% convert output variables to single float

disp('reshape variables')

tgt_ESM = single(round(tgt_ESM,2));

% reshape the variables
obsValue_map = reshape(obsValue_map,length(tgt_lon),length(tgt_lat),[]);
obsValue_test_map = reshape(obsValue_test_map,length(tgt_lon),length(tgt_lat),[]);
tgt_ESM = reshape(tgt_ESM,length(tgt_lon),length(tgt_lat),[]);
time_ESM = reshape(time_ESM,[],1);

%% save the observation maps

% save
var_save = {'lat','lon','metaTable','name_model','name_var','name_experiment',...
    'obsValue_map','obsValue_test_map','tgt_ESM',...
    'time_ESM','unit_var'};
var_desc = {
    'Latitude coordinates of the target grid for downscaling.', ...
    'Longitude coordinates of the target grid for downscaling.',...
    'Metadata table',...
    'Name of the model that is being downscaled',...
    'Name of the variable that is being downscaled',...
    'Name of the PMIP experiment that is being downscaled',...
    'Observation maps using the whole dataset (lon x lat x month)',...
    'Observation maps using only the testing dataset (lon x lat x month)',...
    'Original ESM maps obtained at the target downscaled resolution resolution using GPR (lon x lat x month)',...
    'Time arrayin months since 0 CE',...
    'Units for the output (e.g. °C, or mm/day)'
    };
var_unit = {
    'degree (°)', ...
    'degree (°)', ...
    '/',...
    '/',...
    '/',...
    '/',...
    unit_var,...
    unit_var,...
    unit_var,...
    'Months since 0 CE',...
    '/'
    };

main_description = ['This file contains the observation maps data for ' name_var '.'];
info_var = struct();
for i_var = 1:length(var_save)
    info_var(i_var).variable_name = var_save{i_var};
    info_var(i_var).description = var_desc{i_var};
    info_var(i_var).unit = var_unit{i_var};
end

save_name = fullfile(path_output, ...
    [name_var '_observations' suffix]);
    
% Save the variables
save(save_name, var_save{:}, 'info_var', 'main_description','-v7.3');

% Save the dataset downscaled with PCR to NetCDF

% File name for NetCDF
ncfile_obs = fullfile(path_output, ...
    [name_var '_observations' suffix '.nc']);

% Create the NetCDF file
ncid = netcdf.create(ncfile_obs, 'NETCDF4'); % Create a new NetCDF file

% Main description as global attribute
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'description', ...
    ['This file contains the observation maps data for ' name_var '.']);

% Define dimensions
dim_lon = netcdf.defDim(ncid, 'lon', size(tgt_ESM, 1));
dim_lat = netcdf.defDim(ncid, 'lat', size(tgt_ESM, 2));
dim_time = netcdf.defDim(ncid, 'time', size(tgt_ESM, 3));
dim_site = netcdf.defDim(ncid, 'site', height(metaTable));
dim_bound = netcdf.defDim(ncid, 'bound', 2);
dim_scalar = netcdf.defDim(ncid, 'scalar', 1);

% Define and write other variables to NetCDF
for i_var = 1:length(var_save)
    var_name = var_save{i_var};
    unit_text = var_unit{i_var};
    % skip opt_gpr
    if strcmp(var_name,'opt_gpr') || strcmp(var_name,'unit_var') || strcmp(var_name,'metaTable')
    else
        var_data = eval(var_name); % Evaluate the variable
        
        % Switch based on the data class for correct NetCDF data type
        switch class(var_data)
            case 'single'
                NC_CLASS = 'NC_FLOAT';
            case 'double'
                NC_CLASS = 'NC_DOUBLE';
            case 'char'
                NC_CLASS = 'NC_CHAR';
                % Main description as global attribute
                netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), var_name,var_data);
        end
        
        % Define the variable in the NetCDF file
        if strcmp(NC_CLASS,'NC_CHAR')
        else
            switch ndims(var_data)
                case 2
                    if length(var_data)==1
                        varid = netcdf.defVar(ncid, var_name, NC_CLASS, dim_scalar);
                        netcdf.putVar(ncid, varid, var_data);
                        netcdf.putAtt(ncid, varid, 'description', var_desc{i_var});
                        netcdf.putAtt(ncid, varid, 'unit', unit_text);
                    else
                        if size(var_data,1)>1 && size(var_data,2)>1
                            varid = netcdf.defVar(ncid, var_name, NC_CLASS, [dim_site,dim_time]);
                            netcdf.defVarDeflate(ncid, varid, true, true, 3);
                            netcdf.putVar(ncid, varid, var_data);
                            netcdf.putAtt(ncid, varid, 'description', var_desc{i_var});
                            netcdf.putAtt(ncid, varid, 'unit', unit_text);
                        else
                            switch length(var_data)
                                case length(lat)
                                    varid = netcdf.defVar(ncid, var_name, NC_CLASS, dim_lat);
                                    netcdf.putVar(ncid, varid, var_data);
                                    netcdf.putAtt(ncid, varid, 'description', var_desc{i_var});
                                    netcdf.putAtt(ncid, varid, 'unit', unit_text);
                                case length(lon)
                                    varid = netcdf.defVar(ncid, var_name, NC_CLASS, dim_lon);
                                    netcdf.putVar(ncid, varid, var_data);
                                    netcdf.putAtt(ncid, varid, 'description', var_desc{i_var});
                                    netcdf.putAtt(ncid, varid, 'unit', unit_text);
                                case length(time_ESM)
                                    varid = netcdf.defVar(ncid, var_name, NC_CLASS, dim_time);
                                    netcdf.putVar(ncid, varid, var_data);
                                    netcdf.putAtt(ncid, varid, 'description', var_desc{i_var});
                                    netcdf.putAtt(ncid, varid, 'unit', unit_text);
                            end
                        end
                    end
                case 3
                    varid = netcdf.defVar(ncid, var_name, NC_CLASS, [dim_lon,dim_lat,dim_time]);
                    netcdf.defVarDeflate(ncid, varid, true, true, 5);
                    netcdf.putVar(ncid, varid, var_data);
                    netcdf.putAtt(ncid, varid, 'description', var_desc{i_var});
                    netcdf.putAtt(ncid, varid, 'unit', unit_text);
                case 4
                    varid = netcdf.defVar(ncid, var_name, NC_CLASS, [dim_lon,dim_lat,dim_time,dim_bound]);
                    netcdf.defVarDeflate(ncid, varid, true, true, 5);
                    netcdf.putVar(ncid, varid, var_data);
                    netcdf.putAtt(ncid, varid, 'description', var_desc{i_var});
                    netcdf.putAtt(ncid, varid, 'unit', unit_text);
            end
        end
    end
end

% Close the NetCDF file
netcdf.close(ncid);

% save the metadata table
save_name = fullfile(path_output, ...
    [name_var '_metadata_table' suffix '.csv']);
writetable(metaTable, save_name, 'WriteRowNames', true, 'WriteVariableNames', true);

disp('finished')
