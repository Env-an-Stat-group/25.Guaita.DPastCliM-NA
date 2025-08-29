% LR downscaling
%
% PR GUAITA - Jan 2025

clear
close all

%% parameters (that you should most likely change)

disp('Starting downscaling...')
rng(812)

% parameters
path_main = '/data/pguaita/downscaling/';
addpath(genpath(path_main));
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
path_shp_file = fullfile(path_main,'/matlab_code/visualization/world_borders/ne_10m_admin_0_countries.shp'); 
suffix = '_NA_020';

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

% load metaTable
i_mth = mth_interval(1,:);    
i_mth_txt = [int2str(mth_interval(1,1)) '-' int2str(mth_interval(1,end))];
path_loadfile_tmp = fullfile(path_downmodel,[name_var '_PCR_original_data_mth=' i_mth_txt suffix]);
load(path_loadfile_tmp,'metaTable')

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

%% create downscaled maps (LR)

disp('create downscaled maps and observation maps (LR)')
tic

save_name = fullfile(path_output, ...
    [name_var '_' name_model '_' name_experiment '_raw_LRdownscaled' suffix '.mat']);
load(save_name)
disp('Loaded raw (ungridded) downscaled data')

% get the predicted values by adding u_mat to the predicted average
dsValue_mat_lr = dsEValue_mat_lr + u_mat_lr;

% number of iteration per cycle (for the parfor loop) - every loop saves
% one file
n_iter = 1200;

% filter on calibration
flag_cal = metaTable.flag_cal;
dsValue_mat_lr = dsValue_mat_lr(flag_cal,:);
PI_mat_lr = PI_mat_lr(flag_cal,:,:);
metaTable = metaTable(flag_cal,:);

% define coordinates of interpolated stations
x_lr = metaTable.lon;
y_lr = metaTable.lat;

% aggregate stations that are within the same node grid
n_stations = length(y_lr);
coords = [y_lr, x_lr];

% Flatten grid for comparison
grid_coords = [tgt_latgrid(:), tgt_longrid(:)];

% Preallocate
assigned_grid_idx = nan(n_stations,1);

% Assign each station to nearest grid point
for i = 1:n_stations
    [~, idx] = min((tgt_latgrid(:) - y_lr(i)).^2 + (tgt_longrid(:) - x_lr(i)).^2);
    assigned_grid_idx(i) = idx;
end

% Find unique grid nodes that received at least one station
unique_grid_idx = unique(assigned_grid_idx);

% Preallocate outputs
n_unique = numel(unique_grid_idx);
dsValue_mat_filtered = nan(n_unique, size(dsValue_mat_lr,2));
PI_mat_filtered = nan(n_unique, size(PI_mat_lr,2),2);
lat_filtered = nan(n_unique,1);
lon_filtered = nan(n_unique,1);
metaTable_filtered = table();

% Aggregate station data by grid node
for k = 1:n_unique
    idx = find(assigned_grid_idx == unique_grid_idx(k));
    
    % Mean of dsValue_mat_lr for all stations in this grid node
    dsValue_mat_filtered(k,:) = mean(dsValue_mat_lr(idx,:), 1, 'omitnan');
    PI_mat_filtered(k,:,1) = mean(PI_mat_lr(idx,:,1), 1, 'omitnan');
    PI_mat_filtered(k,:,2) = mean(PI_mat_lr(idx,:,2), 1, 'omitnan');
    
    % Assign average coordinates (optional: you could also use grid_coords(unique_grid_idx(k),:))
    lat_filtered(k) = mean(y_lr(idx));
    lon_filtered(k) = mean(x_lr(idx));
    
    % combine metadata (keep only the first station)
    metaTable_filtered = [metaTable_filtered; metaTable(idx(1),:)]; 
end

metaTable_filtered.lat = lat_filtered;
metaTable_filtered.lon = lon_filtered;
metaTable = metaTable_filtered;

dsValue_mat_lr = dsValue_mat_filtered;
PI_mat_lr = PI_mat_filtered;

disp('Aggregated stations that were in the same grid (and filtered them through the calibration dataset)')
disp(['---> had ' int2str(length(x_lr)) ' calibration stations'])
x_lr = metaTable.lon;
y_lr = metaTable.lat;
disp(['---> now have ' int2str(length(x_lr)) ' stations'])

[y_proj, x_proj] = projfwd(proj, y_lr, x_lr);  % first lat, then lon in degrees
[ygrid_proj, xgrid_proj] = projfwd(proj, tgt_latgrid(:), tgt_longrid(:)); 

for i_iter = 1:ceil(numel(time_ESM)/n_iter)
    % define time indexing for slicing
    idx_ini = n_iter*(i_iter-1)+1;
    idx_fin = min(numel(time_ESM),n_iter*i_iter);
    idx_iter= idx_ini:idx_fin;
    dsValue_mat_slice   = dsValue_mat_lr(:,idx_iter);
    PI_mat_slice        = PI_mat_lr(:,idx_iter,:);

    % initialize variable to store final maps
    dsValue_map_lr      = nan(length(tgt_lon)*length(tgt_lat),numel(idx_iter));
    PI_map_lr           = nan([length(tgt_lon)*length(tgt_lat),numel(idx_iter), 2]);
    PI_map_inf_lr       = nan([length(tgt_lon)*length(tgt_lat),numel(idx_iter)]);
    PI_map_sup_lr       = nan([length(tgt_lon)*length(tgt_lat),numel(idx_iter)]);

    parfor i_time = 1:length(idx_iter)
        % pick z value for the maps
        z_lr = dsValue_mat_slice(:,i_time);
        z_lr_PI_inf = PI_mat_slice(:,i_time,1);
        z_lr_PI_sup = PI_mat_slice(:,i_time,2);
            
        % Interpolate
        F = scatteredInterpolant(x_proj, y_proj, z_lr, 'natural', 'nearest');
        dsValue_map_lr(:,i_time) = F(xgrid_proj, ygrid_proj);
        
        F_inf = scatteredInterpolant(x_proj, y_proj, z_lr, 'natural', 'nearest');
        PI_map_inf_lr(:,i_time) = F_inf(xgrid_proj, ygrid_proj);
        
        F_sup = scatteredInterpolant(x_proj, y_proj, z_lr_PI_sup, 'natural', 'nearest');
        PI_map_sup_lr(:,i_time) = F_sup(xgrid_proj, ygrid_proj);        
    end
    dsValue_map_lr(not(flag_land(:)),:)=nan;
    PI_map_inf_lr(not(flag_land(:)),:)=nan;
    PI_map_sup_lr(not(flag_land(:)),:)=nan;
    disp(['iteration ' int2str(i_iter) ' #station(' int2str(idx_ini) '-' int2str(idx_fin) ') done'])
    disp(['---> Done. Been running for ' num2str(round(toc/60,1)) ' mins'])

    PI_map_lr(:,:,1) = PI_map_inf_lr;
    PI_map_lr(:,:,2) = PI_map_sup_lr;
    clear PI_map_inf_lr PI_map_sup_lr
    
    % convert output variables to single float
    disp('converting output to single float')
    dsValue_map_lr = single(round(dsValue_map_lr,2));
    PI_map_lr = single(round(PI_map_lr,2));

    %% Adjust the gridded product to match the scaling of the target ESM data
    
    disp('adjusting downscaling gridded products to standard deviation of the ESM data')
    n_year_mov = 30;
    
    tgt_ESM = reshape(tgt_ESM,length(tgt_lon)*length(tgt_lat),[]);
    std_correction = movstd(tgt_ESM(:,idx_iter)-movmean(tgt_ESM(:,idx_iter), n_year_mov, 2, 'omitnan'), n_year_mov, 0, 2, 'omitnan') ./...
        movstd(dsValue_map_lr - movmean(dsValue_map_lr, n_year_mov, 2, 'omitnan'), n_year_mov, 0, 2, 'omitnan');
    dsValue_map_lr = movmean(dsValue_map_lr, n_year_mov, 2, 'omitnan') + ...
        (dsValue_map_lr - movmean(dsValue_map_lr, n_year_mov, 2, 'omitnan')) .* std_correction ;

    %% convert output variables to single float

    disp('converting output to single float and reshape variables')
        
    % reshape the variables
    dsValue_map_lr = reshape(dsValue_map_lr,length(tgt_lon),length(tgt_lat),[]);
    PI_map_lr = reshape(PI_map_lr,length(tgt_lon),length(tgt_lat),[],2);
    time_ESM_array = reshape(time_ESM,[],1);
    time_array = time_ESM(idx_iter);

%% save the dataset downscaled with LR
    
    % save
    var_save = {'lat','lon','name_model','name_var','name_experiment',...
        'dsValue_map_lr',...
        'time_array','unit_var'};
    var_desc = {
        'Latitude coordinates of the target grid for downscaling.', ...
        'Longitude coordinates of the target grid for downscaling.',...
        'Name of the model that is being downscaled',...
        'Name of the variable that is being downscaled',...
        'Name of the PMIP experiment that is being downscaled',...
        'Downscaled variable map (realizations) from LR (lon x lat x month)',...
        'Time arrayin months since 0 CE (month x year, corresponds to the dimension of the maps)',...
        'Units for the output (e.g. °C, or mm/day)'
        };
    var_unit = {
        'degree (°)', ...
        'degree (°)', ...
        '/',...
        '/',...
        '/',...
        unit_var,...
        'Months since 0 CE',...
        '/'
        };
    
    main_description = ['This file contains the ' name_model ' data for ' name_var ' downscaled with LR.'];
    info_var = struct();
    for i_var = 1:length(var_save)
        info_var(i_var).variable_name = var_save{i_var};
        info_var(i_var).description = var_desc{i_var};
        info_var(i_var).unit = var_unit{i_var};
    end
    
    save_name = fullfile(path_output, ...
        [name_var '_' name_model '_' name_experiment '_LRdownscaled' suffix '_' int2str(idx_iter(1)) '-' int2str(idx_iter(end))]);
        
    % Save the variables
    save(save_name, var_save{:}, 'info_var', 'main_description','-v7.3');
    
    % Save the dataset downscaled with LR to NetCDF
    
    % File name for NetCDF
    ncfile_lr = fullfile(path_output, ...
        [name_var '_' name_model '_' name_experiment '_LRdownscaled' suffix '_' int2str(idx_iter(1)) '-' int2str(idx_iter(end)) '.nc']);
    
    % Create the NetCDF file
    ncid = netcdf.create(ncfile_lr, 'NETCDF4'); % Create a new NetCDF file
    
    % Main description as global attribute
    netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'description', ...
        ['This file contains the ' name_model ' data for ' name_var ' downscaled with LR.']);
    
    % Define dimensions
    dim_lon = netcdf.defDim(ncid, 'lon', size(dsValue_map_lr, 1));
    dim_lat = netcdf.defDim(ncid, 'lat', size(dsValue_map_lr, 2));
    dim_time = netcdf.defDim(ncid, 'time', size(dsValue_map_lr, 3));
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
                                netcdf.defVarDeflate(ncid, varid, 1, 1, 3);
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
                                    case length(time_array)
                                        varid = netcdf.defVar(ncid, var_name, NC_CLASS, dim_time);
                                        netcdf.putVar(ncid, varid, var_data);
                                        netcdf.putAtt(ncid, varid, 'description', var_desc{i_var});
                                        netcdf.putAtt(ncid, varid, 'unit', unit_text);
                                end
                            end
                        end
                    case 3
                        varid = netcdf.defVar(ncid, var_name, NC_CLASS, [dim_lon,dim_lat,dim_time]);
                        netcdf.defVarDeflate(ncid, varid, 1, 1, 5);
                        netcdf.putVar(ncid, varid, var_data);
                        netcdf.putAtt(ncid, varid, 'description', var_desc{i_var});
                        netcdf.putAtt(ncid, varid, 'unit', unit_text);
                    case 4
                        varid = netcdf.defVar(ncid, var_name, NC_CLASS, [dim_lon,dim_lat,dim_time,dim_bound]);
                        netcdf.defVarDeflate(ncid, varid, 1, 1, 5);
                        netcdf.putVar(ncid, varid, var_data);
                        netcdf.putAtt(ncid, varid, 'description', var_desc{i_var});
                        netcdf.putAtt(ncid, varid, 'unit', unit_text);
                end
            end
        end
    end
    
    % Close the NetCDF file
    netcdf.close(ncid);
    
    %% save the mean predicted downscaling dataset with LR (and prediction intervals)
    
    % save
    var_save = {'lat','lon','name_model','name_var','name_experiment',...
        'PI_map_lr',...
        'time_array','unit_var'};
    var_desc = {
        'Latitude coordinates of the target grid for downscaling.', ...
        'Longitude coordinates of the target grid for downscaling.',...
        'Name of the model that is being downscaled',...
        'Name of the variable that is being downscaled',...
        'Name of the PMIP experiment that is being downscaled',...
        'Prediction interval (mean prediction) from LR (lon x lat x month x boundary)',...
        'Time arrayin months since 0 CE (month x year, corresponds to the dimension of the maps)',...
        'Units for the output (e.g. °C, or mm/day)'
        };
    var_unit = {
        'degree (°)', ...
        'degree (°)', ...
        '/',...
        '/',...
        '/',...
        unit_var,...
        'Months since 0 CE',...
        '/'
        };
    
    main_description = ['This file contains the ' name_model ' data for ' name_var ' downscaled with LR.'];
    info_var = struct();
    for i_var = 1:length(var_save)
        info_var(i_var).variable_name = var_save{i_var};
        info_var(i_var).description = var_desc{i_var};
        info_var(i_var).unit = var_unit{i_var};
    end
    
    save_name = fullfile(path_output, ...
        [name_var '_' name_model '_' name_experiment '_PI_LRdownscaled' suffix '_' int2str(idx_iter(1)) '-' int2str(idx_iter(end)) ]);
        
    % Save the variables
    save(save_name, var_save{:}, 'info_var', 'main_description','-v7.3');
    
    % Save the dataset downscaled with LR to NetCDF
    
    % File name for NetCDF
    ncfile_lr = fullfile(path_output, ...
        [name_var '_' name_model '_' name_experiment '_PI_LRdownscaled' suffix '_' int2str(idx_iter(1)) '-' int2str(idx_iter(end)) '.nc']);
    
    % Create the NetCDF file
    ncid = netcdf.create(ncfile_lr, 'NETCDF4'); % Create a new NetCDF file
    
    % Main description as global attribute
    netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'description', ...
        ['This file contains the ' name_model ' data for ' name_var ' downscaled with LR.']);
    
    % Define dimensions
    dim_lon = netcdf.defDim(ncid, 'lon', size(dsValue_map_lr, 1));
    dim_lat = netcdf.defDim(ncid, 'lat', size(dsValue_map_lr, 2));
    dim_time = netcdf.defDim(ncid, 'time', size(dsValue_map_lr, 3));
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
                                netcdf.defVarDeflate(ncid, varid, 1, 1, 3);
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
                                    case length(time_array)
                                        varid = netcdf.defVar(ncid, var_name, NC_CLASS, dim_time);
                                        netcdf.putVar(ncid, varid, var_data);
                                        netcdf.putAtt(ncid, varid, 'description', var_desc{i_var});
                                        netcdf.putAtt(ncid, varid, 'unit', unit_text);
                                end
                            end
                        end
                    case 3
                        varid = netcdf.defVar(ncid, var_name, NC_CLASS, [dim_lon,dim_lat,dim_time]);
                        netcdf.defVarDeflate(ncid, varid, 1, 1, 5);
                        netcdf.putVar(ncid, varid, var_data);
                        netcdf.putAtt(ncid, varid, 'description', var_desc{i_var});
                        netcdf.putAtt(ncid, varid, 'unit', unit_text);
                    case 4
                        varid = netcdf.defVar(ncid, var_name, NC_CLASS, [dim_lon,dim_lat,dim_time,dim_bound]);
                        netcdf.defVarDeflate(ncid, varid, 1, 1, 5);
                        netcdf.putVar(ncid, varid, var_data);
                        netcdf.putAtt(ncid, varid, 'description', var_desc{i_var});
                        netcdf.putAtt(ncid, varid, 'unit', unit_text);
                end
            end
        end
    end
    
    % Close the NetCDF file
    netcdf.close(ncid);
    
end
    
disp(['Done in ' num2str(round(toc/60,1)) ' mins'])
disp('finished')

    
