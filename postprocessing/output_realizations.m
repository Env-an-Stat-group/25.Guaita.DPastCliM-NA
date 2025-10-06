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
            range_limit = [0 15];
            range_step = 3;
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
            range_limit = [0 15];
            range_step = 3;
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

% load metaTable
metaTable=readtable(fullfile(path_file,[name_var '_metadata_table' suffix '.csv']));

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

clear EOds_hat_mat EOds_hat_mat_lr u_mat u_mat_lr epsilon 

%% save the predicted downscaling dataset (realizations) with PCR (and prediction intervals)

% save
var_save = {'lat','lon','name_model','name_var','name_experiment',...
    'Ods_hat_mat','PI_mat',...
    'time_ESM','unit_var'};
var_desc = {
    'Latitude coordinates of the target grid for downscaling.', ...
    'Longitude coordinates of the target grid for downscaling.',...
    'Name of the model that is being downscaled',...
    'Name of the variable that is being downscaled',...
    'Name of the PMIP experiment that is being downscaled',...
    'Downscaled variable matrix (realizations) from PCR (station x time)',...
    'Prediction intervals from PCR (station x time x 2 [1=lower bound, 2=upper bound])',...
    'Time array in months since 0 CE',...
    'Units for the output (e.g. °C, or mm/day)'
    };
var_unit = {
    'degree (°)', ...
    'degree (°)', ...
    '/',...
    '/',...
    '/',...
    unit_var,...
    unit_var,...
    unit_var,...
    'Months since 0 CE',...
    '/'
    };

main_description = ['This file contains the ' name_model ' data for ' name_var ' downscaled with PCR.'];
info_var = struct();
for i_var = 1:length(var_save)
    info_var(i_var).variable_name = var_save{i_var};
    info_var(i_var).description = var_desc{i_var};
    info_var(i_var).unit = var_unit{i_var};
end

save_name = fullfile(path_file, ...
    [name_var '_' name_model '_' name_experiment '_statreal_PCRdownscaled' suffix]);
    
% Save the variables
save(save_name, var_save{:}, 'info_var', 'main_description','-v7.3');

% Save the dataset downscaled with PCR to NetCDF

% File name for NetCDF
ncfile = fullfile(path_file, ...
    [name_var '_' name_model '_' name_experiment '_statreal_PCRdownscaled' suffix '.nc']);

% Create the NetCDF file
ncid = netcdf.create(ncfile, 'NETCDF4'); % Create a new NetCDF file

% Main description as global attribute
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'description', ...
    ['This file contains the ' name_model ' data for ' name_var ' downscaled with PCR.']);

% Define dimensions
dim_lon = netcdf.defDim(ncid, 'lon', length(lon));
dim_lat = netcdf.defDim(ncid, 'lat', length(lat));
dim_time = netcdf.defDim(ncid, 'time', length(time_ESM));
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
                            netcdf.defVarDeflate(ncid, varid, true, true, 4);
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
                    varid = netcdf.defVar(ncid, var_name, NC_CLASS, [dim_site,dim_time,dim_bound]);
                    netcdf.defVarDeflate(ncid, varid, true, true, 4);
                    netcdf.putVar(ncid, varid, var_data);
                    netcdf.putAtt(ncid, varid, 'description', var_desc{i_var});
                    netcdf.putAtt(ncid, varid, 'unit', unit_text);
            end
        end
    end
end

% Close the NetCDF file
netcdf.close(ncid);

disp('finished')

%% save the predicted downscaling dataset (realizations) with LR (and prediction intervals)

% save
var_save = {'lat','lon','name_model','name_var','name_experiment',...
    'Ods_hat_mat_lr','PI_mat_lr',...
    'time_ESM','unit_var'};
var_desc = {
    'Latitude coordinates of the target grid for downscaling.', ...
    'Longitude coordinates of the target grid for downscaling.',...
    'Name of the model that is being downscaled',...
    'Name of the variable that is being downscaled',...
    'Name of the PMIP experiment that is being downscaled',...
    'Downscaled variable matrix (realizations) from LR (station x time)',...
    'Prediction intervals from LR (station x time x 2 [1=lower bound, 2=upper bound])',...
    'Time array in months since 0 CE',...
    'Units for the output (e.g. °C, or mm/day)'
    };
var_unit = {
    'degree (°)', ...
    'degree (°)', ...
    '/',...
    '/',...
    '/',...
    unit_var,...
    unit_var,...
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

save_name = fullfile(path_file, ...
    [name_var '_' name_model '_' name_experiment '_statreal_LRdownscaled' suffix]);
    
% Save the variables
save(save_name, var_save{:}, 'info_var', 'main_description','-v7.3');

% Save the dataset downscaled with LR to NetCDF

% File name for NetCDF
ncfile = fullfile(path_file, ...
    [name_var '_' name_model '_' name_experiment '_statreal_LRdownscaled' suffix '.nc']);

% Create the NetCDF file
ncid = netcdf.create(ncfile, 'NETCDF4'); % Create a new NetCDF file

% Main description as global attribute
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'description', ...
    ['This file contains the ' name_model ' data for ' name_var ' downscaled with LR.']);

% Define dimensions
dim_lon = netcdf.defDim(ncid, 'lon', length(lon));
dim_lat = netcdf.defDim(ncid, 'lat', length(lat));
dim_time = netcdf.defDim(ncid, 'time', length(time_ESM));
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
                            netcdf.defVarDeflate(ncid, varid, true, true, 4);
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
                    varid = netcdf.defVar(ncid, var_name, NC_CLASS, [dim_site,dim_time,dim_bound]);
                    netcdf.defVarDeflate(ncid, varid, true, true, 4);
                    netcdf.putVar(ncid, varid, var_data);
                    netcdf.putAtt(ncid, varid, 'description', var_desc{i_var});
                    netcdf.putAtt(ncid, varid, 'unit', unit_text);
            end
        end
    end
end

% Close the NetCDF file
netcdf.close(ncid);

disp('finished')