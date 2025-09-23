% LR downscaling
%
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
name_var = 'tas'; % variable name
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

disp('LR')
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

%% get the residual (and fill missing values of residuals on calibration stations)

disp('get the residuals...')
tic

time_ESM_array = time_ESM(:);
residual_mat = nan(height(metaTable),length(time_ESM_array));

switch name_var
    case 'tas'
        for i_ID = 1:height(metaTable)            
            if mod(i_ID,100)==0
                disp(['station number ' int2str(i_ID)])
            end
            if metaTable.flag_cal(i_ID)
                flag_ID = ismember(obsTable.ID,metaTable.ID(i_ID));
                obsTable_tmp = obsTable(flag_ID,:);
                [flag_time,idx_time] = ismember(obsTable_tmp.month_since_0CE,time_ESM_array);
                residual_mat(i_ID,idx_time)=obsTable_tmp.u_lr;
            end
        end
    case 'pr'
        for i_ID = 1:height(metaTable)            
            if mod(i_ID,100)==0
                disp(['station number ' int2str(i_ID)])
            end
            if metaTable.flag_cal(i_ID)
                flag_ID = ismember(obsTable.ID,metaTable.ID(i_ID));
                obsTable_tmp = obsTable(flag_ID,:);
                [flag_time,idx_time] = ismember(obsTable_tmp.month_since_0CE,time_ESM_array);
                residual_mat(i_ID,idx_time)=obsTable_tmp.u_lr;
            end
        end
end

time_residual = (numel(time_ESM)-(12*100-1)):numel(time_ESM); % selecting the period 1915-2014
residual_mat = residual_mat(:,ismember(time_ESM_array,time_residual));

% interpolate missing residuals of the calibration stations from the same
% residuals at the corresponding timestep

warning('off','all')
for i_time = 1:length(time_residual)
    if mod(i_time,100)==0
        disp(['timestep ' int2str(i_time)])
        disp(['Been running for ' num2str(round(toc/60,1)) ' mins'])
    end

    current_time = time_residual(i_time);

    x_obs = metaTable.lon;
    y_obs = metaTable.lat;
    z_Value = residual_mat(:,i_time);
    flag_nan = isnan(z_Value);
    flag_cal = metaTable.flag_cal;

    [y_proj, x_proj] = projfwd(proj, y_obs, x_obs);   % first lat, then lon in degrees
    f = fit([y_proj(not(flag_nan) & flag_cal), x_proj(not(flag_nan) & flag_cal)], z_Value(not(flag_nan) & flag_cal), 'thinplateinterp');
    residual_tmp = f(y_proj, x_proj);
    residual_mat(flag_cal,i_time) = residual_tmp(flag_cal);
end
warning('on','all')

disp(['Done in ' num2str(round(toc/60,1)) ' mins'])

%% calibrate ARMA(p,q) model on the gridded residuals

disp('calibrating arma model')
tic

disp('transforming the residuals through Spatial Error Model')

% Define seasonal flags
%flag_mam = false(size(residual_mat,2),1);
%flag_jja = false(size(residual_mat,2),1);
%flag_son = false(size(residual_mat,2),1);
%flag_djf = false(size(residual_mat,2),1);
%flag_mam([3:12:end 4:12:end 5:12:end]) = true;
%flag_jja([6:12:end 7:12:end 8:12:end]) = true;
%flag_son([9:12:end 10:12:end 11:12:end]) = true;
%flag_djf([12:12:end 1:12:end 2:12:end]) = true;

% Initialize accumulator for dR matrices and a counter
dR_sum = zeros(size(residual_mat,1));
count = 0;

eps_mat = nan(size(residual_mat));
%[lambda, ~, eps_mat(:,flag_mam), W_mam, threshold_best_mam] = fit_SEM_MLE_fmincon(residual_mat(:,flag_mam), metaTable, proj);
%disp(['MAM - Lambda = ' num2str(lambda) '; thresh = ' num2str(threshold_best_mam/1000) ' km']);
%[lambda, ~, eps_mat(:,flag_jja), W_jja, threshold_best_jja] = fit_SEM_MLE_fmincon(residual_mat(:,flag_jja), metaTable, proj);
%disp(['JJA - Lambda = ' num2str(lambda) '; thresh = ' num2str(threshold_best_jja/1000) ' km']);
%[lambda, ~, eps_mat(:,flag_son), W_son, threshold_best_son] = fit_SEM_MLE_fmincon(residual_mat(:,flag_son), metaTable, proj);
%disp(['SON - Lambda = ' num2str(lambda) '; thresh = ' num2str(threshold_best_son/1000) ' km']);
%[lambda, ~, eps_mat(:,flag_djf), W_djf, threshold_best_djf] = fit_SEM_MLE_fmincon(residual_mat(:,flag_djf), metaTable, proj);
%disp(['DJF - Lambda = ' num2str(lambda) '; thresh = ' num2str(threshold_best_djf/1000) ' km']);
[lambda, ~, eps_mat, W, threshold_best] = fit_SEM_MLE_fmincon(residual_mat, metaTable, proj);
disp(['Lambda = ' num2str(lambda) '; thresh = ' num2str(threshold_best/1000) ' km']);

%%

n_lag = 10;
pcorr = zeros(size(residual_mat,1),n_lag+1);
qcorr = zeros(size(residual_mat,1),n_lag+1);
for i_ID = 1:size(residual_mat,1)
    if metaTable.flag_cal(i_ID)
        timeseries = residual_mat(i_ID,:)';
        pcorr(i_ID,:)=parcorr(timeseries,NumLags=n_lag);
        qcorr(i_ID,:)=autocorr(timeseries,NumLags=n_lag);
    end
end

pchange=diff(prctile(pcorr,50,1)>0.05);
qchange=diff(prctile(qcorr,50,1)>0.05);
p = find(pchange(2:end)==-1,1,'first');
q = find(qchange(2:end)==-1,1,'first');

if isempty(p)
    p=1;
end

if isempty(q)
    q=1;
end

disp(['p=' int2str(p) ', q=' int2str(q)])
model = arima(p, 0, q);

disp('for PCR...')
list_ARMA = cell(size(eps_mat,1),1);
flag_cal = metaTable.flag_cal;
parfor i_ID = 1:size(eps_mat,1)
    if flag_cal(i_ID)
        timeseries = eps_mat(i_ID,:)';
        % Estimate ARMA model
        list_ARMA{i_ID} = estimate(model, timeseries,'Display','off');
    end
end

disp(['Done in ' num2str(round(toc/60,1)) ' mins'])

%% for every month, downscale the data in every station and aggregate in the same matrix

disp('downscale data for every month and aggregate in same matrix')
tic

% initialize matrix to store the predicted mean of the downscaled product
EOds_hat_mat_lr = nan(height(metaTable),12,size(time_ESM,2));
PI_mat_lr = nan(height(metaTable),12,size(time_ESM,2),2);

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

    % downscale ESM maps

    % calculate averages in every station
    mu_gO_local = varfun(@mean, obsTable_mth, 'InputVariables', 'mu_gO', 'GroupingVariables', 'ID');
    O_t_local = varfun(@mean, obsTable_mth, 'InputVariables', 'O_t', 'GroupingVariables', 'ID');
    mu_gO_local = mu_gO_local(:,[1 3]);
    O_t_local = O_t_local(:,[1 3]);
    mu_gO_local.Properties.VariableNames(2) = "mu_gO";
    O_t_local.Properties.VariableNames(2) = "O_t";

    % from LR
    [EOds_hat_mat_mth_lr, PI_mat_mth_lr] = ...
        ds_ESM_mat_lr_v1(tgt_ESM_mth, tgt_lon, tgt_lat, time_ESM_mth, lm_list_lr, mu_gO_local, O_t_local, flag_cal,...
        metaTable, name_var);

    % aggregate in same matrix
    EOds_hat_mat_lr(:,i_mth,:) = EOds_hat_mat_mth_lr;
    PI_mat_lr(:,i_mth,:,:) = PI_mat_mth_lr;
    disp(['Been running for ' num2str(round(toc/60,1)) ' mins'])
end

% reshape matrix so that it's station x time
EOds_hat_mat_lr = reshape(EOds_hat_mat_lr,size(EOds_hat_mat_lr,1),[]);
PI_mat_lr = reshape(PI_mat_lr,size(PI_mat_lr,1),[],2);

disp(['Done in ' num2str(round(toc/60,1)) ' mins'])

%% generate residuals from the ARMA model

disp('generate residuals from ARMA and add them to the predicted average')
tic

eps_ARMA = zeros(size(EOds_hat_mat_lr));

% Forecast one step ahead
start_time = max(p,q)+1;

n_iter = 48; % subgroup size for iterations
for i_iter = 1:ceil(height(metaTable)/n_iter)
    idx_ini = n_iter*(i_iter-1)+1;
    idx_fin = min(height(metaTable),n_iter*i_iter);
    idx_iter= idx_ini:idx_fin;

    eps_ARMA_tmp = compute_eps_arma(metaTable(idx_iter,:), list_ARMA(idx_iter), time_ESM, start_time, p, q);
    eps_ARMA(idx_iter,:)=eps_ARMA_tmp;
    disp(['iteration ' int2str(i_iter) ' #station(' int2str(idx_ini) '-' int2str(idx_fin) ') done'])
    disp(['---> Done. Been running for ' num2str(round(toc/60,1)) ' mins'])
end

disp('transform back ARMA residuals through SEM')

%[epsilon] = inverse_SEM(eps_ARMA, W, lambda);
%[u_mat_lr] = inverse_SEM_season(eps_ARMA, W_mam, W_jja, W_son, W_djf, lambda);
[u_mat_lr] = inverse_SEM(eps_ARMA, W, lambda);

disp(['Done in ' num2str(round(toc/60,1)) ' mins'])

%% convert output variables to single float

disp('converting output to single float and reshape variables')

% reshape the variables
time_ESM = reshape(time_ESM,[],1);

%% save the mean predicted downscaling dataset with LR (and prediction intervals)

% save
var_save = {'lat','lon','name_model','name_var','name_experiment',...
    'EOds_hat_mat_lr','u_mat_lr','PI_mat_lr',...
    'time_ESM','unit_var'};
var_desc = {
    'Latitude coordinates of the target grid for downscaling.', ...
    'Longitude coordinates of the target grid for downscaling.',...
    'Name of the model that is being downscaled',...
    'Name of the variable that is being downscaled',...
    'Name of the PMIP experiment that is being downscaled',...
    'Downscaled variable matrix (mean prediction) from LR (station x time)',...
    'Residuals calculated with ARMA model and transformed back through SEM',...
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

save_name = fullfile(path_output, ...
    [name_var '_' name_model '_' name_experiment '_raw_LRdownscaled' suffix]);
    
% Save the variables
save(save_name, var_save{:}, 'info_var', 'main_description','-v7.3');

% Save the dataset downscaled with LR to NetCDF

% File name for NetCDF
ncfile_lr = fullfile(path_output, ...
    [name_var '_' name_model '_' name_experiment '_raw_LRdownscaled' suffix '.nc']);

% Create the NetCDF file
ncid = netcdf.create(ncfile_lr, 'NETCDF4'); % Create a new NetCDF file

% Main description as global attribute
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'description', ...
    ['This file contains the ' name_model ' data for ' name_var ' downscaled with LR.']);

% Define dimensions
dim_lon = netcdf.defDim(ncid, 'lon', length(tgt_lon));
dim_lat = netcdf.defDim(ncid, 'lat', length(tgt_lat));
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
