% this performs the pc regression using the patterns obtained through pca
%
% PR GUAITA - Jan 2025

clear
close all

%% parameters (that you should most likely change)

rng(812)

% parameters
path_main = '/data/pguaita/downscaling/';
addpath(genpath(path_main));
name_model = 'MPI-ESM1-2-LR'; % model name
name_var = 'tas'; % variable name
path_ESM = fullfile(path_main,'CMIP_data');
path_fig = fullfile(path_main,['downscaling_models_' name_model],'figures_PCR');
path_shp_file = fullfile(path_main,'/matlab_code/visualization/world_borders/ne_10m_admin_0_countries.shp'); 
suffix = '_NA_020';

%% parameters that you might or might not change

% number of iterations to calculate RMSE of the realizations
n_err_iter_lr = 50;            % for the linear regression
n_err_iter_selection = 10;      % for the model selection process
n_err_iter_final = 50;         % for the final PCR calibration

%% parameters that you most likely should not change

% number of principal components by month
switch name_var
    case 'tas'
        frac_threshold = 0.001; % [%] the threshold of explained variance above which a PC is considered valid for regression
    case 'pr'
        frac_threshold = 0.001; % [%] the threshold of explained variance above which a PC is considered valid for regression
end

switch name_var
    case 'tas'
        mth_interval = (1:12)';
        unit_var = '°C';
    case 'pr'
        mth_interval = (1:12)';
        unit_var = 'mm/day';
end

% starting and ending year for calibrating the downscaling model/selecting
% the observational data
year_start   = 1875;
year_end = 2014;
year_base_ESM = 1850;

%number of years for validation and testing (the remaining are for
%acalibration)
n_year_data  = year_end-year_start+1; % number of years that you want to consider starting from the indicated starting year
n_year_val   = 50; % amount of random years that you are holding out for validation/model selection
n_year_test   = 15; % amount of random years that you are holding out for validation/model selection

%% Figures parameters

% figure colorbar limits
switch name_var
        case 'pr'
            var_title = 'pr';
            % absolute value parameters
            abs_limit = [1.75 4.25];
            abs_step = 0.5;
            n_color_abs = 10;
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
            abs_limit = [-10 25];
            abs_step = 5;
            n_color_abs = 14;
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
time_bound = [1875 1910 1945 1980 2015];

res_fig = '-r500';
res_fig_low = '-r100';
res_plot= '-r300';

%% setup variables

if not(exist(path_fig,'dir'))
    mkdir(path_fig);
end

% number of years for training/calibration
n_year_train = n_year_data-n_year_val-n_year_test;

% load target grid resolution
load(fullfile(path_main,['static_maps/downscaling_grid' suffix '.mat']));

% filters for the domain 
lim_lat = [min(lat),max(lat)];
lim_lon = [min(lon),max(lon)];
filter_lat = [min(lat),max(lat)];
filter_lon = [min(lon),max(lon)];

tgt_lon = lon;
tgt_lat = lat;
clear lat lon

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

%% load observations

% load the observation table
load(fullfile(path_main,'obs_data',['GHCNm_' name_var suffix '.mat']))

%filter on time
flag_time_obs = (year_start*12+1)<=obsTable.month_since_0CE & ...
    obsTable.month_since_0CE<=((year_end+1)*12);
obsTable = obsTable(flag_time_obs,:);
% update meta table
flag_ID = ismember(metaTable.ID, obsTable.ID);
metaTable = metaTable(flag_ID,:);

%% aggregate input ESM files
% load nc files
file_list_ESM = dir(fullfile(path_ESM,[name_var '_Amon_' name_model '_historical_*.nc']));

[tgt_ESM, time_ESM] = ...
    process_ESM_Data_v1(file_list_ESM, path_ESM, name_model, name_var, ...
        year_start, year_end, year_base_ESM, filter_lat, filter_lon, tgt_lat, tgt_lon);

% define sizes
tgt_size = [length(tgt_lon) length(tgt_lat)];

% for each timestep interpolate ESM data on the station coordinates and
% save the value
for i_time = 1:length(time_ESM)    
    flag_time = obsTable.month_since_0CE==time_ESM(i_time);
    obsTable_tmp = obsTable(flag_time,:);
    % interpolate coarse_map on the target grid and on the station
    % coordinates
    f_int = griddedInterpolant(tgt_longrid,tgt_latgrid,tgt_ESM(:,:,i_time),'nearest','nearest');
    obsTable.ValueESM(flag_time) = f_int(obsTable_tmp.lon,obsTable_tmp.lat);
end

% correct values lower than 0 for precipitation, in case there are any
switch name_var
    case 'pr'
        obsTable.ValueESM = max(0,obsTable.ValueESM);
end

%% filter and separate for training and validation/model selection
% organize the matrix monthly (4D matrix: lon x lat x mth x year) 
tgt_ESM = reshape(tgt_ESM,size(tgt_ESM,1),size(tgt_ESM,2),12,[]);

ind_val  = randperm(n_year_data);
flag_val = false(1,n_year_data);
flag_test= false(1,n_year_data);
flag_val(ind_val(1:n_year_val)) = true;
flag_test(ind_val((n_year_val+1):(n_year_val+n_year_test))) = true;
tgt_ESM_val = tgt_ESM(:,:,:,flag_val);
flag_cal = not(flag_val | flag_test);
tgt_ESM_cal = tgt_ESM(:,:,:,flag_cal);

%% transform observations to Y, calculate mu_Y and mu_ESM

for i_mth_interval = 1:size(mth_interval,1)
    % filter by specified months
    i_mth = mth_interval(i_mth_interval,:);
    time_mth = 12*year_start+...
        reshape(repmat(i_mth,n_year_data,1)'+12*((1:n_year_data)-1),[],1);    
    flag_mth = ismember(obsTable.month_since_0CE,time_mth);
    obsTable_tmp = obsTable(flag_mth,:);

    % transform observed variable
    switch name_var
        case 'tas'
            Y_t = 0;
            obsTable_tmp.Y_t(:) = Y_t;
            obsTable_tmp.Y = obsTable_tmp.Value+Y_t;
        case 'pr'
            for i_ID = 1:height(metaTable)
                flag_ID = ismember(obsTable_tmp.ID,metaTable.ID(i_ID));
                obsTable_tmp_single = obsTable_tmp(flag_ID,:);
                Y_t = 1.1+min(obsTable_tmp_single.Value);
                obsTable_tmp.Y_t(flag_ID) = Y_t;
                obsTable_tmp.Y(flag_ID) = log(obsTable_tmp_single.Value+Y_t);
            end
    end    
    
    mu_ESM = mean(tgt_ESM_cal(:,:,i_mth,:),[3 4],'omitnan');
    f_int = griddedInterpolant(tgt_longrid,tgt_latgrid,reshape(mu_ESM,tgt_size),'nearest','nearest');
    
    % calculate mu_Y
    for i_ID = 1:height(metaTable)
        flag_ID = ismember(obsTable_tmp.ID,metaTable.ID(i_ID));
        obsTable_tmp.mu_Y(flag_ID) = mean(obsTable_tmp.Y(flag_ID),'omitnan');
        % calculate mu_ESM with interpolation
        obsTable_tmp.mu_ESM(flag_ID) = reshape(f_int(metaTable.lon(i_ID),metaTable.lat(i_ID)),[],1);
    end

    % add the transformed variable and mu_Y to the main table
    obsTable.Y(flag_mth) = obsTable_tmp.Y; % Y is the transformed variable
    obsTable.Y_t(flag_mth) = obsTable_tmp.Y_t; 
    obsTable.mu_Y(flag_mth) = obsTable_tmp.mu_Y;
    obsTable.mu_ESM(flag_mth) = obsTable_tmp.mu_ESM;
end

clear obsTable_tmp
%% divide in cal/val tables
% array with the timeslice (months since 0 CE) corresponding to the PC
% tgt_ESM_cal is ordered by month, and then by year
time_pcr_cal = 12*year_start+...
    reshape(repmat(1:12,n_year_train,1)'+(12*(find(flag_cal)-1)),[],1);
time_pcr_val = 12*year_start+...
    reshape(repmat(1:12,n_year_val,1)'+(12*(find(flag_val)-1)),[],1);
time_pcr_test = 12*year_start+...
    reshape(repmat(1:12,n_year_test,1)'+(12*(find(flag_test)-1)),[],1);

% define obsTable for calibration
flag_tmp=ismember(obsTable.month_since_0CE,time_pcr_cal);
obsTable.flag_cal(flag_tmp) = true;
% update metaTable
flag_ID = ismember(metaTable.ID,obsTable.ID(obsTable.flag_cal));
metaTable.flag_cal(flag_ID) = true;

% define obsTable for validation
flag_time_val=ismember(obsTable.month_since_0CE,time_pcr_val);
obsTable.flag_val(flag_time_val) = true;
% update metaTable
flag_ID = ismember(metaTable.ID,obsTable.ID(obsTable.flag_val));
metaTable.flag_val(flag_ID) = true;

% define obsTable for test
flag_time_test=ismember(obsTable.month_since_0CE,time_pcr_test);
obsTable.flag_test(flag_time_test) = true;
% update metaTable
flag_ID = ismember(metaTable.ID,obsTable.ID(obsTable.flag_test));
metaTable.flag_test(flag_ID)=true;

% screen calibration dataset for stations with less than either 10 years
% (after 1925), or than 20 years (before 1925) or 30 years (after 1925) 
% of observations in the calibration years.
% they end up in the test datasets
n_removed = 0;
for i_ID = height(metaTable):-1:1
    flag_ID = ismember(obsTable.ID,metaTable.ID(i_ID));
    obsTable_tmp = obsTable(flag_ID,:);
    length_caltime_year = floor(sum(obsTable_tmp.flag_cal)/12);
    start_timeseries_year = min(obsTable_tmp.month_since_0CE)/12;

    inflim = 20; % inferior limit (to remove station from dataset)
    % filtering of the data (if this is true, the station goes in
    % the testing dataset)
    bool_filter = (length_caltime_year<30);

    if bool_filter

        % if the series is too short just remove it
        if length_caltime_year<inflim
            n_removed = n_removed +1;
            % remove from dataset
            obsTable(flag_ID,:)=[];
            metaTable(i_ID,:)=[];
        else
            % add to testing
            obsTable.flag_test(flag_ID)=true;
            metaTable.flag_test(i_ID)=true;
            % remove from validation
            obsTable.flag_val(flag_ID)=false;
            metaTable.flag_val(i_ID)=false;
            % remove from the calibration dataset
            obsTable.flag_cal(flag_ID)=false;
            metaTable.flag_cal(i_ID)=false;
        end

    end
end

disp([int2str(n_removed) ' timeseries shorter than ' int2str(inflim) ' years. they were removed from the dataset']);
% visualize dataset splitting
visualize_splitting_v1(obsTable, metaTable, year_start, n_year_data, path_fig, name_var, suffix, res_fig)


%% PC regression

for i_mth_interval = 1:size(mth_interval,1)

    i_mth = mth_interval(i_mth_interval,:);
    i_mth_txt = [int2str(mth_interval(i_mth_interval,1)) '-' int2str(mth_interval(i_mth_interval,end))];
    disp(['calibrating model for month(s) ' i_mth_txt])
    tgt_ESM_mth = reshape(squeeze(tgt_ESM(:,:,i_mth,:)),[],size(tgt_ESM,4)*length(i_mth));
    tgt_ESM_cal_mth = reshape(squeeze(tgt_ESM_cal(:,:,i_mth,:)),[],n_year_train*length(i_mth));

    % PCA
    disp('PCA...')
    [eof_all, pc_all, pc_var_all, tsquared_all, frac_var_pc_all, mu_ESM] = ...
        pca(tgt_ESM_cal_mth');
    mu_ESM = mu_ESM';
    n_pc = size(eof_all,2);

    %filter dataset on considered months
    flag_mth = any((mod(obsTable.month_since_0CE-1,12)+1)==i_mth,2);
    obsTable_mth = obsTable(flag_mth,:);
    time_mth = 12*year_start+...
        reshape(repmat(i_mth,n_year_data,1)'+12*((1:n_year_data)-1),[],1);    
    time_mthcal = intersect(time_mth,time_pcr_cal);

    %% simple linear regression for comparison
    disp('calibrate LR...')
    [obsTable_mth,rmse_array_lr,r2_array_lr,lm_list_lr] =  ...
        recal_predval_lr_v2(obsTable_mth,metaTable,time_mth,name_var,n_err_iter_lr);

    rmse_lr = mean(rmse_array_lr);
    disp(['lr RMSE ' num2str(rmse_lr)])
    r2_lr = mean(r2_array_lr);
    disp(['lr R2adj ' num2str(r2_lr)])
    
    %% calibrate model with only the prescripted EOFs and downscale the whole dataset

    disp('initial PCR...')
    % define the initial PCs that are used as regressors (as indicated in
    % the parameters at the beginning)
    % initialize the array that lists the used principal components
    n_pc_min = find(frac_var_pc_all>10,1,'last');
    max_pc = 16-n_pc_min;
    n_pc_array = 1:n_pc_min;

    % get regression overthe whole dataset (using the minimum number of
    % PCs, call it 'the best' selected PCs for the moment: this will change
    % during the model selection process)
    [obsTable_tmpbest,rmse_array,r2_array,lm_list_tmpbest] =  ...
        recal_predval_v2(n_pc_array,pc_all,eof_all,...
            tgt_ESM_mth, mu_ESM, ...
            obsTable_mth,metaTable,time_mthcal,time_mth,name_var,true,false,n_err_iter_selection);
    
    %calculate the RMSE and store it for starting model selection
    rmse_list = [];
    R2_list = [];
    RMSE_tmpbest = mean(rmse_array);
    disp(['initial RMSE ' num2str(RMSE_tmpbest)])
    R2adj_tmpbest =  mean(r2_array);
    disp(['initial R2adj ' num2str(R2adj_tmpbest)])
    
    %% model selection

    disp('PCR model selection...')
    rmse_list(length(rmse_list)+1)=RMSE_tmpbest;
    R2_list(length(R2_list)+1)=R2adj_tmpbest;
    
    % search for PCs that improve RMSE if included among regressors
    % *tmpbest are the ones that are updated
    % *best are the variables that are saved at the end of each cycle
    n_pc_array_best = n_pc_array;
    n_pc_array_tmpbest = n_pc_array;
    no_new_PC = false;
    for i_iter = 1:max_pc
        if not(no_new_PC)
        disp(['Iteration ' int2str(i_iter) '/' int2str(max_pc)])
        added_PC = false;
        % find the best pc among the ones that were not selected yet
        for i_pc = setdiff(1:(n_pc_min+20),n_pc_array_best)
            % regress again using adding one new PC
            n_pc_array_new = [n_pc_array_best i_pc];
            frac_var_ESM_new = frac_var_pc_all(n_pc_array_new);
            if frac_var_ESM_new(end)>frac_threshold
            disp(['PC #' int2str(i_pc)])
            % regress on the test period
            [obsTable_tmp,rmse_array,r2_array,lm_list_tmpbest] =  ...
                recal_predval_v2(n_pc_array_new,pc_all,eof_all,...
                    tgt_ESM_mth, mu_ESM, ...
                    obsTable_mth,metaTable,time_mthcal,time_mth,name_var,false,false,n_err_iter_selection);

            RMSE_new = mean(rmse_array);
            R2adj_new =  mean(r2_array);
            
            % when adding the current PC improves RMSE, add the current PC to
            % the list and update the RMSE
            switch name_var 
                case 'tas'
                    add_PC_flag = RMSE_new<0.999*RMSE_tmpbest;% ...
                        %&& R2adj_new>R2adj_tmpbest;
                case 'pr'
                    add_PC_flag = RMSE_new<0.999*RMSE_tmpbest;
            end
            if add_PC_flag
                RMSE_tmpbest = RMSE_new;
                R2adj_tmpbest = R2adj_new;
                n_pc_array_tmpbest = n_pc_array_new;
                added_PC = true;
                disp('-> improvement')
            end
            end
        end
        % save the list of PCs that provide improvements
        n_pc_array_best = n_pc_array_tmpbest;
        if added_PC
            rmse_list(length(rmse_list)+1)=RMSE_tmpbest;
            R2_list(length(R2_list)+1)=R2adj_tmpbest;
        else
            no_new_PC = true;
        end
        end
    end

    %% final PCR
    
    disp('final PCR...')
    n_pc_array = n_pc_array_best;
    [obsTable_mth,rmse_array,r2_array,lm_list] =  ...
        recal_predval_v2(n_pc_array,pc_all,eof_all,...
            tgt_ESM_mth, mu_ESM, ...
            obsTable_mth,metaTable,time_mthcal,time_mth,name_var,true,true,n_err_iter_final);

    RMSE_tmpbest = mean(rmse_array);
    disp(['final RMSE ' num2str(RMSE_tmpbest)])
    R2adj_tmpbest =  mean(r2_array);
    disp(['final R2adj ' num2str(R2adj_tmpbest)])

    % calculate averages in every station
    mu_Y_local = varfun(@mean, obsTable_mth, 'InputVariables', 'mu_Y', 'GroupingVariables', 'ID');
    Y_t_local = varfun(@mean, obsTable_mth, 'InputVariables', 'Y_t', 'GroupingVariables', 'ID');
    mu_Y_local = mu_Y_local(:,[1 3]);
    Y_t_local = Y_t_local(:,[1 3]);
    mu_Y_local.Properties.VariableNames(2) = "mu_Y";
    Y_t_local.Properties.VariableNames(2) = "Y_t";
      
    %% EOFs

    plot_eof_maps_four(eof_all, n_pc_array, tgt_size, flag_land, tgt_lat, tgt_lon, eof_limit, ...
        lim_lat, lim_lon, eof_step, palette_eof(256), path_shp_file, res_fig, name_var, i_mth_txt, unit_var, path_fig, suffix)

    plot_eof_maps_full(eof_all, n_pc_array, tgt_size, flag_land, tgt_lat, tgt_lon, eof_limit, ...
        lim_lat, lim_lon, eof_step, palette_eof(256), path_shp_file, res_fig, name_var, i_mth_txt, unit_var, path_fig, suffix)
    
    close all

    %% save models
    disp('save...')
    lat = tgt_lat;
    lon = tgt_lon;
    var_save = {'eof_all','frac_threshold','i_mth','flag_cal','flag_test','flag_val'...
        'lm_list','lm_list_lr','lat','lon','mu_Y_local','Y_t_local',...
        'n_pc_array','pc_all','r2_array','r2_array_lr','rmse_array','rmse_array_lr','rmse_list'};
    var_desc = {
        'Array of EOFs (empirical orthogonal functions for projection).', ...
        'Fraction of explained variance threshold used for selecting EOFs in the model.', ...
        'The months being processed in the this downscaling iteration.', ...
        'flag for the calibration years',...
        'flag for the test years',...
        'flag for the model selection/validation years',...
        'List of PCR models for each station, trained on the ESM data.',...
        'List of standard linear regression models for each station, trained on the ESM data',...
        'Latitude coordinates of the target grid for downscaling.',...
        'Longitude coordinates of the target grid for downscaling.',...
        'Local mean values for observed variables (in each station). Used for inverse transformation in the downscaling process.',...
        'Translation of the data. Used for inverse transformation in the downscaling process.',...
        'Array of principal component indices selected for the PCR.',...
        'The set of all principal components used in the EOF decomposition.',...
        'Array of R-squared values for the PCR',...
        'Array of R-squared values for the standard linear regression .',...
        'Array of Root Mean Squared Errors (RMSE) for the PCR.',...
        'Array of RMSE for the standard linear regression model.',...
        'List of RMSE values showing rmse improvements at each selecte PC'...
    };    
    main_description = [];
    reg_formula = "predictand = SX' * EOF, with SX the time anomalies, and EOF the empirical orthogonal functions";
    info_var = struct();
    for i_var = 1:length(var_save)
        info_var(i_var).variable_name = var_save{i_var};
        info_var(i_var).description = var_desc{i_var};
    end
    
    save_name = fullfile(path_main, ['downscaling_models_' name_model], ...
        [name_var '_PCR_models_mth=' i_mth_txt suffix]);
    
    % Save the variables
    save(save_name, var_save{:}, 'info_var', 'main_description', 'reg_formula');

    %% save original dataset with splitting
    var_save = {'lat','lon','n_pc_array','i_mth', 'obsTable_mth','metaTable'};
    var_desc = {
        'Latitude coordinates of the target grid for downscaling.', ...
        'Longitude coordinates of the target grid for downscaling.',...
        'Array of selected principal components used for EOF decomposition.',...
        'Array of EOFs (empirical orthogonal functions for projection).', ...
        'The months being processed in the this downscaling iteration.', ...
        'Table containing the observed data, including the stations and their respective variables.',...
        'Metadata table'
    };
    info_var = struct();
    for i_var = 1:length(var_save)
        info_var(i_var).variable_name = var_save{i_var};
        info_var(i_var).description = var_desc{i_var};
    end

    save_name = fullfile(path_main, ['downscaling_models_' name_model], ...
        [name_var '_PCR_original_data_mth=' i_mth_txt suffix]);
    
    % Save the variables
    save(save_name, var_save{:}, 'info_var', 'main_description', 'reg_formula');

end


