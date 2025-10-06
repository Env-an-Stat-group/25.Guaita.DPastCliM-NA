% retrieve downscaled data and output specific stations timeseries

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
n_min_yr = 75; % minimum number of years for stations
path_fig = fullfile(path_main,['downscaling_output_' name_model],'figures_PCR');
path_file = fullfile(path_main,['downscaling_output_' name_model]);
path_obs = fullfile(path_main,'obs_data');
path_downmodel = fullfile(path_main,['downscaling_models_' name_model]);
path_shp_file = fullfile(path_main,'/matlab_code_git/visualization/world_borders/ne_10m_admin_0_countries.shp'); 
suffix = '_NA_020';

% define stations
switch name_var
    case 'tas'
        list_station = {'USW00014740','USW00093193'};
    case 'pr'
        list_station = {'USC00190666','USC00049452'};
end

%% load raw data and metadata
disp('load data')

load(fullfile(path_file,[name_var '_' name_model '_' name_experiment '_raw_PCRdownscaled' suffix '.mat']))

%metadata
metatable = readtable(fullfile(path_file,[name_var '_metadata_table' suffix '.csv']));

%% get realizations
disp('get realizations')

switch name_var
    case 'tas'
        Ods_hat_mat = EOds_hat_mat + u_mat;
    case 'pr'
        Ods_hat_mat = pr_realizations(EOds_hat_mat,u_mat,...
            path_main,name_var,suffix,name_model);
end

%% extract data

disp('extracting data and csv output ')

list_var = {'month_since_0_CE','realization','mean_prediction','PI95_low','PI95_up'};
list_vartype = repmat({'double'},1,numel(list_var));
for i_stat = 1:length(list_station)
    flag_station = ismember(metatable.ID,list_station{i_stat});
    metatable_station = metatable(flag_station,:);
    lon = metatable_station.lon;
    lat = metatable_station.lat;
    elev = metatable_station.elev;
    ID = metatable_station.ID;
    T = table( ...
        (1:size(Ods_hat_mat,2))', ...
        Ods_hat_mat(flag_station,:)', ...
        EOds_hat_mat(flag_station,:)', ...
        PI_mat(flag_station,:,1)', ...
        PI_mat(flag_station,:,2)', ...
        'VariableNames', list_var);
    outfile = sprintf('station_%s_%s.csv', name_var, ID{1});  % filename per station

    % Open file and write metadata
    fid = fopen(outfile,'w');
    fprintf(fid, 'lon,%f\n', lon);
    fprintf(fid, 'lat,%f\n', lat);
    fprintf(fid, 'elev,%f\n', elev);
    fprintf(fid, 'ID,%s\n', ID{1});
    fclose(fid);
    
    % Now write table (with header row)
    writetable(T, outfile, 'WriteMode', 'append', 'WriteVariableNames', true);
end

disp('finished')
