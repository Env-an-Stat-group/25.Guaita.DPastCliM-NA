clear  % Clear the workspace, removing all variables

% Define the main path where the data and code are located
path_main = '/data/pguaita/downscaling';%
addpath(genpath(path_main));  % Add all subdirectories of the main path to the MATLAB search path

% Define the path for the GHCN dataset and specify the variable name
path_GHCM = fullfile(path_main,'obs_data/GHCNm');  % Path to the GHCN data folder
GHCNm_name_var = 'pr';  % Specify the variable name (e.g., 'tas' for temperature or 'pr' for precipitation)
suffix = '_NA_020';  % Suffix to be used in file names

% Define the path for the figures and shapefile (world borders)
path_fig = fullfile(path_main,'obs_data','figures');  % Path for saving figures
path_shp_file = fullfile(path_main,'matlab_code/visualization/world_borders/ne_10m_admin_0_countries.shp');  % Path to shapefile for world borders
res_fig = '-r500';  % Resolution for figure output
res_plot = '-r300';  % Resolution for plot output

% gpr parameters
        
switch GHCNm_name_var
    case 'pr'
        sigmanoise = 0.1;
        reg_fact = 0.25;
        kernelFunction = "exponential";  % Exponential kernel for precipitation
    case 'tas'
        sigmanoise = 0.1;
        reg_fact = 0.1;
        kernelFunction = "exponential";  % Exponential kernel for temperature
end

opt_gpr = {'KernelFunction', kernelFunction, ...
           'FitMethod','fic',...
           'PredictMethod','fic',...
           'Standardize', true};


% Switch to set file paths based on the selected variable name (tas or pr)
switch GHCNm_name_var
    case 'tas'  % If the variable is temperature
        name_file_dat = 'ghcnm.tavg.v4.0.1.20240916.qcf.dat';  % Data file for temperature
        name_file_inv = 'ghcnm.tavg.v4.0.1.20240916.qcf.inv';  % Inverse file for temperature
        % Define full paths to the temperature data files
        path_file_dat = fullfile(path_GHCM, name_file_dat);  
        path_file_inv = fullfile(path_GHCM, name_file_inv);  
    case 'pr'  % If the variable is precipitation
        name_folder_file = 'ghcn-m_pr';  % Folder for precipitation data
end

% Check if the figures directory exists, and if not, create it
if not(exist(path_fig, 'dir'))
    mkdir(path_fig);  % Create the directory for figures
end

% Define the limits for the time period (historical data range)
% lim_lat = [7 75];  % Latitude limits for the region of interest (commented out, not used in this section)
% lim_lon = [-180 -50];  % Longitude limits for the region of interest (commented out, not used in this section)
lim_year = [1850 2014];  % Define the year range (1850 to 2014)

% Load the grid for downscaling (this is a predefined grid for data processing)
load(fullfile(path_main, 'static_maps', ['downscaling_grid' suffix '.mat']));  % Load downscaling grid based on suffix

% Define latitude and longitude limits for the domain based on the grid
lim_lat = [min(lat), max(lat)];  % Set latitude limits from the grid
lim_lon = [min(lon), max(lon)];  % Set longitude limits from the grid

% Switch to set the units for the GHCN data based on the selected variable
switch GHCNm_name_var
    case 'tas'  % If the variable is temperature
        GHCNm_unit_var = '°C';  % Unit is degrees Celsius
    case 'pr'  % If the variable is precipitation
        GHCNm_unit_var = 'mm/day';  % Unit is millimeters per day
end

%% load files

% Switch based on the type of variable (temperature 'tas' or precipitation 'pr')
switch GHCNm_name_var
    case 'tas'
        % Read the inventory (meta) file for temperature data
        
        fid = fopen(path_file_inv, 'r');  % Open the inventory file for reading
        
        % Define the format of the inventory file (columns for ID, lat, lon, elev, location)
        formatSpec = '%11s %f %f %f %s';
        
        % Read the inventory data
        inv_file = textscan(fid, formatSpec, 'Delimiter', '', 'Whitespace', '', 'MultipleDelimsAsOne', true);
        
        fclose(fid);  % Close the file after reading
        
        % Convert the data into a table for easier manipulation
        metaTable = table(inv_file{1}, inv_file{2}, inv_file{3}, inv_file{4}, inv_file{5}, ...
                          'VariableNames', {'ID', 'lat', 'lon', 'elev', 'location'});
        
        clear inv_file  % Clear the inventory data from memory

        % Filter the data based on latitude and longitude limits (domain)
        flag_domain = lim_lat(1) <= metaTable.lat & metaTable.lat <= lim_lat(2) & ...
                      lim_lon(1) <= metaTable.lon & metaTable.lon <= lim_lon(2);
        metaTable = metaTable(flag_domain, :);  % Keep only rows within the specified region
        
        % Read the data file containing temperature data
        fid = fopen(path_file_dat, 'r');  % Open the data file for reading
        
        % Define the format of the input data file (columns for ID, Year, and temperature data)
        formatSpec = '%11s %4d %4s';
        
        % Add formats for the 12 months of data
        for i = 1:12
            formatSpec = [formatSpec ' %5f %*3s'];  % Each month has a temperature value and a flag
        end
        
        % Read the data from the file
        obs_file = textscan(fid, formatSpec, 'Delimiter', '', 'Whitespace', '', 'MultipleDelimsAsOne', true);
        
        fclose(fid);  % Close the data file after reading
        
        % Create a table to hold the temperature data for each station and month
        obsTable_tmp = table(obs_file{1}, obs_file{2}, obs_file{4}/100, obs_file{5}/100, obs_file{6}/100, obs_file{7}/100, ...
                             obs_file{8}/100, obs_file{9}/100, obs_file{10}/100, obs_file{11}/100, obs_file{12}/100, ...
                             obs_file{13}/100, obs_file{14}/100, obs_file{15}/100, ...
                             'VariableNames', {'ID', 'Year', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'});

        % Make the observed data table correspond to the meta table by matching station IDs
        flag_station = ismember(string(obsTable_tmp.ID), string(metaTable.ID));  % Find matching IDs
        obsTable_tmp = obsTable_tmp(flag_station, :);  % Keep only the matching rows
        flag_station = ismember(string(metaTable.ID), string(obsTable_tmp.ID));  % Find matching IDs in the meta table
        metaTable = metaTable(flag_station, :);  % Keep only the matching rows

        clear obs_file  % Clear the observation data from memory

    case 'pr'
        % Initialize variables for precipitation data (used for performance tracking)
        count = 0;
        t_read = 0;
        t_tabl = 0;
        t_mont = 0;
        t_atta = 0;
        t_clea = 0;
        tot_height_sum = 0;
        
        % Read the inventory (meta) file for precipitation data
        fid = fopen(fullfile(path_main, 'obs_data', 'GHCNm', 'ghcn-m_v4_prcp_inventory.txt'), 'r');
        
        % Define the format of the inventory file (columns for ID, lat, lon, elev, year_start, etc.)
        formatSpec = '%11s %f %f %f %40s %d %d %d';
        
        % Read the inventory data
        inv_file = textscan(fid, formatSpec, 'Delimiter', '', 'Whitespace', '', 'MultipleDelimsAsOne', true);
        
        fclose(fid);  % Close the file after reading
        
        % Convert the data into a table for easier manipulation
        metaTable = table(inv_file{1}, inv_file{2}, inv_file{3}, inv_file{4}, inv_file{7}, ...
                          'VariableNames', {'ID', 'lat', 'lon', 'elev', 'year_start'});
        
        % Filter the data based on the time range (lim_year)
        flag_time = metaTable.year_start <= lim_year(2);  % Keep only data within the specified year range
        metaTable = metaTable(flag_time, 1:4);  % Keep only relevant columns (ID, lat, lon, elev)
        
        clear inv_file  % Clear the inventory data from memory

        % Filter the data based on latitude and longitude limits (domain)
        flag_domain = lim_lat(1) <= metaTable.lat & metaTable.lat <= lim_lat(2) & ...
                      lim_lon(1) <= metaTable.lon & metaTable.lon <= lim_lon(2);
        metaTable = metaTable(flag_domain, :);  % Keep only rows within the specified region

        % Get the list of data files (CSV files) for precipitation
        file_list = dir(fullfile(path_main, 'obs_data', 'GHCNm', name_folder_file, '*.csv'));

        % Remove lat and lon duplicates
        % Identify unique locations (lat, lon) and the first occurrence ID
        [unique_locs, first_idx] = unique(metaTable(:, {'lat', 'lon'}), 'rows', 'stable');
        metaTable = metaTable(first_idx,:);

        % Filter the file list to match the stations in the metaTable
        file_list_tmp = char({file_list.name});
        file_list_tmp = string(file_list_tmp(:, 1:11));  % Extract station IDs from file names
        flag_file = ismember(file_list_tmp, string(metaTable.ID));  % Find matching station IDs
        file_list = file_list(flag_file);  % Keep only the matching files
        flag_file = ismember(string(metaTable.ID), file_list_tmp);  % Find matching file names
        metaTable = metaTable(flag_file, :);  % Keep only the matching rows in the metaTable
        
        % Define the table to store precipitation data (3 columns: ID, month, value)
        obsTable = table('Size', [300 * height(metaTable), 3], ...
                        'VariableTypes', {'string', 'int64', 'single'}, ...
                        'VariableNames', {'ID', 'month_since_0CE', 'Value'});
        input_row_table = 1;  % Keep track of the position in the table

        % Define the format for reading the precipitation data from each file
        formatSpec = '%11s %q %f %f %f %4d %2d %f %c %c %c %s';
        
        % Loop through each file and read the precipitation data
        for i_file = 1:length(file_list)
            % Open the file for reading
            fid = fopen(fullfile(file_list(i_file).folder, file_list(i_file).name), 'r');
            
            % Read the data from the file
            obs_file = textscan(fid, formatSpec, 'Delimiter', ',', 'Whitespace', '');
            
            fclose(fid);  % Close the file after reading
            
            % Convert the precipitation data to mm/day
            if not(isempty(obs_file{3}))
                % Create a temporary table for the data
                ID_array = obs_file{1};
                year_array = obs_file{6};
                month_array = obs_file{7};
                value_array = obs_file{8};
                n_min = min([length(ID_array),length(year_array),length(month_array),length(value_array)]);
                ID_array = ID_array(1:n_min);
                year_array = year_array(1:n_min);
                month_array = month_array(1:n_min);
                value_array = value_array(1:n_min);
                tmp_table = table(ID_array, year_array * 12 + month_array, value_array / 10, ...
                                  'VariableNames', {'ID', 'month_since_0CE', 'Value'});
                
                % Adjust the precipitation value based on the number of days in the month
                tmp_table.Value = tmp_table.Value ./ double(eomday(obs_file{6}, obs_file{7}));
                
                % Append the data to the main observation table
                obsTable(input_row_table:(input_row_table + height(tmp_table) - 1), :) = tmp_table;
                input_row_table = input_row_table + height(tmp_table);  % Update the table position
            else
                count = count + 1;  % Count missing or empty data files
            end

            clear obs_file  % Clear the observation data from memory
        end
end

%% Filter by time

switch GHCNm_name_var
    case 'tas'
        % Filter temperature data (obsTable_tmp) based on the year range (lim_year)
        flag_year = lim_year(1) <= obsTable_tmp.Year & obsTable_tmp.Year <= lim_year(2);
        obsTable_tmp = obsTable_tmp(flag_year, :);  % Keep only the rows within the year range
        
        % Filter the metadata table (metaTable) based on the filtered temperature data
        list_ID = categorical(unique(obsTable_tmp.ID));  % Get the unique station IDs from the filtered data
        % Find rows in the metaTable that match the station IDs
        flag_ID = ismember(metaTable.ID, list_ID);
        metaTable = metaTable(flag_ID, :);  % Keep only the rows in metaTable with matching station IDs
        
    case 'pr'
        % For precipitation data, first remove rows with missing IDs
        obsTable = obsTable(~ismissing(obsTable.ID), :);  % Remove rows with missing IDs
        
        % Filter precipitation data based on the year range (lim_year)
        flag_year = (lim_year(1) * 12 + 1) <= obsTable.month_since_0CE & ...
                    obsTable.month_since_0CE <= ((lim_year(2) + 1) * 12);  % Convert years to months and filter
        obsTable = obsTable(flag_year, :);  % Keep only the rows within the specified month range
end

%% Reorganize table and fix time so that we have month number starting from 0 CE

switch GHCNm_name_var
    case 'tas'
        % Create an empty table with 3 columns: 'ID', 'month_since_0CE', and 'Value'
        obsTable = table('Size', [0 3], 'VariableTypes', {'string', 'int64', 'single'}, ...
                         'VariableNames', {'ID', 'month_since_0CE', 'Value'});
        
        % Loop through each month (1 to 12) to reorganize data by months
        for i_mth = 1:12
            % Extract columns for ID, Year, and the month data from the obsTable_tmp
            tmp_table = obsTable_tmp(:, [1 2 2 + i_mth]);  % Select the relevant columns (ID, Year, and month data)
            
            % Convert 'Year' to 'month_since_0CE' (months since 0 CE)
            tmp_table.Year = tmp_table.Year * 12 + i_mth;  % Convert year to the corresponding month number
            tmp_table.Properties.VariableNames{'Year'} = 'month_since_0CE';  % Rename the 'Year' column
            tmp_table.Properties.VariableNames{3} = 'Value';  % Rename the third column to 'Value'
            
            % If it's the first month, initialize the obsTable with this data
            if i_mth == 1
                obsTable = tmp_table;
            else
                % Otherwise, append the current month's data to the obsTable
                obsTable = vertcat(obsTable, tmp_table);
            end
        end
        
        % Remove erroneous values (-99.99) from the table
        flag_error = obsTable.Value == -99.99;  % Identify rows with erroneous values
        obsTable = obsTable(~flag_error, :);  % Remove rows with the error value
        
    case 'pr'
        % For precipitation, replace erroneous values (less than 0) with NaN
        flag_error = obsTable.Value < 0;  % Identify rows with erroneous values (negative values)
        obsTable.Value(flag_error) = NaN;  % Replace these values with NaN
end

%% Get final table with latitude, longitude, and elevation

% Match the IDs in obsTable with those in metaTable and get the corresponding metadata
[~, IDloc_inv] = ismember(categorical(obsTable.ID), categorical(metaTable.ID));  % Find matching IDs
obsTable.lat = metaTable.lat(IDloc_inv);  % Add latitude to the obsTable
obsTable.lon = metaTable.lon(IDloc_inv);  % Add longitude to the obsTable
obsTable.elev = metaTable.elev(IDloc_inv);  % Add elevation to the obsTable

%% filter out timeseries shorter than 20 years
for i_ID = height(metaTable):-1:1
    flag_ID = ismember(obsTable.ID,metaTable.ID{i_ID});
    if (sum(flag_ID)/12)<20
        % remove from stations and metatable
        obsTable(flag_ID,:)=[];
        metaTable(i_ID,:)=[];
    end
end

%% Save the final table

% Create a structure to store metadata about the dataset
note = struct();
note.time = 'time is expressed in number of months starting from January of 0 BC';  % Time information
note.data_origin = 'origin = this gridded dataset is from the Global Historical Climate Network at monthly resolution (GHCNm)';  % Data origin description

% Define the path to save the final dataset
path_save = fullfile(path_main, 'obs_data', ['GHCNm_' GHCNm_name_var suffix '.mat']);
note.variable.name = GHCNm_name_var;  % Store the variable name (temperature or precipitation)
note.variable.unit = GHCNm_unit_var;  % Store the unit of the variable (°C or mm/day)

% Save the tables (obsTable and metaTable) and the note structure to a .mat file
save(path_save, 'obsTable', 'metaTable', 'note');

%% Grid observations

% Define a time array for the target matrix, ensuring unique months are sorted.
obs_time = sort(unique(obsTable.month_since_0CE));  
obs_map = nan(length(lon), length(lat), length(obs_time));  % Create an empty map for gridded data

% Filter the observation dataset to include only coordinates within the specified domain (latitude and longitude limits).
flag_coordinates = lim_lon(1) <= obsTable.lon & obsTable.lon <= lim_lon(2) & lim_lat(1) <= obsTable.lat & obsTable.lat <= lim_lat(2);
obsTable = obsTable(flag_coordinates, :);  % Keep rows within the domain

% Get query grids for target longitude and latitude using ndgrid to create grid arrays.
[tgt_longrid, tgt_latgrid] = ndgrid(lon, lat);

% Load land-sea mask data to determine the land and sea areas.
landseamask = ncread(fullfile(path_main, 'static_maps', 'IMERG_land_sea_mask.nc'), 'landseamask');
flag_land = not(landseamask == 100);  % Create a flag for land cells (where land is 100, and sea is 0)
lat_land = ncread(fullfile(path_main, 'static_maps', 'IMERG_land_sea_mask.nc'), 'lat');
lon_land = ncread(fullfile(path_main, 'static_maps', 'IMERG_land_sea_mask.nc'), 'lon');
[flag_land, lon_land, lat_land] = cmipcoord_2_map(flag_land, lon_land, lat_land);  % Adjust the land-sea mask coordinates
% Regrid the land-sea mask onto the target grid (lat, lon).
[lonland_grid, latland_grid] = ndgrid(lon_land, lat_land);
f_int_land = griddedInterpolant(lonland_grid, latland_grid, double(flag_land), 'nearest');  % Interpolation function for land mask
flag_land = logical(f_int_land(tgt_longrid, tgt_latgrid));  % Apply the land-sea mask on the target grid

% Regrid for each observation time
for i_time = 1:length(obs_time)
    i_month0BC = obs_time(i_time);  % Get the current month since 0 BCE
    flag_time = obsTable.month_since_0CE == i_month0BC;  % Filter for observations in the current month
    obs_tmp = obsTable(flag_time, :);  % Get the filtered observations for the current month
    
    % Remove any NaN values from the observations
    obs_tmp(isnan(obs_tmp.Value), :) = [];
    
    % If there are enough observations (more than 2), perform GPR.
    if height(obs_tmp) > 2
        x = obs_tmp.lat;  % Latitude
        y = obs_tmp.lon;  % Longitude
        z = obs_tmp.Value;  % Observation values
        
        % Select the appropriate kernel function based on the variable type (temperature or precipitation)
        
        % Train a GPR model using the latitude and longitude as features.
        gprMdl = fitrgp([x, y], z, opt_gpr{:}, 'SigmaLowerBound', sigmanoise * std(z, 'omitnan'),'Regularization',reg_fact * std(z, 'omitnan'));
        
        % Predict the values on the target grid using the trained model and reshape the predictions.
        Y_hat_tmp = reshape(predict(gprMdl, [reshape(tgt_latgrid, [], 1), reshape(tgt_longrid, [], 1)]), size(tgt_longrid));
        obs_map(:, :, i_time) = Y_hat_tmp;  % Store the predicted values in the obs_map
    else
        % If there are not enough observations, use the mean of the available values for the entire grid.
        obs_map(:, :, i_time) = mean(obs_tmp.Value);
        % Apply the land-sea mask to set the sea areas to NaN.
        tmp_map = obs_map(:, :, i_time);
        tmp_map(~flag_land) = nan;
        obs_map(:, :, i_time) = tmp_map;  % Store the land-sea-masked map
    end
end

% Convert the final obs_map to single precision for storage efficiency.
obs_map = single(obs_map);

%% Save

% Create a metadata structure with information about the dataset.
note = struct();
note.time = 'time is expressed in number of months starting from January of 0 BC';  % Description of the time format
note.data_origin = 'origin = this gridded dataset is from the Global Historical Climate Network at monthly resolution (GHCNm)';  % Data origin

% Define the path to save the final gridded dataset.
path_save = fullfile(path_main, 'obs_data', ['gridobs_' GHCNm_name_var suffix '.mat']);
note.variable.name = GHCNm_name_var;  % Store the name of the variable (temperature or precipitation)
note.variable.unit = GHCNm_unit_var;  % Store the unit of the variable (°C or mm/day)

% Save the gridded observation data, latitude, longitude, time, and metadata into a .mat file.
save(path_save, 'obs_map', 'lat', 'lon', 'obs_time', 'note', '-v7.3');

%% plot

close all  % Close all existing figures before creating a new one

% Define figure colorbar limits and settings based on the variable (precipitation or temperature)
switch GHCNm_name_var
    case 'pr'
        % Precipitation settings
        val_limit = [0.75 6.25];  % Set the range of values for the colorbar (precipitation)
        val_step = 0.5;           % Step size for the colorbar ticks
        n_color_abs = 10;         % Number of colors in the color map
        dot_color = 'b';          % Color of the dots representing station locations (blue)
        palette_abs = @(m) Tol_seq_iridescent(m);  % Color palette function for precipitation
        var_title = 'pr';         % Title for the variable (precipitation)
        
    case 'tas'
        % Temperature settings
        val_limit = [-15 35];     % Set the range of values for the colorbar (temperature)
        val_step = 10;             % Step size for the colorbar ticks
        n_color_abs = 10;         % Number of colors in the color map
        dot_color = 'r';          % Color of the dots representing station locations (red)
        palette_abs = @(m) Tol_seq_YlOrBr(m);  % Color palette function for temperature
        var_title = 'tas';        % Title for the variable (temperature)
end

% Define the number of panels (subplots) for the figure
num_panels = 16;  % There will be 16 panels in total for the subplots

% Define time windows for the maps (in years)
time_bound = [1875 1910 1945 1980 2015];  % The start year for each period
time_bound_mth = time_bound * 12;         % Convert years to months (multiply by 12)

% Reshape the obs_map to include a 12-month dimension, ensuring it's ready for monthly plotting
obs_map = reshape(obs_map, size(obs_map, 1), size(obs_map, 2), 12, []);

%% Loop to generate subplots

close all  % Close any previously opened figures

% Create a new figure with a specified size for plotting
f = figure('Name', ['gridded observed ' GHCNm_name_var], 'Position', [50 50 1000 1200]);

for i_panel = 1:4
    % Loop through time periods defined in the 'time_bound' array (e.g., 1875-1910, 1910-1945, etc.)
    % Determine the station locations for each time period by filtering based on months
    flag_time = time_bound_mth(i_panel) <= obsTable.month_since_0CE & obsTable.month_since_0CE <= time_bound_mth(i_panel + 1);
    obs_unique = unique([obsTable.lat(flag_time), obsTable.lon(flag_time)], 'rows');  % Get unique lat/lon pairs for the stations

    % Create title text for the current time period
    title_text = [int2str(time_bound(i_panel)) '-' int2str(time_bound(i_panel + 1) - 1)];

    % Create subplot (num_panels rows, 1 column, current panel index)
    subplot(4, 4, i_panel);
        
    % Plot the station locations using the geosurfm_sub_dot function
    % The function should plot the station dots on a map with specified parameters
    geosurfm_sub_dot(obs_unique(:,1), obs_unique(:,2), 4, dot_color, ...
        title_text, lim_lat, lim_lon, ...
        path_shp_file, res_fig);
end

% Plot maps for the calibration period (years from 1945-1980, 1980-2015, etc.)
for i_panel = 5:8
    % Calculate the index for the subplot (modulo operation for wrapping)
    i_panel_mod = mod(i_panel - 1, 4) + 1;
    
    % Compute the average map for the time period by averaging over months
    plot_map = mean(obs_map(:,:,:,(time_bound(i_panel_mod) - lim_year(1) + 1):(time_bound(i_panel_mod + 1) - lim_year(1))), [3 4], 'omitnan');

    % Create title text for the calibration period (e.g., 1945-1980 year-round)
    title_text = [int2str(time_bound(i_panel_mod)) '-' int2str(time_bound(i_panel_mod + 1) - 1) ' (year-round)'];

    % Create subplot (num_panels rows, 1 column, current panel index)
    subplot(4, 4, i_panel);
        
    % Plot the map using the geosurfm_sub function (this function will plot a map with specific settings)
    geosurfm_sub(plot_map, flag_land, ...
        lat, lon, val_limit, ...
        [], [], title_text, ...
        lim_lat, lim_lon, ...
        val_step, palette_abs(n_color_abs), path_shp_file, res_fig);
end

% Plot seasonal maps (JJA - June, July, August)
for i_panel = 9:12
    % Similar to the previous loop, but for the JJA (summer) months
    i_panel_mod = mod(i_panel - 1, 4) + 1;
    plot_map = mean(obs_map(:,:,6:8,(time_bound(i_panel_mod) - lim_year(1) + 1):(time_bound(i_panel_mod + 1) - lim_year(1))), [3 4], 'omitnan');

    % Title text for JJA (summer) months
    title_text = [int2str(time_bound(i_panel_mod)) '-' int2str(time_bound(i_panel_mod + 1) - 1) ' (JJA)'];

    % Create subplot (num_panels rows, 1 column, current panel index)
    subplot(4, 4, i_panel);
        
    % Plot using the geosurfm_sub function (seasonal map)
    geosurfm_sub(plot_map, flag_land, ...
        lat, lon, val_limit, ...
        [], [], title_text, ...
        lim_lat, lim_lon, ...
        val_step, palette_abs(n_color_abs), path_shp_file, res_fig);
end

% Plot seasonal maps (DJF - December, January, February)
for i_panel = 13:16
    % Similar to the previous loop, but for DJF (winter) months
    i_panel_mod = mod(i_panel - 1, 4) + 1;
    plot_map = mean(obs_map(:,:,[12 1 2],(time_bound(i_panel_mod) - lim_year(1) + 1):(time_bound(i_panel_mod + 1) - lim_year(1))), [3 4], 'omitnan');

    % Title text for DJF (winter) months
    title_text = [int2str(time_bound(i_panel_mod)) '-' int2str(time_bound(i_panel_mod + 1) - 1) ' (DJF)'];

    % Create subplot (num_panels rows, 1 column, current panel index)
    subplot(4, 4, i_panel);
        
    % Plot using the geosurfm_sub function (seasonal map)
    geosurfm_sub(plot_map, flag_land, ...
        lat, lon, val_limit, ...
        [], [], title_text, ...
        lim_lat, lim_lon, ...
        val_step, palette_abs(n_color_abs), path_shp_file, res_fig);
end

% Set the title for the entire figure
title_text = [GHCNm_name_var ' station position and mean'];
sgtitle(title_text, 'fontweight', 'bold');  % Super title for the figure

% Add a colorbar for the entire figure
% The colorbar is positioned outside the subplots at the bottom
c = colorbar('Location', 'southoutside');  
c.Position = [0.1 0.065 0.8 0.023];  % Position and size of the colorbar
c.FontSize = 12;                % Set font size for ticks
c.FontWeight = 'bold';          % Set font weight for ticks
c.Ticks = (round(val_limit(1) / val_step) * val_step):val_step:(round(val_limit(2) / val_step) * val_step);  % Set ticks on colorbar
c.Label.String = GHCNm_unit_var;  % Set the label for the colorbar (unit of the variable)

% Save the final figure with all the subplots
save_name = [GHCNm_name_var ' station position and mean' suffix];
save_path = fullfile(path_fig, save_name);
print(f, [save_path '.tiff'], '-dtiff', res_plot);  % Save the figure as a TIFF file

