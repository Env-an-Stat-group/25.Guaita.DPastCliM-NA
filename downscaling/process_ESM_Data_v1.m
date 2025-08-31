function [tgt_ESM, time_ESM] = process_ESM_Data_v1(file_list_ESM, path_ESM, name_model, name_var, year_start, year_end, year_base_ESM, filter_lat, filter_lon, tgt_lat, tgt_lon)
    % processESMData - Processes ESM data, interpolates it on the target grid and calculates values for stations.
    %
    % Inputs:
    %   file_list_ESM - List of ESM files to load
    %   path_ESM - Path where ESM files are stored
    %   name_var - Variable name to process ('pr' or 'tas')
    %   year_start - Start year for the time filter
    %   year_end - End year for the time filter
    %   filter_lat - Latitude filter range [min_lat, max_lat]
    %   filter_lon - Longitude filter range [min_lon, max_lon]
    %   tgt_lat - Latitude values for target grid
    %   tgt_lon - Longitude values for target grid
    %
    % Outputs:
    %   tgt_ESM - Processed and interpolated ESM data on the target grid
    %   time_ESM - Processed time values in months since 0 CE

    % Load total precipitation and convert from [m] to [mm], concatenate every variable file along the third dimension
    for i_file = 1:length(file_list_ESM)
        path_load_ESM = fullfile(path_ESM, file_list_ESM(i_file).name);
        if i_file == 1
            ESM_lat = double(ncread(path_load_ESM, 'lat'));
            ESM_lon = double(ncread(path_load_ESM, 'lon'));
            switch name_var
                case 'pr'
                    conv_fac = 86400;  % Conversion factor for precipitation
                    ESM_map = ncread(path_load_ESM, name_var) * conv_fac;
                case 'tas'
                    conv_fac = 273.15;  % Conversion factor for temperature
                    ESM_map = ncread(path_load_ESM, name_var) - conv_fac;
            end
            time_ESM = ncread(path_load_ESM, 'time');
        else
            switch name_var
                case 'pr'
                    tmp_var = ncread(path_load_ESM, name_var) * conv_fac;
                case 'tas'
                    tmp_var = ncread(path_load_ESM, name_var) - conv_fac;
            end
            tmp_time = ncread(path_load_ESM, 'time');
            ESM_map = cat(3, ESM_map, tmp_var);
            time_ESM = cat(1, time_ESM, tmp_time);
        end
    end

    % Filter on time the whole loaded CMIP series
    [start_time, end_time] = timegreg_2_timeCMIP(datetime(year_start, 1, 1), datetime(year_end + 1, 1, 1), datetime(1850, 1, 1), name_model);
    flag_time = time_ESM >= start_time & time_ESM <= end_time;
    ESM_map = ESM_map(:, :, flag_time);
    
    % Change time_ESM (days since 1850) to months since 0 CE
    time_ESM = year_base_ESM * 12 + (1:length(time_ESM));
    
    % Filter the time
    time_ESM = time_ESM(flag_time);

    % Fix longitude-latitude of 3D matrix
    [ESM_map, c_lon, c_lat] = cmipcoord_2_map(ESM_map, ESM_lon, ESM_lat);

    % Filter the coarse lat and lon, such that they contain the fine grid
    ind_c_lat_min = find(c_lat >= filter_lat(1), 1, 'first') - 1;
    ind_c_lat_max = find(c_lat <= filter_lat(2), 1, 'last') + 1;
    ind_c_lon_min = find(c_lon >= filter_lon(1), 1, 'first') - 1;
    ind_c_lon_max = find(c_lon <= filter_lon(2), 1, 'last') + 1;
    c_lat = c_lat(ind_c_lat_min:ind_c_lat_max);
    c_lon = c_lon(ind_c_lon_min:ind_c_lon_max);

    % Filter the domain
    cESM_map = ESM_map(ind_c_lon_min:ind_c_lon_max, ind_c_lat_min:ind_c_lat_max, :);

    % Define lon and lat grid
    [c_longrid, c_latgrid] = ndgrid(c_lon, c_lat);

    % Define sizes
    tgt_size = [length(tgt_lon), length(tgt_lat)];

    tgt_ESM = nan([tgt_size, size(ESM_map, 3)]);

    % define lon and lat grid
    [tgt_longrid,tgt_latgrid] = ndgrid(tgt_lon,tgt_lat);
    
    % Interpolate on station coordinates for each time step
    for i_time = 1:length(time_ESM)
        % Interpolate coarse_map on the target grid and on the station coordinates
        f_int = griddedInterpolant(c_longrid, c_latgrid, cESM_map(:,:,i_time), 'makima', 'makima');
        tgt_ESM_tmp = f_int(tgt_longrid, tgt_latgrid);
        [tgt_ESM(:,:,i_time), tgt_lon, tgt_lat] = cmipcoord_2_map(tgt_ESM_tmp, tgt_lon, tgt_lat);
    end

    % Transpose time array
    time_ESM = time_ESM';

    switch name_var
        case 'pr'
            tgt_ESM = max(0,tgt_ESM);
    end

end
