function plot_eof_maps_four(eof_all, n_pc_array, tgt_size, flag_land, tgt_lat, tgt_lon, eof_limit, ...
    lim_lat, lim_lon, eof_step, palette_eof, path_shp_file, res_fig, name_var, i_mth_txt, unit_var, path_fig, suffix)
    % This function creates a set of EOF map plots (up to 4) with shared colorbar
    % Input:
    %   eof_all        - Matrix of EOFs
    %   n_pc_array     - Array of principal component indices to plot
    %   tgt_size       - Size for reshaping the EOF map
    %   flag_land      - Land flag data (for the map)
    %   tgt_lat, tgt_lon - Latitude and Longitude of the target map
    %   eof_limit      - Color limits for the EOF plots
    %   lim_lat, lim_lon - Latitude and longitude limits for the plot
    %   eof_step       - Step size for colorbar ticks
    %   palette_eof    - Color palette for EOF
    %   path_shp_file  - Path to shapefile
    %   res_fig        - Resolution for saving the figure
    %   name_var       - Variable name for the title
    %   i_mth_txt      - Month or time period identifier
    %   unit_var       - Unit for the variable (for colorbar label)
    %   path_fig       - Path to save the figure
    %   suffix         - Suffix to append to the saved figure name

    % Select the EOFs from the provided array
    eof_select = eof_all(:, n_pc_array);

    % Create a new figure with specified size
    f = figure('Position', [50, 0, 1300, 300]);

    % Initialize a variable to collect all subplot axes
    axes_handles = gobjects(4, 1);

    % Loop through panels (maximum 4 EOFs to plot)
    for i_panel = 1:min(length(n_pc_array), 4)
        % Determine the indices for the current time window
        i_panel_mod = mod(i_panel - 1, 4) + 1;
        plot_map = eof_select(:, i_panel);
        plot_map = reshape(plot_map, tgt_size);  % Reshape the map

        % Set title for the subplot
        title_text = ['EOF ' int2str(n_pc_array(i_panel))];

        % Create subplot (1 row, 4 columns)
        ax = subplot(1, 4, i_panel); % Assign to ax for customization
        axes_handles(i_panel) = ax; % Store the axes handle for shared colorbar

        % Plot using the geosurfm_sub function (presumably custom function for geospatial plotting)
        geosurfm_sub(plot_map, flag_land, ...
            tgt_lat, tgt_lon, eof_limit, ...
            [], [], title_text, ...
            lim_lat, lim_lon, ...
            eof_step, palette_eof, path_shp_file, res_fig);
    end

    % Set a global title for all subplots
    sgtitle(['First 4 EOFs for ' name_var ' (mth=' i_mth_txt ')'], 'FontSize', 16, 'FontWeight', 'bold');

    % Add a shared colorbar
    cb = colorbar('Location', 'eastoutside'); % Position the colorbar to the right of all subplots
    cb.Label.String = unit_var; % Add a label to the colorbar
    cb.FontSize = 10; % Adjust font size for clarity

    % Set ticks centered around 0 with step size of cstep
    c_lim = eof_limit;
    cstep = eof_step;
    tick_min = floor(c_lim(1) / cstep) * cstep;  % Round down to nearest multiple of cstep
    tick_max = ceil(c_lim(2) / cstep) * cstep;  % Round up to nearest multiple of cstep
    ticks = tick_min:cstep:tick_max;  % Generate the tick values

    % Set the ticks for the colorbar explicitly
    cb.Ticks = ticks;
    cb.Limits = c_lim;  % Ensure colorbar limits match c_lim

    % Adjust the position of the colorbar to span all subplots
    set(cb, 'Position', [0.93, 0.1, 0.02, 0.73]); % [left, bottom, width, height]

    % Define figure save path
    name_figuresave = fullfile(path_fig, [name_var '_eof_' i_mth_txt suffix]);

    % Save the figure
    print(f, name_figuresave, '-dtiff', res_fig);
end
