function plot_abs_field(dsValue_map, time_bound, year_start, flag_land, ...
    lat, lon, abs_limit, abs_step, palette_abs, path_shp_file, res_fig, ...
    name_var, i_mth_txt, unit_var, lim_lat, lim_lon, n_color_abs, path_fig, suffix)

% reshape maps to space x time if needed
if ndims(dsValue_map)==3
    dsValue_map = reshape(dsValue_map,[],size(dsValue_map,3));
end

tgt_size = [length(lon) length(lat)];
    % Create a new figure with specified size
    f = figure('Position', [50, 0, 1000, 420]);
    
    % Initialize a variable to collect all subplot axes
    axes_handles = gobjects(size(time_bound,1), 1);
    
    % Loop through panels
    for i_panel = 1:size(time_bound,1)
        % Determine the indices for the current time window
        i_panel_mod = mod(i_panel - 1, size(time_bound,1)) + 1;
        plot_map = mean(dsValue_map(:, (time_bound(i_panel_mod,1) - year_start + 1):(time_bound(i_panel_mod,2) - year_start  + 1)), 2, 'omitnan');
        plot_map = reshape(plot_map, tgt_size);
        
        % Set title for the subplot
        title_text = [int2str(time_bound(i_panel_mod,1)) '-' int2str(time_bound(i_panel_mod,2))];
        
        % Create subplot 
        ax = subplot(1, size(time_bound,1), i_panel); % Assign to ax for customization
        axes_handles(i_panel) = ax; % Store the axes handle for shared colorbar
    
        % After plotting first map
        pos1 = get(ax, 'Position');
        pos1(1) = pos1(1) * 0.98; 
        pos1(2) = pos1(2) - 0.06; 
        set(ax, 'Position', pos1);

        % Plot using the geosurfm function
        geosurfm_sub(plot_map, flag_land, ...
            lat, lon, abs_limit, ...
            [], [], title_text, ...
            lim_lat, lim_lon, ...
            abs_step, palette_abs(n_color_abs), path_shp_file, res_fig);

        % Compute mean and standard deviation
        mean_val = mean(plot_map(:), 'omitnan');
        std_val = std(plot_map(:), 'omitnan');
        
        % Add mean ± SD below the subplot
        text(1, -0.03, sprintf(['%.2f ± %.2f ' unit_var], mean_val, std_val), ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 14);
    end

    % Set the overall title for the figure
    sgtitle([i_mth_txt ' mean PCR-downscaled ' name_var], 'FontSize', 16, 'FontWeight', 'bold');

    % Add a shared colorbar
    cb = colorbar('Location', 'eastoutside'); % Position the colorbar to the right of all subplots
    cb.Label.String = unit_var; % Add a label to the colorbar
    cb.FontSize = 14; % Adjust font size for clarity
    
    % Set ticks centered around 0 with step size of abs_step
    c_lim = abs_limit;
    cstep = abs_step;
    tick_min = floor(c_lim(1) / cstep) * cstep;  % Round down to nearest multiple of cstep
    tick_max = ceil(c_lim(2) / cstep) * cstep;  % Round up to nearest multiple of cstep
    ticks = tick_min:cstep:tick_max;  % Generate the tick values

    % Set the ticks for the colorbar explicitly
    cb.Ticks = ticks;
    cb.Limits = c_lim;  % Ensure colorbar limits match c_lim

    % Adjust the position of the colorbar to span all subplots
    set(cb, 'Position', [0.93, 0.1, 0.02, 0.73]); % [left, bottom, width, height]
    
    % Define figure save path
    name_figuresave = fullfile(path_fig, [name_var '_abs field_' i_mth_txt suffix]);
    
    % Save the figure
    print(f, name_figuresave, '-dtiff', res_fig);
end
