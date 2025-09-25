function plot_anomaly_ESM(dsValue_anom_map, time_bound, year_start, flag_land, ...
    lat, lon, bias_limit, bias_step, palette_bias, path_shp_file, res_fig, ...
    name_var, i_mth_txt, unit_var, lim_lat, lim_lon, n_color_bias, path_fig, suffix)

    % reshape maps to space x time if needed
    if ndims(dsValue_anom_map)==3
        dsValue_anom_map = reshape(dsValue_anom_map,[],size(dsValue_anom_map,3));
    end
    if ndims(dsValue_anom_map)==4
        dsValue_anom_map = reshape(dsValue_anom_map,[],size(dsValue_anom_map,3)*size(dsValue_anom_map,4));
    end

    tgt_size = [length(lon) length(lat)];
    % Create a new figure with specified size
    f = figure('Position', [50, 0, 1000, 500]);
    
    % Initialize a variable to collect all subplot axes
    %axes_handles = gobjects(size(time_bound,1), 1);
    
    % Loop through panels
    for i_panel = 1:size(time_bound,1)
        % Determine the indices for the current time window
        i_panel_mod = mod(i_panel - 1, size(time_bound,1)) + 1;
        plot_map = mean(dsValue_anom_map(:, (time_bound(i_panel_mod,1) - year_start + 1):(time_bound(i_panel_mod,2) - year_start  + 1)), 2, 'omitnan');
        plot_map = reshape(plot_map, tgt_size);
        
        % Set title for the subplot
        title_text = [int2str(time_bound(i_panel_mod,1)) '-' int2str(time_bound(i_panel_mod,2))];
        
        % Create subplot 
        ax = subplot(1, size(time_bound,1), i_panel); % Assign to ax for customization    

        % Plot using the geosurfm function
        geosurfm_sub(plot_map, flag_land & not(isnan(plot_map)), ...
            lat, lon, bias_limit, ...
            [], [], title_text, ...
            lim_lat, lim_lon, ...
            bias_step, palette_bias(n_color_bias), path_shp_file, res_fig);

        % Compute mean and standard deviation
        mean_val = mean(plot_map(:), 'omitnan');
        std_val = std(plot_map(:), 'omitnan');
        
        % Place text in normalized figure coordinates (bottom-left)
        text(0.02, 0.08, sprintf('%.2f ± %.2f', mean_val, std_val), ...
            'Units', 'normalized', 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'bottom', 'FontSize', 15, 'FontWeight', 'bold', 'Color', 'k');
    end

    % Set the overall title for the figure
    sgtitle([i_mth_txt ' differences in mean ESM ' name_var], 'FontSize', 18, 'FontWeight', 'bold');

    % Add a shared colorbar at the bottom
    cb = colorbar('Location', 'southoutside');
    cb.Label.String = unit_var;
    cb.FontSize =   15;
    cb.FontWeight = 'bold';
    
    % Set ticks centered around 0 with step size
    c_lim = bias_limit;   % or abs_limit depending on function
    cstep = bias_step;    % or abs_step depending on function
    tick_min = floor(c_lim(1) / cstep) * cstep;
    tick_max = ceil(c_lim(2) / cstep) * cstep;
    ticks = tick_min:cstep:tick_max;
    cb.Ticks = ticks;
    cb.Limits = c_lim;
    
    % Position colorbar nicely below all subplots
    cb.Position = [0.25, 0.2, 0.6, 0.03];  % [x, y, width, height]
        
    % Define figure save path
    name_figuresave = fullfile(path_fig, [name_var '_ESM_anomaly field_' i_mth_txt suffix]);
    
    % Save the figure
    print(f, name_figuresave, '-dtiff', res_fig);
end
