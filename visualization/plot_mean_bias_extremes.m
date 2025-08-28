function plot_mean_bias_extremes(dsValue_map, ESM_map, ...
    obsValue_test_map, ...
    lat, lon, flag_land, lim_lat, lim_lon, bias_limit, bias_step, ...
    palette_bias, n_color_bias, path_shp_file, res_fig, name_var, unit_var, ...
    path_fig, i_mth_txt, suffix)

% Reshape maps to space x time if needed
if ndims(dsValue_map) == 3
    dsValue_map = reshape(dsValue_map, [], size(dsValue_map, 3));
    ESM_map = reshape(ESM_map, [], size(ESM_map, 3));
    obsValue_test_map = reshape(obsValue_test_map, [], size(obsValue_test_map, 3));
end

tgt_size = [length(lon) length(lat)];

% Define the function to create figures for the two bias extremes (lower and upper)
    function plot_bias_for_extremes(data_map, data_map_esm, percentile_1, percentile_2, title_text_prefix, save_suffix)
        % Create a new figure with specified size
        f = figure('Position', [50, 0, 800, 700]);
        
        % Initialize a variable to collect all subplot axes
        axes_handles = gobjects(2, 2);  % 2x2 grid for both percentiles for PCR and ESM
        
        % Calculate bias over the entire period (using all data) for the first percentile
        plot_map_value_1 = prctile(data_map, percentile_1, 2) - prctile(obsValue_test_map, percentile_1, 2);
        plot_map_esm_1 = prctile(data_map_esm, percentile_1, 2) - prctile(obsValue_test_map, percentile_1, 2);

        % Calculate bias over the entire period (using all data) for the second percentile
        plot_map_value_2 = prctile(data_map, percentile_2, 2) - prctile(obsValue_test_map, percentile_2, 2);
        plot_map_esm_2 = prctile(data_map_esm, percentile_2, 2) - prctile(obsValue_test_map, percentile_2, 2);

        % Reshape the result to match the target size for both percentiles
        plot_map_value_1 = reshape(plot_map_value_1, tgt_size);
        plot_map_esm_1 = reshape(plot_map_esm_1, tgt_size);
        plot_map_value_2 = reshape(plot_map_value_2, tgt_size);
        plot_map_esm_2 = reshape(plot_map_esm_2, tgt_size);
        
        % Create a subplot for PCR (dsValue_map) bias at the first percentile (10th)
        ax1 = subplot(2, 2, 1);  % 2x2 grid, first plot (10th percentile PCR)
        axes_handles(1) = ax1;
        geosurfm_sub(plot_map_value_1, flag_land, ...
            lat, lon, bias_limit, ...
            [], [], [ 'PCR - ' num2str(percentile_1) 'p'], ...
            lim_lat, lim_lon, ...
            bias_step, palette_bias(n_color_bias), path_shp_file, res_fig);

        % adjust position
        posit_man = get(ax1, 'Position');
        posit_man(2) = posit_man(2) - 0.03; 
        posit_man(1) = posit_man(1) + 0.03; % move left slightly
        set(ax1, 'Position', posit_man);

        % Compute mean and standard deviation for PCR (dsValue_map) at the first percentile
        mean_val = mean(plot_map_value_1(flag_land), 'omitnan');
        std_val = std(plot_map_value_1(flag_land), 'omitnan');
        
        % Add mean ± SD below the subplot
        text(1, -0.03, sprintf(['%.2f ± %.2f ' unit_var], mean_val, std_val), ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 12);

        % Create a subplot for ESM bias at the first percentile (10th)
        ax2 = subplot(2, 2, 2);  % 2x2 grid, second plot (10th percentile ESM)
        axes_handles(2) = ax2;
        geosurfm_sub(plot_map_esm_1, flag_land, ...
            lat, lon, bias_limit, ...
            [], [], [ 'ESM - ' num2str(percentile_1) 'p'], ...
            lim_lat, lim_lon, ...
            bias_step, palette_bias(n_color_bias), path_shp_file, res_fig);
        % adjust position
        posit_man = get(ax2, 'Position');
        posit_man(2) = posit_man(2) - 0.03; 
        posit_man(1) = posit_man(1) - 0.03; 
        set(ax2, 'Position', posit_man);
        
        % Compute mean and standard deviation for ESM at the first percentile
        mean_val = mean(plot_map_esm_1(flag_land), 'omitnan');
        std_val = std(plot_map_esm_1(flag_land), 'omitnan');
        
        % Add mean ± SD below the subplot
        text(1, -0.03, sprintf(['%.2f ± %.2f ' unit_var], mean_val, std_val), ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 12);

        % Create a subplot for PCR (dsValue_map) bias at the second percentile (90th)
        ax3 = subplot(2, 2, 3);  % 2x2 grid, third plot (90th percentile PCR)
        axes_handles(3) = ax3;
        geosurfm_sub(plot_map_value_2, flag_land, ...
            lat, lon, bias_limit, ...
            [], [], [ 'PCR - ' num2str(percentile_2) 'p'], ...
            lim_lat, lim_lon, ...
            bias_step, palette_bias(n_color_bias), path_shp_file, res_fig);
        % adjust position
        posit_man = get(ax3, 'Position');
        posit_man(1) = posit_man(1) + 0.03; 
        set(ax3, 'Position', posit_man);
        
        % Compute mean and standard deviation for PCR (dsValue_map) at the second percentile
        mean_val = mean(plot_map_value_2(flag_land), 'omitnan');
        std_val = std(plot_map_value_2(flag_land), 'omitnan');
        
        % Add mean ± SD below the subplot
        text(1, -0.03, sprintf(['%.2f ± %.2f ' unit_var], mean_val, std_val), ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 12);

        % Create a subplot for ESM bias at the second percentile (90th)
        ax4 = subplot(2, 2, 4);  % 2x2 grid, fourth plot (90th percentile ESM)
        axes_handles(4) = ax4;
        geosurfm_sub(plot_map_esm_2, flag_land, ...
            lat, lon, bias_limit, ...
            [], [], [ 'ESM - ' num2str(percentile_2) 'p'], ...
            lim_lat, lim_lon, ...
            bias_step, palette_bias(n_color_bias), path_shp_file, res_fig);
        % adjust position
        posit_man = get(ax4, 'Position');
        posit_man(1) = posit_man(1) - 0.03; 
        set(ax4, 'Position', posit_man);

        % Compute mean and standard deviation for ESM at the second percentile
        mean_val = mean(plot_map_esm_2(flag_land), 'omitnan');
        std_val = std(plot_map_esm_2(flag_land), 'omitnan');
        
        % Add mean ± SD below the subplot
        text(1, -0.03, sprintf(['%.2f ± %.2f ' unit_var], mean_val, std_val), ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 12);

        % Set the title for the overall figure
        sgtitle([i_mth_txt ' ' title_text_prefix ' for ' name_var], 'FontSize', 16, 'FontWeight', 'bold');

        % Add a shared colorbar
        cb = colorbar('Location', 'eastoutside'); % Position the colorbar to the right of all subplots
        cb.Label.String = unit_var; % Add a label to the colorbar
        cb.FontSize = 14; % Adjust font size for clarity
        
        % Set ticks centered around 0 with step size of cstep
        c_lim = bias_limit;
        cstep = bias_step;
        tick_min = floor(c_lim(1) / cstep) * cstep;  % Round down to nearest multiple of cstep
        tick_max = ceil(c_lim(2) / cstep) * cstep;  % Round up to nearest multiple of cstep
        ticks = tick_min:cstep:tick_max;  % Generate the tick values
        
        % Set the ticks for the colorbar explicitly
        cb.Ticks = ticks;
        cb.Limits = c_lim;  % Ensure colorbar limits match c_lim

        % Adjust the position of the colorbar to span all subplots
        set(cb, 'Position', [0.93, 0.1, 0.024, 0.8]); % [left, bottom, width, height]
        
        % Define figure save path
        name_figuresave = fullfile(path_fig, [name_var '_pdiff ' save_suffix ' ' i_mth_txt suffix]);
        
        % Save the figure
        print(f, name_figuresave, '-dtiff', res_fig);
    end

% Call the function for both percentiles (10th and 90th)
plot_bias_for_extremes(dsValue_map, ESM_map, 10, 90, 'differences in percentiles', 'extremes_PCR');
end
