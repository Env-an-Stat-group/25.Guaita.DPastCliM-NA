function plot_single(plot_map, tgt_size, flag_land, lat, lon, bias_limit, bias_step, ...
    palette, path_shp_file, res_fig, text_title, sub_title, unit_var, ...
    lim_lat, lim_lon, path_fig, suffix)

plot_map = reshape(plot_map, tgt_size);

% Create figure
f = figure('Position', [100, 50, 800, 600]);

% Plot using geosurfm_sub
geosurfm_sub(plot_map, flag_land, ...
    lat, lon, bias_limit, [], [], sub_title, ...
    lim_lat, lim_lon, bias_step, palette, path_shp_file, res_fig);

% Compute mean ± SD
mean_val = mean(plot_map(flag_land), 'omitnan');
std_val = std(plot_map(flag_land), 'omitnan');
text(1, -0.03, sprintf('%.2f ± %.2f %s', mean_val, std_val, unit_var), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'Units', 'normalized', 'FontSize', 12);

% Overall title
sgtitle([text_title ' - ' sub_title], 'FontSize', 18, 'FontWeight', 'bold');

% Add colorbar
cb = colorbar('eastoutside');
cb.Label.String = unit_var;
cb.FontSize = 12;
tick_min = floor(bias_limit(1)/bias_step)*bias_step;
tick_max = ceil(bias_limit(2)/bias_step)*bias_step;
cb.Ticks = tick_min:bias_step:tick_max;
cb.Limits = bias_limit;

% Save figure
name_fig = fullfile(path_fig, [text_title '_' sub_title suffix]);
print(f, name_fig, '-dtiff', res_fig);

end