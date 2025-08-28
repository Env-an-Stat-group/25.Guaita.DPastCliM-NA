function plot_bias_testing_two_plot(dsValue_map, obsValue_test_map, ESM_map, ...
    flag_land, lat, lon, bias_limit, bias_step, palette_bias, path_shp_file, res_fig, ...
    name_var, i_mth_txt, unit_var, lim_lat, lim_lon, n_color_bias, path_fig, suffix)

% Reshape maps to space x time if needed
if ndims(dsValue_map) == 3
    dsValue_map = reshape(dsValue_map, [], size(dsValue_map, 3));
    obsValue_test_map = reshape(obsValue_test_map, [], size(obsValue_test_map, 3));
    ESM_map = reshape(ESM_map, [], size(ESM_map, 3));
end
tgt_size = [length(lon) length(lat)];

% Create a new figure with specified size
f = figure('Position', [50, 0, 1300, 650]); % Adjusted for two subplots

%% PCR bias
% Calculate PCR bias
plot_map = mean(dsValue_map - obsValue_test_map, 2, 'omitnan');
plot_map = reshape(plot_map, tgt_size);

% Set title for the subplot
title_text = 'PCR';

% Create subplot (1 row, 2 columns, current panel index)
ax1 = subplot(1, 2, 1); % First subplot for PCR

% Plot using the geosurfm function
geosurfm_sub(plot_map, flag_land, ...
    lat, lon, bias_limit, ...
    [], [], title_text, ...
    lim_lat, lim_lon, ...
    bias_step, palette_bias(n_color_bias), path_shp_file, res_fig);

% After plotting first map
pos1 = get(ax1, 'Position');
pos1(1) = pos1(1) + 0.03; % move right slightly
set(ax1, 'Position', pos1);

% Compute mean and standard deviation
mean_val = mean(plot_map(flag_land), 'omitnan');
std_val = std(plot_map(flag_land), 'omitnan');

% Add mean ± SD below the subplot
text(1, -0.03, sprintf(['%.2f ± %.2f ' unit_var], mean_val, std_val), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 12);

%% ESM bias
% Calculate ESM bias
plot_map = mean(ESM_map - obsValue_test_map, 2, 'omitnan');
plot_map = reshape(plot_map, tgt_size);

% Set title for the subplot
title_text = 'ESM';

% Create subplot (1 row, 2 columns, current panel index)
ax2 = subplot(1, 2, 2); % Second subplot for ESM

% Plot using the geosurfm function
geosurfm_sub(plot_map, flag_land, ...
    lat, lon, bias_limit, ...
    [], [], title_text, ...
    lim_lat, lim_lon, ...
    bias_step, palette_bias(n_color_bias), path_shp_file, res_fig);

% After plotting second map
pos2 = get(ax2, 'Position');
pos2(1) = pos2(1) - 0.03; % move left slightly
set(ax2, 'Position', pos2);


% Compute mean and standard deviation
mean_val = mean(plot_map(flag_land), 'omitnan');
std_val = std(plot_map(flag_land), 'omitnan');

% Add mean ± SD below the subplot
text(1, -0.03, sprintf(['%.2f ± %.2f ' unit_var], mean_val, std_val), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 14);

%% Set the overall title for the figure
sgtitle([i_mth_txt ' MB of ' name_var ], 'FontSize', 20, 'FontWeight', 'bold');

%% Add a shared colorbar
cb = colorbar('Location', 'eastoutside'); % Position the colorbar to the right of all subplots
cb.Label.String = unit_var; % Add a label to the colorbar
cb.FontSize = 14; % Adjust font size for clarity

% Set ticks centered around 0 with step size of bias_step
c_lim = bias_limit;
cstep = bias_step;
tick_min = floor(c_lim(1) / cstep) * cstep;  % Round down to nearest multiple of cstep
tick_max = ceil(c_lim(2) / cstep) * cstep;  % Round up to nearest multiple of cstep
ticks = tick_min:cstep:tick_max;  % Generate the tick values

% Set the ticks for the colorbar explicitly
cb.Ticks = ticks;
cb.Limits = c_lim;  % Ensure colorbar limits match c_lim

% Adjust the position of the colorbar to span all subplots
set(cb, 'Position', [0.93, 0.1, 0.015, 0.8]); % [left, bottom, width, height]

%% Define figure save path
name_figuresave = fullfile(path_fig, [name_var '_mean bias_test period_' i_mth_txt suffix]);

% Save the figure
print(f, name_figuresave, '-dtiff', res_fig);

end
