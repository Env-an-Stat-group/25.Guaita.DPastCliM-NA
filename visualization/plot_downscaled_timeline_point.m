function plot_downscaled_timeline_point(name_var, unit_var, time_ESM, time_obs_tmp, ...
    y_obs, y_esm, y_ds, y_inf, y_sup, ...
    path_fig, suffix, name_model, res_plot, x_limits, opt_print_2kmean)
    % Plot a timeline of observations, ESM, and downscaled data
    % with yearly means and shaded prediction intervals.
    %
    % Inputs:
    %   time_ESM     - time vector for ESM (months since 0 CE)
    %   time_obs_tmp - time vector for observations (months since 0 CE)
    %   y_obs        - observed values (monthly)
    %   y_esm        - ESM values (monthly)
    %   y_ds         - downscaled values (monthly)
    %   y_inf        - lower PI bound (monthly)
    %   y_sup        - upper PI bound (monthly)

    % --- Create figure ---
    f = figure('Position', [50, 100, 700, 300]);
    hold on;

    % --- Improved colorblind-friendly palette ---
    colors = [
        0.2, 0.2, 0.2;    % Dark Gray for observations
        0.90, 0.45, 0.30; % Vermilion for ESM
        0.35, 0.70, 0.90; % Sky Blue for downscaled PCR
        0.7, 0.7, 0.7; % Grey for PI
    ];

    % --- Convert months to years ---
    years_obs_tmp = (double(time_obs_tmp) / 12);
    years_obs = (min(double(time_obs_tmp)):max(double(time_obs_tmp))) / 12;
    years_esm = (time_ESM / 12);

    % correct y_obs with nans
    y_obs_tmp = y_obs;
    y_obs = nan(length(years_obs),1);
    y_obs(ismember(years_obs,years_obs_tmp)) = y_obs_tmp;

    % --- Shaded PI area (same color as downscaled PCR, transparent) ---
    fill([years_esm; flipud(years_esm)], ...
        [y_inf'; flipud(y_sup')], colors(4,:), ...
         'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    % --- Plot lines ---
    plot(years_esm, y_esm, 'LineWidth', 0.8, 'Color', colors(2,:), 'DisplayName', name_model);
    plot(years_obs, y_obs, 'LineWidth', 0.8, 'Color', colors(1,:), 'DisplayName', 'GHCN-m');
    plot(years_esm, y_ds,  'LineWidth', 1, 'Color', colors(3,:), 'DisplayName', 'PCR downscaled');

    % --- Calculate mean over the time window
    window_centre = round(((x_limits(1)*12+1)+(x_limits(2)*12))/2);
    mean_window = mean(y_ds((window_centre-15*12+1):(window_centre+15*12)));
    mean_all = mean(y_ds);
    anomaly = mean_window-mean_all;

    % --- Add horizontal line for mean over selected window ---
    plot([x_limits(1)-3, x_limits(2)+3], [mean_window, mean_window], ...
        '--', 'Color', colors(3,:), 'LineWidth', 1);    
    
    % --- Title and labels ---
    title([name_var ' timeline in ' suffix ' ' int2str(x_limits(1)) '-' int2str(x_limits(2))], 'FontSize', 16, 'FontWeight', 'bold');
    xlabel('Year', 'FontSize', 18);
    ylabel(unit_var, 'FontSize', 18);

    % --- Grid and axis formatting ---
    grid on;
    set(gca, 'FontSize', 18);

    % Define ticks every 100 years
    %maxYear = max([max(years_obs), max(years_esm)]);  % full span
    %yearTicks = 0:100:maxYear;
    xticks(x_limits(1):5:x_limits(2));
    %xticklabels(yearTicks);
    xlim([x_limits(1)-0.25 x_limits(2)+0.25]);

    % --- Variable-specific y-limits ---
    switch name_var
        case 'pr'
            ylim([0 ceil(prctile([y_sup'; y_ds'; y_esm; y_obs],99.5,'all'))])
        case 'tas'
            ylim([floor(min([y_inf'; y_ds'; y_esm; y_obs])) ceil(max([y_sup'; y_ds'; y_esm; y_obs]))])
    end

    % --- Add anomaly text in bottom-left corner ---
    yl = ylim;
    xl = xlim;
    text(xl(1) + 0.02*(xl(2)-xl(1)), yl(2) - 0.05*(yl(2)-yl(1)), ...
        sprintf('Anomaly (30-yr): %.2f %s', anomaly, unit_var), ...
        'FontSize', 14, 'Color', colors(3,:), 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

    if opt_print_2kmean
        text(xl(2) - 0.02*(xl(2)-xl(1)), yl(2) - 0.05*(yl(2)-yl(1)), ...
            sprintf('2k mean: %.2f %s', mean_all, unit_var), ...
            'FontSize', 14, 'Color', colors(3,:), 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    end

    % --- Save the figure ---
    name_figuresave = fullfile(path_fig, [name_var '_timeline_' suffix '_' int2str(x_limits(1)) '-' int2str(x_limits(2))]);
    print(f, [name_figuresave '.png'], '-dpng', res_plot);
end
