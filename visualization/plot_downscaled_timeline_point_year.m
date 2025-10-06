function plot_downscaled_timeline_point_year(name_var, unit_var, time_ESM, time_obs, ...
    y_obs, y_esm, y_ds, ...
    path_fig, suffix, name_model, res_plot)
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
    f = figure('Position', [50, 100, 1200, 400]);
    hold on;

    % --- Improved colorblind-friendly palette ---
    colors = [
        0.2, 0.2, 0.2;    % Dark Gray for observations
        0.90, 0.45, 0.30; % Vermilion for ESM
        0.35, 0.70, 0.90; % Sky Blue for downscaled PCR
        0.7, 0.7, 0.7; % Grey for PI
    ];

    % Adjust obs array so that it matches the esm and downscaling arrays

    % --- Convert months to years ---
    
    years_esm = unique(floor(time_ESM / 12));
    years_esm = years_esm(1:(end-1));

    % correct y_obs with nans
    y_ds = mean(reshape(y_ds,12,[]),1,'omitnan');
    y_esm = mean(reshape(y_esm,12,[]),1,'omitnan');
    [~,idx_obs] = ismember(time_obs,time_ESM);
    y_obs_all = nan(size(time_ESM));
    y_obs_all(idx_obs) = y_obs;
    y_obs = mean(reshape(y_obs_all,12,[]),1);
    clear y_obs_all
    
    % --- Plot lines ---
    plot(years_esm, y_esm, 'LineWidth', 0.3, 'Color', colors(2,:), 'DisplayName', name_model);
    plot(years_esm, y_obs, 'LineWidth', 0.3, 'Color', colors(1,:), 'DisplayName', 'GHCN-m');
    plot(years_esm, y_ds,  'LineWidth', 0.3, 'Color', colors(3,:), 'DisplayName', 'PCR downscaled');

    % --- Plot 30-yr movmean---
    plot(years_esm, movmean(y_esm,30,'omitnan'), 'LineWidth', 4, 'Color', 0.8*colors(2,:), 'DisplayName', name_model);
    plot(years_esm, movmean(y_obs,30,'omitnan'), 'LineWidth', 4, 'Color', 0.8*colors(1,:), 'DisplayName', 'GHCN-m');
    plot(years_esm, movmean(y_ds,30,'omitnan'),  'LineWidth', 4, 'Color', 0.8*colors(3,:), 'DisplayName', 'PCR downscaled');

    % --- Calculate mean over the time window
    mean_all = mean(y_ds);

    % --- Add horizontal line for mean over selected window ---
    plot([years_esm(1), years_esm(end)], [mean_all, mean_all], ...
        '--', 'Color', 0.5*colors(3,:), 'LineWidth', 1);    
    
    % --- Title and labels ---
    title([name_var ' timeline in ' suffix ], 'FontSize', 16, 'FontWeight', 'bold');
    xlabel('Year', 'FontSize', 18);
    ylabel(unit_var, 'FontSize', 18);

    % --- Grid and axis formatting ---
    grid on;
    set(gca, 'FontSize', 18);

    % Define ticks every 100 years
    %maxYear = max([max(years_obs), max(years_esm)]);  % full span
    %yearTicks = 0:100:maxYear;
    xticks(years_esm(1):100:years_esm(end));
    %xticklabels(yearTicks);
    xlim([years_esm(1) years_esm(end)]);

    % --- Variable-specific y-limits ---
    switch name_var
        case 'pr'
            ylim([0 ceil(max([y_ds, y_esm, y_obs]))])
        case 'tas'
            ylim([floor(min([y_ds y_esm, y_obs])) ceil(max([y_ds, y_esm, y_obs]))])
    end

    % --- Add anomaly text in bottom-left corner ---
    yl = ylim;
    xl = xlim;
    text(xl(2) - 0.02*(xl(2)-xl(1)), yl(1) + 0.05*(yl(2)-yl(1)), ...
            sprintf('2k mean: %.2f %s', mean_all, unit_var), ...
            'FontSize', 14, 'Color', colors(3,:), 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

    % --- Save the figure ---
    name_figuresave = fullfile(path_fig, [name_var '_timeline_' suffix '_year']);
    print(f, [name_figuresave '.png'], '-dpng', res_plot);
end
