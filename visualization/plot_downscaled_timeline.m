function plot_downscaled_timeline(name_var, unit_var, i_mth_txt, obsValue_test_map, ESM_map, ...
    dsValue_map, dsPIinf_map, dsPIsup_map, path_fig, suffix, name_model, res_plot)
    % This function plots a timeline for various data variables with improved
    % color selection for accessibility (colorblind-friendly).
    
    if isempty(dsValue_map)
        return
    end


    % Reshape maps to space x time if needed
    if ndims(dsValue_map) == 3
        dsValue_map = reshape(dsValue_map, [], size(dsValue_map, 3));
        obsValue_test_map = reshape(obsValue_test_map, [], size(obsValue_test_map, 3));
        ESM_map = reshape(ESM_map, [], size(ESM_map, 3));
        dsPIinf_map = reshape(dsPIinf_map, [], size(dsPIinf_map, 3));
        dsPIsup_map = reshape(dsPIsup_map, [], size(dsPIsup_map, 3));
    end

    % Create a figure for plotting
    f = figure('Position', [50, 0, 1500, 300]);
    hold on;

    % Define a more colorblind-friendly palette (modified for higher contrast)
    colors = [
        0.2, 0.2, 0.2;   % Dark Gray for 'Original Data'
        0.1, 0.6, 0.6;   % Teal for 'Original ESM'
        0.0, 0.4, 0.8;   % Dark Blue for 'Downscaled PCR'
        0.9, 0.4, 0.1;   % Orange for 'PCR low 95% PI'
        0.95, 0.85, 0.15 % Yellow for 'PCR up 95% PI'
    ];

    % Plot original data
    plot(mean(obsValue_test_map, 1, 'omitnan'), 'LineWidth', 0.7, 'DisplayName', 'GHCN-m', 'Color', colors(1, :));

    % Plot target ESM monthly values
    plot(mean(ESM_map, 1, 'omitnan'), 'LineWidth', 0.7, 'DisplayName', name_model, 'Color', colors(2, :));

    % Plot downscaled PCR timeline
    plot(mean(dsValue_map, 1, 'omitnan'), 'LineWidth', 0.7, 'DisplayName', 'PCR Downscaled', 'Color', colors(3, :));
    
    % 30-y moving averages
    plot(movmean(mean(obsValue_test_map, 1, 'omitnan'), 30), 'LineWidth', 2, 'Color', colors(1, :), 'HandleVisibility', 'off');
    plot(movmean(mean(ESM_map, 1, 'omitnan'), 30), 'LineWidth', 2, 'Color', colors(2, :), 'HandleVisibility', 'off');
    plot(movmean(mean(dsValue_map, 1, 'omitnan'), 30), 'LineWidth', 2, 'Color', colors(3, :), 'HandleVisibility', 'off');
    
    % Plot downscaled PCR 95% PI (Lower Bound)
    plot(mean(dsPIinf_map, 1, 'omitnan'), 'LineWidth', 0.7, 'DisplayName', 'PCR low 95% PI', 'Color', colors(4, :));

    % Plot downscaled PCR 95% PI (Upper Bound)
    plot(mean(dsPIsup_map, 1, 'omitnan'), 'LineWidth', 0.7, 'DisplayName', 'PCR up 95% PI', 'Color', colors(5, :));

    % Add title and labels
    title([i_mth_txt ' ' name_var ' ' name_var ' timeline'], 'FontSize', 16, 'FontWeight', 'bold');
    xlabel('Year', 'FontSize', 14);
    ylabel(unit_var, 'FontSize', 14);

    % Adjust grid and axis
    grid on;
    set(gca, 'FontSize', 12);

    % Set X-axis limits and labels (ensure data fits within the plot)
    xlim([1, size(obsValue_test_map, 2)]);  % Adjust according to your data dimensions
    xticks(1:250:size(obsValue_test_map, 2)); 
    xticklabels(0:250:size(obsValue_test_map, 2)); % Customize if needed
    
    switch name_var
        case 'pr'
            ylim([0 ceil(prctile(mean(dsPIsup_map, 1, 'omitnan'),99))])
        case 'tas'
            ylim([floor(prctile(mean(dsPIinf_map, 1, 'omitnan'),1)) ceil(prctile(mean(dsPIsup_map, 1, 'omitnan'),99))])
    end

    % Add a legend
    legend('show', 'Location', 'northeastoutside', 'FontSize', 12);

    % Define figure save path
    name_figuresave = fullfile(path_fig, [name_var '_timeline_' i_mth_txt suffix]);

    % Save the figure
    print(f, name_figuresave, '-dtiff', res_plot);
end
