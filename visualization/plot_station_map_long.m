function plot_station_map_long(meta, lim_lat, lim_lon, path_shp_file, fig_path, plot_title)

    figure('Color','w','Position',[100 100 1000 600]);

    % Load world borders
    world = shaperead(path_shp_file,'UseGeoCoords',true);

    % Colorblind-friendly colors
    color_cal = [213 94 0]/255;  % orange-red
    
    % Count stations using height of meta table
    n_stations = height(meta);

    % KDE for calibration and test stations
    lat_smooth = linspace(lim_lat(1), lim_lat(2), 200);
    [f_cal, ~] = ksdensity(meta.lat, lat_smooth, 'Bandwidth', 2);

    % --- Map subplot (left) ---
    ax_map = subplot(1,5,1:4);
    hold on
    for k = 1:length(world)
        plot(world(k).Lon, world(k).Lat, 'k');
    end

    % Plot stations with smaller markers using new colors
    plot(meta.lon, meta.lat, 'o', ...
         'MarkerFaceColor', color_cal, 'MarkerEdgeColor', color_cal, 'MarkerSize', 3.5);

    % Add total number of stations in the bottom-left corner
    text(lim_lon(1)+0.03*range(lim_lon), lim_lat(1)+0.04*range(lim_lat), ...
        sprintf('# stations: %d', n_stations), ...
        'FontSize', 12, 'Color', 'k', 'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'left', 'BackgroundColor', 'w', 'EdgeColor', 'k');

    xlim(lim_lon)
    ylim(lim_lat)
    xlabel('Longitude', 'FontSize', 14)
    ylabel('Latitude', 'FontSize', 14)
    title(plot_title, 'FontSize', 16)
    set(ax_map, 'FontSize', 12)
    box on

    % --- Density subplot (right) ---
    ax_density = subplot(1,5,5);
    hold on
    plot(f_cal, lat_smooth, 'Color', color_cal, 'LineWidth', 2)
    
    % Set axes limits and ticks
    xlim([0, max(f_cal) * 1.1]);  % Add 10% padding
    ylim(lim_lat)
    set(gca, 'YAxisLocation', 'right', ...
             'FontSize', 12)
    title('KDE', 'FontSize', 12)
    xlabel('Density', 'FontSize', 12)  % x-axis label
    box on



    % Save figure
    saveas(gcf, fig_path);
end
