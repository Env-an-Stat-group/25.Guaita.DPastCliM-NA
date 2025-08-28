function geosurfm_sub_dot(lat, lon, dot_size, color, title_text, latlim_map, lonlim_map, path_shp_file, res)
% function to plot dots on a map
% INPUT
% lat = array of latitudes for the dots
% lon = array of longitudes for the dots
% dot_size = size of the dots to plot
% color = color of the dots
% title_text = title of the plot
% latlim_map, lonlim_map = latitude and longitude limits to cut the map
% shp_file_path = full path to the .shp file (optional)
% save_name = name of the file to save the plot (optional)

if not(exist('res','var'))
    res = '-r300';
end

% plot world map borders
ax = worldmap(latlim_map, lonlim_map);
setm(ax);
hold on

% Plot the dots
geoshow(lat, lon, 'DisplayType', 'point', 'Marker', '.', 'Color', color, 'MarkerSize', dot_size);
title(title_text,'fontsize',10,'fontweight','bold')
% fix view options
% Hide latitude and longitude labels and ticks
setm(ax, 'MLabelLocation', [], 'PLabelLocation', [], 'MLineLocation', [], 'PLineLocation', []);

% Plot shapefile if provided
if exist('path_shp_file', 'var') && ~isempty(path_shp_file)
    geoshow(path_shp_file, 'DefaultFaceColor', 'none', 'DefaultEdgeColor', 'black', 'LineWidth', 0.9);
end


end
