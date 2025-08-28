clear
close all

% Define the path to the shapefile
path_shp = 'ne_10m_admin_0_countries.shp';

% Read the shapefile
S = shaperead(path_shp, 'UseGeoCoords', true);

% Plot it
figure;
axesm('eqdcylin', 'Frame', 'on', 'Grid', 'on'); % Equidistant Cylindrical
geoshow(S, 'FaceColor', 'none', 'EdgeColor', 'k');

% Optional: title or inspect attributes
title('Admin-1 States/Provinces (Natural Earth 110m)');
