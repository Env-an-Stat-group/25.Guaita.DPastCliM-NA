function geosurfm_meansdlabel(matrix_map, filter,...
    lat, lon, c_lim, ...
    save_name, clabel, title_text,...
    latlim_map, lonlim_map, ...
    cstep, palette, shp_file_path,res)
% function to plot on the domain with borders and annotate mean ± SD

matrix_map(~filter) = nan;

if ~exist('res','var')
    res = '-r300';
end

% size of the map
[lat_mat, lon_mat] = meshgrid(lat, lon);

if ~isempty(save_name)
    f = figure('Name', save_name, 'Position', [100 100 1000 250]);
end

% plot world map borders
ax = worldmap(latlim_map, lonlim_map);
setm(ax);

% turn off latitude and longitude labels
setm(ax, 'MeridianLabel', 'off', 'ParallelLabel', 'off', ...
    'MLineLocation', NaN, 'PLineLocation', NaN, 'MapProjection', 'eckert4');

hold on

%% plot map
surfacem(lat_mat, lon_mat, matrix_map);
clim(c_lim)
colormap(palette)
c = colorbar;
c.Label.String = clabel;
c.Ticks = (round(c_lim(1)/cstep)*cstep):cstep:(round(c_lim(2)/cstep)*cstep);
title(title_text)

% fix view options
setm(gca,'fontsize',18,'fontweight','bold','glinewidth',1);
set(gca,'FontWeight','bold','FontSize',18,'Tickdir','out','linewidth',1.5,'box','on')


% position colorbar south outside
set(c, 'Location', 'southoutside');
c.FontSize =   18;

% map shp file
if exist('shp_file_path','var')
    geoshow(shp_file_path,'DefaultFaceColor', 'none', 'DefaultEdgeColor', 'black','LineWidth', 0.9)
end

%% Add mean ± SD annotation
mean_val = mean(matrix_map(filter), 'omitnan');
std_val  = std(matrix_map(filter), 'omitnan');

% Place text in normalized figure coordinates (bottom-left)
text(0.02, 0.08, sprintf('%.2f ± %.2f', mean_val, std_val), ...
    'Units', 'normalized', 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');

if ~isempty(save_name)
    print(f, [save_name '.png'], '-dpng', res);
end

end
