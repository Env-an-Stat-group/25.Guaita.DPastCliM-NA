function geosurfm(matrix_map, filter,...
    lat, lon, c_lim, ...
    save_name, clabel, title_text,...
    latlim_map, lonlim_map, ...
    cstep, palette, shp_file_path,res)
% function to plot on the domain with borders
% INPUT
% matrix_map = matrix with the data to plot. 
%              latitude on rows, longitude on columns
% filter = boolean matrix to filter the data (same dimension as matrix_map)
% lat, lon = latitude and longitude arrays
% c_lim = colorbar limits
% latlim_map, lonlim_map = latitude longitude limits to cut the map
% differently from the original coordinates
% cstep = the step you want to use for the colorbar
% palette = the colormap palette you want to use e.g. flipud(hot)
% shp_file_path = full path to the .shp file

matrix_map(not(filter))=nan;

if not(exist('res','var'))
    res = '-r300';
end

% size of the map, i.e. number of nodes
[lat_mat, lon_mat] = meshgrid(lat,lon);

if not(isempty(save_name))
    f = figure('Name',save_name,'Position',[100 100 1000 500]);
end

% plot world map borders
ax = worldmap(latlim_map, lonlim_map);
setm(ax);
hold on

%%plot map
surfacem(lat_mat,lon_mat,matrix_map);
clim(c_lim)
colormap(palette)
c=colorbar;
c.Label.String = clabel;
c.Ticks=(round(c_lim(1)/cstep)*cstep):cstep:(round(c_lim(2)/cstep)*cstep);
title(title_text)
% fix view options
setm(gca,'fontsize',10,'fontweight','bold','glinewidth',1);
set(gca,'FontWeight','bold','FontSize',15,'Tickdir','out','linewidth',1.5,'box','on')

% map shp file
if exist('shp_file_path','var')
    geoshow(shp_file_path,'DefaultFaceColor', 'none', 'DefaultEdgeColor', 'black','LineWidth', 0.9)
end

if not(isempty(save_name))
    print(f,save_name,'-dtiff',res);
end


end