function [newmap, newlon, newlat] = cmipcoord_2_map(oldmap, oldlon, oldlat)
% this function is used for changing a map coordinate system from using a
% lon range of (0, 360) to a range of (-180, 180). Additionally sorts the
% latitude and longitude arrays in ascending order, fixing the map
% accordingly

% fix old longitude values
oldlon(oldlon>180)=oldlon(oldlon>180)-360;

% sort by longitude both the longitude vector and the map
[newlon,sortind_lon]=sort(oldlon);
newmap = oldmap(sortind_lon,:,:);
% do the same with latitude, updating the new map
[newlat,sortind_lat]=sort(oldlat);
newmap = newmap(:,sortind_lat,:);

end