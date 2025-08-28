function newmap = br_gr_y(m)
% color map with 12 colors for absolute values.

if ~exist('m','var') 
    m = 12;
end


color_list = [    0.2081    0.1663    0.5292   %%% dark blue 
    0.1403    0.3147    0.8168 
    0.0410    0.4502    0.8685 
    0.0734    0.5410    0.8257 
    0.0232    0.6407    0.7925 
    0.1024    0.6984    0.6934 
    0.3187    0.7395    0.5625   %%% green 
    0.5745    0.7484    0.4479 
    0.7798    0.7361    0.3658 
    0.9613    0.7281    0.2774 
    0.9763    0.8328    0.1590 
    0.9763    0.9831    0.0538   %%% yellow   
     ];

newmap = interp1(linspace(0,1,size(color_list,1)),...
    color_list,...
    linspace(0,1,m));


end