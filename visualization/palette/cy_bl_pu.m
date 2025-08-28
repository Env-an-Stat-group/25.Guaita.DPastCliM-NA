function newmap = cy_bl_pu(m)
% color map for TOAR special issue guidelines 2023/24.

if ~exist('m','var') 
    m = 10;
end


color_list = [...
   230   245   252
   175   220   243
    77   182   234
   113   152   223
   149    93   213
   180    19   201
    ]/255;

newmap = interp1(linspace(0,1,size(color_list,1)),...
    color_list,...
    linspace(0,1,m));


end