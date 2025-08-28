function newmap = Tol_div_PRGn(m)
% Iridescent colorbar from Tol sequential color palettes

if ~exist('m','var') 
    m = 10;
end


color_list = [...
    118,42,131
    153,112,171
    194,165,207
    231,212,232
    247,247,247
    247,247,247
    217,240,211
    172,211,158
    90,174,97
    27,120,55
    ]/255;

newmap = interp1(linspace(0,1,size(color_list,1)),...
    color_list,...
    linspace(0,1,m));


end