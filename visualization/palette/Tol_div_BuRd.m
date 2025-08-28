function newmap = Tol_div_BuRd(m)
% Iridescent colorbar from Tol sequential color palettes

if ~exist('m','var') 
    m = 10;
end


color_list = [...
    33,102,172
    67,147,195
    146,197,222
    209,229,240
    247,247,247
    253,219,199
    244,165,130
    214,96,77
    178,24,43]/255;

newmap = interp1(linspace(0,1,size(color_list,1)),...
    color_list,...
    linspace(0,1,m));


end