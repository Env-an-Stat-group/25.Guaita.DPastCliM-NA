function newmap = Tol_div_sunset(m)
% Iridescent colorbar from Tol sequential color palettes

if ~exist('m','var') 
    m = 11;
end


color_list = [...
    54,75,154
    74,123,183
    110,166,205
    152,202,225
    194,228,239
    234,236,204
    254,218,139
    253,179,102
    246,126,75
    221,61,45
    165,0,38
    ]/255;

newmap = interp1(linspace(0,1,size(color_list,1)),...
    color_list,...
    linspace(0,1,m));


end