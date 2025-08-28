function newmap = Tol_seq_iridescent(m)
% Iridescent colorbar from Tol sequential color palettes

if ~exist('m','var') 
    m = 8;
end


color_list = [...
    254,251,233
    234,240,181
    194,227,210
    155,210,225
    123,188,231
    147,152,210
    154,112,158
    104,73,87
    ]/255;

newmap = interp1(linspace(0,1,size(color_list,1)),...
    color_list,...
    linspace(0,1,m));


end