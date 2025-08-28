function newmap = Tol_seq_YlOrBr(m)
% Iridescent colorbar from Tol sequential color palettes

if ~exist('m','var') 
    m = 9;
end


color_list = [...
    255,255,229
    255,247,188
    254,227,145
    254,196,79
    251,154,41
    236,112,20
    204,76,2
    153,52,4
    102,37,6
    ]/255;

newmap = interp1(linspace(0,1,size(color_list,1)),...
    color_list,...
    linspace(0,1,m));


end