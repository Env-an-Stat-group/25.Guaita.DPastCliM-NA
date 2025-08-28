function newmap = color10_abs(m)
% color map for TOAR special issue guidelines 2023/24.

if ~exist('m','var') 
    m = 10;
end


color_list = [...
    0.2081    0.1663    0.5292   %%% dark blue            1 
    0.3961    0.3176    0.8000   %%% purple               2 
    0.0123    0.4213    0.8802   %%% blue                 3 
    0.4941    0.7647    0.8980   %%% light blue           4 
    0.1157    0.7022    0.6843   %%% blue-green           5 
    0.5216    0.6980    0.1725   %%% muted green          6 
    0.9968    0.7513    0.2325   %%% yellow-orange        7 
    1.0000    0.4863    0.0000   %%% dark orange          8 
    0.8000    0.3176    0.3176   %%% light reddish brown  9 
    0.6980    0.1725    0.1725   %%% reddish brown        10
    ];

newmap = interp1(linspace(0,1,size(color_list,1)),...
    color_list,...
    linspace(0,1,m));


end