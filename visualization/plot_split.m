% this code is for visualizing the dataset split in calibration/test/model
% selection datasets


path_main = '/data/pguaita/downscaling/';

name_var = 'tas';

year_start = 1875;
n_year_data = 140;

name_model = 'MPI-ESM1-2-LR';
path_modeldata = fullfile(path_main,['downscaling_models_' name_model]);
path_fig = fullfile(path_modeldata,'figures_PCR');
suffix = '_NA_020';

disp('defined parameters')

%% load metaTable and observations (monthly) and aggregate them

for i_mth = 1:12
    txt_mth = [int2str(i_mth) '-' int2str(i_mth)];

    load(fullfile(path_modeldata,[name_var '_PCR_original_data_mth=' txt_mth suffix '.mat']))
    if i_mth==1
        obsTable = obsTable_mth;
    else
        obsTable = vertcat(obsTable,obsTable_mth);
    end
end

disp('loaded data')

%% plot splitting
visualize_splitting_v1(obsTable, metaTable, year_start, ...
    n_year_data, path_fig, name_var, suffix, '-r500')

disp('done')
