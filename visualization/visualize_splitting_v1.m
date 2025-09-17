function visualize_splitting_v1(obsTable, metaTable, year_start, n_year_data, path_fig, name_var, suffix, res)
    % VISUALIZE_SPLITTING: Visualize the data splitting for calibration, 
    % model selection, and testing over time and across stations.
    %
    % Parameters:
    % - obsTable: Table with observation data, including flags for calibration, validation, and testing
    % - metaTable: Table with metadata for stations
    % - year_start: Starting year for the dataset
    % - n_year_data: Number of years in the dataset
    % - path_fig: Path to save the output figure
    % - name_var: Variable name for the output figure
    % - suffix: Suffix for the figure filename
    % - res_plot: Resolution for saving the figure (e.g., 300 for 300 dpi)

    % Initialize split matrix with (#timesteps) x (#stations)
    split_mat = zeros(n_year_data * 12, height(metaTable));
    for i_ID = 1:height(metaTable)
        flag_ID = ismember(obsTable.ID, metaTable.ID(i_ID));
        obsTable_tmp = obsTable(flag_ID, :);
        time_ind = (obsTable_tmp.month_since_0CE - year_start * 12);
        
        % Assign codes to the split matrix
        split_mat(time_ind(obsTable_tmp.flag_cal), i_ID) = 1;
        split_mat(time_ind(obsTable_tmp.flag_val), i_ID) = 2;
        split_mat(time_ind(obsTable_tmp.flag_test), i_ID) = 3;
    end

    % Calculate percentages for calibration, model selection, and testing
    total_cells = numel(split_mat); % Total number of cells
    no_data_count = sum(split_mat(:) == 0);
    cal_count = sum(split_mat(:) == 1);
    val_count = sum(split_mat(:) == 2);
    test_count = sum(split_mat(:) == 3);
    valid_cells = total_cells - no_data_count; % Exclude "No data"

    cal_percentage = (cal_count / valid_cells) * 100;
    val_percentage = (val_count / valid_cells) * 100;
    test_percentage = (test_count / valid_cells) * 100;

    % Define the color palette
    palette = [0.9, 0.9, 0.9;         % No data
               240/256, 165/256, 5/256; % Calibration
               0, 105/256, 162/256;     % Model selection
               0, 193/256, 141/256];    % Testing

    % Create the figure
    f = figure('Position', [50, 50, 1400, 500]);

    % Assign colors to the matrix based on the palette
    colormap(palette);
    surf(split_mat, 'LineStyle', 'none');
    view(0, 90); % Overhead view

    % Adjust axes limits to fit the data perfectly
    axis tight; % Ensures the axes fit tightly to the data

    % Increase tick label font size
    ax = gca; % get current axes
    ax.FontSize = 14; % set tick label font size

    % Set axis labels with larger font size
    xlabel('Station', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('Year', 'FontSize', 16, 'FontWeight', 'bold');

    % Determine total years and define tick step
    years = year_start:(year_start + size(split_mat, 1)/12 - 1);
    y_tick_step = 10; % Tick every 10 years
    
    % Calculate corresponding indices for y-ticks
    ytick_indices = (1:12*y_tick_step:size(split_mat, 1)); % Every 10 years
    ytick_years = years(1:y_tick_step:end); % Extract years corresponding to tick indices
    
    % Apply ticks and labels
    yticks(ytick_indices); % Set y-ticks
    yticklabels(cellstr(num2str(ytick_years'))); % Set y-tick labels

    % Define the label text with percentages
    label_text = {
        sprintf('No data'), ...
        sprintf('Calibration (%.1f%%)', cal_percentage), ...
        sprintf('Model selection (%.1f%%)', val_percentage), ...
        sprintf('Testing (%.1f%%)', test_percentage)};

    % Create a custom legend with specified labels and colors
    hold on;
    legend_patches = zeros(1, numel(label_text));
    for i = 1:numel(label_text)
        % Create an invisible patch for each label with the corresponding color
        legend_patches(i) = patch([0 1 1 0], [0 0 1 1], palette(i, :), ...
                                  'EdgeColor', 'none', ...
                                  'FaceColor', palette(i, :), ...
                                  'DisplayName', label_text{i});
    end
    hold off;

    % Add the legend
    legend(legend_patches, label_text, 'Location', 'northeastoutside');

    % Generate the output filename
    name_figure_save = fullfile(path_fig, [name_var '_split_dataset_PCR' suffix]);

    % Save the figure
    print(f, name_figure_save, '-dtiff', res);

    %% mikado plot with only first 100 stations     
        % Initialize split matrix with (#timesteps) x (#stations)
    split_mat = zeros(n_year_data * 12, height(metaTable));
    for i_ID = 1:height(metaTable)
        flag_ID = ismember(obsTable.ID, metaTable.ID(i_ID));
        obsTable_tmp = obsTable(flag_ID, :);
        time_ind = (obsTable_tmp.month_since_0CE - year_start * 12);
        
        % Assign codes to the split matrix
        split_mat(time_ind(obsTable_tmp.flag_cal), i_ID) = 1;
        split_mat(time_ind(obsTable_tmp.flag_val), i_ID) = 2;
        split_mat(time_ind(obsTable_tmp.flag_test), i_ID) = 3;
    end

    % Calculate percentages for calibration, model selection, and testing
    total_cells = numel(split_mat); % Total number of cells
    no_data_count = sum(split_mat(:) == 0);
    cal_count = sum(split_mat(:) == 1);
    val_count = sum(split_mat(:) == 2);
    test_count = sum(split_mat(:) == 3);
    valid_cells = total_cells - no_data_count; % Exclude "No data"

    cal_percentage = (cal_count / valid_cells) * 100;
    val_percentage = (val_count / valid_cells) * 100;
    test_percentage = (test_count / valid_cells) * 100;

    % Define the color palette
    palette = [0.9, 0.9, 0.9;         % No data
               240/256, 165/256, 5/256; % Calibration
               0, 105/256, 162/256;     % Model selection
               0, 193/256, 141/256];    % Testing

    % Create the figure
    f = figure('Position', [50, 50, 1400, 500]);

    % Assign colors to the matrix based on the palette
    colormap(palette);
    surf(split_mat(:,1:100), 'LineStyle', 'none');
    view(0, 90); % Overhead view

    % Adjust axes limits to fit the data perfectly
    axis tight; % Ensures the axes fit tightly to the data

    % Increase tick label font size
    ax = gca; % get current axes
    ax.FontSize = 14; % set tick label font size

    % Set axis labels with larger font size
    xlabel('Station', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('Year', 'FontSize', 16, 'FontWeight', 'bold');

    % Determine total years and define tick step
    years = year_start:(year_start + size(split_mat, 1)/12 - 1);
    y_tick_step = 10; % Tick every 10 years
    
    % Calculate corresponding indices for y-ticks
    ytick_indices = (1:12*y_tick_step:size(split_mat, 1)); % Every 10 years
    ytick_years = years(1:y_tick_step:end); % Extract years corresponding to tick indices
    
    % Apply ticks and labels
    yticks(ytick_indices); % Set y-ticks
    yticklabels(cellstr(num2str(ytick_years'))); % Set y-tick labels

    % Define the label text with percentages
    label_text = {
        sprintf('No data'), ...
        sprintf('Calibration (%.1f%%)', cal_percentage), ...
        sprintf('Model selection (%.1f%%)', val_percentage), ...
        sprintf('Testing (%.1f%%)', test_percentage)};

    % Create a custom legend with specified labels and colors
    hold on;
    legend_patches = zeros(1, numel(label_text));
    for i = 1:numel(label_text)
        % Create an invisible patch for each label with the corresponding color
        legend_patches(i) = patch([0 1 1 0], [0 0 1 1], palette(i, :), ...
                                  'EdgeColor', 'none', ...
                                  'FaceColor', palette(i, :), ...
                                  'DisplayName', label_text{i});
    end
    hold off;

    % Add the legend
    legend(legend_patches, label_text, 'Location', 'northeastoutside');

    % Generate the output filename
    name_figure_save = fullfile(path_fig, [name_var '_split_dataset_100_PCR' suffix '.png']);

    % Save the figure
    print(f, name_figure_save, '-dpng', res);
end
