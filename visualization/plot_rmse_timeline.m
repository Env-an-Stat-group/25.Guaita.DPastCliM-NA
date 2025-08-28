function plot_rmse_timeline(rmse_timeline_PCR, rmse_timeline_LR, rmse_timeline_ESM, flag_season, time_array,...
            name_var, i_mth_txt, suffix, path_fig, res_plot,year_start_hist,year_end)

    % Define colorblind-friendly colors
    color_pcr = [0.0, 0.45, 0.70]; % Blue for PCR
    color_lr = [0.30, 0.74, 0.93]; % Light Blue for LR
    color_esm = [0.93, 0.69, 0.13]; % Gold for ESM

    rmse_timeline_PCR(not(flag_season))=nan;
    rmse_timeline_LR(not(flag_season))=nan;
    rmse_timeline_ESM(not(flag_season))=nan;

    rmse_timeline_PCR = mean(reshape(rmse_timeline_PCR,12,[]),1,'omitnan');
    rmse_timeline_LR = mean(reshape(rmse_timeline_LR,12,[]),1,'omitnan');
    rmse_timeline_ESM = mean(reshape(rmse_timeline_ESM,12,[]),1,'omitnan');
    
    time_array = reshape(time_array,12,[]);
    time_array = time_array(12,:)/12;

    % Perform linear regression for RMSE (PCR)
    X_pcr = [ones(length(time_array), 1), time_array'];
    [b_pcr, ~, ~, ~, stats_pcr] = regress(rmse_timeline_PCR(:), X_pcr);
    slope_pcr = b_pcr(2);
    se_pcr = std(rmse_timeline_PCR' - (X_pcr * b_pcr),'omitnan') / sqrt(sum(~isnan(rmse_timeline_PCR)));

    % Perform linear regression for RMSE (LR)
    X_lr = [ones(length(time_array), 1), time_array'];
    [b_lr, ~, ~, ~, stats_lr] = regress(rmse_timeline_LR(:), X_lr);
    slope_lr = b_lr(2);
    se_lr = std(rmse_timeline_LR' - (X_lr * b_lr),'omitnan') / sqrt(sum(~isnan(rmse_timeline_LR)));

    % Perform linear regression for RMSE (ESM)
    X_esm = [ones(length(time_array), 1), time_array'];
    [b_esm, ~, ~, ~, stats_esm] = regress(rmse_timeline_ESM(:), X_esm);
    slope_esm = b_esm(2);
    se_esm = std(rmse_timeline_ESM' - (X_esm * b_esm),'omitnan') / sqrt(sum(~isnan(rmse_timeline_ESM)));

    % Compare slopes (PCR vs ESM and LR vs ESM)
    z_pcr_vs_esm = (slope_pcr - slope_esm) / sqrt(se_pcr^2 + se_esm^2);
    pval_pcr_vs_esm = 2 * (1 - normcdf(abs(z_pcr_vs_esm))); % Two-tailed test

    z_lr_vs_esm = (slope_lr - slope_esm) / sqrt(se_lr^2 + se_esm^2);
    pval_lr_vs_esm = 2 * (1 - normcdf(abs(z_lr_vs_esm))); % Two-tailed test

    % Create a new figure
    f = figure('Position', [50, 0, 1500, 400]);
    hold on;

    % Plot RMSE for the ESM
    plot(time_array, rmse_timeline_ESM, 'LineWidth', 1.5, 'DisplayName', 'ESM', 'Color', color_esm, 'LineStyle', '-'); 

    % Plot RMSE for the downscaled PCR
    scatter(time_array, rmse_timeline_PCR, 20, color_pcr, 'filled', 'o', 'DisplayName', 'PCR'); 

    % Plot RMSE for the downscaled LR
    scatter(time_array, rmse_timeline_LR, 20, color_lr, 'filled', '^' , 'DisplayName', 'LR'); 

    % Plot linear regression fits
    plot(time_array, X_pcr * b_pcr, 'LineWidth', 1.2, 'LineStyle', '-', 'Color', 'k', 'DisplayName', sprintf('PCR'));
    plot(time_array, X_lr * b_lr, 'LineWidth', 1.2, 'LineStyle', '--', 'Color', 'k', 'DisplayName', sprintf('LR'));
    plot(time_array, X_esm * b_esm, 'LineWidth', 1.2, 'LineStyle', ':', 'Color', 'k', 'DisplayName', sprintf('ESM'));

    % Title and labels
    title(['RMSE Timeline of ' name_var ' (' i_mth_txt ')'], 'FontSize', 16, 'FontWeight', 'bold');
    xlabel('Time (Years)', 'FontSize', 14);
    ylabel('Mean RMSE', 'FontSize', 14);

    % Add legend
    legend('Location', 'northeastoutside', 'FontSize', 12);

    % Adjust grid and axis properties
    grid on;
    set(gca, 'FontSize', 12);
    xlim([year_start_hist year_end]);
    xticks(year_start_hist:5:year_end);
    xtickangle(45);

    % Display p-values of slope comparisons
    annotation('textbox', [0.15, 0.82, 0.3, 0.1], 'String', ...
        sprintf('PCR vs ESM p=%.4f\nLR vs ESM p=%.4f', pval_pcr_vs_esm, pval_lr_vs_esm), ...
        'FontSize', 10, 'EdgeColor', 'none');

    % Save the figure
    outputFilename = fullfile(path_fig, [name_var '_RMSE_timeline_' i_mth_txt suffix '.tif']);
    print(f, outputFilename, '-dtiff', res_plot);

end
