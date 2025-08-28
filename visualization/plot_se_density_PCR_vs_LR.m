function plot_se_density_PCR_vs_LR(se_PCR_mat, se_LR_mat, se_ESM_mat, flag_season, ...
    name_var, i_mth_txt, suffix, path_fig, res_plot)

%% Plot histogram for RMSE for linear regression vs PCR vs ESM, including means

    % Define colorblind-friendly colors
    color_lr = [0.35, 0.7, 0.9]; % Light blue for LR
    color_pcr = [0.9, 0.6, 0];   % Orange for PCR
    color_esm = [0.5, 0.5, 0.5]; % Gray for ESM

    % RMSE histograms
    se_pcr = se_PCR_mat(:,flag_season); % PCR data
    se_lr  = se_LR_mat(:,flag_season); % LR data
    se_esm = se_ESM_mat(:,flag_season); % ESM data

    se_pcr = se_pcr(:);
    se_lr = se_lr(:);
    se_esm = se_esm(:);
    
    % Calculate means
    mean_pcr = mean(se_pcr,'omitnan');
    mean_lr  = mean(se_lr,'omitnan');
    mean_esm = mean(se_esm,'omitnan');

    % Define the number of bins and bin edges
    switch name_var
        case 'tas'
            binstep=0.01;
        case 'pr'
            binstep=0.01;
    end
    binEdges = (floor(min([se_pcr, se_lr, se_esm])/binstep)*binstep):binstep:(ceil(2*max([mean_pcr, mean_lr, mean_esm])/binstep)*binstep); % Common bin edges
    
    % calculate Kyllback-Leibler divergence from histogram distributions
    % with the given bins
    h_LR = histcounts(se_lr,binEdges);
    h_PCR= histcounts(se_pcr,binEdges);
    h_ESM= histcounts(se_esm,binEdges);


    % Perform Kolmogorov-Smirnov test
    [h, p_ks, ks_stat] = kstest2(se_pcr, se_lr);

    % Normalize the histograms to represent probability distributions
    probdist_LR     = h_LR  / sum(h_LR);
    probdist_PCR    = h_PCR / sum(h_PCR);
    probdist_ESM    = h_ESM / sum(h_ESM);
    
    % Avoid dividing by zero or taking log of zero
    probdist_LR(probdist_LR == 0) = eps;
    probdist_PCR(probdist_PCR == 0) = eps;
    probdist_ESM(probdist_ESM == 0) = eps;
    
    % Calculate (symmetric Kullback-Leibler divergence)
    kldiv_PCR_LR = sum(probdist_PCR .* log(probdist_PCR ./ probdist_LR));
    kldiv_LR_PCR = sum(probdist_LR .* log(probdist_LR ./ probdist_PCR));
    kldiv_PCR_ESM = sum(probdist_PCR .* log(probdist_PCR ./ probdist_ESM));
    kldiv_LR_ESM  = sum(probdist_LR .* log(probdist_LR ./ probdist_ESM));
    kldiv_ESM_PCR = sum(probdist_ESM .* log(probdist_ESM ./ probdist_PCR));
    kldiv_ESM_LR  = sum(probdist_ESM .* log(probdist_ESM ./ probdist_PCR));


    % Create a new figure
    f = figure('Position', [50 50 1500 500]);
    %axes1 = axes('Parent',f);
    hold on;

    % Plot histograms with common bins
    %histogram(se_lr, binEdges, 'Normalization', 'probability', 'FaceColor', color_lr, 'EdgeColor', 'none');
    %histogram(se_pcr, binEdges, 'Normalization', 'probability', 'FaceColor', color_pcr, 'EdgeColor', 'none');
    %histogram(se_esm, binEdges, 'Normalization', 'probability', 'FaceColor', color_esm, 'EdgeColor', 'none');

    % Plot smooth density curves using ksdensity
    [f_lr, x_lr] = ksdensity(se_lr,binEdges);
    [f_pcr, x_pcr] = ksdensity(se_pcr,binEdges);
    [f_esm, x_esm] = ksdensity(se_esm,binEdges);
    
    plot(x_lr, f_lr, '-', 'Color', color_lr, 'LineWidth', 2);   % LR density
    plot(x_pcr, f_pcr, '-', 'Color', color_pcr, 'LineWidth', 2); % PCR density
    plot(x_esm, f_esm, '-', 'Color', color_esm, 'LineWidth', 2); % ESM density

    % Dashed lines for the means
    plot([mean_lr, mean_lr], ylim, '--', 'Color', color_lr, 'LineWidth', 2); % Dashed line for LR mean
    plot([mean_pcr, mean_pcr], ylim, '--', 'Color', color_pcr, 'LineWidth', 2); % Dashed line for PCR mean
    plot([mean_esm, mean_esm], ylim, '--', 'Color', color_esm, 'LineWidth', 2); % Dashed line for ESM mean

    % Add labels, title, and legend
    ylabel('Density', 'FontSize', 14);
    xlabel('SE', 'FontSize', 14);
    title([name_var ' SE Histograms (Testing dataset over cal. stat. only; PCR vs LR vs ESM; Month = ' i_mth_txt ')'], 'FontSize', 16);
    legend({'Linear Regression', 'PCR', 'ESM'}, 'Location', 'northeastoutside', 'FontSize', 12);

    switch name_var
        case 'tas'
            xlim([0 max(binEdges)])
        case 'pr'
            xlim([0 max(binEdges)])
    end
    grid on;
    
    % Display Jeffreys divergence (symmetric Kullback-Leibler divergence)
    annotation('textbox', [0.15, 0.82, 0.3, 0.1], 'String', ...
        sprintf(['Two-Sample Kolmogorov-Smirnov (PCR vs LR) p=%.4f\n' ...
        'Jeffreys divergence (PCR vs LR)=%.4f\n' ...
        'Jeffreys divergence (LR vs ESM)=%.4f\n' ...
        'Jeffreys divergence (PCR vs ESM)=%.4f'], ...
        p_ks, ...
        kldiv_PCR_LR+kldiv_LR_PCR, ...
        kldiv_ESM_LR+kldiv_LR_ESM, ...
        kldiv_ESM_PCR+kldiv_PCR_ESM), ...
        'FontSize', 10, 'EdgeColor', 'none');


    % Save the histogram as a TIFF file with specified resolution
    outputFilename = [name_var '_SE_density_PCR_vs_LR_vs_ESM_mth=' i_mth_txt suffix '.tif'];
    print(f, fullfile(path_fig, outputFilename), '-dtiff', res_plot);

end
