function plot_ssim_density_PCR_vs_LR(ssim_PCR_mat, ssim_LR_mat, ssim_ESM_mat, ...
    flag_season, name_var, i_mth_txt, suffix, path_fig, res_plot)

%% Plot SSIM distribution (KDE) for PCR vs LR vs ESM, including means and KL divergence

    % Define colorblind-friendly colors
    color_lr = [0.35, 0.7, 0.9]; % Light blue for LR
    color_pcr = [0.9, 0.6, 0];   % Orange for PCR
    color_esm = [0.5, 0.5, 0.5]; % Gray for ESM

    % Filter season and flatten
    ssim_pcr = ssim_PCR_mat(flag_season); 
    ssim_lr  = ssim_LR_mat(flag_season);  
    ssim_esm = ssim_ESM_mat(flag_season); 

    % Remove NaNs
    ssim_pcr = ssim_pcr(~isnan(ssim_pcr));
    ssim_lr  = ssim_lr(~isnan(ssim_lr));
    ssim_esm = ssim_esm(~isnan(ssim_esm));

    % Means
    mean_pcr = mean(ssim_pcr);
    mean_lr  = mean(ssim_lr);
    mean_esm = mean(ssim_esm);

    % Bins for histograms
    binstep = 0.01;
    binEdges = -1:binstep:1;

    % Histogram counts
    h_LR  = histcounts(ssim_lr, binEdges);
    h_PCR = histcounts(ssim_pcr, binEdges);
    h_ESM = histcounts(ssim_esm, binEdges);

    % Normalize to probability distributions
    prob_LR  = h_LR  / sum(h_LR);
    prob_PCR = h_PCR / sum(h_PCR);
    prob_ESM = h_ESM / sum(h_ESM);

    % Avoid zero division/log(0)
    prob_LR(prob_LR == 0)   = eps;
    prob_PCR(prob_PCR == 0) = eps;
    prob_ESM(prob_ESM == 0) = eps;

    % Compute symmetric KL divergence (Jeffreys divergence)
    kldiv_PCR_LR = sum(prob_PCR .* log(prob_PCR ./ prob_LR));
    kldiv_LR_PCR = sum(prob_LR .* log(prob_LR ./ prob_PCR));
    kldiv_PCR_ESM = sum(prob_PCR .* log(prob_PCR ./ prob_ESM));
    kldiv_LR_ESM  = sum(prob_LR  .* log(prob_LR  ./ prob_ESM));
    kldiv_ESM_PCR = sum(prob_ESM .* log(prob_ESM ./ prob_PCR));
    kldiv_ESM_LR  = sum(prob_ESM .* log(prob_ESM ./ prob_LR));

    % KS test
    [~, p_ks, ~] = kstest2(ssim_pcr, ssim_lr);

    % KDE
    [f_lr, x_lr]     = ksdensity(ssim_lr, binEdges);
    [f_pcr, x_pcr]   = ksdensity(ssim_pcr, binEdges);
    [f_esm, x_esm]   = ksdensity(ssim_esm, binEdges);

    % Plot
    f = figure('Position', [50 50 1500 500]);
    hold on;

    plot(x_lr, f_lr, '-', 'Color', color_lr, 'LineWidth', 2);
    plot(x_pcr, f_pcr, '-', 'Color', color_pcr, 'LineWidth', 2);
    plot(x_esm, f_esm, '-', 'Color', color_esm, 'LineWidth', 2);

    % Means
    plot([mean_lr, mean_lr], ylim, '--', 'Color', color_lr, 'LineWidth', 2);
    plot([mean_pcr, mean_pcr], ylim, '--', 'Color', color_pcr, 'LineWidth', 2);
    plot([mean_esm, mean_esm], ylim, '--', 'Color', color_esm, 'LineWidth', 2);

    % Labels and legend
    ylabel('Density', 'FontSize', 14);
    xlabel('SSIM', 'FontSize', 14);
    title([name_var ' SSIM KDEs (Testing set; PCR vs LR vs ESM; Month = ' i_mth_txt ')'], 'FontSize', 16);
    legend({'Linear Regression', 'PCR', 'ESM'}, 'Location', 'northeastoutside', 'FontSize', 12);

    xlim([-1 1]); grid on;

    % Divergence annotation
    annotation('textbox', [0.15, 0.82, 0.3, 0.1], 'String', ...
        sprintf(['KS test (PCR vs LR) p=%.4f\n' ...
        'Jeffreys div (PCR vs LR)=%.4f\n' ...
        'Jeffreys div (LR vs ESM)=%.4f\n' ...
        'Jeffreys div (PCR vs ESM)=%.4f'], ...
        p_ks, ...
        kldiv_PCR_LR + kldiv_LR_PCR, ...
        kldiv_ESM_LR + kldiv_LR_ESM, ...
        kldiv_PCR_ESM + kldiv_ESM_PCR), ...
        'FontSize', 10, 'EdgeColor', 'none');

    % Save figure
    outputFilename = [name_var '_SSIM_density_PCR_vs_LR_vs_ESM_mth=' i_mth_txt suffix '.tif'];
    print(f, fullfile(path_fig, outputFilename), '-dtiff', res_plot);

end
