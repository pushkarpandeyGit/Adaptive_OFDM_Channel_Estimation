function plot_day2_mse()
% =========================================================================
% PLOT_DAY2_MSE: Publication Log-Scale MSE vs SNR Plot
% =========================================================================

thisFile = mfilename('fullpath');
here = fileparts(thisFile);
matFile = fullfile(here, 'Novelty_Day2_MSE_results.mat');

if ~exist(matFile, 'file')
    error('Novelty_Day2_MSE_results.mat not found.');
end

data = load(matFile);
SNR_dB = data.SNR_dB;

fig = figure('Name', 'Novelty Day 2 - Log-Scale MSE Performance', ...
    'Color', 'w', 'Position', [100 100 1150 500]);

% --- Subplot 1: SUI-2 ---
subplot(1, 2, 1);
semilogy(SNR_dB, data.MSE_LMMSE_SUI2, ':k', 'LineWidth', 2.0); hold on;
semilogy(SNR_dB, data.MSE_FIXED_SUI2, '-xg', 'LineWidth', 1.8, 'MarkerSize', 8);
semilogy(SNR_dB, data.MSE_PARK_SUI2,  '-sb', 'LineWidth', 1.8, 'MarkerFaceColor', 'b', 'MarkerSize', 6);
semilogy(SNR_dB, data.MSE_ML_SUI2,    '-or', 'LineWidth', 2.2, 'MarkerFaceColor', 'r', 'MarkerSize', 6);
grid on; box on;
set(gca, 'YScale', 'log', 'YMinorGrid', 'on', 'FontSize', 10);
xlim([0 30]); ylim([1e-4 1]);
xlabel('SNR (dB)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Mean Squared Error (MSE)', 'FontSize', 11, 'FontWeight', 'bold');
title('(a) SUI-2 Channel Model', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Ideal LMMSE', 'Fixed SVD (L=64)', 'Park 2022 Adaptive SVD', 'PROPOSED NOVEL ML-SVD'}, ...
    'Location', 'southwest', 'FontSize', 9);

% --- Subplot 2: SUI-5 ---
subplot(1, 2, 2);
semilogy(SNR_dB, data.MSE_LMMSE_SUI5, ':k', 'LineWidth', 2.0); hold on;
semilogy(SNR_dB, data.MSE_FIXED_SUI5, '-xg', 'LineWidth', 1.8, 'MarkerSize', 8);
semilogy(SNR_dB, data.MSE_PARK_SUI5,  '-sb', 'LineWidth', 1.8, 'MarkerFaceColor', 'b', 'MarkerSize', 6);
semilogy(SNR_dB, data.MSE_ML_SUI5,    '-or', 'LineWidth', 2.2, 'MarkerFaceColor', 'r', 'MarkerSize', 6);
grid on; box on;
set(gca, 'YScale', 'log', 'YMinorGrid', 'on', 'FontSize', 10);
xlim([0 30]); ylim([5e-4 1]);
xlabel('SNR (dB)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Mean Squared Error (MSE)', 'FontSize', 11, 'FontWeight', 'bold');
title('(b) SUI-5 Channel Model', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Ideal LMMSE', 'Fixed SVD (L=64)', 'Park 2022 Adaptive SVD', 'PROPOSED NOVEL ML-SVD'}, ...
    'Location', 'southwest', 'FontSize', 9);

sgtitle('Proposed ML-Adaptive SVD vs State-of-the-Art Park 2022 Channel Estimation', ...
    'FontSize', 13, 'FontWeight', 'bold');

outPng = fullfile(here, 'Novelty_Day2_MSE_vs_SNR.png');
saveas(fig, outPng);
fprintf('Successfully updated log-scale plot: %s\n', outPng);

end
