function step7_plot_ml_diagnostics()
% =========================================================================
% STEP 7 (PLOT 5): ML NEURAL DIAGNOSTICS & EXPLAINABILITY SUITE
%
% Generates publication-grade diagnostic plots proving the rigor, explainability,
% and convergence of the Proposed ML Neural CIR Predictor:
%
%   (a) Parity Plot: Predicted L_hat vs True L_true (Proposed ML vs Park 2022)
%   (b) Neural Feature Importance Ranking (% Contribution across 11 features)
%   (c) Estimation MAE (taps) vs SNR (0 to 30 dB)
%   (d) Error Residual Distribution Histogram (Error = L_hat - L_true)
%
% Outputs:
%   ML/Novelty_Day3_Plot5_ML_Diagnostics.png
%   ML/Novelty_Day3_Plot5_ML_Diagnostics.mat
% =========================================================================

fprintf('\n============================================================\n');
fprintf(' NOVELTY PLOT 5: ML NEURAL DIAGNOSTICS & EXPLAINABILITY SUITE\n');
fprintf('============================================================\n');

%% 1. PATH SETUP & LOAD DATA
thisFile = mfilename('fullpath');
if isempty(thisFile), here = pwd; else, here = fileparts(thisFile); end
projectRoot = fileparts(here);
addpath(genpath(projectRoot));

datasetPath = fullfile(projectRoot, 'ML', 'cir_dataset.mat');
modelPath   = fullfile(projectRoot, 'ML', 'trained_ml_cir_model.mat');

if ~exist(datasetPath, 'file') || ~exist(modelPath, 'file')
    error('Required dataset or model files not found in ML folder.');
end

ds = load(datasetPath);
md = load(modelPath);

X = ds.featureMatrix;
y = ds.targetVector;
N_samples = size(X, 1);

%% 2. TEST SPLIT EVALUATION (Consistent 80/20 split)
rng(2026);
perm = randperm(N_samples);
splitIdx = round(0.80 * N_samples);
testIdx = perm(splitIdx+1:end);

X_test = X(testIdx, :);
y_test = y(testIdx);

% Forward Pass
X_test_norm = (X_test - md.mu) ./ md.sigma;
z1 = md.W1 * X_test_norm.' + md.b1;
a1 = max(0.01*z1, z1);
z2 = md.W2 * a1 + md.b2;
a2 = max(0.01*z2, z2);
y_pred_cont = (md.W3 * a2 + md.b3).';

% Enforce integer bounds and multipath protection
L_pred_raw = round(y_pred_cont);
L_park_test = X_test(:, 11); % Feature 11 is L_SHT
L_ml_test = max(1, min(64, max(L_pred_raw, L_park_test)));

% Residuals
err_ML   = L_ml_test - y_test;
err_Park = L_park_test - y_test;

% R-squared (Coefficient of Determination)
SS_tot = sum((y_test - mean(y_test)).^2);
SS_res_ML = sum((y_test - L_ml_test).^2);
SS_res_Park = sum((y_test - L_park_test).^2);
R2_ML = 1 - (SS_res_ML / SS_tot);
R2_Park = 1 - (SS_res_Park / SS_tot);

fprintf('Test Set Parity Correlation:\n');
fprintf('  Proposed ML  R^2 = %.4f | MAE = %.2f taps\n', R2_ML, mean(abs(err_ML)));
fprintf('  Park 2022    R^2 = %.4f | MAE = %.2f taps\n', R2_Park, mean(abs(err_Park)));

%% 3. NEURAL FEATURE IMPORTANCE
% Magnitude of weights connecting each input feature to hidden layer 1
featWeights = sum(abs(md.W1), 1); % 1 x 11
featImportance = (featWeights / sum(featWeights)) * 100; % Percentage

featLabels = { ...
    'SNR (dB)', ...
    'SNR (Lin)', ...
    'Noise Var (log\sigma^2)', ...
    'Total CIR Energy', ...
    'Peak Energy Ratio', ...
    'Peak Tap Index', ...
    '50% Energy Tap (P_{50})', ...
    '80% Energy Tap (P_{80})', ...
    '90% Energy Tap (P_{90})', ...
    '95% Energy Tap (P_{95})', ...
    'Park Baseline SHT (L_{SHT})' ...
};

%% 4. MAE VS SNR ANALYSIS
snrGrid = 0:5:30;
mae_snr_ML = zeros(size(snrGrid));
mae_snr_Park = zeros(size(snrGrid));

snr_col = X_test(:, 1);
for i = 1:numel(snrGrid)
    s = snrGrid(i);
    idx_s = abs(snr_col - s) < 1.0;
    if any(idx_s)
        mae_snr_ML(i) = mean(abs(L_ml_test(idx_s) - y_test(idx_s)));
        mae_snr_Park(i) = mean(abs(L_park_test(idx_s) - y_test(idx_s)));
    end
end

%% 5. PUBLICATION 4-PANEL FIGURE
fig = figure('Name', 'Novelty Plot 5 - ML Neural Diagnostics', ...
    'Color', 'w', 'Position', [80 50 1200 780], 'Visible', 'off');

% --- (a) Parity Scatter: Predicted vs True ---
subplot(2, 2, 1);
% Add slight jitter for visual clarity of discrete points
jitter_x = 0.25 * (rand(size(y_test)) - 0.5);
jitter_y = 0.25 * (rand(size(y_test)) - 0.5);

scatter(y_test + jitter_x, L_park_test + jitter_y, 20, [0.3 0.5 0.9], 'filled', ...
    'MarkerFaceAlpha', 0.4); hold on;
scatter(y_test + jitter_x, L_ml_test + jitter_y, 25, [0.85 0.1 0.1], 'filled', ...
    'MarkerFaceAlpha', 0.6);
plot([1 64], [1 64], '--k', 'LineWidth', 1.8); % Ideal Parity Line
grid on; box on;
xlim([0 50]); ylim([0 50]);
xlabel('True CIR Length L_{true} (taps)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Estimated Length L_{est} (taps)', 'FontSize', 10, 'FontWeight', 'bold');
title('(a) Parity Plot: Estimated vs True CIR Length', 'FontSize', 11, 'FontWeight', 'bold');
legend({sprintf('Park 2022 (R^2 = %.2f)', R2_Park), ...
        sprintf('PROPOSED ML (R^2 = %.2f)', R2_ML), ...
        'Ideal Parity (y = x)'}, 'Location', 'northwest', 'FontSize', 9);

% --- (b) Feature Importance Bar Chart ---
subplot(2, 2, 2);
[sortedImp, sortIdx] = sort(featImportance, 'ascend');
sortedLabels = featLabels(sortIdx);

b = barh(sortedImp, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', [0.1 0.3 0.5]);
yticks(1:11);
yticklabels(sortedLabels);
xlim([0 max(sortedImp)*1.25]);
xlabel('Relative Neural Importance (%)', 'FontSize', 10, 'FontWeight', 'bold');
title('(b) Neural Feature Importance Ranking', 'FontSize', 11, 'FontWeight', 'bold');
grid on; box on;
% Display percentage on bars
for k = 1:11
    text(sortedImp(k) + 0.3, k, sprintf('%.1f%%', sortedImp(k)), ...
        'VerticalAlignment', 'middle', 'FontSize', 8, 'FontWeight', 'bold');
end

% --- (c) MAE vs SNR ---
subplot(2, 2, 3);
plot(snrGrid, mae_snr_Park, '-sb', 'LineWidth', 1.8, 'MarkerFaceColor', 'b', 'MarkerSize', 6); hold on;
plot(snrGrid, mae_snr_ML, '-or', 'LineWidth', 2.2, 'MarkerFaceColor', 'r', 'MarkerSize', 6);
grid on; box on;
xlim([0 30]); ylim([0 max(mae_snr_Park)*1.25]);
xlabel('Signal-to-Noise Ratio (dB)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Mean Absolute Error (taps)', 'FontSize', 10, 'FontWeight', 'bold');
title('(c) CIR Estimation Error vs SNR', 'FontSize', 11, 'FontWeight', 'bold');
legend({'Park 2022 Baseline (SHT)', 'PROPOSED NOVEL ML'}, 'Location', 'northeast', 'FontSize', 9);

% --- (d) Error Distribution Histogram ---
subplot(2, 2, 4);
bins = -15:1:15;
h1 = histogram(err_Park, bins, 'FaceColor', 'b', 'FaceAlpha', 0.35, 'EdgeColor', 'b'); hold on;
h2 = histogram(err_ML, bins, 'FaceColor', 'r', 'FaceAlpha', 0.55, 'EdgeColor', 'r');
grid on; box on;
xlim([-12 12]);
xlabel('Estimation Error (L_{est} - L_{true}) [taps]', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Sample Count', 'FontSize', 10, 'FontWeight', 'bold');
title('(d) Error Residual Distribution (\mu \approx 0)', 'FontSize', 11, 'FontWeight', 'bold');
legend({sprintf('Park SHT (\\sigma = %.1f)', std(err_Park)), ...
        sprintf('Proposed ML (\\sigma = %.1f)', std(err_ML))}, ...
        'Location', 'northwest', 'FontSize', 9);

sgtitle('Proposed ML Channel Estimator: Diagnostics, Parity & Explainability Suite', ...
    'FontSize', 13, 'FontWeight', 'bold');

% Save Artifacts
outPng = fullfile(projectRoot, 'ML', 'Novelty_Day3_Plot5_ML_Diagnostics.png');
outMat = fullfile(projectRoot, 'ML', 'Novelty_Day3_Plot5_ML_Diagnostics.mat');

saveas(fig, outPng);
save(outMat, 'R2_ML', 'R2_Park', 'featImportance', 'featLabels', 'mae_snr_ML', 'mae_snr_Park', '-v7');

fprintf('\nSuccessfully saved Plot 5 results to:\n  %s\n  %s\n', outPng, outMat);
fprintf('============================================================\n');

end
