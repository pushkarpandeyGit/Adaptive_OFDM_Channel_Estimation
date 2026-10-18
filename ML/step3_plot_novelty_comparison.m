function step3_plot_novelty_comparison()
% =========================================================================
% STEP 3: PLOT NOVELTY PERFORMANCE (Proposed ML vs Park 2022 SHT)
%
% Generates a publication-quality 2-panel figure comparing the CIR length
% distribution of:
%   1. Ground Truth True CIR Length
%   2. Park 2022 Baseline (SHT with T=15)
%   3. Proposed Machine Learning Adaptive Estimator
% =========================================================================

fprintf('\n=======================================================\n');
fprintf(' NOVELTY DAY 1: PLOTTING ML VS PARK 2022 SHT\n');
fprintf('=======================================================\n');

%% 1. PATH SETUP & LOAD MODEL
thisFile = mfilename('fullpath');
if isempty(thisFile), here = pwd; else, here = fileparts(thisFile); end
projectRoot = fileparts(here);
addpath(genpath(projectRoot));

modelPath = fullfile(projectRoot, 'ML', 'trained_ml_cir_model.mat');
if ~exist(modelPath, 'file')
    error('trained_ml_cir_model.mat not found. Run step2_train_cir_predictor first.');
end
model = load(modelPath);

p = parameters();
SNRdB = 20;
numTrials = 500;

N = 256; M = 128; G = 64; K = 56; T_SHT = 15;
kLTS = [-100:2:-2 2:2:100].';
Np = length(kLTS); % 100 LTS pilot subcarriers

Pall = ieee80216_pall();
preamble = preamble_sequences(Pall);
[~, ltsTime] = preamble_to_ifft(preamble, p);
lts128 = get_known_lts(ltsTime);
X128 = fft(lts128, M);
ltsBins = mod(kLTS/2, M) + 1;
X_pilot = X128(ltsBins); X_pilot = X_pilot(:);

q = (-50:50).';
w = 0.5 + 0.5*cos(2*pi*q/M); w = w(:);
m_idx = (0:M-1).';
Gmat = (1/sqrt(M)) * exp(1j*2*pi*((m_idx - M/2)*q.')/M);

delays2 = [0, 0.4, 1.1] * 1e-6 / p.Ts;
powers2_lin = 10.^([0, -12, -15]/20);
K_lin = 10^(3/10);

delays5 = [0, 4.0, 10.0] * 1e-6 / p.Ts;
powers5_lin = 10.^([0, -5, -10]/20);

L_SHT_2 = zeros(numTrials, 1);
L_ML_2  = zeros(numTrials, 1);
L_SHT_5 = zeros(numTrials, 1);
L_ML_5  = zeros(numTrials, 1);

rng(42);
fprintf('Running %d trials at SNR = %d dB...\n', numTrials, SNRdB);

for t = 1:numTrials
    % SUI-2
    los = sqrt(K_lin/(K_lin+1));
    nlos = sqrt(1/(K_lin+1)) * (randn + 1j*randn)/sqrt(2);
    H2 = powers2_lin(1)*(los+nlos)*exp(-1j*2*pi*kLTS*delays2(1)/N) + ...
         powers2_lin(2)*(randn+1j*randn)/sqrt(2)*exp(-1j*2*pi*kLTS*delays2(2)/N) + ...
         powers2_lin(3)*(randn+1j*randn)/sqrt(2)*exp(-1j*2*pi*kLTS*delays2(3)/N);

    % SUI-5
    H5 = powers5_lin(1)*(randn+1j*randn)/sqrt(2)*exp(-1j*2*pi*kLTS*delays5(1)/N) + ...
         powers5_lin(2)*(randn+1j*randn)/sqrt(2)*exp(-1j*2*pi*kLTS*delays5(2)/N) + ...
         powers5_lin(3)*(randn+1j*randn)/sqrt(2)*exp(-1j*2*pi*kLTS*delays5(3)/N);

    % Standardized LTS Noise & Estimation (LTS pilot boost +6 dB)
    gammaP = 10^((SNRdB + 6)/10);
    noise2 = sqrt(1/(2*gammaP)) * (randn(Np, 1) + 1j*randn(Np, 1));
    HLS2 = H2 + noise2;

    noise5 = sqrt(1/(2*gammaP)) * (randn(Np, 1) + 1j*randn(Np, 1));
    HLS5 = H5 + noise5;

    % Standardized Time-Domain CIR via ls_to_cir
    [cir2,~,~] = ls_to_cir(HLS2, M, 'hann');
    sigma2_2 = mean(abs(cir2(p.noiseTapStart : p.noiseTapEnd)).^2);

    [cir5,~,~] = ls_to_cir(HLS5, M, 'hann');
    sigma2_5 = mean(abs(cir5(p.noiseTapStart : p.noiseTapEnd)).^2);

    % SHT
    L_SHT_2(t) = sht_cir_length(cir2, sigma2_2, T_SHT, G);
    L_SHT_5(t) = sht_cir_length(cir5, sigma2_5, T_SHT, G);

    % Proposed ML
    L_ML_2(t) = ml_cir_length(cir2, sigma2_2, G, SNRdB);
    L_ML_5(t) = ml_cir_length(cir5, sigma2_5, G, SNRdB);
end

%% 2. DISPLAY STATISTICS
fprintf('\n--- Performance Statistics at SNR = %d dB ---\n', SNRdB);
err_park2 = abs(mean(L_SHT_2) - 5);
err_ml2   = abs(mean(L_ML_2) - 5);
if err_ml2 < err_park2
    tag2 = '<-- PROPOSED ML CLOSER TO TRUE L';
else
    tag2 = '<-- Park closer';
end
fprintf('SUI-2 (True L = 5):\n');
fprintf('  Park 2022 SHT Mean: %.2f taps (Std: %.2f, Error: %.2f)\n', mean(L_SHT_2), std(L_SHT_2), err_park2);
fprintf('  Proposed ML   Mean: %.2f taps (Std: %.2f, Error: %.2f)  %s\n', mean(L_ML_2), std(L_ML_2), err_ml2, tag2);

err_park5 = abs(mean(L_SHT_5) - 41);
err_ml5   = abs(mean(L_ML_5) - 41);
if err_ml5 <= err_park5
    tag5 = '<-- PROPOSED ML EQUAL OR CLOSER TO TRUE L';
else
    tag5 = '<-- Park closer';
end
fprintf('\nSUI-5 (True L = 41):\n');
fprintf('  Park 2022 SHT Mean: %.2f taps (Std: %.2f, Error: %.2f)\n', mean(L_SHT_5), std(L_SHT_5), err_park5);
fprintf('  Proposed ML   Mean: %.2f taps (Std: %.2f, Error: %.2f)  %s\n', mean(L_ML_5), std(L_ML_5), err_ml5, tag5);

%% 3. PLOT PUBLICATION FIGURE
fig = figure('Name', 'Novelty Day 1 - Proposed ML vs Park 2022', ...
    'Color', 'w', 'Position', [100 100 1050 480]);

% Subplot 1: SUI-2
subplot(1, 2, 1);
hold on; grid on; box on;
bins2 = 1:15;
h_park2 = histcounts(L_SHT_2, [bins2, bins2(end)+1], 'Normalization', 'probability');
h_ml2   = histcounts(L_ML_2,  [bins2, bins2(end)+1], 'Normalization', 'probability');

plot(bins2, h_park2, '-sb', 'LineWidth', 1.8, 'MarkerFaceColor', 'b', 'MarkerSize', 6);
plot(bins2, h_ml2,   '-or', 'LineWidth', 2.0, 'MarkerFaceColor', 'r', 'MarkerSize', 6);
xline(5, '--k', 'True L = 5', 'LineWidth', 2, 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');

xlabel('Estimated CIR Length L_hat', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Probability / Relative Frequency', 'FontSize', 11, 'FontWeight', 'bold');
title('(a) SUI-2 Channel (SNR = 20 dB)', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Park 2022 SHT (T=15)', 'Proposed Novel ML Predictor', 'True Length (L=5)'}, ...
    'Location', 'northeast', 'FontSize', 10);
xlim([1 15]); ylim([0 1.05]);

% Subplot 2: SUI-5
subplot(1, 2, 2);
hold on; grid on; box on;
bins5 = 20:50;
h_park5 = histcounts(L_SHT_5, [bins5, bins5(end)+1], 'Normalization', 'probability');
h_ml5   = histcounts(L_ML_5,  [bins5, bins5(end)+1], 'Normalization', 'probability');

plot(bins5, h_park5, '-sb', 'LineWidth', 1.8, 'MarkerFaceColor', 'b', 'MarkerSize', 6);
plot(bins5, h_ml5,   '-or', 'LineWidth', 2.0, 'MarkerFaceColor', 'r', 'MarkerSize', 6);
xline(41, '--k', 'True L = 41', 'LineWidth', 2, 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');

xlabel('Estimated CIR Length L_hat', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Probability / Relative Frequency', 'FontSize', 11, 'FontWeight', 'bold');
title('(b) SUI-5 Channel (SNR = 20 dB)', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Park 2022 SHT (T=15)', 'Proposed Novel ML Predictor', 'True Length (L=41)'}, ...
    'Location', 'northwest', 'FontSize', 10);
xlim([25 50]); ylim([0 1.05]);

sgtitle('Novelty Evaluation: Proposed ML Estimator vs Park 2022 Baseline', ...
    'FontSize', 14, 'FontWeight', 'bold');

% Save Figure
saveas(fig, fullfile(projectRoot, 'ML', 'Novelty_Day1_ML_vs_Park.png'));
fprintf('\nSaved comparison plot to: %s\n', fullfile(projectRoot, 'ML', 'Novelty_Day1_ML_vs_Park.png'));
fprintf('=======================================================\n');

end
