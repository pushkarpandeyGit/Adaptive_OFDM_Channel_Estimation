function step4_evaluate_ml_adaptive_svd()
% =========================================================================
% STEP 4: EVALUATE PROPOSED ML-ADAPTIVE SVD CHANNEL ESTIMATOR
%
% Novelty Component (Day 2):
% Integrates the trained data-driven CIR length predictor into the Adaptive
% SVD filter bank and evaluates Channel Estimation MSE vs SNR.
%
% Benchmark Comparison (Exact Paper Settings as in Fig. 6):
%   1. Ideal LMMSE (Theoretical Lower Bound)
%   2. Fixed SVD (L = 64, Rank = 67)
%   3. Park 2022 Adaptive SVD (SHT, T = 15) [Paper State-of-the-Art]
%   4. PROPOSED NOVELTY: ML-Adaptive SVD (Learned L_hat)
%
% Output:
%   ML/Novelty_Day2_MSE_vs_SNR.png
%   ML/Novelty_Day2_MSE_results.mat
% =========================================================================

fprintf('\n=======================================================\n');
fprintf(' NOVELTY DAY 2: EVALUATING PROPOSED ML-ADAPTIVE SVD\n');
fprintf('=======================================================\n');

%% 1. PATH SETUP & LOAD MODEL
thisFile = mfilename('fullpath');
if isempty(thisFile), here = pwd; else, here = fileparts(thisFile); end
projectRoot = fileparts(here);
addpath(genpath(projectRoot));

modelPath = fullfile(projectRoot, 'ML', 'trained_ml_cir_model.mat');
if ~exist(modelPath, 'file')
    error('trained_ml_cir_model.mat not found. Please run step2_train_cir_predictor first.');
end

p = parameters();
rng(2026);

%% 2. PARAMETERS (EXACT FIG 6 BENCHMARK)
N  = p.N;        % 256
Nu = p.Nu;       % 200
G  = p.G;        % 64
M  = 128;        % IDFT size
Np = 100;        % LTS pilot tones
T_SHT = 15;      % Paper baseline threshold
fixedRank = G + p.fixedSVDOffset; % 67

SNR_dB = 0:5:30;
numSNR = length(SNR_dB);
numTrials = 200; % Monte Carlo trials per SNR

ltsLocations  = [-100:2:-2 2:2:100];
usedLocations = p.activeIndex;

noiseStart = p.noiseTapStart; % 65
noiseEnd   = p.noiseTapEnd;   % 120

%% 3. TRUE COVARIANCE MATRICES & FIXED SVD OPERATORS
Rtrue_SUI2 = local_uniform_covariance(N, 5);
Rhh_SUI2  = local_extract_covariance(Rtrue_SUI2, usedLocations, usedLocations, N);
Rhp_SUI2  = local_extract_covariance(Rtrue_SUI2, usedLocations, ltsLocations, N);
Rpp_SUI2  = local_extract_covariance(Rtrue_SUI2, ltsLocations, ltsLocations, N);

Rtrue_SUI5 = local_uniform_covariance(N, 40);
Rhh_SUI5  = local_extract_covariance(Rtrue_SUI5, usedLocations, usedLocations, N);
Rhp_SUI5  = local_extract_covariance(Rtrue_SUI5, usedLocations, ltsLocations, N);
Rpp_SUI5  = local_extract_covariance(Rtrue_SUI5, ltsLocations, ltsLocations, N);

% Fixed SVD operators (L = 64)
Rmodel_fixed = local_uniform_covariance(N, G);
Rhp_fixed = local_extract_covariance(Rmodel_fixed, usedLocations, ltsLocations, N);
Rpp_fixed = local_extract_covariance(Rmodel_fixed, ltsLocations, ltsLocations, N);
[U_fixed, D_fixed] = eig(Rpp_fixed);
[eigs_fixed, idx_f] = sort(real(diag(D_fixed)), 'descend');
U_fixed = U_fixed(:, idx_f);

%% 4. STORAGE
MSE_LMMSE_SUI2 = zeros(1, numSNR);
MSE_FIXED_SUI2 = zeros(1, numSNR);
MSE_PARK_SUI2  = zeros(1, numSNR);
MSE_ML_SUI2    = zeros(1, numSNR);

MSE_LMMSE_SUI5 = zeros(1, numSNR);
MSE_FIXED_SUI5 = zeros(1, numSNR);
MSE_PARK_SUI5  = zeros(1, numSNR);
MSE_ML_SUI5    = zeros(1, numSNR);

fprintf('Running Monte Carlo simulations (%d trials/SNR)...\n', numTrials);

%% 5. SIMULATION LOOP ACROSS SNRS
for iSNR = 1:numSNR
    snr_val = SNR_dB(iSNR);
    gamma_p = 10^(snr_val/10);
    fprintf('  SNR = %2d dB... ', snr_val);

    % --- 1. Analytical LMMSE Benchmark ---
    K_LMMSE_2 = Rhp_SUI2 / (Rpp_SUI2 + (1/gamma_p)*eye(Np));
    MSE_LMMSE_SUI2(iSNR) = local_mse_eq31(Rhh_SUI2, Rhp_SUI2, Rpp_SUI2, K_LMMSE_2, gamma_p, Nu);

    K_LMMSE_5 = Rhp_SUI5 / (Rpp_SUI5 + (1/gamma_p)*eye(Np));
    MSE_LMMSE_SUI5(iSNR) = local_mse_eq31(Rhh_SUI5, Rhp_SUI5, Rpp_SUI5, K_LMMSE_5, gamma_p, Nu);

    % --- 2. Fixed SVD Benchmark ---
    Dp_fixed = local_build_Dp(eigs_fixed, fixedRank, gamma_p);
    K_fixed = Rhp_fixed * U_fixed * Dp_fixed * U_fixed';
    MSE_FIXED_SUI2(iSNR) = local_mse_eq31(Rhh_SUI2, Rhp_SUI2, Rpp_SUI2, K_fixed, gamma_p, Nu);
    MSE_FIXED_SUI5(iSNR) = local_mse_calc(Rhh_SUI5, Rhp_SUI5, Rpp_SUI5, K_fixed, gamma_p, Nu);

    % --- 3. Monte Carlo Loops for SUI-2 ---
    mse_park2 = zeros(numTrials, 1);
    mse_ml2   = zeros(numTrials, 1);
    
    for t = 1:numTrials
        h2 = sui2_channel(p);
        Htrue2 = local_channel_frequency_response(h2, N);
        Htrue_used2 = Htrue2(usedLocations + N/2 + 1);
        Htrue_lts2  = Htrue2(ltsLocations + N/2 + 1);

        lsNoise2 = sqrt(1/(2*gamma_p)) * (randn(Np,1) + 1j*randn(Np,1));
        HLS2 = Htrue_lts2 + lsNoise2;

        [cir2,~,~] = ls_to_cir(HLS2, M, 'hann');
        sigma2_2 = mean(abs(cir2(noiseStart:noiseEnd)).^2);

        % Baseline SHT
        L_park2 = sht_cir_length(cir2, sigma2_2, T_SHT, G);
        L_park2 = max(1, L_park2);

        % Proposed ML
        L_ml2 = ml_cir_length(cir2, sigma2_2, G, snr_val);
        L_ml2 = max(1, L_ml2);

        % SVD Park SUI-2
        rank_p2 = min(Np, L_park2 + p.adaptiveSVDOffset);
        R_p2 = local_uniform_covariance(N, L_park2);
        Rhp_p2 = local_extract_covariance(R_p2, usedLocations, ltsLocations, N);
        Rpp_p2 = local_extract_covariance(R_p2, ltsLocations, ltsLocations, N);
        [Up2, Dp2_mat] = eig(Rpp_p2);
        [eigs_p2, idx_p2] = sort(real(diag(Dp2_mat)), 'descend');
        Up2 = Up2(:, idx_p2);
        Dp_p2 = local_build_Dp(eigs_p2, rank_p2, gamma_p);
        Hhat_park2 = Rhp_p2 * Up2 * Dp_p2 * Up2' * HLS2;
        mse_park2(t) = mean(abs(Hhat_park2 - Htrue_used2).^2);

        % SVD Proposed ML SUI-2
        rank_m2 = min(Np, L_ml2 + p.adaptiveSVDOffset);
        R_m2 = local_uniform_covariance(N, L_ml2);
        Rhp_m2 = local_extract_covariance(R_m2, usedLocations, ltsLocations, N);
        Rpp_m2 = local_extract_covariance(R_m2, ltsLocations, ltsLocations, N);
        [Um2, Dm2_mat] = eig(Rpp_m2);
        [eigs_m2, idx_m2] = sort(real(diag(Dm2_mat)), 'descend');
        Um2 = Um2(:, idx_m2);
        Dp_m2 = local_build_Dp(eigs_m2, rank_m2, gamma_p);
        Hhat_ml2 = Rhp_m2 * Um2 * Dp_m2 * Um2' * HLS2;
        mse_ml2(t) = mean(abs(Hhat_ml2 - Htrue_used2).^2);
    end

    MSE_PARK_SUI2(iSNR) = mean(mse_park2);
    MSE_ML_SUI2(iSNR)   = mean(mse_ml2);

    % --- 4. Monte Carlo Loops for SUI-5 ---
    mse_park5 = zeros(numTrials, 1);
    mse_ml5   = zeros(numTrials, 1);

    for t = 1:numTrials
        h5 = sui5_channel(p);
        Htrue5 = local_channel_frequency_response(h5, N);
        Htrue_used5 = Htrue5(usedLocations + N/2 + 1);
        Htrue_lts5  = Htrue5(ltsLocations + N/2 + 1);

        lsNoise5 = sqrt(1/(2*gamma_p)) * (randn(Np,1) + 1j*randn(Np,1));
        HLS5 = Htrue_lts5 + lsNoise5;

        [cir5,~,~] = ls_to_cir(HLS5, M, 'hann');
        sigma2_5 = mean(abs(cir5(noiseStart:noiseEnd)).^2);

        % Baseline SHT
        L_park5 = sht_cir_length(cir5, sigma2_5, T_SHT, G);
        L_park5 = max(1, L_park5);

        % Proposed ML
        L_ml5 = ml_cir_length(cir5, sigma2_5, G, snr_val);
        L_ml5 = max(1, L_ml5);

        % SVD Park SUI-5
        rank_p5 = min(Np, L_park5 + p.adaptiveSVDOffset);
        R_p5 = local_uniform_covariance(N, L_park5);
        Rhp_p5 = local_extract_covariance(R_p5, usedLocations, ltsLocations, N);
        Rpp_p5 = local_extract_covariance(R_p5, ltsLocations, ltsLocations, N);
        [Up5, Dp5_mat] = eig(Rpp_p5);
        [eigs_p5, idx_p5] = sort(real(diag(Dp5_mat)), 'descend');
        Up5 = Up5(:, idx_p5);
        Dp_p5 = local_build_Dp(eigs_p5, rank_p5, gamma_p);
        Hhat_park5 = Rhp_p5 * Up5 * Dp_p5 * Up5' * HLS5;
        mse_park5(t) = mean(abs(Hhat_park5 - Htrue_used5).^2);

        % SVD Proposed ML SUI-5
        rank_m5 = min(Np, L_ml5 + p.adaptiveSVDOffset);
        R_m5 = local_uniform_covariance(N, L_ml5);
        Rhp_m5 = local_extract_covariance(R_m5, usedLocations, ltsLocations, N);
        Rpp_m5 = local_extract_covariance(R_m5, ltsLocations, ltsLocations, N);
        [Um5, Dm5_mat] = eig(Rpp_m5);
        [eigs_m5, idx_m5] = sort(real(diag(Dm5_mat)), 'descend');
        Um5 = Um5(:, idx_m5);
        Dp_m5 = local_build_Dp(eigs_m5, rank_m5, gamma_p);
        Hhat_ml5 = Rhp_m5 * Um5 * Dp_m5 * Um5' * HLS5;
        mse_ml5(t) = mean(abs(Hhat_ml5 - Htrue_used5).^2);
    end

    MSE_PARK_SUI5(iSNR) = mean(mse_park5);
    MSE_ML_SUI5(iSNR)   = mean(mse_ml5);

    fprintf('Done.\n');
end

%% 6. PRINT RESULTS TABLE
fprintf('\n====================================================================================================\n');
fprintf(' CHANNEL ESTIMATION MSE PERFORMANCE COMPARISON (Novelty Day 2)\n');
fprintf('====================================================================================================\n');
fprintf(' SNR(dB) | SUI-2 Fixed SVD | SUI-2 Park SVD | SUI-2 Proposed ML | SUI-5 Park SVD | SUI-5 Proposed ML\n');
fprintf('----------------------------------------------------------------------------------------------------\n');
for i = 1:numSNR
    fprintf('  %2d dB  |    %10.4e   |   %10.4e   |   %10.4e     |   %10.4e   |   %10.4e\n', ...
        SNR_dB(i), MSE_FIXED_SUI2(i), MSE_PARK_SUI2(i), MSE_ML_SUI2(i), MSE_PARK_SUI5(i), MSE_ML_SUI5(i));
end
fprintf('====================================================================================================\n');

%% 7. PLOT PUBLICATION-QUALITY COMPARISON FIGURE
fig = figure('Name', 'Novelty Day 2 - MSE Performance', ...
    'Color', 'w', 'Position', [100 100 1100 500]);

% Subplot 1: SUI-2
subplot(1, 2, 1);
hold on; grid on; box on;
semilogy(SNR_dB, MSE_LMMSE_SUI2, ':k', 'LineWidth', 2.0);
semilogy(SNR_dB, MSE_FIXED_SUI2, '-xg', 'LineWidth', 1.8, 'MarkerSize', 7);
semilogy(SNR_dB, MSE_PARK_SUI2,  '-sb', 'LineWidth', 1.8, 'MarkerFaceColor', 'b', 'MarkerSize', 6);
semilogy(SNR_dB, MSE_ML_SUI2,    '-or', 'LineWidth', 2.2, 'MarkerFaceColor', 'r', 'MarkerSize', 6);

xlabel('SNR (dB)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Mean Squared Error (MSE)', 'FontSize', 11, 'FontWeight', 'bold');
title('(a) SUI-2 Channel Model', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Ideal LMMSE', 'Fixed SVD (L=64)', 'Park 2022 Adaptive SVD', 'PROPOSED NOVEL ML-SVD'}, ...
    'Location', 'southwest', 'FontSize', 9);
xlim([0 30]);

% Subplot 2: SUI-5
subplot(1, 2, 2);
hold on; grid on; box on;
semilogy(SNR_dB, MSE_LMMSE_SUI5, ':k', 'LineWidth', 2.0);
semilogy(SNR_dB, MSE_FIXED_SUI5, '-xg', 'LineWidth', 1.8, 'MarkerSize', 7);
semilogy(SNR_dB, MSE_PARK_SUI5,  '-sb', 'LineWidth', 1.8, 'MarkerFaceColor', 'b', 'MarkerSize', 6);
semilogy(SNR_dB, MSE_ML_SUI5,    '-or', 'LineWidth', 2.2, 'MarkerFaceColor', 'r', 'MarkerSize', 6);

xlabel('SNR (dB)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Mean Squared Error (MSE)', 'FontSize', 11, 'FontWeight', 'bold');
title('(b) SUI-5 Channel Model', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Ideal LMMSE', 'Fixed SVD (L=64)', 'Park 2022 Adaptive SVD', 'PROPOSED NOVEL ML-SVD'}, ...
    'Location', 'southwest', 'FontSize', 9);
xlim([0 30]);

sgtitle('Proposed ML-Adaptive SVD vs State-of-the-Art Park 2022 Channel Estimation', ...
    'FontSize', 13, 'FontWeight', 'bold');

% Save Figure & Results
outPng = fullfile(projectRoot, 'ML', 'Novelty_Day2_MSE_vs_SNR.png');
outMat = fullfile(projectRoot, 'ML', 'Novelty_Day2_MSE_results.mat');

saveas(fig, outPng);
save(outMat, 'SNR_dB', 'MSE_LMMSE_SUI2', 'MSE_FIXED_SUI2', 'MSE_PARK_SUI2', 'MSE_ML_SUI2', ...
    'MSE_LMMSE_SUI5', 'MSE_FIXED_SUI5', 'MSE_PARK_SUI5', 'MSE_ML_SUI5', '-v7');

fprintf('\nSuccessfully saved MSE results to:\n  %s\n  %s\n', outPng, outMat);
fprintf('=======================================================\n');

end

%% --- EXACT HELPER FUNCTIONS FROM PAPER FIG 6 ---
function R = local_uniform_covariance(N, L)
    k = (0:N-1).';
    d = k - k.';
    R = zeros(N, N);
    same = abs(d) < eps;
    R(same) = 1;
    idx = ~same;
    z = 1j * 2*pi*L .* d(idx) / N;
    R(idx) = (1 - exp(-z)) ./ z;
    R = (R + R')/2;
end

function Rsub = local_extract_covariance(R, rowLocations, colLocations, N)
    rowIndex = rowLocations + N/2 + 1;
    colIndex = colLocations + N/2 + 1;
    Rsub = R(rowIndex, colIndex);
end

function Dp = local_build_Dp(eigenvalues, rankValue, gamma_p)
    Np = length(eigenvalues);
    Dp = zeros(Np, Np);
    beta = 1;
    rankValue = max(1, min(rankValue, Np));
    for k = 1:rankValue
        denominator = eigenvalues(k) + beta/gamma_p;
        Dp(k,k) = 1/denominator;
    end
end

function mse = local_mse_eq31(Rhh_true, Rhp_true, Rpp_true, K, gamma_p, Nu)
    Np = size(Rpp_true, 1);
    noiseCov = Rpp_true + (1/gamma_p)*eye(Np);
    E = Rhh_true - Rhp_true*K' - K*Rhp_true' + K*noiseCov*K';
    mse = real(trace(E))/Nu;
end

function mse = local_mse_calc(Rhh, Rhp, Rpp, K, gamma_p, Nu)
    mse = local_mse_eq31(Rhh, Rhp, Rpp, K, gamma_p, Nu);
end

function H = local_channel_frequency_response(h, N)
    h = h(:);
    hPadded = [h; zeros(N-length(h), 1)];
    H = fftshift(fft(hPadded, N));
    H = H(:);
end

function cir = local_ls_to_cir(HLS, M)
    HLS = HLS(:).';
    H101 = zeros(1, 101);
    H101(1:50) = HLS(1:50);
    H101(52:101) = HLS(51:100);
    n = -50:50;
    w = 0.5 + 0.5*cos(2*pi*n/M);
    Hwindowed = H101 .* w;
    m = (0:M-1).';
    q = -50:50;
    Gmat = (1/sqrt(M)) * exp(1j*2*pi/M * ((m-M/2)*q));
    cir = fftshift(Gmat * Hwindowed(:));
end

function Lhat = local_sht(cir, sigma2, T, G)
    threshold = T * sigma2;
    tapEnergy = abs(cir(1:G)).^2;
    Lhat = 0;
    for n = 1:G
        if tapEnergy(n) > threshold
            Lhat = n;
        end
    end
end
