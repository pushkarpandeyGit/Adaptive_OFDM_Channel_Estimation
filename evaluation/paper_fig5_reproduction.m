%% ================================================================
% paper_fig5_reproduction.m
%
% FIGURE 5 REPRODUCTION
%
% "Comparison of the average values of the CIR length estimates
% between SHT and SOT"
%
% (a) SUI-2
% (b) SUI-5
%
% Methods:
%   SHT : T = 15
%   SHT : T = 25
%   SOT
%
% SNR = 0:5:40 dB
%
% IMPORTANT:
%   SHT uses the Hann-windowed CIR.
%   SOT uses the UNWINDOWED CIR so that channel-energy leakage
%   remains present, as required for reproducing the behavior
%   discussed in the paper.
%
% ================================================================

% clear;
clc;
close all;

%% ================================================================
% 1. FIND PROJECT ROOT
% ================================================================

thisDir = fileparts(mfilename('fullpath'));

projectRoot = thisDir;

while ~exist(fullfile(projectRoot, 'config', 'parameters.m'), 'file')

    parentDir = fileparts(projectRoot);

    if strcmp(parentDir, projectRoot)
        error('Could not find project root containing config/parameters.m');
    end

    projectRoot = parentDir;
end

addpath(genpath(projectRoot));

fprintf('\n');
fprintf('====================================================\n');
fprintf('FIGURE 5 REPRODUCTION\n');
fprintf('====================================================\n');

fprintf('Project root:\n%s\n\n', projectRoot);

%% ================================================================
% 2. LOAD PARAMETERS
% ================================================================

p = parameters();

rng(p.randomSeed);

%% Simulation settings

SNR_values = 0:5:40;

numTrials = 1000;

%% SHT parameters

G = p.G;                 % 64
K = p.K;                 % 56

noiseStart = G + 1;      % 65
noiseEnd   = G + K;      % 120

T_SHT_15 = 15;
T_SHT_25 = 25;

%% SOT parameters

% Number of subcarriers used in the SOT theoretical threshold (100 LTS tones)
Np_SOT = 100;

% Estimated number of significant taps
Nt_hat = 3;

%% CIR parameters

M = 128;
Np_active = 100;

%% ================================================================
% 3. BUILD IEEE 802.16 PREAMBLE
% ================================================================

fprintf('Building IEEE 802.16 preamble...\n');

Pall = ieee80216_pall();

preamble = preamble_sequences(Pall);

[stsTime, ltsTime] = ...
    preamble_to_ifft(preamble, p);

txStruct = ...
    build_80216_preamble(stsTime, ltsTime, p);

tx = txStruct.complete;

fprintf('Transmit preamble length = %d samples\n', length(tx));

%% ================================================================
% 4. KNOWN LTS
% ================================================================

knownLTS = get_known_lts(ltsTime);

knownLTS_fft_full = known_lts_fft(knownLTS);

%% ================================================================
% 5. ACTIVE LTS SUBCARRIERS
% ================================================================

kLTS = [-100:2:-2 2:2:100];

q = mod(kLTS/2, 128);

ltsBins = q + 1;

X100 = knownLTS_fft_full(ltsBins);

fprintf('Active LTS tones = %d\n', length(kLTS));

%% ================================================================
% 6. VALIDATE SOT THRESHOLD
% ================================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf('SOT THRESHOLD VALIDATION\n');
fprintf('====================================================\n');

rho20 = 10^(20/10);

xi20 = sqrt( ...
    log(((Np_SOT - Nt_hat) * rho20 * Np_SOT) ...
    / Nt_hat^2) ...
    / ...
    (rho20 * Np_SOT - Nt_hat) );

T_SOT_20 = rho20 * Np_SOT * xi20^2;

fprintf('Np_SOT              = %d\n', Np_SOT);
fprintf('Nt_hat              = %d\n', Nt_hat);
fprintf('SNR                 = 20 dB\n');
fprintf('SOT multiplier      = %.4f\n', T_SOT_20);
fprintf('Expected paper val. = approximately 12.97\n');

%% ================================================================
% 7. STORAGE
% ================================================================

numSNR = length(SNR_values);

meanSHT15_SUI2 = zeros(1,numSNR);
meanSHT25_SUI2 = zeros(1,numSNR);
meanSOT_SUI2   = zeros(1,numSNR);

meanSHT15_SUI5 = zeros(1,numSNR);
meanSHT25_SUI5 = zeros(1,numSNR);
meanSOT_SUI5   = zeros(1,numSNR);

%% ================================================================
% 8. MONTE-CARLO SIMULATION
% ================================================================

for iSNR = 1:numSNR

    SNR = SNR_values(iSNR);

    fprintf('\n');
    fprintf('====================================================\n');
    fprintf('SNR = %d dB\n', SNR);
    fprintf('====================================================\n');

    %% ============================================================
    % SOT THEORETICAL THRESHOLD
    % =============================================================

    rho = 10^(SNR/10);

    numerator = ...
        log(((Np_SOT - Nt_hat) * rho * Np_SOT) ...
        / Nt_hat^2);

    denominator = ...
        rho * Np_SOT - Nt_hat;

    xi_SOT = sqrt(numerator / denominator);

    % Convert normalized amplitude threshold to
    % energy threshold relative to noise variance.

    T_SOT = rho * Np_SOT * xi_SOT^2;

    fprintf('SOT threshold multiplier = %.4f\n', T_SOT);

    %% ============================================================
    % STORAGE FOR THIS SNR
    % =============================================================

    L_SHT15_SUI2 = zeros(numTrials,1);
    L_SHT25_SUI2 = zeros(numTrials,1);
    L_SOT_SUI2   = zeros(numTrials,1);

    L_SHT15_SUI5 = zeros(numTrials,1);
    L_SHT25_SUI5 = zeros(numTrials,1);
    L_SOT_SUI5   = zeros(numTrials,1);

    %% ============================================================
    % SUI-2
    % =============================================================

    fprintf('\nRunning SUI-2...\n');

    for trial = 1:numTrials

        %% --------------------------------------------------------
        % Generate SUI-2 channel
        % ---------------------------------------------------------

        h = sui2_channel(p);

        %% --------------------------------------------------------
        % Pass complete preamble through channel
        % ---------------------------------------------------------

        rxClean = apply_channel_to_preamble(tx,h);

        %% --------------------------------------------------------
        % AWGN
        %
        % Use transmit signal power as SNR reference.
        % ---------------------------------------------------------

        txPower = mean(abs(tx).^2);

        noisePower = ...
            txPower / 10^(SNR/10);

        noise = sqrt(noisePower/2) .* ...
            (randn(size(rxClean)) + ...
             1j*randn(size(rxClean)));

        rx = rxClean + noise;

        %% --------------------------------------------------------
        % Extract second LTS
        % ---------------------------------------------------------

        [~, lts2] = extract_lts(rx,p);

        %% --------------------------------------------------------
        % FFT of received LTS
        % ---------------------------------------------------------

        Y_LTS = lts_fft(lts2);

        %% --------------------------------------------------------
        % Extract 100 active LTS tones
        % ---------------------------------------------------------

        Y100 = Y_LTS(ltsBins);

        %% --------------------------------------------------------
        % LS channel estimate
        % ---------------------------------------------------------

        HLS100 = Y100 ./ X100;

        %% ========================================================
        % SHT PATH
        % ========================================================

        % Hann-windowed CIR

        cirSHT = paper_cir_transform(HLS100);

        %% Noise variance estimated from taps 65:120

        sigma2_SHT = ...
            mean(abs(cirSHT(noiseStart:noiseEnd)).^2);

        %% Threshold for T = 15

        eta15 = T_SHT_15 * sigma2_SHT;

        %% Threshold for T = 25

        eta25 = T_SHT_25 * sigma2_SHT;

        %% CIR length estimates

        L_SHT15_SUI2(trial) = ...
            local_sht_length(cirSHT,G,eta15);

        L_SHT25_SUI2(trial) = ...
            local_sht_length(cirSHT,G,eta25);

        %% ========================================================
        % SOT PATH
        % ========================================================

        % IMPORTANT:
        %
        % SOT does NOT use the Hann-windowed CIR.
        %
        % We deliberately retain the leakage.

        cirSOT = ...
            paper_cir_transform_unwindowed(HLS100);

        %% --------------------------------------------------------
        % Theoretical LS/CIR noise variance
        %
        % FFT:
        %     noise variance becomes M * noisePower
        %
        % LS:
        %     divide by |X|^2
        %
        % CIR transform:
        %     1/sqrt(M)
        %
        % Therefore:
        %
        % sigma_CIR^2 =
        %     sum( M*sigma_n^2/|X_k|^2 ) / M
        %
        % ---------------------------------------------------------

        sigma2_LS = mean( ...
            (M * noisePower) ./ abs(X100).^2 );

        sigma2_SOT = ...
            (Np_active/M) * sigma2_LS;

        %% --------------------------------------------------------
        % SOT threshold
        % ---------------------------------------------------------

        etaSOT = T_SOT * sigma2_SOT;

        %% --------------------------------------------------------
        % SOT CIR length
        % ---------------------------------------------------------

        L_SOT_SUI2(trial) = ...
            local_sht_length(cirSOT,G,etaSOT);

    end

    %% ============================================================
    % SUI-5
    % =============================================================

    fprintf('\nRunning SUI-5...\n');

    for trial = 1:numTrials

        %% --------------------------------------------------------
        % Generate SUI-5 channel
        % ---------------------------------------------------------

        h = sui5_channel(p);

        %% --------------------------------------------------------
        % Pass complete preamble through channel
        % ---------------------------------------------------------

        rxClean = apply_channel_to_preamble(tx,h);

        %% --------------------------------------------------------
        % AWGN
        % ---------------------------------------------------------

        txPower = mean(abs(tx).^2);

        noisePower = ...
            txPower / 10^(SNR/10);

        noise = sqrt(noisePower/2) .* ...
            (randn(size(rxClean)) + ...
             1j*randn(size(rxClean)));

        rx = rxClean + noise;

        %% --------------------------------------------------------
        % Extract second LTS
        % ---------------------------------------------------------

        [~, lts2] = extract_lts(rx,p);

        %% --------------------------------------------------------
        % FFT
        % ---------------------------------------------------------

        Y_LTS = lts_fft(lts2);

        %% --------------------------------------------------------
        % Extract active LTS tones
        % ---------------------------------------------------------

        Y100 = Y_LTS(ltsBins);

        %% --------------------------------------------------------
        % LS channel estimate
        % ---------------------------------------------------------

        HLS100 = Y100 ./ X100;

        %% ========================================================
        % SHT PATH
        % ========================================================

        cirSHT = ...
            paper_cir_transform(HLS100);

        %% Noise variance

        sigma2_SHT = ...
            mean(abs(cirSHT(noiseStart:noiseEnd)).^2);

        %% Thresholds

        eta15 = T_SHT_15 * sigma2_SHT;

        eta25 = T_SHT_25 * sigma2_SHT;

        %% CIR length estimates

        L_SHT15_SUI5(trial) = ...
            local_sht_length(cirSHT,G,eta15);

        L_SHT25_SUI5(trial) = ...
            local_sht_length(cirSHT,G,eta25);

        %% ========================================================
        % SOT PATH
        % ========================================================

        cirSOT = ...
            paper_cir_transform_unwindowed(HLS100);

        %% --------------------------------------------------------
        % Theoretical SOT noise variance
        % --------------------------------------------------------

        sigma2_LS = mean( ...
            (M * noisePower) ./ abs(X100).^2 );

        sigma2_SOT = ...
            (Np_active/M) * sigma2_LS;

        %% SOT threshold

        etaSOT = T_SOT * sigma2_SOT;

        %% SOT CIR length

        L_SOT_SUI5(trial) = ...
            local_sht_length(cirSOT,G,etaSOT);

    end

    %% ============================================================
    % AVERAGES
    % =============================================================

    meanSHT15_SUI2(iSNR) = ...
        mean(L_SHT15_SUI2);

    meanSHT25_SUI2(iSNR) = ...
        mean(L_SHT25_SUI2);

    meanSOT_SUI2(iSNR) = ...
        mean(L_SOT_SUI2);

    meanSHT15_SUI5(iSNR) = ...
        mean(L_SHT15_SUI5);

    meanSHT25_SUI5(iSNR) = ...
        mean(L_SHT25_SUI5);

    meanSOT_SUI5(iSNR) = ...
        mean(L_SOT_SUI5);

    %% ============================================================
    % PRINT CURRENT RESULTS
    % =============================================================

    fprintf('\n');

    fprintf('SUI-2:\n');
    fprintf('  SHT T=15 : %.3f\n', ...
        meanSHT15_SUI2(iSNR));

    fprintf('  SHT T=25 : %.3f\n', ...
        meanSHT25_SUI2(iSNR));

    fprintf('  SOT      : %.3f\n', ...
        meanSOT_SUI2(iSNR));

    fprintf('\n');

    fprintf('SUI-5:\n');
    fprintf('  SHT T=15 : %.3f\n', ...
        meanSHT15_SUI5(iSNR));

    fprintf('  SHT T=25 : %.3f\n', ...
        meanSHT25_SUI5(iSNR));

    fprintf('  SOT      : %.3f\n', ...
        meanSOT_SUI5(iSNR));

end

%% ================================================================
% 9. PRINT COMPLETE RESULT TABLE
% ================================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf('COMPLETE FIG. 5 RESULTS\n');
fprintf('====================================================\n');

fprintf('\n');
fprintf('SUI-2\n');
fprintf('----------------------------------------------------\n');
fprintf(' SNR      SHT15       SHT25        SOT\n');
fprintf('----------------------------------------------------\n');

for i = 1:numSNR

    fprintf('%3d      %7.3f     %7.3f      %7.3f\n', ...
        SNR_values(i), ...
        meanSHT15_SUI2(i), ...
        meanSHT25_SUI2(i), ...
        meanSOT_SUI2(i));

end

fprintf('\n');
fprintf('SUI-5\n');
fprintf('----------------------------------------------------\n');
fprintf(' SNR      SHT15       SHT25        SOT\n');
fprintf('----------------------------------------------------\n');

for i = 1:numSNR

    fprintf('%3d      %7.3f     %7.3f      %7.3f\n', ...
        SNR_values(i), ...
        meanSHT15_SUI5(i), ...
        meanSHT25_SUI5(i), ...
        meanSOT_SUI5(i));

end

%% ================================================================
% 10. PAPER-STYLE FIGURE
% ================================================================

figure( ...
    'Name','Figure 5 - Paper Style', ...
    'Color','w', ...
    'Position',[100 50 700 850]);

%% ================================================================
% FIGURE 5(a) — SUI-2
% ================================================================

subplot(2,1,1);

%% SHT T = 15

plot( ...
    SNR_values, ...
    meanSHT15_SUI2, ...
    '-k*', ...
    'LineWidth',1.0, ...
    'MarkerSize',7);

hold on;

%% SHT T = 25

plot( ...
    SNR_values, ...
    meanSHT25_SUI2, ...
    '-bs', ...
    'LineWidth',1.0, ...
    'MarkerSize',5);

%% SOT

plot( ...
    SNR_values, ...
    meanSOT_SUI2, ...
    '-ro', ...
    'LineWidth',1.0, ...
    'MarkerSize',6);

%% True CIR length

yline( ...
    4.4, ...
    '--k', ...
    'LineWidth',1.0);

%% Formatting

grid on;
box on;

xlim([0 40]);
ylim([0 25]);

xticks(0:5:40);
yticks(0:5:25);

xlabel('SNR (dB)');

ylabel( ...
    'Average Values of CIR Length Estimates');

legend( ...
    'SHT (T=15)', ...
    'SHT (T=25)', ...
    'SOT', ...
    'True CIR Length', ...
    'Location','northeast');

%% Panel label

text( ...
    0.5, ...
    -0.25, ...
    '(a)', ...
    'Units','normalized', ...
    'HorizontalAlignment','center', ...
    'FontSize',11);

%% ================================================================
% FIGURE 5(b) — SUI-5
% ================================================================

subplot(2,1,2);

%% SHT T = 15

plot( ...
    SNR_values, ...
    meanSHT15_SUI5, ...
    '-k*', ...
    'LineWidth',1.0, ...
    'MarkerSize',7);

hold on;

%% SHT T = 25

plot( ...
    SNR_values, ...
    meanSHT25_SUI5, ...
    '-bs', ...
    'LineWidth',1.0, ...
    'MarkerSize',5);

%% SOT

plot( ...
    SNR_values, ...
    meanSOT_SUI5, ...
    '-ro', ...
    'LineWidth',1.0, ...
    'MarkerSize',6);

%% True CIR length

yline( ...
    40, ...
    '--k', ...
    'LineWidth',1.0);

%% Formatting

grid on;
box on;

xlim([0 40]);
ylim([10 64]);

xticks(0:5:40);
yticks(10:5:65);

xlabel('SNR (dB)');

ylabel( ...
    'Average Values of CIR Length Estimates');

legend( ...
    'SHT (T=15)', ...
    'SHT (T=25)', ...
    'SOT', ...
    'True CIR Length', ...
    'Location','southeast');

%% Panel label

text( ...
    0.5, ...
    -0.25, ...
    '(b)', ...
    'Units','normalized', ...
    'HorizontalAlignment','center', ...
    'FontSize',11);

%% ================================================================
% 11. FIGURE FONT
% ================================================================

set( ...
    findall(gcf,'Type','axes'), ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'LineWidth',0.8);

%% ================================================================
% 12. SAVE FIGURE
% ================================================================

figFile = fullfile( ...
    projectRoot, ...
    'evaluation', ...
    'paper_fig5_reproduction.png');

exportgraphics( ...
    gcf, ...
    figFile, ...
    'Resolution',300);

%% ================================================================
% 13. SAVE NUMERICAL RESULTS
% ================================================================

resultFile = fullfile( ...
    projectRoot, ...
    'evaluation', ...
    'paper_fig5_results.mat');

save( ...
    resultFile, ...
    'SNR_values', ...
    'meanSHT15_SUI2', ...
    'meanSHT25_SUI2', ...
    'meanSOT_SUI2', ...
    'meanSHT15_SUI5', ...
    'meanSHT25_SUI5', ...
    'meanSOT_SUI5');

fprintf('\n');
fprintf('====================================================\n');
fprintf('FIGURE 5 COMPLETE\n');
fprintf('====================================================\n');

fprintf('Figure saved to:\n%s\n\n',figFile);

fprintf('Results saved to:\n%s\n',resultFile);


%% =================================================================
% LOCAL FUNCTION 1
% SHT/SOT sequential CIR-length detector
% =================================================================

function Lhat = local_sht_length(cir,G,threshold)

    %% Energy of first G CIR samples

    tapEnergy = abs(cir(1:G)).^2;

    %% Initial estimate

    Lhat = 0;

    %% Sequential hypothesis test

    for n = 1:G

        if tapEnergy(n) > threshold

            Lhat = n;

        end

    end

end


%% =================================================================
% LOCAL FUNCTION 2
% Hann-windowed CIR transformation
%
% Used by SHT.
% =================================================================

function cir = paper_cir_transform(HLS100)

    %% Parameters

    M = 128;
    Np = 100;

    %% -------------------------------------------------------------
    % Reconstruct 101 frequency samples
    %
    % Negative frequencies:
    %
    %   -50 ... -1
    %
    % DC:
    %
    %    0
    %
    % Positive frequencies:
    %
    %   +1 ... +50
    %
    % HLS100 contains:
    %
    %   -100,-98,...,-2,+2,...,+98,+100
    %
    % After the frequency decimation by 2:
    %
    %   -50,-49,...,-1,+1,...,+49,+50
    %
    % -------------------------------------------------------------

    H101 = zeros(101,1);

    H101(1:50) = HLS100(1:50);
    H101(51)   = 0.5 * (HLS100(50) + HLS100(51)); % DC interpolation
    H101(52:101) = HLS100(51:100);

    %% -------------------------------------------------------------
    % Hann window
    %
    % w_n = 0.5 + 0.5*cos(2*pi*n/M)
    %
    % n = -50,...,+50
    %
    % -------------------------------------------------------------

    n = (-Np/2):(Np/2);

    w = ...
        0.5 + 0.5*cos(2*pi*n/M);

    w = w(:);

    %% Apply window

    Hw = H101 .* w;

    %% -------------------------------------------------------------
    % Paper Eq. (23)
    %
    % G(m,n) =
    % 1/sqrt(M)
    % exp(j*2*pi/M*(m-M/2)*(n-Np/2))
    %
    % Here the centered frequency index is
    %
    % q = -50,...,+50
    %
    % -------------------------------------------------------------

    m = (0:M-1).';

    q = (-Np/2):(Np/2);

    Gmat = ...
        (1/sqrt(M)) * ...
        exp( ...
            1j * 2*pi/M * ...
            ((m-M/2) * q));

    %% Transform

    gCentered = Gmat * Hw;

    %% Shift to normal CIR ordering

    cir = fftshift(gCentered);

end


%% =================================================================
% LOCAL FUNCTION 3
% UNWINDOWED CIR transformation
%
% Used by SOT.
%
% IMPORTANT:
% No Hann window is applied here.
%
% This allows channel-energy leakage to remain in the CIR.
% =================================================================

function cir = paper_cir_transform_unwindowed(HLS100)

    %% Parameters

    M = 128;
    Np = 100;

    %% -------------------------------------------------------------
    % Reconstruct 101 frequency samples
    % -------------------------------------------------------------

    H101 = zeros(101,1);

    H101(1:50) = HLS100(1:50);
    H101(51)   = 0.5 * (HLS100(50) + HLS100(51)); % DC interpolation
    H101(52:101) = HLS100(51:100);

    %% -------------------------------------------------------------
    % NO WINDOW
    %
    % This is the important difference from SHT.
    % -------------------------------------------------------------

    Hw = H101;

    %% -------------------------------------------------------------
    % Paper Eq. (23) transform
    % -------------------------------------------------------------

    m = (0:M-1).';

    q = (-Np/2):(Np/2);

    Gmat = ...
        (1/sqrt(M)) * ...
        exp( ...
            1j * 2*pi/M * ...
            ((m-M/2) * q));

    %% Transform

    gCentered = Gmat * Hw;

    %% CIR ordering

    cir = fftshift(gCentered);

end