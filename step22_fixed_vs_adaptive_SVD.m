%% STEP 22 - FIXED SVD VS ADAPTIVE SVD
%
% Purpose:
% Fair Monte-Carlo comparison between:
%
%   1. Fixed-rank SVD
%      p = G + 3 = 67
%
%   2. Adaptive-rank SVD
%      p = L_hat + 3
%
% IMPORTANT:
% Both estimators use the SAME:
%   - channel realization
%   - transmitted preamble
%   - AWGN realization
%   - LS channel estimate
%   - covariance matrix
%
% Therefore the comparison is fair.

clear;
clc;

%% ============================================================
% STEP 0: PARAMETERS
% =============================================================

p = parameters();

rng(p.randomSeed);

SNRdB = 20;
M = 128;
windowType = p.paperWindow;

numTrials = 100;

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' STEP 22 - FIXED SVD VS ADAPTIVE SVD\n');
fprintf('=============================================================\n');

fprintf('SNR                  = %d dB\n',SNRdB);
fprintf('Number of trials     = %d\n',numTrials);
fprintf('CIR/IDFT length M    = %d\n',M);
fprintf('Window               = %s\n',windowType);

fprintf('\n');
fprintf('Fixed rank           = G + %d = %d\n', ...
    p.fixedSVDOffset, ...
    p.G + p.fixedSVDOffset);

fprintf('Adaptive rank        = L_hat + %d\n', ...
    p.adaptiveSVDOffset);

%% ============================================================
% STEP 1: GENERATE KNOWN PREAMBLE
% =============================================================

fprintf('\n--- STEP 1: Generate known preamble ---\n');

Pall = ieee80216_pall();

preambleSeq = preamble_sequences(Pall);

[stsTime,ltsTime] = ...
    preamble_to_ifft(preambleSeq,p);

preamble = ...
    build_80216_preamble(stsTime,ltsTime,p);

txPreamble = preamble.complete;

fprintf('Preamble length = %d samples\n', ...
    length(txPreamble));

%% ============================================================
% STEP 2: KNOWN LTS
% =============================================================

fprintf('\n--- STEP 2: Prepare known LTS reference ---\n');

knownLTS = get_known_lts(ltsTime);

X_LTS = known_lts_fft(knownLTS);

fprintf('Known LTS length     = %d\n',length(knownLTS));
fprintf('Known LTS FFT length = %d\n',length(X_LTS));

%% ============================================================
% STEP 3: STORAGE
% =============================================================

Lhat_all = zeros(numTrials,1);

adaptiveRank_all = zeros(numTrials,1);

fixedRank_all = ...
    (p.G + p.fixedSVDOffset) * ones(numTrials,1);

mseFixed_all = zeros(numTrials,1);

mseAdaptive_all = zeros(numTrials,1);

noiseVariance_all = zeros(numTrials,1);

%% ============================================================
% STEP 4: MONTE-CARLO LOOP
% =============================================================

fprintf('\n--- STEP 3: Monte-Carlo comparison ---\n');

for trial = 1:numTrials

    %% --------------------------------------------------------
    % A. Generate stochastic SUI-2 channel
    % ---------------------------------------------------------

    h = sui2_channel(p);

    %% --------------------------------------------------------
    % B. Pass SAME preamble through SAME channel
    % ---------------------------------------------------------

    rxClean = ...
        apply_channel(txPreamble,h);

    %% --------------------------------------------------------
    % C. Add SAME AWGN realization
    % ---------------------------------------------------------

    rx = ...
        add_awgn(rxClean,SNRdB);

    %% --------------------------------------------------------
    % D. Extract second LTS
    % ---------------------------------------------------------

    [~,rxLTS2] = ...
        extract_lts(rx,p);

    %% --------------------------------------------------------
    % E. FFT
    % ---------------------------------------------------------

    Y_LTS = ...
        lts_fft(rxLTS2);

    %% --------------------------------------------------------
    % F. Extract 100 LTS training tones
    % ---------------------------------------------------------

    [Yp,Xp,ltsLocations,~] = ...
        extract_lts_100(Y_LTS,X_LTS);

    %% --------------------------------------------------------
    % G. LS channel estimate
    % ---------------------------------------------------------

    H_LS_100 = ...
        ls_channel_estimate(Yp,Xp);

    H_LS_100 = H_LS_100(:);

    %% --------------------------------------------------------
    % H. LS -> CIR
    % ---------------------------------------------------------

    [cir,~,~] = ...
        ls_to_cir( ...
        H_LS_100, ...
        M, ...
        windowType);

    %% --------------------------------------------------------
    % I. Noise variance
    % ---------------------------------------------------------

    noiseVariance = ...
    estimate_noise_variance( ...
    cir, ...
    p.noiseTapStart, ...
    p.noiseTapEnd);

    noiseVariance_all(trial) = ...
        noiseVariance;

    %% --------------------------------------------------------
    % J. SHT CIR length
    % ---------------------------------------------------------

    L_hat = ...
        sht_cir_length( ...
        cir, ...
        noiseVariance, ...
        p.T, ...
        p.G);

    Lhat_all(trial) = ...
        L_hat;

    %% --------------------------------------------------------
    % K. Adaptive rank
    % ---------------------------------------------------------

    pRankAdaptive = ...
        L_hat + p.adaptiveSVDOffset;

    pRankAdaptive = ...
        max( ...
        1, ...
        min(pRankAdaptive,length(ltsLocations)));

    adaptiveRank_all(trial) = ...
        pRankAdaptive;

    %% --------------------------------------------------------
    % L. Fixed rank
    % ---------------------------------------------------------

    pRankFixed = ...
        p.G + p.fixedSVDOffset;

    %% --------------------------------------------------------
    % M. Build covariance using estimated CIR length
    % ---------------------------------------------------------

    R = ...
        build_uniform_covariance( ...
        p.N,L_hat);

    %% --------------------------------------------------------
    % N. Rpp
    % ---------------------------------------------------------

    Rpp = ...
        build_Rpp( ...
        R, ...
        ltsLocations, ...
        p.N);

    %% --------------------------------------------------------
    % O. Rhp
    % ---------------------------------------------------------

    Rhp = ...
        build_Rhp( ...
        R, ...
        p.activeIndex, ...
        ltsLocations, ...
        p.N);

    %% --------------------------------------------------------
    % P. Eigendecomposition
    % ---------------------------------------------------------

    [U,D] = ...
        decompose_Rpp(Rpp);

    %% --------------------------------------------------------
    % Q. Eigenvalues
    % ---------------------------------------------------------

    eigenvalues = ...
        diag(D);

    eigenvalues = ...
        eigenvalues(:);

    %% --------------------------------------------------------
    % R. Regularization
    % ---------------------------------------------------------

    beta = 1;

    gamma_p = ...
        10^(SNRdB/10);

    %% --------------------------------------------------------
    % S. FIXED SVD
    % ---------------------------------------------------------

    DpFixed = ...
        build_Dp( ...
        eigenvalues, ...
        pRankFixed, ...
        beta, ...
        gamma_p);

    H_SVD_fixed = ...
        svd_channel_estimate( ...
        Rhp, ...
        U, ...
        DpFixed, ...
        H_LS_100);

    H_SVD_fixed = ...
        H_SVD_fixed(:);

    %% --------------------------------------------------------
    % T. ADAPTIVE SVD
    % ---------------------------------------------------------

    DpAdaptive = ...
        build_Dp( ...
        eigenvalues, ...
        pRankAdaptive, ...
        beta, ...
        gamma_p);

    H_SVD_adaptive = ...
        svd_channel_estimate( ...
        Rhp, ...
        U, ...
        DpAdaptive, ...
        H_LS_100);

    H_SVD_adaptive = ...
        H_SVD_adaptive(:);

    %% --------------------------------------------------------
    % U. TRUE CHANNEL FREQUENCY RESPONSE
    % ---------------------------------------------------------

    hPadded = [ ...
        h(:).' ...
        zeros(1,p.N-length(h)) ...
        ];

    H_true_full = ...
        fftshift( ...
        fft(hPadded,p.N));

    usedBins = ...
        p.activeIndex + p.N/2 + 1;

    H_true_used = ...
        H_true_full(usedBins);

    H_true_used = ...
        H_true_used(:);

    %% --------------------------------------------------------
    % V. MSE - FIXED SVD
    % ---------------------------------------------------------

    mseFixed_all(trial) = ...
        mean( ...
        abs( ...
        H_SVD_fixed - H_true_used ...
        ).^2);

    %% --------------------------------------------------------
    % W. MSE - ADAPTIVE SVD
    % ---------------------------------------------------------

    mseAdaptive_all(trial) = ...
        mean( ...
        abs( ...
        H_SVD_adaptive - H_true_used ...
        ).^2);

    %% --------------------------------------------------------
    % X. Progress
    % ---------------------------------------------------------

    if mod(trial,10) == 0 || trial == 1

        fprintf( ...
            ['Trial %3d / %3d | ' ...
             'L_hat = %2d | ' ...
             'Fixed p = %2d | ' ...
             'Adaptive p = %2d | ' ...
             'Fixed MSE = %.6e | ' ...
             'Adaptive MSE = %.6e\n'], ...
            trial, ...
            numTrials, ...
            L_hat, ...
            pRankFixed, ...
            pRankAdaptive, ...
            mseFixed_all(trial), ...
            mseAdaptive_all(trial));
    end
end

%% ============================================================
% STEP 5: STATISTICS
% =============================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' MONTE-CARLO RESULTS\n');
fprintf('=============================================================\n');

%% CIR statistics

meanLhat = ...
    mean(Lhat_all);

stdLhat = ...
    std(Lhat_all);

fprintf('\n');
fprintf('CIR LENGTH STATISTICS\n');
fprintf('---------------------\n');

fprintf('Mean L_hat       = %.4f\n',meanLhat);
fprintf('Std L_hat        = %.4f\n',stdLhat);
fprintf('Minimum L_hat    = %d\n',min(Lhat_all));
fprintf('Maximum L_hat    = %d\n',max(Lhat_all));

%% Rank statistics

meanAdaptiveRank = ...
    mean(adaptiveRank_all);

stdAdaptiveRank = ...
    std(adaptiveRank_all);

fprintf('\n');
fprintf('RANK STATISTICS\n');
fprintf('---------------\n');

fprintf('Fixed rank       = %d\n',p.G+p.fixedSVDOffset);

fprintf('Mean adaptive rank = %.4f\n', ...
    meanAdaptiveRank);

fprintf('Std adaptive rank  = %.4f\n', ...
    stdAdaptiveRank);

fprintf('Minimum adaptive rank = %d\n', ...
    min(adaptiveRank_all));

fprintf('Maximum adaptive rank = %d\n', ...
    max(adaptiveRank_all));

%% MSE statistics

meanMSEFixed = ...
    mean(mseFixed_all);

stdMSEFixed = ...
    std(mseFixed_all);

meanMSEAdaptive = ...
    mean(mseAdaptive_all);

stdMSEAdaptive = ...
    std(mseAdaptive_all);

fprintf('\n');
fprintf('MSE STATISTICS\n');
fprintf('--------------\n');

fprintf('Fixed SVD mean MSE      = %.12e\n', ...
    meanMSEFixed);

fprintf('Fixed SVD std MSE       = %.12e\n', ...
    stdMSEFixed);

fprintf('Adaptive SVD mean MSE   = %.12e\n', ...
    meanMSEAdaptive);

fprintf('Adaptive SVD std MSE    = %.12e\n', ...
    stdMSEAdaptive);

%% ============================================================
% STEP 6: MSE COMPARISON
% =============================================================

mseImprovement = ...
    100 * ...
    (meanMSEFixed - meanMSEAdaptive) ...
    / meanMSEFixed;

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' MSE COMPARISON\n');
fprintf('=============================================================\n');

fprintf('Fixed SVD MSE       = %.12e\n', ...
    meanMSEFixed);

fprintf('Adaptive SVD MSE    = %.12e\n', ...
    meanMSEAdaptive);

fprintf('MSE improvement     = %.4f %%\n', ...
    mseImprovement);

%% ============================================================
% STEP 7: RANK REDUCTION
% =============================================================

fixedRank = ...
    p.G + p.fixedSVDOffset;

rankReduction = ...
    100 * ...
    (fixedRank - meanAdaptiveRank) ...
    / fixedRank;

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' RANK / COMPLEXITY COMPARISON\n');
fprintf('=============================================================\n');

fprintf('Fixed rank             = %d\n', ...
    fixedRank);

fprintf('Mean adaptive rank     = %.4f\n', ...
    meanAdaptiveRank);

fprintf('Average rank reduction = %.4f %%\n', ...
    rankReduction);

fprintf('Adaptive uses          = %.4f %% of fixed rank\n', ...
    100*meanAdaptiveRank/fixedRank);

%% ============================================================
% STEP 8: L_HAT DISTRIBUTION
% =============================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' L_HAT DISTRIBUTION\n');
fprintf('=============================================================\n');

uniqueL = ...
    unique(Lhat_all);

for k = 1:length(uniqueL)

    currentL = ...
        uniqueL(k);

    count = ...
        sum(Lhat_all == currentL);

    percentage = ...
        100*count/numTrials;

    fprintf( ...
        'L_hat = %2d : %3d trials (%6.2f%%)\n', ...
        currentL, ...
        count, ...
        percentage);
end

%% ============================================================
% STEP 9: FINAL SUMMARY
% =============================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' STEP 22 COMPLETE\n');
fprintf('=============================================================\n');

fprintf('Trials                    : %d\n',numTrials);
fprintf('SNR                       : %d dB\n',SNRdB);

fprintf('Mean CIR length           : %.4f\n', ...
    meanLhat);

fprintf('Fixed rank                : %d\n', ...
    fixedRank);

fprintf('Mean adaptive rank        : %.4f\n', ...
    meanAdaptiveRank);

fprintf('Fixed SVD MSE             : %.12e\n', ...
    meanMSEFixed);

fprintf('Adaptive SVD MSE          : %.12e\n', ...
    meanMSEAdaptive);

fprintf('MSE improvement           : %.4f %%\n', ...
    mseImprovement);

fprintf('Average rank reduction    : %.4f %%\n', ...
    rankReduction);

fprintf('=============================================================\n');