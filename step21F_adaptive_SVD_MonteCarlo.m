%% STEP 21F - ADAPTIVE SVD MONTE-CARLO VALIDATION
%
% Purpose:
% Validate the paper's adaptive SVD estimator statistically.
%
% For every trial:
%
% Preamble
%    ↓
% SUI-2 channel
%    ↓
% AWGN
%    ↓
% Extract second LTS
%    ↓
% FFT
%    ↓
% Extract 100 LTS tones
%    ↓
% LS channel estimate
%    ↓
% Hann window + IDFT
%    ↓
% Noise variance
%    ↓
% SHT
%    ↓
% L_hat
%    ↓
% Adaptive rank = L_hat + 3
%    ↓
% Covariance matrix
%    ↓
% SVD
%    ↓
% Adaptive SVD estimate
%    ↓
% MSE
%
% This script does NOT modify the existing Step 21E implementation.

clear;
clc;

%% ============================================================
%  STEP 0: PARAMETERS
% =============================================================

p = parameters();

rng(p.randomSeed);

SNRdB = 20;
M = 128;
windowType = p.paperWindow;

numTrials = 100;

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' STEP 21F - ADAPTIVE SVD MONTE-CARLO VALIDATION\n');
fprintf('=============================================================\n');

fprintf('SNR                  = %d dB\n',SNRdB);
fprintf('Number of trials     = %d\n',numTrials);
fprintf('CIR/IDFT length M    = %d\n',M);
fprintf('Window               = %s\n',windowType);
fprintf('SHT threshold T      = %d\n',p.T);
fprintf('Maximum CIR length G= %d\n',p.G);
fprintf('Adaptive rank        = L_hat + %d\n', ...
    p.adaptiveSVDOffset);


%% ============================================================
%  STEP 1: GENERATE KNOWN PREAMBLE
% =============================================================

fprintf('\n--- STEP 1: Generate known preamble ---\n');

Pall = ieee80216_pall();

preambleSeq = preamble_sequences(Pall);

[stsTime,ltsTime] = preamble_to_ifft(preambleSeq,p);

preamble = build_80216_preamble(stsTime,ltsTime,p);

txPreamble = preamble.complete;

fprintf('Preamble length = %d samples\n', ...
    length(txPreamble));


%% ============================================================
%  STEP 2: KNOWN LTS REFERENCE
% =============================================================

fprintf('\n--- STEP 2: Prepare known LTS reference ---\n');

knownLTS = get_known_lts(ltsTime);

X_LTS = known_lts_fft(knownLTS);

fprintf('Known LTS length = %d\n',length(knownLTS));
fprintf('Known LTS FFT length = %d\n',length(X_LTS));


%% ============================================================
%  STEP 3: PREPARE STORAGE
% =============================================================

Lhat_all = zeros(numTrials,1);

rank_all = zeros(numTrials,1);

mse_all = zeros(numTrials,1);

noiseVariance_all = zeros(numTrials,1);

channelLength_all = zeros(numTrials,1);


%% ============================================================
%  STEP 4: MONTE-CARLO LOOP
% =============================================================

fprintf('\n--- STEP 3: Monte-Carlo simulation ---\n');

for trial = 1:numTrials

    %% --------------------------------------------------------
    % Generate new stochastic SUI-2 channel
    % ---------------------------------------------------------

    h = sui2_channel(p);

    channelLength_all(trial) = length(h);


    %% --------------------------------------------------------
    % Pass known preamble through channel
    % ---------------------------------------------------------

    rxClean = apply_channel(txPreamble,h);


    %% --------------------------------------------------------
    % Add AWGN
    % ---------------------------------------------------------

    rx = add_awgn(rxClean,SNRdB);


    %% --------------------------------------------------------
    % Extract second LTS
    % ---------------------------------------------------------

    [~,rxLTS2] = extract_lts(rx,p);


    %% --------------------------------------------------------
    % FFT of received LTS
    % ---------------------------------------------------------

    Y_LTS = lts_fft(rxLTS2);


    %% --------------------------------------------------------
    % Extract 100 LTS training tones
    % ---------------------------------------------------------

    [Yp,Xp,ltsLocations,~] = ...
        extract_lts_100(Y_LTS,X_LTS);


    %% --------------------------------------------------------
    % LS channel estimate
    % ---------------------------------------------------------

    H_LS_100 = ls_channel_estimate(Yp,Xp);

    H_LS_100 = H_LS_100(:);


    %% --------------------------------------------------------
    % LS -> CIR
    % ---------------------------------------------------------

    [cir,~,~] = ...
        ls_to_cir(H_LS_100,M,windowType);


    %% --------------------------------------------------------
    % Noise variance estimation
    % ---------------------------------------------------------

    noiseVariance = estimate_noise_variance( ...
        cir, ...
        p.noiseTapStart, ...
        p.noiseTapEnd);

    noiseVariance_all(trial) = noiseVariance;


    %% --------------------------------------------------------
    % SHT CIR length estimation
    % ---------------------------------------------------------

    L_hat = sht_cir_length( ...
        cir, ...
        noiseVariance, ...
        p.T, ...
        p.G);

    Lhat_all(trial) = L_hat;


    %% --------------------------------------------------------
    % Adaptive SVD rank
    % ---------------------------------------------------------

    pRankAdaptive = ...
        L_hat + p.adaptiveSVDOffset;

    pRankAdaptive = max( ...
        1, ...
        min(pRankAdaptive,length(ltsLocations)));

    rank_all(trial) = pRankAdaptive;


    %% --------------------------------------------------------
    % Build covariance matrix
    % ---------------------------------------------------------

    R = build_uniform_covariance( ...
        p.N,L_hat);


    %% --------------------------------------------------------
    % Rpp
    % ---------------------------------------------------------

    Rpp = build_Rpp( ...
        R, ...
        ltsLocations, ...
        p.N);


    %% --------------------------------------------------------
    % Rhp
    % ---------------------------------------------------------

    Rhp = build_Rhp( ...
        R, ...
        p.activeIndex, ...
        ltsLocations, ...
        p.N);


    %% --------------------------------------------------------
    % Eigendecomposition
    % ---------------------------------------------------------

    [U,D] = decompose_Rpp(Rpp);


    %% --------------------------------------------------------
    % Eigenvalues
    % ---------------------------------------------------------

    eigenvalues = diag(D);

    eigenvalues = eigenvalues(:);


    %% --------------------------------------------------------
    % SVD filtering matrix
    % ---------------------------------------------------------

    beta = 1;

    gamma_p = 10^(SNRdB/10);

    Dp = build_Dp( ...
        eigenvalues, ...
        pRankAdaptive, ...
        beta, ...
        gamma_p);


    %% --------------------------------------------------------
    % Adaptive SVD channel estimate
    % ---------------------------------------------------------

    H_SVD_used = svd_channel_estimate( ...
        Rhp, ...
        U, ...
        Dp, ...
        H_LS_100);

    H_SVD_used = H_SVD_used(:);


    %% --------------------------------------------------------
    % True channel frequency response
    % ---------------------------------------------------------

    hPadded = [ ...
        h(:).' ...
        zeros(1,p.N-length(h)) ...
        ];

    H_true_full = fftshift( ...
        fft(hPadded,p.N));


    usedBins = ...
        p.activeIndex + p.N/2 + 1;


    H_true_used = ...
        H_true_full(usedBins);

    H_true_used = H_true_used(:);


    %% --------------------------------------------------------
    % MSE
    % ---------------------------------------------------------

    mse_all(trial) = ...
        mean(abs(H_SVD_used-H_true_used).^2);


    %% --------------------------------------------------------
    % Progress display
    % ---------------------------------------------------------

    if mod(trial,10) == 0 || trial == 1

        fprintf( ...
            'Trial %3d / %3d | L_hat = %2d | Rank = %2d | MSE = %.6e\n', ...
            trial, ...
            numTrials, ...
            L_hat, ...
            pRankAdaptive, ...
            mse_all(trial));

    end

end


%% ============================================================
%  STEP 5: STATISTICS
% =============================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' MONTE-CARLO RESULTS\n');
fprintf('=============================================================\n');


meanLhat = mean(Lhat_all);

stdLhat = std(Lhat_all);


meanRank = mean(rank_all);

stdRank = std(rank_all);


meanMSE = mean(mse_all);

stdMSE = std(mse_all);


fprintf('\n');
fprintf('CIR LENGTH STATISTICS\n');
fprintf('---------------------\n');

fprintf('Mean L_hat       = %.4f\n',meanLhat);

fprintf('Std L_hat        = %.4f\n',stdLhat);

fprintf('Minimum L_hat    = %d\n',min(Lhat_all));

fprintf('Maximum L_hat    = %d\n',max(Lhat_all));


fprintf('\n');
fprintf('ADAPTIVE RANK STATISTICS\n');
fprintf('------------------------\n');

fprintf('Mean rank        = %.4f\n',meanRank);

fprintf('Std rank         = %.4f\n',stdRank);

fprintf('Minimum rank     = %d\n',min(rank_all));

fprintf('Maximum rank     = %d\n',max(rank_all));


fprintf('\n');
fprintf('MSE STATISTICS\n');
fprintf('--------------\n');

fprintf('Mean adaptive SVD MSE = %.12e\n',meanMSE);

fprintf('Std adaptive SVD MSE  = %.12e\n',stdMSE);

fprintf('Minimum MSE           = %.12e\n',min(mse_all));

fprintf('Maximum MSE           = %.12e\n',max(mse_all));


%% ============================================================
%  STEP 6: L_HAT DISTRIBUTION
% =============================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' L_HAT DISTRIBUTION\n');
fprintf('=============================================================\n');

uniqueL = unique(Lhat_all);

for k = 1:length(uniqueL)

    currentL = uniqueL(k);

    count = sum(Lhat_all == currentL);

    percentage = 100*count/numTrials;

    fprintf( ...
        'L_hat = %2d : %3d trials (%6.2f%%)\n', ...
        currentL, ...
        count, ...
        percentage);

end


%% ============================================================
%  STEP 7: RANK DISTRIBUTION
% =============================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' ADAPTIVE RANK DISTRIBUTION\n');
fprintf('=============================================================\n');

uniqueRank = unique(rank_all);

for k = 1:length(uniqueRank)

    currentRank = uniqueRank(k);

    count = sum(rank_all == currentRank);

    percentage = 100*count/numTrials;

    fprintf( ...
        'Rank = %2d : %3d trials (%6.2f%%)\n', ...
        currentRank, ...
        count, ...
        percentage);

end


%% ============================================================
%  STEP 8: FINAL SUMMARY
% =============================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' STEP 21F COMPLETE\n');
fprintf('=============================================================\n');

fprintf('Number of trials       : %d\n',numTrials);

fprintf('SNR                    : %d dB\n',SNRdB);

fprintf('Mean CIR length        : %.4f\n',meanLhat);

fprintf('Std CIR length         : %.4f\n',stdLhat);

fprintf('Mean adaptive rank     : %.4f\n',meanRank);

fprintf('Mean adaptive SVD MSE  : %.12e\n',meanMSE);

fprintf('Std adaptive SVD MSE   : %.12e\n',stdMSE);

fprintf('=============================================================\n');