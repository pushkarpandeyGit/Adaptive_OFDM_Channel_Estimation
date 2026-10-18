%% STEP 21E - PAPER ADAPTIVE SVD
% Adaptive rank = L_hat + 3
%
% Flow:
% Preamble
%    ↓
% SUI-2 channel
%    ↓
% AWGN
%    ↓
% Extract second LTS
%    ↓
% 128-point FFT
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
% Estimated CIR length L_hat
%    ↓
% Adaptive rank p = L_hat + 3
%    ↓
% Build covariance matrix
%    ↓
% SVD
%    ↓
% Adaptive SVD channel estimate
%
% IMPORTANT:
% This is the paper's quasi-static SVD path.
% It uses the 100 LTS training tones, NOT the 8 data pilots.

clear;
clc;

%% ================================================================
% 1. PARAMETERS
% ================================================================

p = parameters();
rng(p.randomSeed);

SNRdB = 20;
M = 128;
windowType = p.paperWindow;

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' STEP 21E - PAPER ADAPTIVE SVD\n');
fprintf('=============================================================\n');

fprintf('SNR                 = %d dB\n', SNRdB);
fprintf('CIR/IDFT length M   = %d\n', M);
fprintf('Window               = %s\n', windowType);
fprintf('SHT threshold T      = %d\n', p.T);
fprintf('Maximum CIR length G = %d\n', p.G);

%% ================================================================
% 2. GENERATE IEEE 802.16 PREAMBLE
% ================================================================

fprintf('\n--- STEP 1: Generate preamble ---\n');

Pall = ieee80216_pall();

preambleSeq = preamble_sequences(Pall);

[stsTime,ltsTime] = preamble_to_ifft(preambleSeq,p);

preamble = build_80216_preamble(stsTime,ltsTime,p);

% IMPORTANT:
% build_80216_preamble returns a STRUCT.
% The actual waveform is .complete
txPreamble = preamble.complete;

fprintf('Preamble length = %d samples\n',length(txPreamble));

%% ================================================================
% 3. GENERATE SUI-2 CHANNEL
% ================================================================

fprintf('\n--- STEP 2: Generate SUI-2 channel ---\n');

h = sui2_channel(p);

fprintf('Channel length = %d samples\n',length(h));
fprintf('Channel taps:\n');

for k = 1:length(h)
    fprintf('  h(%d) = %+.6e %+.6ej\n', ...
        k,real(h(k)),imag(h(k)));
end

%% ================================================================
% 4. PASS PREAMBLE THROUGH CHANNEL
% ================================================================

fprintf('\n--- STEP 3: Apply channel ---\n');

rxClean = apply_channel(txPreamble,h);

fprintf('Received clean length = %d samples\n',length(rxClean));

%% ================================================================
% 5. ADD AWGN
% ================================================================

fprintf('\n--- STEP 4: Add AWGN ---\n');

rx = add_awgn(rxClean,SNRdB);

fprintf('SNR = %d dB\n',SNRdB);
fprintf('Received signal length = %d samples\n',length(rx));

%% ================================================================
% 6. EXTRACT LTS
% ================================================================

fprintf('\n--- STEP 5: Extract LTS ---\n');

[rxLTS1,rxLTS2] = extract_lts(rx,p);

fprintf('LTS1 length = %d\n',length(rxLTS1));
fprintf('LTS2 length = %d\n',length(rxLTS2));

%% ================================================================
% 7. FFT OF RECEIVED AND KNOWN LTS
% ================================================================

fprintf('\n--- STEP 6: FFT ---\n');

Y_LTS = lts_fft(rxLTS2);

knownLTS = get_known_lts(ltsTime);

X_LTS = known_lts_fft(knownLTS);

fprintf('Y_LTS length = %d\n',length(Y_LTS));
fprintf('X_LTS length = %d\n',length(X_LTS));

%% ================================================================
% 8. EXTRACT 100 LTS TRAINING TONES
% ================================================================

fprintf('\n--- STEP 7: Extract 100 LTS tones ---\n');

[Yp,Xp,ltsLocations,ltsBins] = ...
    extract_lts_100(Y_LTS,X_LTS);

fprintf('Number of LTS tones = %d\n',length(ltsLocations));
fprintf('Yp size = %d x %d\n',size(Yp,1),size(Yp,2));
fprintf('Xp size = %d x %d\n',size(Xp,1),size(Xp,2));

%% ================================================================
% 9. LS CHANNEL ESTIMATE
% ================================================================

fprintf('\n--- STEP 8: LS channel estimation ---\n');

H_LS_100 = ls_channel_estimate(Yp,Xp);

H_LS_100 = H_LS_100(:);

fprintf('H_LS_100 size = %d x %d\n', ...
    size(H_LS_100,1),size(H_LS_100,2));

%% ================================================================
% 10. CONVERT LS FREQUENCY RESPONSE TO CIR
% ================================================================

fprintf('\n--- STEP 9: LS -> CIR ---\n');

[cir,windowedLS,window] = ...
    ls_to_cir(H_LS_100,M,windowType);

fprintf('CIR length = %d\n',length(cir));

%% ================================================================
% 11. ESTIMATE NOISE VARIANCE
% ================================================================

fprintf('\n--- STEP 10: Estimate noise variance ---\n');

noiseVariance = estimate_noise_variance(cir,p.noiseTapStart,p.noiseTapEnd);

fprintf('Estimated noise variance = %.12e\n',noiseVariance);

%% ================================================================
% 12. SHT CIR LENGTH ESTIMATION
% ================================================================

fprintf('\n--- STEP 11: SHT CIR length estimation ---\n');

L_hat = sht_cir_length( ...
    cir, ...
    noiseVariance, ...
    p.T, ...
    p.G);

fprintf('Estimated CIR length L_hat = %d\n',L_hat);

%% ================================================================
% 13. ADAPTIVE SVD RANK
% ================================================================

fprintf('\n--- STEP 12: Adaptive SVD rank ---\n');

% Paper:
%
% p_adaptive = L_hat + 3

pRankAdaptive = L_hat + p.adaptiveSVDOffset;

% Make sure rank is valid
pRankAdaptive = max(1,min(pRankAdaptive,length(ltsLocations)));

fprintf('Adaptive SVD rank = %d\n',pRankAdaptive);
fprintf('Formula           = L_hat + %d\n',p.adaptiveSVDOffset);
fprintf('                   = %d + %d\n', ...
    L_hat,p.adaptiveSVDOffset);

%% ================================================================
% 14. BUILD CHANNEL COVARIANCE MATRIX
% ================================================================

fprintf('\n--- STEP 13: Build covariance matrix ---\n');

R = build_uniform_covariance(p.N,L_hat);

fprintf('R size = %d x %d\n',size(R,1),size(R,2));

%% ================================================================
% 15. BUILD Rpp
% ================================================================

fprintf('\n--- STEP 14: Build Rpp ---\n');

Rpp = build_Rpp( ...
    R, ...
    ltsLocations, ...
    p.N);

fprintf('Rpp size = %d x %d\n',size(Rpp,1),size(Rpp,2));

%% ================================================================
% 16. BUILD Rhp
% ================================================================

fprintf('\n--- STEP 15: Build Rhp ---\n');

Rhp = build_Rhp( ...
    R, ...
    p.activeIndex, ...
    ltsLocations, ...
    p.N);

fprintf('Rhp size = %d x %d\n',size(Rhp,1),size(Rhp,2));

%% ================================================================
% 17. EIGENVALUE DECOMPOSITION
% ================================================================

fprintf('\n--- STEP 16: Eigenvalue decomposition ---\n');

[U,D] = decompose_Rpp(Rpp);

fprintf('U size = %d x %d\n',size(U,1),size(U,2));
fprintf('D size = %d x %d\n',size(D,1),size(D,2));

%% ================================================================
% 18. EXTRACT EIGENVALUES CORRECTLY
% ================================================================

fprintf('\n--- STEP 17: Extract eigenvalues ---\n');

% D is a diagonal MATRIX.
%
% We need the diagonal values as a VECTOR.
%
% Therefore:
%
% eigenvalues = diag(D)

eigenvalues = diag(D);

% Force column-vector format
eigenvalues = eigenvalues(:);

fprintf('Eigenvalue vector size = %d x %d\n', ...
    size(eigenvalues,1),size(eigenvalues,2));

fprintf('Number of eigenvalues = %d\n',length(eigenvalues));

%% ================================================================
% 19. CHECK EIGENVALUES
% ================================================================

fprintf('\nFirst 10 eigenvalues:\n');

numToShow = min(10,length(eigenvalues));

for k = 1:numToShow
    fprintf('  lambda(%d) = %.12e\n', ...
        k,eigenvalues(k));
end

%% ================================================================
% 20. SVD REGULARIZATION
% ================================================================

fprintf('\n--- STEP 18: SVD filtering matrix ---\n');

% Paper SVD equation uses:
%
% delta_k = 1 / (lambda_k + beta/gamma_p)
%
% for k = 1,...,p
%
% and zero afterwards.

beta = 1;

gamma_p = 10^(SNRdB/10);

fprintf('beta   = %.6f\n',beta);
fprintf('gamma_p = %.6e\n',gamma_p);

%% ================================================================
% 21. BUILD Dp
% ================================================================

Dp = build_Dp( ...
    eigenvalues, ...
    pRankAdaptive, ...
    beta, ...
    gamma_p);

fprintf('\nDp size = %d x %d\n', ...
    size(Dp,1),size(Dp,2));

%% ================================================================
% 22. CRITICAL DIMENSION CHECK
% ================================================================

fprintf('\n=============================================================\n');
fprintf(' DIMENSION CHECK\n');
fprintf('=============================================================\n');

fprintf('Rhp          : %d x %d\n', ...
    size(Rhp,1),size(Rhp,2));

fprintf('U            : %d x %d\n', ...
    size(U,1),size(U,2));

fprintf('Dp           : %d x %d\n', ...
    size(Dp,1),size(Dp,2));

fprintf('H_LS_100     : %d x %d\n', ...
    size(H_LS_100,1),size(H_LS_100,2));

fprintf('Adaptive rank: %d\n',pRankAdaptive);

%% ================================================================
% 23. VERIFY DIMENSIONS
% ================================================================

if size(Rhp,2) ~= size(U,1)
    error('Rhp and U dimensions are inconsistent.');
end

if size(U,2) ~= size(Dp,1)
    error('U and Dp dimensions are inconsistent.');
end

if size(Dp,2) ~= size(U,1)
    error('Dp and U dimensions are inconsistent.');
end

if size(U,2) ~= length(H_LS_100)
    error('U and H_LS_100 dimensions are inconsistent.');
end

fprintf('\nDimension checks: PASS\n');

%% ================================================================
% 24. ADAPTIVE SVD CHANNEL ESTIMATION
% ================================================================

fprintf('\n--- STEP 19: Adaptive SVD channel estimation ---\n');

H_SVD_used = svd_channel_estimate( ...
    Rhp, ...
    U, ...
    Dp, ...
    H_LS_100);

H_SVD_used = H_SVD_used(:);

fprintf('H_SVD_used size = %d x %d\n', ...
    size(H_SVD_used,1),size(H_SVD_used,2));

%% ================================================================
% 25. TRUE CHANNEL FREQUENCY RESPONSE
% ================================================================

fprintf('\n--- STEP 20: Calculate true channel response ---\n');

% Zero-pad channel to N = 256 samples
hPadded = [h(:).' zeros(1,p.N-length(h))];

% 256-point channel frequency response
H_true_full = fftshift(fft(hPadded,p.N));

% Extract the 200 used subcarriers
usedBins = p.activeIndex + p.N/2 + 1;

H_true_used = H_true_full(usedBins);

H_true_used = H_true_used(:);

fprintf('H_true_used size = %d x %d\n', ...
    size(H_true_used,1),size(H_true_used,2));

%% ================================================================
% 26. CALCULATE MSE
% ================================================================

fprintf('\n--- STEP 21: Calculate MSE ---\n');

mseAdaptiveSVD = mean( ...
    abs(H_SVD_used - H_true_used).^2);

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' ADAPTIVE SVD RESULT\n');
fprintf('=============================================================\n');

fprintf('SNR                  = %d dB\n',SNRdB);
fprintf('Estimated CIR length= %d\n',L_hat);
fprintf('Adaptive SVD rank    = %d\n',pRankAdaptive);
fprintf('Rhp                  = %d x %d\n',size(Rhp,1),size(Rhp,2));
fprintf('U                    = %d x %d\n',size(U,1),size(U,2));
fprintf('Dp                   = %d x %d\n',size(Dp,1),size(Dp,2));
fprintf('H estimate           = %d x %d\n', ...
    size(H_SVD_used,1),size(H_SVD_used,2));
fprintf('Adaptive SVD MSE      = %.12e\n',mseAdaptiveSVD);

fprintf('=============================================================\n');

%% ================================================================
% 27. CHECK NUMBER OF RETAINED MODES
% ================================================================

retainedModes = sum(abs(diag(Dp)) > 0);

fprintf('\nRetained SVD modes = %d\n',retainedModes);

if retainedModes == pRankAdaptive
    fprintf('Retained-mode check: PASS\n');
else
    fprintf('Retained-mode check: WARNING\n');
end

%% ================================================================
% 28. FINAL CHECK
% ================================================================

fprintf('\n=============================================================\n');
fprintf(' STEP 21E COMPLETE\n');
fprintf('=============================================================\n');

fprintf('SHT CIR length       : %d\n',L_hat);
fprintf('Adaptive rank        : %d\n',pRankAdaptive);
fprintf('Adaptive SVD MSE     : %.12e\n',mseAdaptiveSVD);
fprintf('=============================================================\n');