 clc;
clear;

addpath(genpath('/MATLAB Drive/Adaptive_OFDM_Channel_Estimation'));

fprintf('\n');
fprintf('====================================================\n');
fprintf(' STEP 32A: END-TO-END SUI-2 TEST\n');
fprintf(' NO AWGN — SIGNAL PROCESSING ONLY\n');
fprintf('====================================================\n');

%% ====================================================
% 1. Load project parameters
% =====================================================

p = parameters;

%% ====================================================
% 2. Generate IEEE 802.16 LTS
% =====================================================

Pall = ieee80216_pall;

preamble = preamble_sequences(Pall);

[stsTime,ltsTime] = ...
    preamble_to_ifft(preamble,p);

knownLTS = get_known_lts(ltsTime);

X_LTS = known_lts_fft(knownLTS);

%% ====================================================
% 3. Create the SUI-2 channel
% =====================================================

h = sui2_channel(p);

%% ====================================================
% 4. Pass the COMPLETE preamble through channel
% =====================================================

preambleStruct = ...
    build_80216_preamble(stsTime,ltsTime,p);

tx = preambleStruct.complete;

rxClean = apply_channel(tx,h);

%% ====================================================
% 5. Extract received second LTS
% =====================================================

[lts1_rx,lts2_rx] = ...
    extract_lts(rxClean,p);

%% ====================================================
% 6. FFT received LTS
% =====================================================

Y_LTS = lts_fft(lts2_rx);

%% ====================================================
% 7. Extract 100 LTS tones
% =====================================================

[Y_p,X_p,ltsLocations,ltsBins] = ...
    extract_lts_100(Y_LTS,X_LTS);

%% ====================================================
% 8. 100-point LS channel estimate
% =====================================================

hLS_100 = pilot_ls_estimate(Y_p,X_p);

%% ====================================================
% 9. Convert LS estimate to CIR
% =====================================================

[cir,windowedLS,window,Gmatrix] = ...
    ls100_to_cir128(hLS_100);

%% ====================================================
% 10. Estimate noise variance
%
% IMPORTANT:
% There is NO AWGN in this test.
%
% Therefore this is NOT a final noise-performance
% experiment.
%
% We only use the CIR tail to verify the complete
% processing chain.
% =====================================================

[noiseVariance,noiseTaps] = ...
    estimate_noise_variance(cir,65,120);

%% ====================================================
% 11. SHT
% =====================================================

[L_hat,threshold,tapEnergy,decision] = ...
    sht_cir_length(cir,noiseVariance,p.T,p.G);

%% ====================================================
% 12. Display results
% =====================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' STEP 32A RESULTS\n');
fprintf('====================================================\n');

fprintf('SUI-2 channel length       : %d samples\n',length(h));
fprintf('Number of LTS tones        : %d\n',length(ltsLocations));
fprintf('LS estimate length         : %d\n',length(hLS_100));
fprintf('CIR length                 : %d\n',length(cir));
fprintf('Noise variance estimate   : %.6e\n',noiseVariance);
fprintf('SHT threshold              : %.6e\n',threshold);
fprintf('Estimated CIR length L_hat : %d samples\n',L_hat);

fprintf('====================================================\n');

%% ====================================================
% 13. Basic structural checks
% =====================================================

assert(length(tx) == 640);
assert(length(lts2_rx) == 128);
assert(length(Y_LTS) == 128);
assert(length(Y_p) == 100);
assert(length(X_p) == 100);
assert(length(hLS_100) == 100);
assert(length(cir) == 128);
assert(length(noiseTaps) == 56);
assert(length(tapEnergy) == 64);
assert(length(decision) == 64);

fprintf('\nPASS: all signal dimensions are correct.\n');

fprintf('\n');
fprintf('====================================================\n');
fprintf(' STEP 32A STRUCTURAL TEST PASSED\n');
fprintf('====================================================\n');