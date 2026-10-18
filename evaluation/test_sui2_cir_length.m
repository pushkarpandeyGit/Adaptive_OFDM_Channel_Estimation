function test_sui2_cir_length()

clc;
clear;
close all;

addpath(genpath('/MATLAB Drive/Adaptive_OFDM_Channel_Estimation'));

fprintf('\n');
fprintf('============================================================\n');
fprintf(' STEP 41 - SUI-2 CIR LENGTH MONTE CARLO TEST\n');
fprintf('============================================================\n');
%% -----------------------------------------------------------
% Parameters
% ------------------------------------------------------------

p = parameters;

K = p.K;
T = p.T;
G = p.G;

snrDb = 20;

numTrials = 1000;

rng(p.randomSeed);

%% -----------------------------------------------------------
% Generate IEEE 802.16 known LTS
% ------------------------------------------------------------

Pall = ieee80216_pall();

preamble = preamble_sequences(Pall);

[~, ltsTime] = preamble_to_ifft(preamble, p);

knownLTS = get_known_lts(ltsTime);

X_LTS = known_lts_fft(knownLTS);

%% -----------------------------------------------------------
% Generate SUI-2 channel
% ------------------------------------------------------------

h = sui2_channel(p);

%% -----------------------------------------------------------
% Storage for estimated CIR lengths
% ------------------------------------------------------------

Lhat_all = zeros(1, numTrials);

%% -----------------------------------------------------------
% Monte Carlo simulation
% ------------------------------------------------------------

for trial = 1:numTrials

    % -------------------------------------------------------
    % Pass known LTS through SUI-2 channel
    % -------------------------------------------------------

    receivedLTS = apply_channel(knownLTS, h);

    % -------------------------------------------------------
    % Add AWGN at 20 dB SNR
    % -------------------------------------------------------

    receivedLTS = add_awgn(receivedLTS, snrDb);

    % -------------------------------------------------------
    % Keep the received LTS corresponding to the useful
    % 128-sample LTS portion
    % -------------------------------------------------------

    receivedLTS = receivedLTS(1:length(knownLTS));

    % -------------------------------------------------------
    % FFT of received LTS
    % -------------------------------------------------------

    Y_LTS = lts_fft(receivedLTS);

    % -------------------------------------------------------
    % Extract the 100 IEEE 802.16 LTS training tones
    % -------------------------------------------------------

    [Y_p, X_p, ~, ~] = extract_lts_100(Y_LTS, X_LTS);

    % -------------------------------------------------------
    % Least-squares channel estimate
    % -------------------------------------------------------

    H_LS = ls_channel_estimate(Y_p, X_p);

    % -------------------------------------------------------
    % Convert LS frequency response to CIR
    %
    % Paper method:
    %   Hann window
    %   M = 128
    % -------------------------------------------------------

    cir = ls_to_cir(H_LS, 128, p.paperWindow);

    % -------------------------------------------------------
    % Estimate noise variance from CIR taps 65:120
    % -------------------------------------------------------

    noiseVariance = estimate_noise_variance( ...
        cir, ...
        p.noiseTapStart, ...
        p.noiseTapEnd);

    % -------------------------------------------------------
    % Sequential Hypothesis Test
    %
    % Input order:
    % sht_cir_length(cir, noiseVariance, T, G)
    % -------------------------------------------------------

    Lhat = sht_cir_length( ...
        cir, ...
        noiseVariance, ...
        T, ...
        G);

    % -------------------------------------------------------
    % Store estimated CIR length
    % -------------------------------------------------------

    Lhat_all(trial) = Lhat;

end

%% -----------------------------------------------------------
% Calculate Monte Carlo statistics
% ------------------------------------------------------------

meanL = mean(Lhat_all);

stdL = std(Lhat_all);

minL = min(Lhat_all);

maxL = max(Lhat_all);

%% -----------------------------------------------------------
% Display results
% ------------------------------------------------------------

fprintf('\n');
fprintf('============================================================\n');
fprintf(' STEP 41 RESULTS\n');
fprintf('============================================================\n');

fprintf('Channel               = SUI-2\n');
fprintf('SNR                   = %.1f dB\n', snrDb);
fprintf('SHT T                 = %.2f\n', T);
fprintf('G                     = %d\n', G);
fprintf('CIR size M            = %d\n', p.M);
fprintf('Window                = %s\n', p.paperWindow);
fprintf('Number of trials      = %d\n', numTrials);

fprintf('\n');

fprintf('Mean L_hat            = %.4f\n', meanL);
fprintf('Std L_hat             = %.4f\n', stdL);
fprintf('Minimum L_hat         = %d\n', minL);
fprintf('Maximum L_hat         = %d\n', maxL);

fprintf('\n');

fprintf('Paper reference:\n');
fprintf('Mean L_hat            = 5.47\n');
fprintf('Std L_hat             = 0.80\n');

fprintf('\n');

fprintf('Current channel model:\n');
fprintf('SUI-2 implementation is deterministic.\n');
fprintf('The stochastic Ricean K-factor is not yet implemented.\n');

fprintf('============================================================\n');
fprintf('\n');

end