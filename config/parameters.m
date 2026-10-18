function p = parameters()
% ==============================================================
% parameters.m
%
% PROJECT:
% Adaptive OFDM Channel Estimation Based on CIR Length Adaptation
%
% BASE PAPER:
% "A Low-Complexity Channel Estimation for OFDM Systems
%  Based on CIR Length Adaptation"
%
% This file contains:
%   1. IEEE 802.16 OFDM parameters
%   2. Channel parameters
%   3. SHT parameters
%   4. SVD parameters
%   5. Our future novelty parameters
%
% IMPORTANT:
% Paper parameters and our proposed parameters are kept
% separately so that we can clearly distinguish reproduction
% from novelty.
% ==============================================================


%% =============================================================
% 1. IEEE 802.16 OFDM PARAMETERS
% =============================================================

% Total number of subcarriers
p.N = 256;

% Number of used subcarriers
p.Nu = 200;

% Number of data subcarriers
p.Nd = 192;

% Number of pilot subcarriers
p.Np = 8;

% Guard interval / Cyclic Prefix
p.G = 64;


%% =============================================================
% 2. FREQUENCY-DOMAIN SUBCARRIER INDICES
% =============================================================

% IEEE 802.16 uses subcarrier indices:
%
%       -100 ... -1, 0, 1 ... 100
%
% Subcarrier 0 is null.
%
% Total used subcarriers:
% 200
%
% Null subcarriers:
% 56

p.subcarrierIndex = -100:100;

% Remove DC subcarrier
p.activeIndex = p.subcarrierIndex( ...
    p.subcarrierIndex ~= 0);


%% =============================================================
% 3. PILOT SUBCARRIERS
% =============================================================

% Pilot locations specified in the paper

p.pilotLocations = [ ...
    -88 ...
    -63 ...
    -38 ...
    -13 ...
     13 ...
     38 ...
     63 ...
     88 ...
];

% Pilot modulation
p.pilotModulation = 'BPSK';


%% =============================================================
% 4. DATA MODULATION
% =============================================================

% 16-QAM is used for the SER simulations

p.M = 16;

% Number of bits per QAM symbol

p.bitsPerQAMSymbol = log2(p.M);


%% =============================================================
% 5. SYSTEM BANDWIDTH
% =============================================================

% Paper:
% BW = 3.5 MHz

p.BW = 3.5e6;


%% =============================================================
% 6. SAMPLING INTERVAL
% =============================================================

% Paper:
% Ts = 0.25 microseconds

p.Ts = 0.25e-6;


%% =============================================================
% 7. SHT PARAMETERS
% =============================================================

% Guard interval length

p.G = 64;

% Number of samples used to estimate noise variance
%
% Paper uses:
%
%       n = G+1 ... G+K
%
%       65 ... 120
%
% Therefore:
%
%       K = 56

p.K = 56;


% Scaling factor used in the paper's SHT simulations

p.T = 15;


% Candidate T values
%
% Used later for:
%   - studying SHT behavior
%   - ML dataset generation
%   - ML-based adaptive threshold

p.candidateT = 5:25;


%% =============================================================
% 8. SHT TAP RANGES
% =============================================================

% Noise-estimation taps:
%
% G+1 to G+K
%
% = 65 to 120

p.noiseTapStart = p.G + 1;

p.noiseTapEnd = p.G + p.K;


% CIR taps that are tested sequentially

p.testTapStart = 1;

p.testTapEnd = p.G;


%% =============================================================
% 9. SVD PARAMETERS
% =============================================================

% Fixed SVD estimator:
%
%       p = G + 3

p.fixedSVDOffset = 3;

% Adaptive SVD estimator:
%
%       p_hat = L_hat + 3

p.adaptiveSVDOffset = 3;


%% =============================================================
% 10. CHANNEL MODELS
% =============================================================

p.channelModels = { ...
    'SUI-2', ...
    'SUI-5', ...
    'Vehicular-A' ...
};


%% =============================================================
% 11. SUI-2 CHANNEL PARAMETERS
% =============================================================

% Relative path powers in dB

p.SUI2.pathPower_dB = [0 -12 -15];

% Path delays in microseconds

p.SUI2.pathDelay_us = [0 0.4 1.1];

% First path Ricean K-factor

p.SUI2.RiceanK_dB = 3;


%% =============================================================
% 12. SUI-5 CHANNEL PARAMETERS
% =============================================================

% Relative path powers in dB

p.SUI5.pathPower_dB = [0 -5 -10];

% Path delays in microseconds

p.SUI5.pathDelay_us = [0 4.0 10.0];


%% =============================================================
% 13. ITU-R VEHICULAR A PARAMETERS
% =============================================================

% Relative path powers in dB

p.VehA.pathPower_dB = ...
    [0 -1 -9 -10 -15 -20];


% Path delays in microseconds

p.VehA.pathDelay_us = ...
    [0 0.31 0.71 1.09 1.73 2.01];


%% =============================================================
% 14. TIME-VARYING CHANNEL
% =============================================================

% Doppler frequency used in the paper's
% time-varying-channel experiment

p.fd = 30;


% Number of decision-directed iterations
%
% Paper uses two iterations for the
% time-varying channel simulation.

p.numIterations = 2;


%% =============================================================
% 15. SER SIMULATION PARAMETERS
% =============================================================

% Carrier frequency used for the 16-QAM SER simulation

p.fc = 1e9;


% Maximum CFO due to oscillator inaccuracy

p.maxCFOppm = 20;


%% =============================================================
% 16. SIMULATION SNR RANGE
% =============================================================

% General simulation range.
%
% We will change this for individual paper figures
% when reproducing them exactly.

p.SNRdB = 0:5:40;


%% =============================================================
% 17. MONTE CARLO PARAMETERS
% =============================================================

% Temporary value for debugging.
%
% We will increase this later for final results.

p.numFrames = 100;


% Reproducible random seed

p.randomSeed = 42;


%% =============================================================
% 18. ORIGINAL PAPER WINDOW
% =============================================================

% The original paper uses a Hann window before
% the M-point IDFT for CIR estimation.

p.paperWindow = 'hann';


%% =============================================================
% 19. OUR NOVELTY — ADAPTIVE WINDOW
% =============================================================

% These are NOT parameters from the original paper.
%
% These belong to our proposed modification.

p.adaptiveWindows = { ...
    'hann', ...
    'hamming', ...
    'blackman' ...
};


%% =============================================================
% 20. OUR NOVELTY — ML FEATURES
% =============================================================

% Candidate features for predicting T_ML

p.ML.features = { ...
    'SNR', ...
    'NoiseVariance', ...
    'DelaySpread', ...
    'CIRFeature', ...
    'WindowType' ...
};


%% =============================================================
% 21. DISPLAY PARAMETERS
% =============================================================

p.showParameters = true;


%% =============================================================
% DISPLAY
% =============================================================

if p.showParameters

    fprintf('\n');
    fprintf('====================================================\n');
    fprintf(' Adaptive OFDM Channel Estimation Project\n');
    fprintf('====================================================\n');

    fprintf('\n--- IEEE 802.16 OFDM ---\n');

    fprintf('N  (total subcarriers) : %d\n', p.N);
    fprintf('Nu (used subcarriers)  : %d\n', p.Nu);
    fprintf('Nd (data subcarriers)  : %d\n', p.Nd);
    fprintf('Np (pilot subcarriers) : %d\n', p.Np);
    fprintf('G  (guard interval)    : %d\n', p.G);

    fprintf('\n--- Modulation ---\n');

    fprintf('Modulation             : %d-QAM\n', p.M);
    fprintf('Bits/QAM symbol        : %d\n', ...
        p.bitsPerQAMSymbol);

    fprintf('\n--- System ---\n');

    fprintf('Bandwidth              : %.2f MHz\n', ...
        p.BW / 1e6);

    fprintf('Sampling interval      : %.2f us\n', ...
        p.Ts * 1e6);

    fprintf('\n--- SHT ---\n');

    fprintf('K                      : %d\n', p.K);
    fprintf('Noise taps             : %d to %d\n', ...
        p.noiseTapStart, p.noiseTapEnd);

    fprintf('Paper T                : %d\n', p.T);

    fprintf('\n--- SVD ---\n');

    fprintf('Fixed rank offset      : %d\n', ...
        p.fixedSVDOffset);

    fprintf('Adaptive rank offset   : %d\n', ...
        p.adaptiveSVDOffset);

    fprintf('\n--- Time-varying channel ---\n');

    fprintf('Doppler frequency      : %d Hz\n', p.fd);
    fprintf('Iterations             : %d\n', ...
        p.numIterations);

    fprintf('\n--- Novelty ---\n');

    fprintf('Adaptive windows       : Hann / Hamming / Blackman\n');
    fprintf('ML threshold           : Enabled for later stage\n');

    fprintf('\n====================================================\n\n');

end

end