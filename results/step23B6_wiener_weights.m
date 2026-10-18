%% ============================================================
% STEP 23-B6 - WIENER WEIGHT GENERATION
% Adaptive OFDM Channel Estimation Project
%
% Purpose:
% Generate scalar time-domain Wiener weights for the
% WDFT baseline estimator.
%
% Output:
% results/wdft_23B6_wiener_weights.mat
%% ============================================================

clear;
clc;

fprintf('\n============================================================\n');
fprintf(' STEP 23-B6 - WIENER WEIGHT GENERATION\n');
fprintf('============================================================\n\n');

%% ============================================================
% PARAMETERS
%% ============================================================

p = parameters();
rng(p.randomSeed);

N = p.N;
M = 128;
Np = 100;
L = 64;

SNRdB = 45;
numTrials = 2000;

fprintf('N          = %d\n', N);
fprintf('M          = %d\n', M);
fprintf('Np         = %d\n', Np);
fprintf('L          = %d\n', L);
fprintf('SNR        = %d dB\n', SNRdB);
fprintf('Trials     = %d\n', numTrials);

%% ============================================================
% LOAD IEEE 802.16 LTS
%% ============================================================

Pall = ieee80216_pall();

preamble = preamble_sequences(Pall);

[~, ltsTime] = preamble_to_ifft(preamble, p);

ltsReference = get_known_lts(ltsTime);

X = known_lts_fft(ltsReference);
X = X(:);

fprintf('\nKnown LTS prepared.\n');
fprintf('LTS length = %d\n', length(ltsReference));
fprintf('FFT length = %d\n', length(X));

%% ============================================================
% LTS SUBCARRIER LOCATIONS
%
% IEEE 802.16:
%
% active subcarriers = -100 : 100 excluding DC
%
% LTS occupies the even subcarriers:
%
% -100,-98,...,-2,+2,...,+98,+100
%
% Total = 100 tones
%% ============================================================

kAll = -100:100;

kNoDC = kAll(kAll ~= 0);

ltsMask = mod(kNoDC,2) == 0;

kLTS = kNoDC(ltsMask);

if length(kLTS) ~= Np
    error('Expected %d LTS tones, found %d.', ...
        Np, length(kLTS));
end

fprintf('\nLTS tone verification:\n');
fprintf('First tone = %d\n', kLTS(1));
fprintf('Last tone  = %d\n', kLTS(end));
fprintf('Number     = %d\n', length(kLTS));

%% ============================================================
% PRECOMPUTE PAPER EQ. (23) MATRIX
%
% [Gtilde]_(m,n) =
%
% 1/sqrt(M) *
% exp(j*2*pi*(m-M/2)*(n-Np/2)/M)
%
% m = 0,...,M-1
% n = 0,...,Np
%
% Here:
% M  = 128
% Np = 100
%
% Therefore G is 128 x 101.
%% ============================================================

m = (0:M-1).';
n = 0:Np;

Gwdft = (1/sqrt(M)) .* ...
    exp(1j * 2*pi .* ...
    ((m - M/2) * (n - Np/2)) / M);

if ~isequal(size(Gwdft), [M Np+1])
    error('WDFT Eq. (23) matrix has incorrect dimensions.');
end

fprintf('\nWDFT matrix:\n');
fprintf('Rows    = %d\n', size(Gwdft,1));
fprintf('Columns = %d\n', size(Gwdft,2));

%% ============================================================
% HANN WINDOW
%
% Paper:
%
% w_n = 0.5 + 0.5*cos(2*pi*n/M)
%
% n = -Np/2,...,+Np/2
%
% Since Np = 100:
%
% n = -50,...,+50
%
% Total = 101 samples
%% ============================================================

nWindow = (-Np/2):(Np/2);

w = 0.5 + ...
    0.5*cos(2*pi*nWindow/M);

w = w(:);

if length(w) ~= Np+1
    error('Hann window must contain exactly 101 samples.');
end

fprintf('\nHann window length = %d\n', length(w));

%% ============================================================
% WIENER STATISTICS
%
% Scalar Wiener coefficient for tap l:
%
%       w_l =
%
% E[g_l * conj(gHat_l)]
% ----------------------
% E[|gHat_l|^2]
%
% where:
%
% g_l    = true CIR
% gHat_l = raw WDFT CIR estimate
%% ============================================================

Rcross = zeros(L,1);
Rhat   = zeros(L,1);

fprintf('\nGenerating %d Monte-Carlo trials...\n', numTrials);

%% ============================================================
% MONTE-CARLO LOOP
%% ============================================================

for trial = 1:numTrials

    %% --------------------------------------------------------
    % Generate stochastic SUI-2 channel
    %% --------------------------------------------------------

    h = sui2_channel(p);
    h = h(:);

    %% --------------------------------------------------------
    % True N-point CFR
    %% --------------------------------------------------------

    Htrue256 = fft(h, N);

    %% --------------------------------------------------------
    % Extract the 100 LTS tones
    %
    % MATLAB FFT indexing:
    %
    % k >= 0 -> k+1
    % k <  0 -> N+k+1
    %% --------------------------------------------------------

    Htrue = zeros(Np,1);

    for ii = 1:Np

        k = kLTS(ii);

        if k >= 0
            idx = k + 1;
        else
            idx = N + k + 1;
        end

        Htrue(ii) = Htrue256(idx);

    end

    %% --------------------------------------------------------
    % Add AWGN to the LTS channel observations
    %% --------------------------------------------------------

    signalPower = mean(abs(Htrue).^2);

    noisePower = signalPower / ...
        (10^(SNRdB/10));

    noise = sqrt(noisePower/2) .* ...
        (randn(Np,1) + ...
         1j*randn(Np,1));

    Hnoisy = Htrue + noise;

    %% --------------------------------------------------------
    % DC ESTIMATION
    %
    % Paper implementation:
    %
    % H(0) = [H(-2) + H(+2)] / 2
    %% --------------------------------------------------------

    idxMinus2 = find(kLTS == -2,1);
    idxPlus2  = find(kLTS == 2,1);

    if isempty(idxMinus2) || isempty(idxPlus2)
        error('Could not locate -2/+2 LTS tones.');
    end

    Hdc = ...
        (Hnoisy(idxMinus2) + ...
         Hnoisy(idxPlus2)) / 2;

    %% --------------------------------------------------------
    % CONSTRUCT CORRECT 101-POINT VECTOR
    %
    % Correct order:
    %
    % -100 ... -2, DC, +2 ... +100
    %
    % First 50 samples:
    % -100 ... -2
    %
    % Sample 51:
    % DC
    %
    % Last 50 samples:
    % +2 ... +100
    %% --------------------------------------------------------

    H101 = [
        Hnoisy(1:50);
        Hdc;
        Hnoisy(51:100)
    ];

    if length(H101) ~= 101
        error('H101 must contain exactly 101 samples.');
    end

    %% --------------------------------------------------------
    % APPLY HANN WINDOW
    %% --------------------------------------------------------

    Hwindowed = H101 .* w;

    %% --------------------------------------------------------
    % ZERO PADDING
    %
    % 101 original samples
    % 128 desired samples
    %
    % 27 zeros required.
    %
    % Paper implementation:
    %
    % 14 front zeros
    % 13 back zeros
    %
    % 14 + 101 + 13 = 128
    %% --------------------------------------------------------

    numOriginal = Np + 1;
    numZeros = M - numOriginal;

    if numZeros ~= 27
        error('Expected 27 zero-padding samples, found %d.', ...
            numZeros);
    end

    frontZeros = ceil(numZeros/2);
    backZeros  = floor(numZeros/2);

    if frontZeros ~= 14 || backZeros ~= 13
        error('Incorrect WDFT zero-padding arrangement.');
    end

    Hpad = [
        zeros(frontZeros,1);
        Hwindowed;
        zeros(backZeros,1)
    ];

    if length(Hpad) ~= M
        error('Final WDFT vector must have length 128.');
    end

    %% --------------------------------------------------------
    % WDFT TRANSFORM
    %
    % Direct implementation of paper Eq. (23)
    %% --------------------------------------------------------

    gCentered = Gwdft * Hwindowed;

    if length(gCentered) ~= M
        error('WDFT output must contain 128 samples.');
    end

    %% --------------------------------------------------------
    % REORDER CENTERED CIR
    %
    % Eq. (23) produces a centered representation.
    %
    % fftshift() moves the physical CIR region to the beginning
    % so that the channel taps can be compared using:
    %
    % gHat(1),...,gHat(L)
    %% --------------------------------------------------------

    gHatFull = fftshift(gCentered);

    %% --------------------------------------------------------
    % TRUE CIR
    %
    % SUI-2 channel has 5 physical samples.
    %
    % Pad to L=64 for tap-wise Wiener statistics.
    %% --------------------------------------------------------

    gTrue = zeros(L,1);

    copyLength = min(length(h), L);

    gTrue(1:copyLength) = h(1:copyLength);

    %% --------------------------------------------------------
    % ACCUMULATE WIENER STATISTICS
    %% --------------------------------------------------------

    gHat = gHatFull(1:L);

    Rcross = Rcross + ...
        gTrue .* conj(gHat);

    Rhat = Rhat + ...
        abs(gHat).^2;

    %% --------------------------------------------------------
    % PROGRESS
    %% --------------------------------------------------------

    if mod(trial,500) == 0
        fprintf('Completed %d / %d trials\n', ...
            trial, numTrials);
    end

end

%% ============================================================
% AVERAGE MONTE-CARLO STATISTICS
%% ============================================================

Rcross = Rcross / numTrials;
Rhat   = Rhat / numTrials;

%% ============================================================
% SCALAR WIENER WEIGHTS
%% ============================================================

wWiener = zeros(L,1);

valid = Rhat > eps;

wWiener(valid) = ...
    real(Rcross(valid) ./ Rhat(valid));

wWiener(~isfinite(wWiener)) = 0;

%% ============================================================
% DISPLAY RESULTS
%% ============================================================

fprintf('\n============================================\n');
fprintf(' WIENER WEIGHT RESULTS\n');
fprintf('============================================\n');

fprintf('Minimum weight = %.10f\n', ...
    min(wWiener));

fprintf('Maximum weight = %.10f\n', ...
    max(wWiener));

fprintf('Mean weight    = %.10f\n', ...
    mean(wWiener));

fprintf('Number of taps = %d\n', L);

%% ============================================================
% BASIC SANITY CHECK
%% ============================================================

if any(~isfinite(wWiener))
    error('Wiener weights contain NaN or Inf.');
end

if any(wWiener < -1e-6)
    warning('Some Wiener weights are negative. Check WDFT alignment.');
end

%% ============================================================
% SAVE
%% ============================================================

fprintf('\n============================================\n');
fprintf(' SAVING WIENER WEIGHTS\n');
fprintf('============================================\n');

outDir = fullfile(pwd,'results');

if ~exist(outDir,'dir')
    mkdir(outDir);
end

outFile = fullfile( ...
    outDir, ...
    'wdft_23B6_wiener_weights.mat');

save(outFile, ...
    'wWiener', ...
    'Rcross', ...
    'Rhat', ...
    'Gwdft', ...
    'w', ...
    'M', ...
    'L', ...
    'Np', ...
    'SNRdB', ...
    'numTrials', ...
    '-mat');

fprintf('Saved:\n%s\n', outFile);

%% ============================================================
% IMMEDIATE MAT-FILE VALIDATION
%% ============================================================

fprintf('\n============================================\n');
fprintf(' MAT-FILE VALIDATION\n');
fprintf('============================================\n');

test = load(outFile);

fprintf('MAT-file loaded successfully.\n');
fprintf('wWiener size = %s\n', ...
    mat2str(size(test.wWiener)));

fprintf('Minimum weight = %.10f\n', ...
    min(test.wWiener));

fprintf('Maximum weight = %.10f\n', ...
    max(test.wWiener));

fprintf('Mean weight    = %.10f\n', ...
    mean(test.wWiener));

fprintf('\n============================================\n');
fprintf(' 23-B6 COMPLETE\n');
fprintf('============================================\n');