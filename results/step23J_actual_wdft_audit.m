%% STEP 23-J
% ACTUAL WDFT ESTIMATOR AUDIT
%
% Flow:
%   True CFR
%       |
%       v
%   100-tone LS CFR
%       |
%       v
%   DC substitution
%       |
%       v
%   101-point Hann window
%       |
%       v
%   Eq. (23) WDFT transform
%       |
%       v
%   128-point CIR
%       |
%       v
%   Time-domain Wiener weighting
%       |
%       v
%   Reconstructed CFR
%
% This script does NOT use the unstable TSVD/Tikhonov inverse
% from Step 23-H.
%
% The purpose of this step is to audit the ACTUAL WDFT estimator
% behavior before locking Step 23.

clear;
clc;
close all;

%% ============================================================
% 0. PARAMETERS
% =============================================================

p = parameters();
rng(p.randomSeed);

N = p.N;                  % 256
M = 128;                  % WDFT transform size
G = p.G;                  % 64
T = p.T;

SNRdB = p.SNRdB;

numTrain = 2000;
numTest  = 500;

fprintf('\n');
fprintf('============================================================\n');
fprintf(' STEP 23-J: ACTUAL WDFT ESTIMATOR AUDIT\n');
fprintf('============================================================\n');

fprintf('N              = %d\n', N);
fprintf('M              = %d\n', M);
fprintf('G              = %d\n', G);
fprintf('Training trials = %d\n', numTrain);
fprintf('Test trials     = %d\n', numTest);
fprintf('\n');


%% ============================================================
% 1. BUILD EXACT EQ. (23) MATRIX
% =============================================================

% Pilot subcarrier locations:
%
% k = -100:2:100, excluding DC
%
% There are exactly 100 LTS tones.

evenK = -100:2:100;
evenK_noDC = evenK(evenK ~= 0);

assert(length(evenK_noDC) == 100, ...
    'Expected exactly 100 LTS tones.');

% Eq. (23):
%
% [Gtilde]_(m,n)
% =
% 1/sqrt(M) *
% exp(j*2*pi*(m-M/2)*(n-Np/2)/M)
%
% MATLAB indices are shifted by +1.

Np = 100;

Gmat = zeros(M, Np + 1);

for m = 0:M-1
    for n = 0:Np
        Gmat(m+1,n+1) = ...
            (1/sqrt(M)) * ...
            exp(1j*2*pi*(m-M/2)*(n-Np/2)/M);
    end
end

fprintf('Eq. (23) matrix size: %d x %d\n', ...
    size(Gmat,1), size(Gmat,2));


%% ============================================================
% 2. BUILD PHYSICAL 256-POINT CFR -> 100 LTS MAPPING
% =============================================================

% Physical 256-point FFT indexing for signed subcarrier k:
%
% MATLAB index = k + N/2 + 1

idx256 = evenK_noDC + N/2 + 1;

assert(all(idx256 >= 1 & idx256 <= N), ...
    'Invalid physical CFR indices.');


%% ============================================================
% 3. DC LOCATION IN THE 101-ELEMENT VECTOR
% =============================================================

% The 100-tone vector is ordered:
%
% [-100,-98,...,-2, +2,+4,...,+100]
%
% Insert DC between -2 and +2.
%
% Therefore:
%
% H101 =
% [first 50 negative-frequency tones
%  DC estimate
%  last 50 positive-frequency tones]

numNegative = sum(evenK_noDC < 0);

assert(numNegative == 50, ...
    'Expected 50 negative-frequency LTS tones.');

dcPosition = numNegative + 1;

fprintf('DC position in 101-vector = %d\n', dcPosition);


%% ============================================================
% 4. HANN WINDOW
% =============================================================

nWin = -Np/2:Np/2;

hannWindow = 0.5 + ...
             0.5*cos(2*pi*nWin/M);

hannWindow = hannWindow(:);

assert(length(hannWindow) == 101, ...
    'Hann window must contain 101 samples.');

W = diag(hannWindow);


%% ============================================================
% 5. ZERO-PADDING STRUCTURE
% =============================================================

% Paper:
%
% (M-Np)/2 zeros at the front
% (M-Np)/2 - 1 zeros at the back
%
% M = 128
% Np = 100
%
% front = 14
% back  = 13

numFrontZeros = (M-Np)/2;
numBackZeros  = numFrontZeros - 1;

fprintf('Front zeros = %d\n', numFrontZeros);
fprintf('Back zeros  = %d\n', numBackZeros);

assert(numFrontZeros == 14, ...
    'Unexpected front zero count.');

assert(numBackZeros == 13, ...
    'Unexpected back zero count.');


%% ============================================================
% 6. GENERATE KNOWN LTS
% =============================================================

Pall = ieee80216_pall();

preamble = preamble_sequences(Pall);

[~, ltsTime] = preamble_to_ifft(preamble, p);

knownLTS = get_known_lts(ltsTime);

X128 = known_lts_fft(knownLTS);

% Extract the 100 physical LTS tones.

idx128 = mod(evenK_noDC, M) + 1;

X100 = X128(idx128);

assert(length(X100) == 100, ...
    'Expected 100 LTS pilot symbols.');


%% ============================================================
% 7. FUNCTION: BUILD 101-TONE VECTOR
% =============================================================

% We use a local function at the bottom of this script.
%
% DC is estimated as:
%
% H_DC = (H_-2 + H_+2)/2
%
% Then:
%
% H101 =
% [H_-100 ... H_-2 H_DC H_+2 ... H_+100]


%% ============================================================
% 8. FUNCTION: WDFT TRANSFORM
% =============================================================

% The actual WDFT transform used here is:
%
% g_centered = Gmat * W * H101
%
% followed by fftshift to put the physical channel taps
% at the beginning of the CIR representation.
%
% This is consistent with the tap ordering established in
% the previous WDFT audits.

%% ============================================================
% 9. TRAIN WIENER WEIGHTS
% =============================================================

fprintf('\n');
fprintf('------------------------------------------------------------\n');
fprintf(' TRAINING TIME-DOMAIN WIENER WEIGHTS\n');
fprintf('------------------------------------------------------------\n');

% We train one scalar weight for each of the first G taps:
%
% theta_l =
% E[g_true(l) * conj(g_hat(l))]
% --------------------------------
% E[|g_hat(l)|^2]
%
% Then take the real component.

Rcross = zeros(G,1);
Rhat   = zeros(G,1);

for trial = 1:numTrain

    % Generate independent stochastic SUI-2 channel.
    h = sui2_channel(p);

    h = h(:);

    % Pad to 256 samples.
    h256 = zeros(N,1);
    h256(1:min(length(h),N)) = ...
        h(1:min(length(h),N));

    % True CFR.
    Htrue256 = fft(h256);

    % True 100-tone CFR.
    Htrue100 = Htrue256(idx256);

    % Build ideal 101-tone vector.
    Htrue101 = build_H101(Htrue100);

    % Apply Hann window.
    Hwindowed = hannWindow .* Htrue101;

    % WDFT transform.
    gCentered = Gmat * Hwindowed;

    % Put physical channel taps at beginning.
    gFull = fftshift(gCentered);

    % First G taps.
    gHat = gFull(1:G);

    % True CIR, padded to G.
    gTrue = zeros(G,1);

    Ltrue = min(length(h),G);
    gTrue(1:Ltrue) = h(1:Ltrue);

    % Accumulate Wiener statistics.
    Rcross = Rcross + ...
        gTrue .* conj(gHat);

    Rhat = Rhat + ...
        abs(gHat).^2;

end

% Average implicitly cancels because it appears in numerator
% and denominator.

theta = real(Rcross ./ max(Rhat,eps));

% Avoid negative scalar weights.
theta(theta < 0) = 0;

fprintf('Training complete.\n');
fprintf('\n');
fprintf('Wiener weight statistics:\n');
fprintf('  Minimum = %.12g\n', min(theta));
fprintf('  Maximum = %.12g\n', max(theta));
fprintf('  Mean    = %.12g\n', mean(theta));
fprintf('\n');


%% ============================================================
% 10. TEST ACTUAL WDFT ESTIMATOR
% =============================================================

numSNR = length(SNRdB);

mseWDFT  = zeros(numSNR,1);
nmseWDFT = zeros(numSNR,1);

mseRAW  = zeros(numSNR,1);
nmseRAW = zeros(numSNR,1);

cirNMSE_WDFT = zeros(numSNR,1);
cirNMSE_RAW  = zeros(numSNR,1);

for isnr = 1:numSNR

    snrNow = SNRdB(isnr);

    fprintf('Testing SNR = %d dB\n', snrNow);

    errW = 0;
    errR = 0;

    powerTrue = 0;

    cirErrW = 0;
    cirErrR = 0;
    cirPower = 0;

    for trial = 1:numTest

        %% ----------------------------------------------------
        % Generate stochastic SUI-2 channel
        % -----------------------------------------------------

        h = sui2_channel(p);
        h = h(:);

        h256 = zeros(N,1);
        h256(1:min(length(h),N)) = ...
            h(1:min(length(h),N));

        %% ----------------------------------------------------
        % True CFR
        % -----------------------------------------------------

        Htrue256 = fft(h256);

        Htrue100 = Htrue256(idx256);

        %% ----------------------------------------------------
        % Add CFR-domain complex AWGN
        %
        % This isolates the WDFT estimator behavior.
        %
        % Signal power is calculated from the true 100-tone CFR.
        % -----------------------------------------------------

        signalPower = mean(abs(Htrue100).^2);

        noisePower = signalPower / ...
                     (10^(snrNow/10));

        noise = sqrt(noisePower/2) * ...
            (randn(size(Htrue100)) + ...
             1j*randn(size(Htrue100)));

        Hnoisy100 = Htrue100 + noise;

        %% ----------------------------------------------------
        % Build 101-tone vector with DC substitution
        % -----------------------------------------------------

        Hnoisy101 = build_H101(Hnoisy100);

        %% ----------------------------------------------------
        % Hann window
        % -----------------------------------------------------

        Hwindowed = hannWindow .* Hnoisy101;

        %% ----------------------------------------------------
        % WDFT transform
        % -----------------------------------------------------

        gCentered = Gmat * Hwindowed;

        gFull = fftshift(gCentered);

        %% ----------------------------------------------------
        % RAW WDFT CIR
        % -----------------------------------------------------

        gRaw = gFull(1:G);

        %% ----------------------------------------------------
        % Wiener-weighted CIR
        % -----------------------------------------------------

        gWeighted = theta .* gRaw;

        %% ----------------------------------------------------
        % Build true CFR for comparison
        % -----------------------------------------------------

        Htrue100_col = Htrue100(:);

        %% ----------------------------------------------------
        % CFR reconstruction
        %
        % IMPORTANT:
        %
        % We reconstruct only through the same physical
        % frequency-domain basis used by the WDFT transform.
        %
        % This is NOT an arbitrary matrix inverse.
        %
        % The reconstruction here is the conventional
        % time-domain tap -> CFR relation.
        %
        % k values are physical signed subcarrier indices.
        %
        % H(k) = sum_l g(l) exp(-j2*pi*k*l/N)
        % -----------------------------------------------------

        Hraw = cir_to_cfr(gRaw, evenK_noDC, N);

        Hweighted = cir_to_cfr( ...
            gWeighted, evenK_noDC, N);

        %% ----------------------------------------------------
        % CFR error
        % -----------------------------------------------------

        eRaw = Hraw - Htrue100_col;

        eWeighted = ...
            Hweighted - Htrue100_col;

        errR = errR + sum(abs(eRaw).^2);

        errW = errW + sum(abs(eWeighted).^2);

        powerTrue = powerTrue + ...
            sum(abs(Htrue100_col).^2);

        %% ----------------------------------------------------
        % CIR error
        % -----------------------------------------------------

        gTrue = zeros(G,1);

        Ltrue = min(length(h),G);
        gTrue(1:Ltrue) = h(1:Ltrue);

        cirErrR = cirErrR + ...
            sum(abs(gRaw-gTrue).^2);

        cirErrW = cirErrW + ...
            sum(abs(gWeighted-gTrue).^2);

        cirPower = cirPower + ...
            sum(abs(gTrue).^2);

    end

    %% --------------------------------------------------------
    % Average MSE / NMSE
    % ---------------------------------------------------------

    mseRAW(isnr) = ...
        errR/(numTest*100);

    mseWDFT(isnr) = ...
        errW/(numTest*100);

    nmseRAW(isnr) = ...
        errR/max(powerTrue,eps);

    nmseWDFT(isnr) = ...
        errW/max(powerTrue,eps);

    cirNMSE_RAW(isnr) = ...
        cirErrR/max(cirPower,eps);

    cirNMSE_WDFT(isnr) = ...
        cirErrW/max(cirPower,eps);

end


%% ============================================================
% 11. DISPLAY RESULTS
% =============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' WDFT CFR ESTIMATION RESULTS\n');
fprintf('============================================================\n');

fprintf('\n');
fprintf('   SNR(dB)       RAW NMSE        WIENER NMSE\n');
fprintf('------------------------------------------------------------\n');

for i = 1:numSNR
    fprintf('%10d     %12.6e     %12.6e\n', ...
        SNRdB(i), ...
        nmseRAW(i), ...
        nmseWDFT(i));
end


fprintf('\n');
fprintf('============================================================\n');
fprintf(' WDFT CIR ESTIMATION RESULTS\n');
fprintf('============================================================\n');

fprintf('\n');
fprintf('   SNR(dB)       RAW CIR NMSE    WIENER CIR NMSE\n');
fprintf('------------------------------------------------------------\n');

for i = 1:numSNR
    fprintf('%10d     %12.6e     %12.6e\n', ...
        SNRdB(i), ...
        cirNMSE_RAW(i), ...
        cirNMSE_WDFT(i));
end


%% ============================================================
% 12. HIGH-SNR CHECK
% =============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' HIGH-SNR CHECK\n');
fprintf('============================================================\n');

highIdx = find(SNRdB >= 30);

if ~isempty(highIdx)

    rawHigh = nmseRAW(highIdx);
    wHigh   = nmseWDFT(highIdx);

    fprintf('Raw WDFT high-SNR NMSE range:\n');
    fprintf('  min = %.6e\n', min(rawHigh));
    fprintf('  max = %.6e\n', max(rawHigh));

    fprintf('Wiener WDFT high-SNR NMSE range:\n');
    fprintf('  min = %.6e\n', min(wHigh));
    fprintf('  max = %.6e\n', max(wHigh));

end


%% ============================================================
% 13. PLOT CFR NMSE
% =============================================================

figure;

semilogy(SNRdB,nmseRAW,'o-','LineWidth',1.5);
hold on;
semilogy(SNRdB,nmseWDFT,'s-','LineWidth',1.5);

grid on;

xlabel('SNR (dB)');
ylabel('CFR NMSE');

title('Step 23-J: WDFT CFR Estimation');

legend( ...
    'Raw WDFT', ...
    'Wiener-weighted WDFT', ...
    'Location','southwest');


%% ============================================================
% 14. PLOT CIR NMSE
% =============================================================

figure;

semilogy(SNRdB,cirNMSE_RAW,'o-','LineWidth',1.5);
hold on;
semilogy(SNRdB,cirNMSE_WDFT,'s-','LineWidth',1.5);

grid on;

xlabel('SNR (dB)');
ylabel('CIR NMSE');

title('Step 23-J: WDFT CIR Estimation');

legend( ...
    'Raw WDFT', ...
    'Wiener-weighted WDFT', ...
    'Location','southwest');


%% ============================================================
% 15. SAVE RESULTS
% =============================================================

step23J = struct();

step23J.SNRdB = SNRdB;

step23J.theta = theta;

step23J.mseRAW = mseRAW;
step23J.nmseRAW = nmseRAW;

step23J.mseWDFT = mseWDFT;
step23J.nmseWDFT = nmseWDFT;

step23J.cirNMSE_RAW = cirNMSE_RAW;
step23J.cirNMSE_WDFT = cirNMSE_WDFT;

step23J.Gmat = Gmat;
step23J.hannWindow = hannWindow;

step23J.numTrain = numTrain;
step23J.numTest = numTest;

save('step23J_actual_wdft_audit.mat','step23J');

fprintf('\n');
fprintf('Results saved to:\n');
fprintf('step23J_actual_wdft_audit.mat\n');


%% ============================================================
% 16. FINAL STRUCTURAL CHECKS
% =============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' STRUCTURAL CHECKS\n');
fprintf('============================================================\n');

pass = true;

if size(Gmat,1) ~= 128 || size(Gmat,2) ~= 101
    fprintf('FAIL: Eq. (23) matrix dimensions.\n');
    pass = false;
else
    fprintf('PASS: Eq. (23) matrix dimensions = 128 x 101.\n');
end

if length(hannWindow) ~= 101
    fprintf('FAIL: Hann window length.\n');
    pass = false;
else
    fprintf('PASS: Hann window length = 101.\n');
end

if numFrontZeros ~= 14 || numBackZeros ~= 13
    fprintf('FAIL: zero-padding structure.\n');
    pass = false;
else
    fprintf('PASS: zero-padding structure = 14 / 13.\n');
end

if length(theta) ~= G
    fprintf('FAIL: Wiener weight length.\n');
    pass = false;
else
    fprintf('PASS: Wiener weight length = %d.\n',G);
end

if any(~isfinite(theta))
    fprintf('FAIL: non-finite Wiener weights.\n');
    pass = false;
else
    fprintf('PASS: Wiener weights are finite.\n');
end

fprintf('\n');

if pass
    fprintf('STRUCTURAL AUDIT = PASS\n');
else
    fprintf('STRUCTURAL AUDIT = FAIL\n');
end


%% ============================================================
% LOCAL FUNCTION 1
% BUILD 101-TONE VECTOR
% =============================================================

function H101 = build_H101(H100)

    H100 = H100(:);

    if length(H100) ~= 100
        error('build_H101 requires exactly 100 tones.');
    end

    % 50 negative-frequency tones
    Hnegative = H100(1:50);

    % 50 positive-frequency tones
    Hpositive = H100(51:100);

    % DC substitution:
    %
    % average of nearest available tones:
    %
    % H_DC = (H_-2 + H_+2)/2
    %
    Hdc = 0.5*(Hnegative(end) + Hpositive(1));

    H101 = [
        Hnegative;
        Hdc;
        Hpositive
    ];

end


%% ============================================================
% LOCAL FUNCTION 2
% CIR -> 100-TONE CFR
% =============================================================

function H = cir_to_cfr(g, k, N)

    g = g(:);
    k = k(:);

    H = zeros(length(k),1);

    for r = 1:length(k)

        kr = k(r);

        for q = 1:length(g)

            tap = q-1;

            H(r) = H(r) + ...
                g(q) * ...
                exp(-1j*2*pi*kr*tap/N);

        end

    end

end