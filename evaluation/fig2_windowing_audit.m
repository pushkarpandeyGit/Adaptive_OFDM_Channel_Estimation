%% FIGURE 2 - FINAL WINDOWING AUDIT
%
% Purpose:
%   Numerically verify the physical claim behind Fig. 2:
%
%       Rectangular-window leakage  ~ noise variance
%       Hann-window leakage        << noise variance
% %
% This is an INDEPENDENT AUDIT.
%
% It does NOT modify:
%   - estimators/ls_to_cir.m
%   - transmitter/preamble_sequences.m
%   - any SVD/SHT estimator
%
% Paper:
%   H.-R. Park,
%   "A Low-Complexity Channel Estimation for OFDM Systems
%    Based on CIR Length Adaptation,"
%   IEEE Access, 2022.
%
% Conditions:
%   N       = 256
%   M       = 128
%   Np      = 100
%   SNR     = 20 dB
%
% The audit calculates:
%
%   1. LS-domain noise variance
%   2. post-window CIR-domain noise variance
%   3. clean leakage energy
%   4. leakage/noise ratio
%   5. noisy leakage/noise ratio
%
% for:
%
%   SUI-2
%   SUI-5
%
% and:
%
%   Rectangular window
%   Hann window

clear;
clc;
close all;

fprintf('\n');
fprintf('====================================================\n');
fprintf(' FIGURE 2 - FINAL WINDOWING AUDIT\n');
fprintf('====================================================\n');

%% =========================================================
% 1. PARAMETERS
% ==========================================================

p = parameters();

N      = p.N;
M      = 128;
SNR_dB = 20;

fprintf('\nN       = %d\n', N);
fprintf('M       = %d\n', M);
fprintf('SNR     = %.1f dB\n', SNR_dB);

%% =========================================================
% 2. IEEE 802.16 LTS
% ==========================================================

Pall = ieee80216_pall();

preamble = preamble_sequences(Pall);

PLTS = preamble.PLTS(:);

X100 = PLTS(PLTS ~= 0);

if length(X100) ~= 100
    error( ...
        'Expected 100 LTS tones. Found %d.', ...
        length(X100));
end

% LTS subcarrier locations.
kLTS = [-100:2:-2 2:2:100].';

fprintf('\nNumber of LTS tones = %d\n', length(X100));

%% =========================================================
% 3. CHANNELS
% ==========================================================

rng(2026);

% SUI-2
H2 = local_sui2_cfr(p, kLTS, N);

% SUI-5
h5 = sui5_channel(p);

H5 = local_cfr_from_cir(h5, kLTS, N);

%% =========================================================
% 4. NOISY LS ESTIMATES
% ==========================================================

rng(2027);

[H2_noisy, sigmaLS2] = ...
    local_noisy_ls(H2, X100, SNR_dB);

[H5_noisy, sigmaLS5] = ...
    local_noisy_ls(H5, X100, SNR_dB);

fprintf('\nLS-domain noise variance:\n');
fprintf('  SUI-2 = %.6e\n', sigmaLS2);
fprintf('  SUI-5 = %.6e\n', sigmaLS5);

%% =========================================================
% 5. TRANSFORM FOUR CASES
% ==========================================================

% ----------------------------------------------------------
% SUI-2
% ----------------------------------------------------------

g2_rect_clean = ...
    local_fig2_transform(H2, M, 'rectangular');

g2_hann_clean = ...
    local_fig2_transform(H2, M, 'hann');

g2_rect_noisy = ...
    local_fig2_transform(H2_noisy, M, 'rectangular');

g2_hann_noisy = ...
    local_fig2_transform(H2_noisy, M, 'hann');

% ----------------------------------------------------------
% SUI-5
% ----------------------------------------------------------

g5_rect_clean = ...
    local_fig2_transform(H5, M, 'rectangular');

g5_hann_clean = ...
    local_fig2_transform(H5, M, 'hann');

g5_rect_noisy = ...
    local_fig2_transform(H5_noisy, M, 'rectangular');

g5_hann_noisy = ...
    local_fig2_transform(H5_noisy, M, 'hann');

%% =========================================================
% 6. WINDOW NOISE VARIANCE
% ==========================================================
%
% After LS:
%
%       H_LS = H + n/X
%
% If the LS noise variance is sigma_LS^2,
% the windowed CIR-domain noise variance becomes
%
%       sigma_CIR^2 =
%           sigma_LS^2 * ||w||^2 / M
%
% because the transform matrix contains 1/sqrt(M).
%
% This is the quantity we compare leakage against.

q = (-50:50).';

wRect = ones(101,1);

wHann = ...
    0.5 + ...
    0.5*cos(2*pi*q/M);

sigmaCIR2_rect = ...
    sigmaLS2 * sum(abs(wRect).^2) / M;

sigmaCIR2_hann = ...
    sigmaLS2 * sum(abs(wHann).^2) / M;

sigmaCIR5_rect = ...
    sigmaLS5 * sum(abs(wRect).^2) / M;

sigmaCIR5_hann = ...
    sigmaLS5 * sum(abs(wHann).^2) / M;

fprintf('\nPost-window CIR-domain noise variance:\n');

fprintf('\nSUI-2:\n');
fprintf('  Rectangular = %.6e\n', sigmaCIR2_rect);
fprintf('  Hann        = %.6e\n', sigmaCIR2_hann);

fprintf('\nSUI-5:\n');
fprintf('  Rectangular = %.6e\n', sigmaCIR5_rect);
fprintf('  Hann        = %.6e\n', sigmaCIR5_hann);

%% =========================================================
% 7. DEFINE PHYSICAL CIR REGIONS
% ==========================================================
%
% SUI-2:
%
%   delays = 0, 1.6, 4.4 samples
%
% We regard taps above approximately 4.4 samples as
% outside the physical CIR.
%
% Therefore use:
%
%       taps 0:5   = physical-channel region
%       taps 6:127 = leakage region
%
% SUI-5:
%
%       delays = 0, 16, 40 samples
%
% Physical region:
%
%       taps 0:40
%
% Leakage:
%
%       taps 41:127

tap = 0:M-1;

physical2 = tap <= 5;
leakage2  = tap > 5;

physical5 = tap <= 40;
leakage5  = tap > 40;

%% =========================================================
% 8. CLEAN LEAKAGE ENERGY
% ==========================================================

% ----------------------------------------------------------
% SUI-2
% ----------------------------------------------------------

rect2_clean_leak_mean = ...
    mean(abs(g2_rect_clean(leakage2)).^2);

hann2_clean_leak_mean = ...
    mean(abs(g2_hann_clean(leakage2)).^2);

rect2_clean_leak_max = ...
    max(abs(g2_rect_clean(leakage2)).^2);

hann2_clean_leak_max = ...
    max(abs(g2_hann_clean(leakage2)).^2);

% ----------------------------------------------------------
% SUI-5
% ----------------------------------------------------------

rect5_clean_leak_mean = ...
    mean(abs(g5_rect_clean(leakage5)).^2);

hann5_clean_leak_mean = ...
    mean(abs(g5_hann_clean(leakage5)).^2);

rect5_clean_leak_max = ...
    max(abs(g5_rect_clean(leakage5)).^2);

hann5_clean_leak_max = ...
    max(abs(g5_hann_clean(leakage5)).^2);

%% =========================================================
% 9. LEAKAGE / NOISE RATIOS
% ==========================================================

ratio2_rect_mean = ...
    rect2_clean_leak_mean / sigmaCIR2_rect;

ratio2_hann_mean = ...
    hann2_clean_leak_mean / sigmaCIR2_hann;

ratio2_rect_max = ...
    rect2_clean_leak_max / sigmaCIR2_rect;

ratio2_hann_max = ...
    hann2_clean_leak_max / sigmaCIR2_hann;

ratio5_rect_mean = ...
    rect5_clean_leak_mean / sigmaCIR5_rect;

ratio5_hann_mean = ...
    hann5_clean_leak_mean / sigmaCIR5_hann;

ratio5_rect_max = ...
    rect5_clean_leak_max / sigmaCIR5_rect;

ratio5_hann_max = ...
    hann5_clean_leak_max / sigmaCIR5_hann;

%% =========================================================
% 10. PRINT CLEAN LEAKAGE RESULTS
% ==========================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' CLEAN CHANNEL LEAKAGE\n');
fprintf('====================================================\n');

fprintf('\nSUI-2:\n');

fprintf('  Rectangular mean leakage = %.6e\n', ...
    rect2_clean_leak_mean);

fprintf('  Hann mean leakage        = %.6e\n', ...
    hann2_clean_leak_mean);

fprintf('  Rectangular max leakage  = %.6e\n', ...
    rect2_clean_leak_max);

fprintf('  Hann max leakage         = %.6e\n', ...
    hann2_clean_leak_max);

fprintf('\n  Mean leakage / noise variance:\n');

fprintf('    Rectangular = %.4f  (%.2f dB)\n', ...
    ratio2_rect_mean, ...
    10*log10(ratio2_rect_mean));

fprintf('    Hann        = %.4f  (%.2f dB)\n', ...
    ratio2_hann_mean, ...
    10*log10(ratio2_hann_mean));

fprintf('\n  Max leakage / noise variance:\n');

fprintf('    Rectangular = %.4f  (%.2f dB)\n', ...
    ratio2_rect_max, ...
    10*log10(ratio2_rect_max));

fprintf('    Hann        = %.4f  (%.2f dB)\n', ...
    ratio2_hann_max, ...
    10*log10(ratio2_hann_max));

fprintf('\nSUI-5:\n');

fprintf('  Rectangular mean leakage = %.6e\n', ...
    rect5_clean_leak_mean);

fprintf('  Hann mean leakage        = %.6e\n', ...
    hann5_clean_leak_mean);

fprintf('  Rectangular max leakage  = %.6e\n', ...
    rect5_clean_leak_max);

fprintf('  Hann max leakage         = %.6e\n', ...
    hann5_clean_leak_max);

fprintf('\n  Mean leakage / noise variance:\n');

fprintf('    Rectangular = %.4f  (%.2f dB)\n', ...
    ratio5_rect_mean, ...
    10*log10(ratio5_rect_mean));

fprintf('    Hann        = %.4f  (%.2f dB)\n', ...
    ratio5_hann_mean, ...
    10*log10(ratio5_hann_mean));

fprintf('\n  Max leakage / noise variance:\n');

fprintf('    Rectangular = %.4f  (%.2f dB)\n', ...
    ratio5_rect_max, ...
    10*log10(ratio5_rect_max));

fprintf('    Hann        = %.4f  (%.2f dB)\n', ...
    ratio5_hann_max, ...
    10*log10(ratio5_hann_max));

%% =========================================================
% 11. NOISY LEAKAGE
% ==========================================================
%
% This is useful because the paper's Fig. 2 visually compares
% the leakage against the noise floor.
%
% We calculate the average noisy energy in the region where
% no physical path exists.

noisy2_rect = ...
    mean(abs(g2_rect_noisy(leakage2)).^2);

noisy2_hann = ...
    mean(abs(g2_hann_noisy(leakage2)).^2);

noisy5_rect = ...
    mean(abs(g5_rect_noisy(leakage5)).^2);

noisy5_hann = ...
    mean(abs(g5_hann_noisy(leakage5)).^2);

fprintf('\n');
fprintf('====================================================\n');
fprintf(' NOISY LEAKAGE REGION\n');
fprintf('====================================================\n');

fprintf('\nSUI-2:\n');

fprintf('  Rectangular = %.6e\n', noisy2_rect);
fprintf('  Hann        = %.6e\n', noisy2_hann);

fprintf('  Rectangular / noise = %.4f (%.2f dB)\n', ...
    noisy2_rect/sigmaCIR2_rect, ...
    10*log10(noisy2_rect/sigmaCIR2_rect));

fprintf('  Hann / noise        = %.4f (%.2f dB)\n', ...
    noisy2_hann/sigmaCIR2_hann, ...
    10*log10(noisy2_hann/sigmaCIR2_hann));

fprintf('\nSUI-5:\n');

fprintf('  Rectangular = %.6e\n', noisy5_rect);
fprintf('  Hann        = %.6e\n', noisy5_hann);

fprintf('  Rectangular / noise = %.4f (%.2f dB)\n', ...
    noisy5_rect/sigmaCIR5_rect, ...
    10*log10(noisy5_rect/sigmaCIR5_rect));

fprintf('  Hann / noise        = %.4f (%.2f dB)\n', ...
    noisy5_hann/sigmaCIR5_hann, ...
    10*log10(noisy5_hann/sigmaCIR5_hann));

%% =========================================================
% 12. WINDOW ATTENUATION
% ==========================================================

windowReduction2 = ...
    10*log10( ...
    ratio2_rect_mean / ...
    ratio2_hann_mean);

windowReduction5 = ...
    10*log10( ...
    ratio5_rect_mean / ...
    ratio5_hann_mean);

fprintf('\n');
fprintf('====================================================\n');
fprintf(' HANN LEAKAGE REDUCTION\n');
fprintf('====================================================\n');

fprintf('\nSUI-2:\n');
fprintf('  Mean leakage reduction = %.2f dB\n', ...
    windowReduction2);

fprintf('\nSUI-5:\n');
fprintf('  Mean leakage reduction = %.2f dB\n', ...
    windowReduction5);

%% =========================================================
% 13. BASIC INTERPRETATION
% ==========================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' INTERPRETATION\n');
fprintf('====================================================\n');

fprintf('\nA ratio of:\n');
fprintf('  1.0  = leakage equals noise variance\n');
fprintf('  <1   = leakage is below noise variance\n');
fprintf('  >1   = leakage is above noise variance\n');

fprintf('\nThe paper''s qualitative claim is:\n');
fprintf('  Rectangular: leakage comparable to noise\n');
fprintf('  Hann:        leakage significantly below noise\n');

fprintf('\nDO NOT judge this audit from one random realization alone.\n');
fprintf('The important quantity is the leakage/noise relationship.\n');

fprintf('\n====================================================\n');
fprintf(' AUDIT COMPLETE\n');
fprintf('====================================================\n');


%% ========================================================================
% LOCAL FUNCTION
% SUI-2 CFR
% ========================================================================

function H = local_sui2_cfr(p, k, N)

    delay_us = [0 0.4 1.1];

    power_dB = [0 -12 -15];

    delay_samples = ...
        delay_us * 1e-6 / p.Ts;

    K_dB = 3;

    if isfield(p.SUI2,'RiceanK_dB')
        K_dB = p.SUI2.RiceanK_dB;

    elseif isfield(p.SUI2,'K_dB')
        K_dB = p.SUI2.K_dB;
    end

    K = 10^(K_dB/10);

    alpha = zeros(1,3);

    % Ricean first path.

    los = sqrt(K/(K+1));

    scattered = ...
        (randn + 1j*randn)/sqrt(2);

    scattered = ...
        sqrt(1/(K+1))*scattered;

    alpha(1) = ...
        10^(power_dB(1)/20) * ...
        (los + scattered);

    % Rayleigh remaining paths.

    for ii = 2:3

        rayleigh = ...
            (randn + 1j*randn)/sqrt(2);

        alpha(ii) = ...
            10^(power_dB(ii)/20) * ...
            rayleigh;

    end

    H = zeros(size(k));

    for ii = 1:3

        H = H + ...
            alpha(ii) .* ...
            exp( ...
            -1j*2*pi*k*delay_samples(ii)/N);

    end

end


%% ========================================================================
% LOCAL FUNCTION
% CFR FROM CIR
% ========================================================================

function H = local_cfr_from_cir(h, k, N)

    h = h(:).';

    n = 0:length(h)-1;

    H = zeros(size(k));

    for ii = 1:length(k)

        H(ii) = sum( ...
            h .* ...
            exp( ...
            -1j*2*pi*k(ii)*n/N));

    end

end


%% ========================================================================
% LOCAL FUNCTION
% NOISY LS
% ========================================================================

function [H_LS, sigmaLS] = ...
    local_noisy_ls(H, X, SNR_dB)

    H = H(:);
    X = X(:);

    Y_clean = H .* X;

    signalPower = ...
        mean(abs(Y_clean).^2);

    noisePower = ...
        signalPower / ...
        10^(SNR_dB/10);

    noise = ...
        sqrt(noisePower/2) .* ...
        (randn(size(Y_clean)) + ...
         1j*randn(size(Y_clean)));

    Y = Y_clean + noise;

    H_LS = Y ./ X;

    % LS-domain noise variance.
    nLS = noise ./ X;

    sigmaLS = ...
        mean(abs(nLS).^2);

end


%% ========================================================================
% LOCAL FUNCTION
% FIG. 2 TRANSFORM
% ========================================================================

function g = ...
    local_fig2_transform(H_LS, M, windowType)

    H_LS = H_LS(:);

    if length(H_LS) ~= 100
        error( ...
            'Expected 100 LTS channel estimates.');
    end

    % ---------------------------------------------------------
    % 101-point centered vector.
    %
    % 50 negative-frequency tones
    % 1 DC zero
    % 50 positive-frequency tones
    % ---------------------------------------------------------

    p101 = zeros(101,1);

    p101(1:50) = H_LS(1:50);

    p101(51) = 0;

    p101(52:101) = H_LS(51:100);

    % ---------------------------------------------------------
    % Window
    % ---------------------------------------------------------

    n = (-50:50).';

    switch lower(windowType)

        case 'rectangular'

            w = ones(101,1);

        case 'hann'

            w = ...
                0.5 + ...
                0.5*cos(2*pi*n/M);

        otherwise

            error( ...
                'Unknown window: %s', ...
                windowType);

    end

    p_windowed = ...
        p101 .* w;

    % ---------------------------------------------------------
    % Paper transform.
    % ---------------------------------------------------------

    m = (0:M-1).';

    q = (-50:50);

    G = ...
        (1/sqrt(M)) .* ...
        exp( ...
        1j*2*pi .* ...
        ((m-M/2)*q)/M);

    % CIR.

    g = G * p_windowed;

    % Convert centered coordinate to ordinary delay index.
    g = fftshift(g);

    g = g(:).';

end