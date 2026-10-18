%% FIGURE 2 - EFFECT OF WINDOWING ON CIR LENGTH ESTIMATION
%
% Reproduction of Fig. 2 from:
%
% H.-R. Park,
% "A Low-Complexity Channel Estimation for OFDM Systems
%  Based on CIR Length Adaptation,"
% IEEE Access, 2022.
%
% IMPORTANT:
%   This script is intentionally self-contained.
%
%   DO NOT MODIFY:
%       estimators/ls_to_cir.m
%       transmitter/preamble_sequences.m
%       any other estimator/channel file
%
% Paper conditions:
%       SNR = 20 dB
%       SUI-2
%       SUI-5
%       M   = 128
%       Np  = 100 LTS training tones
%
% Fig. 2 curves:
%       1. Rectangular window
%       2. Hann window
%       3. Rectangular window (noise-free)
%       4. Hann window (noise-free)

clear;
clc;
close all;

fprintf('\n');
fprintf('====================================================\n');
fprintf(' FIGURE 2 - EFFECT OF WINDOWING\n');
fprintf('====================================================\n');

%% =========================================================
% 1. PROJECT PARAMETERS
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
%
% Paper:
%
% P_LTS,k = sqrt(2) P_all,k
%
% for even k, with 100 non-zero LTS tones.
%
% The existing project functions are used exactly as they are.

Pall = ieee80216_pall();

preamble = preamble_sequences(Pall);

PLTS_full = preamble.PLTS(:);

% Extract the 100 actual LTS tones.
X100 = PLTS_full(PLTS_full ~= 0);

if length(X100) ~= 100
    error('Expected 100 LTS tones, found %d.', length(X100));
end

% Corresponding subcarrier indices.
kLTS = [-100:2:-2 2:2:100].';

fprintf('\nLTS tones = %d\n', length(X100));

%% =========================================================
% 3. CHANNEL GENERATION
% ==========================================================

% Fixed seeds make this reproduction repeatable.
rng(2026);

% ----------------------------------------------------------
% SUI-2
% ----------------------------------------------------------
%
% Paper:
%   delays = 0, 0.4, 1.1 us
%   powers = 0, -12, -15 dB
%   first path Ricean K = 3 dB
%   other paths Rayleigh
%
H_SUI2 = local_sui2_cfr(p, kLTS, N);

% ----------------------------------------------------------
% SUI-5
% ----------------------------------------------------------
%
% Paper:
%   delays = 0, 4, 10 us
%   powers = 0, -5, -10 dB
%   all paths Rayleigh
%
h_SUI5 = sui5_channel(p);

H_SUI5 = local_cfr_from_cir(h_SUI5, kLTS, N);

%% =========================================================
% 4. GENERATE LS CHANNEL ESTIMATES
% ==========================================================
%
% Received LTS:
%
%       Y = H X + W
%
% LS estimate:
%
%       H_LS = Y/X
%
% We generate:
%
%   - noisy LS estimate
%   - noise-free LS estimate
%
% The paper compares the noisy and noise-free cases.

rng(2027);

[HLS2_noisy, noiseVar2] = ...
    local_noisy_ls(H_SUI2, X100, SNR_dB);

[HLS5_noisy, noiseVar5] = ...
    local_noisy_ls(H_SUI5, X100, SNR_dB);

HLS2_clean = H_SUI2;
HLS5_clean = H_SUI5;

%% =========================================================
% 5. FIGURE 2 TRANSFORMATION
% ==========================================================
%
% Paper Eq. (21):
%
%   w_n = 0.5 + 0.5 cos(2*pi*n/M)
%
% with:
%
%   |n| <= Np/2
%
% For Np = 100:
%
%   n = -50,...,+50
%
% giving 101 positions including the DC position.
%
% Paper Eq. (23):
%
%   G(m,n) =
%       1/sqrt(M) *
%       exp(j*2*pi*(m-M/2)*(n-Np/2)/M)
%
% m = 0,...,M-1
% n = 0,...,Np
%
% Here Np = 100, hence 101 columns.

%% SUI-2

g2_rect_noisy = local_fig2_transform( ...
    HLS2_noisy, M, 'rectangular');

g2_hann_noisy = local_fig2_transform( ...
    HLS2_noisy, M, 'hann');

g2_rect_clean = local_fig2_transform( ...
    HLS2_clean, M, 'rectangular');

g2_hann_clean = local_fig2_transform( ...
    HLS2_clean, M, 'hann');

%% SUI-5

g5_rect_noisy = local_fig2_transform( ...
    HLS5_noisy, M, 'rectangular');

g5_hann_noisy = local_fig2_transform( ...
    HLS5_noisy, M, 'hann');

g5_rect_clean = local_fig2_transform( ...
    HLS5_clean, M, 'rectangular');

g5_hann_clean = local_fig2_transform( ...
    HLS5_clean, M, 'hann');

%% =========================================================
% 6. TAP ENERGY
% ==========================================================

epsPower = 1e-14;

% 10 log10 |g|^2
%
% This is equivalent to:
%
% 20 log10 |g|

E2_rect_noisy = ...
    10*log10( ...
    (abs(g2_rect_noisy).^2 + epsPower) / noiseVar2);

E2_hann_noisy = ...
    10*log10( ...
    (abs(g2_hann_noisy).^2 + epsPower) / noiseVar2);

E2_rect_clean = ...
    10*log10( ...
    (abs(g2_rect_clean).^2 + epsPower) / noiseVar2);

E2_hann_clean = ...
    10*log10( ...
    (abs(g2_hann_clean).^2 + epsPower) / noiseVar2);
E5_rect_noisy = ...
    10*log10( ...
    (abs(g5_rect_noisy).^2 + epsPower) / noiseVar5);

E5_hann_noisy = ...
    10*log10( ...
    (abs(g5_hann_noisy).^2 + epsPower) / noiseVar5);

E5_rect_clean = ...
    10*log10( ...
    (abs(g5_rect_clean).^2 + epsPower) / noiseVar5);

E5_hann_clean = ...
    10*log10( ...
    (abs(g5_hann_clean).^2 + epsPower) / noiseVar5);
tap = 0:M-1;

%% =========================================================
% 7. FIGURE 2(a) - SUI-2
% ==========================================================

figure( ...
    'Name','Fig. 2(a) - SUI-2', ...
    'Color','w');

% Noisy rectangular
plot(tap, E2_rect_noisy, ...
    'LineWidth',1.2);
hold on;

% Noisy Hann
plot(tap, E2_hann_noisy, ...
    'LineWidth',1.2);

% Noise-free rectangular
plot(tap, E2_rect_clean, ...
    '--.', ...
    'LineWidth',1.0, ...
    'MarkerSize',5);

% Noise-free Hann
plot(tap, E2_hann_clean, ...
    '-.', ...
    'LineWidth',1.0);

grid on;

xlabel('Time Index (T_s)');
ylabel('Tap Energy (dB)');

title( ...
    'Effect of Windowing on CIR Length Estimation - SUI-2');

legend( ...
    'Rectangular window', ...
    'Hann window', ...
    'Rectangular window (noise-free)', ...
    'Hann window (noise-free)', ...
    'Location','northeast');

xlim([0 M-1]);

% Paper's SUI-2 figure uses approximately this range.
ylim([-40 40]);

%% Save

plotDir = fileparts(mfilename('fullpath'));

saveas( ...
    gcf, ...
    fullfile(plotDir, ...
    'Fig2_SUI2_windowing_effect.png'));

%% =========================================================
% 8. FIGURE 2(b) - SUI-5
% ==========================================================

figure( ...
    'Name','Fig. 2(b) - SUI-5', ...
    'Color','w');

% Noisy rectangular
plot(tap, E5_rect_noisy, ...
    'LineWidth',1.2);
hold on;

% Noisy Hann
plot(tap, E5_hann_noisy, ...
    'LineWidth',1.2);

% Noise-free rectangular
plot(tap, E5_rect_clean, ...
    '--.', ...
    'LineWidth',1.0, ...
    'MarkerSize',5);

% Noise-free Hann
plot(tap, E5_hann_clean, ...
    '-.', ...
    'LineWidth',1.0);

grid on;

xlabel('Time Index (T_s)');
ylabel('Tap Energy (dB)');

title( ...
    'Effect of Windowing on CIR Length Estimation - SUI-5');

legend( ...
    'Rectangular window', ...
    'Hann window', ...
    'Rectangular window (noise-free)', ...
    'Hann window (noise-free)', ...
    'Location','northeast');

xlim([0 M-1]);

ylim([-50 50]);

%% Save

saveas( ...
    gcf, ...
    fullfile(plotDir, ...
    'Fig2_SUI5_windowing_effect.png'));

%% =========================================================
% 9. IMPORTANT PHYSICAL CHECKS
% ==========================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' FIGURE 2 CHECKS\n');
fprintf('====================================================\n');

fprintf('\nSUI-2:\n');
fprintf('  Delays (us)    : 0, 0.4, 1.1\n');
fprintf('  Delays (Ts)    : 0, 1.6, 4.4\n');
fprintf('  True CIR length: 4.4 samples\n');

fprintf('\nSUI-5:\n');
fprintf('  Delays (us)    : 0, 4, 10\n');
fprintf('  Delays (Ts)    : 0, 16, 40\n');
fprintf('  True CIR length: 40 samples\n');

fprintf('\nCentered representation:\n');
fprintf('  M/2            : %d\n', M/2);
fprintf('  SUI-5 tap 0    -> index 64\n');
fprintf('  SUI-5 tap 16   -> index 80\n');
fprintf('  SUI-5 tap 40   -> index 104\n');

fprintf('\nSNR               : %.1f dB\n', SNR_dB);
fprintf('M                 : %d\n', M);
fprintf('LTS tones         : %d\n', length(X100));

fprintf('\nSaved:\n');
fprintf('  Fig2_SUI2_windowing_effect.png\n');
fprintf('  Fig2_SUI5_windowing_effect.png\n');

fprintf('\n');
fprintf('====================================================\n');
fprintf(' DONE\n');
fprintf('====================================================\n');


%% ========================================================================
% LOCAL FUNCTION
% SUI-2 FREQUENCY RESPONSE
% ========================================================================

function H = local_sui2_cfr(p, k, N)

    % SUI-2 parameters from the paper.
    delay_us = [0 0.4 1.1];
    power_dB = [0 -12 -15];

    % Convert physical delays to samples.
    delay_samples = ...
        delay_us * 1e-6 / p.Ts;

    % Ricean K-factor.
    K_dB = 3;

    if isfield(p.SUI2,'RiceanK_dB')

        K_dB = p.SUI2.RiceanK_dB;

    elseif isfield(p.SUI2,'K_dB')

        K_dB = p.SUI2.K_dB;

    end

    K = 10^(K_dB/10);

    alpha = zeros(1,3);

    % ---------------------------------------------------------
    % First path: Ricean
    % ---------------------------------------------------------

    los = sqrt(K/(K+1));

    scattered = ...
        (randn + 1j*randn)/sqrt(2);

    scattered = ...
        sqrt(1/(K+1))*scattered;

    alpha(1) = ...
        10^(power_dB(1)/20) * ...
        (los + scattered);

    % ---------------------------------------------------------
    % Second and third paths: Rayleigh
    % ---------------------------------------------------------

    for ii = 2:3

        rayleigh = ...
            (randn + 1j*randn)/sqrt(2);

        alpha(ii) = ...
            10^(power_dB(ii)/20) * ...
            rayleigh;

    end

    % ---------------------------------------------------------
    % CFR
    %
    % H(k) =
    %   sum alpha_l exp(-j 2*pi*k*d_l/N)
    %
    % ---------------------------------------------------------

    H = zeros(size(k));

    for ii = 1:3

        H = H + ...
            alpha(ii) .* ...
            exp(-1j*2*pi*k*delay_samples(ii)/N);

    end

end


%% ========================================================================
% LOCAL FUNCTION
% SUI-5 CFR FROM EXISTING CHANNEL
% ========================================================================

function H = local_cfr_from_cir(h, k, N)

    h = h(:).';

    n = 0:length(h)-1;

    H = zeros(size(k));

    for ii = 1:length(k)

        H(ii) = sum( ...
            h .* ...
            exp(-1j*2*pi*k(ii)*n/N) ...
            );

    end

end


%% ========================================================================
% LOCAL FUNCTION
% NOISY LS ESTIMATE
% ========================================================================
function [H_LS, noiseVarianceLS] = local_noisy_ls(H, X, SNR_dB)

H = H(:);
X = X(:);

% Clean received LTS
Y_clean = H .* X;

% Average received signal power
signalPower = mean(abs(Y_clean).^2);

% AWGN variance
noisePower = signalPower / 10^(SNR_dB/10);

% Complex AWGN
noise = sqrt(noisePower/2) .* ...
    (randn(size(Y_clean)) + ...
    1j*randn(size(Y_clean)));

% Received LTS
Y = Y_clean + noise;

% LS estimate
H_LS = Y ./ X;

% Noise variance after LS division.
%
% Since H_LS = H + W/X,
% this is the variance of the LS-domain noise.
noiseVarianceLS = mean(abs(noise ./ X).^2);

end

%% ========================================================================
% LOCAL FUNCTION
% PAPER EQ. (21)-(23)
% ========================================================================

function g = local_fig2_transform(H_LS, M, windowType)

    H_LS = H_LS(:);

    if length(H_LS) ~= 100
        error( ...
            'Fig. 2 transformation requires exactly 100 LTS tones.');
    end

    % ---------------------------------------------------------
    % Construct the 101-point centered vector.
    %
    % 50 negative-frequency LTS tones
    % 1 DC zero
    % 50 positive-frequency LTS tones
    %
    % q = -50,...,+50
    % ---------------------------------------------------------

    p101 = zeros(101,1);

    p101(1:50) = H_LS(1:50);

    % DC:
    p101(51) = 0;

    p101(52:101) = H_LS(51:100);

    % ---------------------------------------------------------
    % Eq. (21)
    %
    % w_n = 0.5 + 0.5*cos(2*pi*n/M)
    %
    % n = -50,...,+50
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
                'Unknown window type: %s', ...
                windowType);

    end

    % Apply window.
    p_windowed = p101 .* w;

    % ---------------------------------------------------------
    % Eq. (23)
    %
    % G(m,n) =
    %
    % 1/sqrt(M) *
    % exp[j*2*pi*(m-M/2)*(n-Np/2)/M]
    %
    % m = 0,...,M-1
    % n = 0,...,100
    %
    % Since Np = 100:
    %
    % n - Np/2 = -50,...,+50
    % ---------------------------------------------------------

    m = (0:M-1).';

    q = (-50:50);

    G = ...
        (1/sqrt(M)) .* ...
        exp( ...
        1j*2*pi .* ...
        ((m - M/2) * q) / M);

    % ---------------------------------------------------------
    % CIR estimate
    % ---------------------------------------------------------

 g = G * p_windowed;

% Eq. (23) uses a centered CIR coordinate system.
% The paper's CIR-length index is the ordinary delay index,
% so circularly shift the centered result.
g = fftshift(g);

g = g(:).';

end