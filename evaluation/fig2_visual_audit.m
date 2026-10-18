%% FIGURE 2 - VISUAL / CIR LEAKAGE AUDIT
% Paper reproduction audit for:
% "A Low-Complexity Channel Estimation for OFDM Systems
%  Based on CIR Length Adaptation"
%
% Purpose:
%   1. Reproduce the paper's Eq. (21) Hann window
%   2. Reproduce the paper's Eq. (23) 128-point transform
%   3. Compare rectangular vs Hann windowing
%   4. Separate CLEAN and NOISY channel behavior
%   5. Focus on the actual CIR leakage shape
%
% IMPORTANT:
%   This script is self-contained.
%   It does NOT modify ls_to_cir.m or any shared estimator.

clear;
clc;
close all;

fprintf('\n');
fprintf('====================================================\n');
fprintf(' FIGURE 2 - VISUAL / CIR LEAKAGE AUDIT\n');
fprintf('====================================================\n\n');

%% ---------------------------------------------------
% 1. PROJECT PARAMETERS
% ----------------------------------------------------
p = parameters();

N   = 256;
M   = 128;
SNR = 20;

fprintf('N       = %d\n', N);
fprintf('M       = %d\n', M);
fprintf('SNR     = %.1f dB\n\n', SNR);

%% ---------------------------------------------------
% 2. IEEE 802.16 LTS
% ----------------------------------------------------
Pall = ieee80216_pall();
preamble = preamble_sequences(Pall);

PLTS = preamble.PLTS(:);

% Extract only the 100 non-zero LTS tones.
X100 = PLTS(PLTS ~= 0);

if length(X100) ~= 100
    error('Expected exactly 100 LTS tones.');
end

fprintf('Number of LTS tones = %d\n', length(X100));

% The 100 LTS tones correspond to:
%
% k = -100,-98,...,-2,+2,...,+100
%
% Therefore q = -50,...,-1,+1,...,+50
%
% The paper's Eq. (23) uses the centered 101-point representation
% including the DC position.

q = -50:50;

if length(q) ~= 101
    error('Internal q-axis error.');
end

%% ---------------------------------------------------
% 3. GENERATE SUI-2 CHANNEL
% ----------------------------------------------------
[h2, delay2_us, power2_dB] = local_sui2_channel(p);

fprintf('\nSUI-2 channel:\n');
fprintf('  Delays (us)    : ');
fprintf('%.4f ', delay2_us);
fprintf('\n');

fprintf('  Delays (samples): ');
fprintf('%g ', delay2_us * 1e-6 / p.Ts);
fprintf('\n');

fprintf('  Powers (dB)    : ');
fprintf('%.2f ', power2_dB);
fprintf('\n');

%% ---------------------------------------------------
% 4. GENERATE SUI-5 CHANNEL
% ----------------------------------------------------
h5 = sui5_channel(p);

delay5_us = p.SUI5.pathDelay_us;
delay5_samples = delay5_us * 1e-6 / p.Ts;

fprintf('\nSUI-5 channel:\n');
fprintf('  Delays (us)     : ');
fprintf('%.4f ', delay5_us);
fprintf('\n');

fprintf('  Delays (samples): ');
fprintf('%g ', delay5_samples);
fprintf('\n');

fprintf('  Powers (dB)     : ');
fprintf('%.2f ', p.SUI5.pathPower_dB);
fprintf('\n');

%% ---------------------------------------------------
% 5. TRUE FREQUENCY RESPONSE
% ----------------------------------------------------
%
% Actual OFDM subcarrier:
%
% H(k) = sum_l h(l) exp(-j 2*pi*k*l/N)
%
% LTS tones:
%
% k = 2q
%
% because the LTS occupies every second subcarrier.

k100 = 2 * q(q ~= 0);

if length(k100) ~= 100
    error('Expected 100 physical LTS subcarrier indices.');
end

H2_100 = local_cfr(h2, k100, N);
H5_100 = local_cfr(h5, k100, N);

%% ---------------------------------------------------
% 6. INSERT DC TO FORM 101-POINT CENTERED VECTOR
% ----------------------------------------------------
%
% q = -50,...,-1,0,+1,...,+50
%
% DC is placed at the center and is zero.
%
% This matches the centered representation used by Eq. (23).

H2_101 = zeros(101,1);
H5_101 = zeros(101,1);

H2_101(1:50)  = H2_100(1:50);
H2_101(52:101) = H2_100(51:100);

H5_101(1:50)  = H5_100(1:50);
H5_101(52:101) = H5_100(51:100);

%% ---------------------------------------------------
% 7. BUILD EQ. (23) TRANSFORM
% ----------------------------------------------------
%
% Paper:
%
% [G]_(m,n) =
% 1/sqrt(M) *
% exp(j*2*pi*(m-M/2)*(n-Np/2)/M)
%
% m = 0,...,M-1
% n = 0,...,Np
%
% Here Np = 100 LTS tones, so the 101-point vector is
% represented using the centered q = -50,...,+50 grid.

G = zeros(M,101);

for m = 0:M-1
    for n = 0:100

        G(m+1,n+1) = ...
            (1/sqrt(M)) * ...
            exp(1j * 2*pi * ...
            (m - M/2) * ...
            (n - 50) / M);

    end
end

%% ---------------------------------------------------
% 8. WINDOWS
% ----------------------------------------------------
%
% Rectangular:
%       w(n) = 1
%
% Hann, paper Eq. (21):
%
%       w(n) = 0.5 + 0.5*cos(2*pi*n/M)
%
% where n = -50,...,+50.

w_rect = ones(101,1);

w_hann = 0.5 + 0.5*cos(2*pi*q(:)/M);

%% ---------------------------------------------------
% 9. TRANSFORM FUNCTION
% ----------------------------------------------------
%
% Apply:
%
%       h_tilde = G * (H .* w)
%
% Then fftshift so the ordinary delay axis is:
%
%       0,1,...,127
%
% rather than the centered transform ordering.

h2_rect_clean = fftshift(G * (H2_101 .* w_rect));
h2_hann_clean = fftshift(G * (H2_101 .* w_hann));

h5_rect_clean = fftshift(G * (H5_101 .* w_rect));
h5_hann_clean = fftshift(G * (H5_101 .* w_hann));

%% ---------------------------------------------------
% 10. ADD 20-dB AWGN IN RECEIVED PILOT DOMAIN
% ----------------------------------------------------
%
% Y = X*H + n
%
% LS estimate:
%
% H_LS = Y/X
%
% Therefore the frequency-domain LS noise is naturally:
%
% n/X
%
% This avoids inserting an artificial constant LS-noise value.

[H2_rect_noisy, H2_hann_noisy] = ...
    local_noisy_estimates(H2_100, X100, SNR, G, w_rect, w_hann);

[H5_rect_noisy, H5_hann_noisy] = ...
    local_noisy_estimates(H5_100, X100, SNR, G, w_rect, w_hann);

%% ---------------------------------------------------
% 11. NORMALIZE EACH CHANNEL TO ITS OWN PEAK
% ----------------------------------------------------
%
% This is only for visual comparison of the leakage shape.
%
% No signal is artificially amplified or suppressed.
%
% Reference:
%       P_dB(n) = 10log10(|h(n)|^2 / max|h(n)|^2)

dB = @(x) 10*log10(max(abs(x).^2,eps) / ...
                    max(abs(x).^2));

%% ---------------------------------------------------
% 12. DELAY AXIS
% ----------------------------------------------------
delay_samples = 0:M-1;
delay_us = delay_samples * p.Ts * 1e6;

%% ===================================================
% 13. SUI-2 CLEAN CHANNEL
% ====================================================

figure('Name','Fig. 2 Audit - SUI-2 Clean');

plot(delay_samples, dB(h2_rect_clean), ...
    'LineWidth',1.2);
hold on;

plot(delay_samples, dB(h2_hann_clean), ...
    'LineWidth',1.2);

grid on;

xlabel('Delay index (samples)');
ylabel('Normalized CIR power (dB)');

title('SUI-2: CLEAN CHANNEL - Windowing Effect');

legend('Rectangular','Hann', ...
    'Location','best');

xlim([0 M-1]);

ylim([-80 5]);

%% ===================================================
% 14. SUI-2 NOISY CHANNEL
% ====================================================

figure('Name','Fig. 2 Audit - SUI-2 Noisy');

plot(delay_samples, dB(H2_rect_noisy), ...
    'LineWidth',1.2);
hold on;

plot(delay_samples, dB(H2_hann_noisy), ...
    'LineWidth',1.2);

grid on;

xlabel('Delay index (samples)');
ylabel('Normalized CIR power (dB)');

title('SUI-2: 20-dB SNR - Windowing Effect');

legend('Rectangular','Hann', ...
    'Location','best');

xlim([0 M-1]);

ylim([-80 5]);

%% ===================================================
% 15. SUI-5 CLEAN CHANNEL
% ====================================================

figure('Name','Fig. 2 Audit - SUI-5 Clean');

plot(delay_samples, dB(h5_rect_clean), ...
    'LineWidth',1.2);
hold on;

plot(delay_samples, dB(h5_hann_clean), ...
    'LineWidth',1.2);

% Mark the actual physical paths.
for d = delay5_samples

    xline(d,'--');

end

grid on;

xlabel('Delay index (samples)');
ylabel('Normalized CIR power (dB)');

title('SUI-5: CLEAN CHANNEL - Windowing Effect');

legend('Rectangular','Hann', ...
    'Location','best');

xlim([0 M-1]);

ylim([-80 5]);

%% ===================================================
% 16. SUI-5 NOISY CHANNEL
% ====================================================

figure('Name','Fig. 2 Audit - SUI-5 Noisy');

plot(delay_samples, dB(H5_rect_noisy), ...
    'LineWidth',1.2);
hold on;

plot(delay_samples, dB(H5_hann_noisy), ...
    'LineWidth',1.2);

for d = delay5_samples

    xline(d,'--');

end

grid on;

xlabel('Delay index (samples)');
ylabel('Normalized CIR power (dB)');

title('SUI-5: 20-dB SNR - Windowing Effect');

legend('Rectangular','Hann', ...
    'Location','best');

xlim([0 M-1]);

ylim([-80 5]);

%% ===================================================
% 17. NUMERICAL SHAPE AUDIT
% ====================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' VISUAL / SHAPE AUDIT NUMERICAL SUMMARY\n');
fprintf('====================================================\n');

% Paper uses taps 65 to 120 as the noise/leakage region.
noiseRegion = 65:120;

sui2_region = false(size(delay_samples));
sui5_region = false(size(delay_samples));

sui2_region(noiseRegion) = true;
sui5_region(noiseRegion) = true;

fprintf('\nSUI-2 physical path region: approximately 0-5 samples\n');
fprintf('SUI-2 leakage region      : samples 65-120\n');

P2_rect = abs(h2_rect_clean).^2;
P2_hann = abs(h2_hann_clean).^2;

fprintf('\nSUI-2 CLEAN:\n');

fprintf('  Rectangular leakage / peak = %.4f (%.2f dB)\n', ...
    mean(P2_rect(sui2_region))/max(P2_rect), ...
    10*log10(mean(P2_rect(sui2_region))/max(P2_rect)));

fprintf('  Hann leakage / peak        = %.4f (%.2f dB)\n', ...
    mean(P2_hann(sui2_region))/max(P2_hann), ...
    10*log10(mean(P2_hann(sui2_region))/max(P2_hann)));

fprintf('\nSUI-5 physical path region: approximately 0-40 samples\n');
fprintf('SUI-5 leakage region      : samples 65-120\n');

P5_rect = abs(h5_rect_clean).^2;
P5_hann = abs(h5_hann_clean).^2;

fprintf('\nSUI-5 CLEAN:\n');

fprintf('  Rectangular leakage / peak = %.4f (%.2f dB)\n', ...
    mean(P5_rect(sui5_region))/max(P5_rect), ...
    10*log10(mean(P5_rect(sui5_region))/max(P5_rect)));

fprintf('  Hann leakage / peak        = %.4f (%.2f dB)\n', ...
    mean(P5_hann(sui5_region))/max(P5_hann), ...
    10*log10(mean(P5_hann(sui5_region))/max(P5_hann)));

%% ---------------------------------------------------

%% ---------------------------------------------------
% 18. PATH PEAK LOCATIONS
% ----------------------------------------------------

[~,idx2_rect] = max(abs(h2_rect_clean).^2);
[~,idx2_hann] = max(abs(h2_hann_clean).^2);

[~,idx5_rect] = max(abs(h5_rect_clean).^2);
[~,idx5_hann] = max(abs(h5_hann_clean).^2);

fprintf('\nPeak locations:\n');

fprintf('  SUI-2 Rectangular = %d samples\n',idx2_rect-1);
fprintf('  SUI-2 Hann        = %d samples\n',idx2_hann-1);

fprintf('  SUI-5 Rectangular = %d samples\n',idx5_rect-1);
fprintf('  SUI-5 Hann        = %d samples\n',idx5_hann-1);

fprintf('\n');
fprintf('====================================================\n');
fprintf(' AUDIT COMPLETE\n');
fprintf('====================================================\n');

%% ===================================================
% LOCAL FUNCTIONS
% ====================================================

function H = local_cfr(h,k,N)

    h = h(:);
    k = k(:);

    l = (0:length(h)-1).';

    H = zeros(length(k),1);

    for ii = 1:length(k)

        H(ii) = sum(h .* ...
            exp(-1j*2*pi*k(ii)*l/N));

    end

end

% ----------------------------------------------------

function [H_rect,H_hann] = ...
    local_noisy_estimates(H_true,X,SNR,G,w_rect,w_hann)

    H_true = H_true(:);
    X = X(:);

    % Received pilot signal.
    Y_clean = X .* H_true;

    % Average received pilot-signal power.
    signalPower = mean(abs(Y_clean).^2);

    % AWGN variance for requested SNR.
    noisePower = signalPower / 10^(SNR/10);

    noise = sqrt(noisePower/2) .* ...
        (randn(size(Y_clean)) + ...
         1j*randn(size(Y_clean)));

    Y_noisy = Y_clean + noise;

    % LS estimate.
    H_LS = Y_noisy ./ X;

    % Insert DC.
    H101 = zeros(101,1);

    H101(1:50) = H_LS(1:50);
    H101(52:101) = H_LS(51:100);

    % Window.
    H_rect = fftshift(G * (H101 .* w_rect));
    H_hann = fftshift(G * (H101 .* w_hann));

end

% ----------------------------------------------------

function [h,delay_us,pathPower_dB] = local_sui2_channel(p)

    delay_us = p.SUI2.pathDelay_us;
    pathPower_dB = p.SUI2.pathPower_dB;

    delay_samples = delay_us * 1e-6 / p.Ts;

    % SUI-2:
    %
    % First path: Ricean, K = 3 dB
    % Remaining paths: Rayleigh

    K_dB = 3;
    K = 10^(K_dB/10);

    L = floor(max(delay_samples)) + 1;

    h = zeros(1,L);

    for ii = 1:length(delay_samples)

        amp = 10^(pathPower_dB(ii)/20);

        if ii == 1

            % Ricean coefficient.
            los = sqrt(K/(K+1));

            scatter = sqrt(1/(K+1)) * ...
                (randn + 1j*randn)/sqrt(2);

            coeff = los + scatter;

        else

            % Rayleigh coefficient.
            coeff = ...
                (randn + 1j*randn)/sqrt(2);

        end

        d = delay_samples(ii);

        % Fractional-delay path.
        %
        % sinc interpolation places the path at the
        % requested physical delay without rounding it.

        n = 0:L-1;

        h = h + amp * coeff .* ...
            sinc(n - d);

    end

end