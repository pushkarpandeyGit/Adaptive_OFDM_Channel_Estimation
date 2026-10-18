function step8_plot_adaptive_window_effect()
% =========================================================================
% STEP 8 (PLOT 6): EFFECT OF ADAPTIVE WINDOWING ON CIR ENERGY COMPACTION
%
% Replicates and significantly upgrades Paper Figure 2 (H.-R. Park 2022)
% by evaluating:
%   1. Rectangular Window (Raw frequency truncation)
%   2. Fixed Hann Window (Park 2022 Baseline)
%   3. PROPOSED ADAPTIVE WINDOW (Hamming for SUI-2, Blackman for SUI-5)
%
% Evaluates across both noise-free channels and SNR = 20 dB noisy channels
% in SUI-2 (short dispersion) and SUI-5 (long dispersion).
%
% Outputs:
%   novelty/Novelty_Day3_Plot6_Adaptive_Windowing_Effect.png
%   novelty/Novelty_Day3_Plot6_Adaptive_Windowing_results.mat
% =========================================================================

fprintf('\n============================================================\n');
fprintf(' NOVELTY PLOT 6: ADAPTIVE WINDOWING ON CIR ENERGY COMPACTION\n');
fprintf('============================================================\n');

%% 1. SYSTEM SETUP
thisFile = mfilename('fullpath');
if isempty(thisFile), here = pwd; else, here = fileparts(thisFile); end
projectRoot = fileparts(here);
addpath(genpath(projectRoot));

p = parameters();
rng(2026);

N      = p.N;   % 256
M      = 128;   % CIR size
SNR_dB = 20;    % Paper condition (Fig. 2)
gammaP = 10^((SNR_dB + 6)/10);

kLTS = [-100:2:-2 2:2:100].';
Np = length(kLTS); % 100 tones

% LTS preamble reference
Pall = ieee80216_pall();
preamble = preamble_sequences(Pall);
PLTS_full = preamble.PLTS(:);
X100 = PLTS_full(PLTS_full ~= 0);

%% 2. CHANNELS
% SUI-2: delays = [0, 0.4, 1.1] us, powers = [0, -12, -15] dB
tau2 = [0, 0.4, 1.1] * 1e-6 / p.Ts;
p2_lin = 10.^([0, -12, -15]/20);
K2_lin = 10^(3/10);
los2 = sqrt(K2_lin/(K2_lin+1));
nlos2 = sqrt(1/(K2_lin+1)) * (randn + 1j*randn)/sqrt(2);
a1 = p2_lin(1) * (los2 + nlos2);
a2 = p2_lin(2) * (randn + 1j*randn)/sqrt(2);
a3 = p2_lin(3) * (randn + 1j*randn)/sqrt(2);
H2_clean = a1 * exp(-1j*2*pi*kLTS*tau2(1)/N) + ...
            a2 * exp(-1j*2*pi*kLTS*tau2(2)/N) + ...
            a3 * exp(-1j*2*pi*kLTS*tau2(3)/N);

% SUI-5: delays = [0, 4.0, 10.0] us, powers = [0, -5, -10] dB
tau5 = [0, 4.0, 10.0] * 1e-6 / p.Ts;
p5_lin = 10.^([0, -5, -10]/20);
b1 = p5_lin(1) * (randn + 1j*randn)/sqrt(2);
b2 = p5_lin(2) * (randn + 1j*randn)/sqrt(2);
b3 = p5_lin(3) * (randn + 1j*randn)/sqrt(2);
H5_clean = b1 * exp(-1j*2*pi*kLTS*tau5(1)/N) + ...
            b2 * exp(-1j*2*pi*kLTS*tau5(2)/N) + ...
            b3 * exp(-1j*2*pi*kLTS*tau5(3)/N);

% Noisy observations at 20 dB
noise2 = sqrt(1/(2*gammaP)) * (randn(Np,1) + 1j*randn(Np,1));
noise5 = sqrt(1/(2*gammaP)) * (randn(Np,1) + 1j*randn(Np,1));

H2_noisy = H2_clean + noise2;
H5_noisy = H5_clean + noise5;

%% 3. WINDOW TRANSFORMATIONS
% Transformation operator matrix (101 tones -> 128 CIR)
n = -50:50;
m_idx = (0:M-1).';
Gmat = (1/sqrt(M)) * exp(1j*2*pi*((m_idx - M/2)*n)/M);

% Windows
w_rect = ones(1, 101);
w_hann = 0.50 + 0.50 * cos(2*pi*n/M);

% Adaptive Window Selector
[winName2, w_adapt2, met2] = adaptive_window_selector(H2_noisy, M, SNR_dB);
[winName5, w_adapt5, met5] = adaptive_window_selector(H5_noisy, M, SNR_dB);

fprintf('SUI-2 Adaptive Selection: %s (P90 = %d, TailRatio = %.3f)\n', ...
    winName2, met2.p90, met2.tailEnergyRatio);
fprintf('SUI-5 Adaptive Selection: %s (P90 = %d, TailRatio = %.3f)\n', ...
    winName5, met5.p90, met5.tailEnergyRatio);

% Construct DC-interpolated 101 grids
H2_101_clean = [H2_clean(1:50); 0.5*(H2_clean(50)+H2_clean(51)); H2_clean(51:100)];
H2_101_noisy = [H2_noisy(1:50); 0.5*(H2_noisy(50)+H2_noisy(51)); H2_noisy(51:100)];

H5_101_clean = [H5_clean(1:50); 0.5*(H5_clean(50)+H5_clean(51)); H5_clean(51:100)];
H5_101_noisy = [H5_noisy(1:50); 0.5*(H5_noisy(50)+H5_noisy(51)); H5_noisy(51:100)];

% CIRs for SUI-2
cir2_rect_clean = fftshift(Gmat * (H2_101_clean .* w_rect.'));
cir2_hann_clean = fftshift(Gmat * (H2_101_clean .* w_hann.'));
cir2_adpt_clean = fftshift(Gmat * (H2_101_clean .* w_adapt2.'));

cir2_rect_noisy = fftshift(Gmat * (H2_101_noisy .* w_rect.'));
cir2_hann_noisy = fftshift(Gmat * (H2_101_noisy .* w_hann.'));
cir2_adpt_noisy = fftshift(Gmat * (H2_101_noisy .* w_adapt2.'));

% CIRs for SUI-5
cir5_rect_clean = fftshift(Gmat * (H5_101_clean .* w_rect.'));
cir5_hann_clean = fftshift(Gmat * (H5_101_clean .* w_hann.'));
cir5_adpt_clean = fftshift(Gmat * (H5_101_clean .* w_adapt5.'));

cir5_rect_noisy = fftshift(Gmat * (H5_101_noisy .* w_rect.'));
cir5_hann_noisy = fftshift(Gmat * (H5_101_noisy .* w_hann.'));
cir5_adpt_noisy = fftshift(Gmat * (H5_101_noisy .* w_adapt5.'));

% Normalization by noise variance as in Paper Fig. 2
noiseVar2 = mean(abs(cir2_hann_noisy(65:120)).^2);
noiseVar5 = mean(abs(cir5_hann_noisy(65:120)).^2);

% Log Tap Energy (dB)
E2_rect_clean = 10*log10(abs(cir2_rect_clean).^2 / noiseVar2);
E2_hann_clean = 10*log10(abs(cir2_hann_clean).^2 / noiseVar2);
E2_adpt_clean = 10*log10(abs(cir2_adpt_clean).^2 / noiseVar2);

E2_rect_noisy = 10*log10(abs(cir2_rect_noisy).^2 / noiseVar2);
E2_hann_noisy = 10*log10(abs(cir2_hann_noisy).^2 / noiseVar2);
E2_adpt_noisy = 10*log10(abs(cir2_adpt_noisy).^2 / noiseVar2);

E5_rect_clean = 10*log10(abs(cir5_rect_clean).^2 / noiseVar5);
E5_hann_clean = 10*log10(abs(cir5_hann_clean).^2 / noiseVar5);
E5_adpt_clean = 10*log10(abs(cir5_adpt_clean).^2 / noiseVar5);

E5_rect_noisy = 10*log10(abs(cir5_rect_noisy).^2 / noiseVar5);
E5_hann_noisy = 10*log10(abs(cir5_hann_noisy).^2 / noiseVar5);
E5_adpt_noisy = 10*log10(abs(cir5_adpt_noisy).^2 / noiseVar5);

taps = 0:M-1;

%% 4. LEAKAGE SUPPRESSION ANALYSIS TABLE
leakage_noise_region_rect_2 = mean(E2_rect_clean(65:120));
leakage_noise_region_hann_2 = mean(E2_hann_clean(65:120));
leakage_noise_region_adpt_2 = mean(E2_adpt_clean(65:120));

leakage_noise_region_rect_5 = mean(E5_rect_clean(65:120));
leakage_noise_region_hann_5 = mean(E5_hann_clean(65:120));
leakage_noise_region_adpt_5 = mean(E5_adpt_clean(65:120));

fprintf('\n========================================================================================\n');
fprintf(' SPECTRAL LEAKAGE IN NOISE OBSERVATION REGION (Taps 65 to 120)\n');
fprintf('========================================================================================\n');
fprintf(' Channel | Rectangular Window | Park 2022 Fixed Hann | Proposed Adaptive Window | Gain (dB)\n');
fprintf('----------------------------------------------------------------------------------------\n');
fprintf('  SUI-2  |     %7.2f dB       |       %7.2f dB      |   %7.2f dB (%s)   |  %5.1f dB\n', ...
    leakage_noise_region_rect_2, leakage_noise_region_hann_2, leakage_noise_region_adpt_2, ...
    winName2, leakage_noise_region_hann_2 - leakage_noise_region_adpt_2);
fprintf('  SUI-5  |     %7.2f dB       |       %7.2f dB      |   %7.2f dB (%s) |  %5.1f dB\n', ...
    leakage_noise_region_rect_5, leakage_noise_region_hann_5, leakage_noise_region_adpt_5, ...
    winName5, leakage_noise_region_hann_5 - leakage_noise_region_adpt_5);
fprintf('========================================================================================\n');

%% 5. PUBLICATION FIGURE
fig = figure('Name', 'Novelty Plot 6 - Adaptive Windowing Effect', ...
    'Color', 'w', 'Position', [80 60 1200 540], 'Visible', 'off');

% (a) SUI-2 Channel
subplot(1, 2, 1);
plot(taps, E2_rect_clean, ':k', 'LineWidth', 1.4); hold on;
plot(taps, E2_hann_clean, '--b', 'LineWidth', 1.6);
plot(taps, E2_adpt_clean, '-r', 'LineWidth', 2.0);
plot(taps, E2_adpt_noisy, '-m', 'LineWidth', 1.0);
grid on; box on;
xlim([0 127]); ylim([-35 40]);
xlabel('Time Index (samples)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Tap Energy (dB)', 'FontSize', 11, 'FontWeight', 'bold');
title('(a) SUI-2 Channel (Short Delay Spread, P_{95} \leq 8)', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Rectangular (clean)', 'Fixed Hann (clean) [Park 2022]', ...
        ['Proposed ' upper(winName2) ' (clean)'], ['Proposed ' upper(winName2) ' (noisy 20dB)']}, ...
        'Location', 'northeast', 'FontSize', 9);
% Annotate Noise Region
patch([64 120 120 64], [-35 -35 40 40], [0.9 0.9 0.9], 'FaceAlpha', 0.2, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
text(92, -28, 'Noise Region', 'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', [0.3 0.3 0.3]);

% (b) SUI-5 Channel
subplot(1, 2, 2);
plot(taps, E5_rect_clean, ':k', 'LineWidth', 1.4); hold on;
plot(taps, E5_hann_clean, '--b', 'LineWidth', 1.6);
plot(taps, E5_adpt_clean, '-r', 'LineWidth', 2.0);
plot(taps, E5_adpt_noisy, '-m', 'LineWidth', 1.0);
grid on; box on;
xlim([0 127]); ylim([-40 40]);
xlabel('Time Index (samples)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Tap Energy (dB)', 'FontSize', 11, 'FontWeight', 'bold');
title('(b) SUI-5 Channel (High Delay Spread, P_{95} > 20)', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Rectangular (clean)', 'Fixed Hann (clean) [Park 2022]', ...
        ['Proposed ' upper(winName5) ' (clean)'], ['Proposed ' upper(winName5) ' (noisy 20dB)']}, ...
        'Location', 'northeast', 'FontSize', 9);
% Annotate Noise Region
patch([64 120 120 64], [-40 -40 40 40], [0.9 0.9 0.9], 'FaceAlpha', 0.2, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
text(92, -32, 'Noise Region', 'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', [0.3 0.3 0.3]);

sgtitle('Effect of Proposed Adaptive Windowing on CIR Energy Compaction vs Park 2022 Baseline', ...
    'FontSize', 13, 'FontWeight', 'bold');

% Save Artifacts
outPng = fullfile(projectRoot, 'novelty', 'Novelty_Day3_Plot6_Adaptive_Windowing_Effect.png');
outMat = fullfile(projectRoot, 'novelty', 'Novelty_Day3_Plot6_Adaptive_Windowing_results.mat');

saveas(fig, outPng);
save(outMat, 'taps', 'E2_rect_clean', 'E2_hann_clean', 'E2_adpt_clean', ...
     'E5_rect_clean', 'E5_hann_clean', 'E5_adpt_clean', ...
     'leakage_noise_region_hann_5', 'leakage_noise_region_adpt_5', '-v7');

fprintf('\nSuccessfully saved Plot 6 results to:\n  %s\n  %s\n', outPng, outMat);
fprintf('============================================================\n');

end
