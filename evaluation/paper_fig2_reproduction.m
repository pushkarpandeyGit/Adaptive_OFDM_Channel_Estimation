%% ========================================================================
% paper_fig2_reproduction.m
%
% FIGURE 2 REPRODUCTION - EFFECT OF WINDOWING ON CIR LENGTH ESTIMATION
% Paper: H.-R. Park, IEEE Access, 2022
%
% (a) SUI-2 model
% (b) SUI-5 model
%
% Conditions:
%   SNR = 20 dB
%   M   = 128
%   Np  = 100 active LTS tones
%
% Curves:
%   1. Rectangular window (solid blue)
%   2. Hann window (solid black)
%   3. Rectangular window (noise-free) (dash-dot blue)
%   4. Hann window (noise-free) (dash-dot black)
% ========================================================================

function paper_fig2_reproduction()

fprintf('\n============================================================\n');
fprintf(' FIGURE 2 REPRODUCTION - EFFECT OF WINDOWING\n');
fprintf('============================================================\n');

%% 1. PARAMETERS & PATHS
thisFile = mfilename('fullpath');
if isempty(thisFile), here = pwd; else, here = fileparts(thisFile); end

if exist(fullfile(here,'config','parameters.m'),'file')
    projectRoot = here;
elseif exist(fullfile(fileparts(here),'config','parameters.m'),'file')
    projectRoot = fileparts(here);
else
    projectRoot = pwd;
end
addpath(genpath(projectRoot));

p = parameters();
N      = p.N;     % 256
M      = 128;     % CIR length
SNR_dB = 20;

%% 2. LTS ACTIVE TONES & PREAMBLE
Pall = ieee80216_pall();
preamble = preamble_sequences(Pall);
PLTS_full = preamble.PLTS(:);
X100 = PLTS_full(PLTS_full ~= 0);
kLTS = [-100:2:-2 2:2:100].';

%% 3. CHANNEL GENERATION
rng(2026);
H_SUI2 = local_sui2_cfr(p, kLTS, N);
h_SUI5 = sui5_channel(p);
H_SUI5 = local_cfr_from_cir(h_SUI5, kLTS, N);

%% 4. LS ESTIMATES
rng(2027);
[HLS2_noisy, noiseVar2] = local_noisy_ls(H_SUI2, X100, SNR_dB);
[HLS5_noisy, noiseVar5] = local_noisy_ls(H_SUI5, X100, SNR_dB);

HLS2_clean = H_SUI2;
HLS5_clean = H_SUI5;

%% 5. FIGURE 2 TRANSFORMATION (Eq. 21 - 23)
% SUI-2
g2_rect_noisy = local_fig2_transform(HLS2_noisy, M, 'rectangular');
g2_hann_noisy = local_fig2_transform(HLS2_noisy, M, 'hann');
g2_rect_clean = local_fig2_transform(HLS2_clean, M, 'rectangular');
g2_hann_clean = local_fig2_transform(HLS2_clean, M, 'hann');

% SUI-5
g5_rect_noisy = local_fig2_transform(HLS5_noisy, M, 'rectangular');
g5_hann_noisy = local_fig2_transform(HLS5_noisy, M, 'hann');
g5_rect_clean = local_fig2_transform(HLS5_clean, M, 'rectangular');
g5_hann_clean = local_fig2_transform(HLS5_clean, M, 'hann');

%% 6. TAP ENERGY (dB)
epsPower = 1e-14;

E2_rect_noisy = 10*log10((abs(g2_rect_noisy).^2 + epsPower) / noiseVar2);
E2_hann_noisy = 10*log10((abs(g2_hann_noisy).^2 + epsPower) / noiseVar2);
E2_rect_clean = 10*log10((abs(g2_rect_clean).^2 + epsPower) / noiseVar2);
E2_hann_clean = 10*log10((abs(g2_hann_clean).^2 + epsPower) / noiseVar2);

E5_rect_noisy = 10*log10((abs(g5_rect_noisy).^2 + epsPower) / noiseVar5);
E5_hann_noisy = 10*log10((abs(g5_hann_noisy).^2 + epsPower) / noiseVar5);
E5_rect_clean = 10*log10((abs(g5_rect_clean).^2 + epsPower) / noiseVar5);
E5_hann_clean = 10*log10((abs(g5_hann_clean).^2 + epsPower) / noiseVar5);

tap = 0:M-1;

%% 7. PLOT FIGURE 2 (MATCHING IEEE ACCESS 2022 1:1)
fig2 = figure('Name','Figure 2 - Effect of Windowing','Color','w', ...
              'Units','pixels','Position',[100 50 650 800]);

% (a) SUI-2 model
ax1 = subplot(2,1,1);
plot(tap, E2_rect_noisy, '-b', 'LineWidth', 1.0); hold on;
plot(tap, E2_hann_noisy, '-k', 'LineWidth', 1.0);
plot(tap, E2_rect_clean, '-.b', 'LineWidth', 0.9);
plot(tap, E2_hann_clean, '-.k', 'LineWidth', 0.9);
grid on; box on;
xlim([0 128]); ylim([-40 50]);
xlabel('Time Index (T_s)', 'Color', 'k', 'FontSize', 10);
ylabel('Tap Energy (dB)', 'Color', 'k', 'FontSize', 10);
lgd1 = legend('Rectangular window', 'Hann window', ...
       'Rectangular window (noise-free)', 'Hann window (noise-free)', ...
       'Location', 'northeast', 'FontSize', 8, 'TextColor', 'k', ...
       'Color', 'w', 'EdgeColor', [0.7 0.7 0.7]);
text(0.5, -0.22, '(a)', 'Units', 'normalized', ...
     'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'normal', 'Color', 'k');
set(ax1, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
         'GridColor', [0.8 0.8 0.8], 'GridAlpha', 0.7, ...
         'FontName', 'Times New Roman', 'FontSize', 10);

% (b) SUI-5 model
ax2 = subplot(2,1,2);
plot(tap, E5_rect_noisy, '-b', 'LineWidth', 1.0); hold on;
plot(tap, E5_hann_noisy, '-k', 'LineWidth', 1.0);
plot(tap, E5_rect_clean, '-.b', 'LineWidth', 0.9);
plot(tap, E5_hann_clean, '-.k', 'LineWidth', 0.9);
grid on; box on;
xlim([0 128]); ylim([-40 50]);
xlabel('Time Index (T_s)', 'Color', 'k', 'FontSize', 10);
ylabel('Tap Energy (dB)', 'Color', 'k', 'FontSize', 10);
lgd2 = legend('Rectangular window', 'Hann window', ...
       'Rectangular window (noise-free)', 'Hann window (noise-free)', ...
       'Location', 'northeast', 'FontSize', 8, 'TextColor', 'k', ...
       'Color', 'w', 'EdgeColor', [0.7 0.7 0.7]);
text(0.5, -0.22, '(b)', 'Units', 'normalized', ...
     'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'normal', 'Color', 'k');
set(ax2, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
         'GridColor', [0.8 0.8 0.8], 'GridAlpha', 0.7, ...
         'FontName', 'Times New Roman', 'FontSize', 10);

%% 8. EXPORT HIGH-RES GRAPHICS
evalDir = fullfile(projectRoot, 'evaluation');
if ~exist(evalDir, 'dir'), mkdir(evalDir); end

exportgraphics(fig2, fullfile(evalDir, 'paper_fig2_reproduction.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(here, 'paper_fig2_reproduction.png'), 'Resolution', 300);

save(fullfile(evalDir, 'paper_fig2_results.mat'), 'tap', ...
     'E2_rect_noisy', 'E2_hann_noisy', 'E2_rect_clean', 'E2_hann_clean', ...
     'E5_rect_noisy', 'E5_hann_noisy', 'E5_rect_clean', 'E5_hann_clean');
save(fullfile(here, 'paper_fig2_results.mat'), 'tap', ...
     'E2_rect_noisy', 'E2_hann_noisy', 'E2_rect_clean', 'E2_hann_clean', ...
     'E5_rect_noisy', 'E5_hann_noisy', 'E5_rect_clean', 'E5_hann_clean');

fprintf('Figure 2 saved successfully to:\n');
fprintf('  %s\n', fullfile(evalDir, 'paper_fig2_reproduction.png'));

end


%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function H = local_sui2_cfr(p, k, N)
    delay_us = [0 0.4 1.1];
    power_dB = [0 -12 -15];
    delay_samples = delay_us * 1e-6 / p.Ts;
    K_dB = 3;
    if isfield(p.SUI2,'RiceanK_dB')
        K_dB = p.SUI2.RiceanK_dB;
    elseif isfield(p.SUI2,'K_dB')
        K_dB = p.SUI2.K_dB;
    end
    K = 10^(K_dB/10);
    alpha = zeros(1,3);
    los = sqrt(K/(K+1));
    scattered = sqrt(1/(K+1)) * (randn + 1j*randn)/sqrt(2);
    alpha(1) = 10^(power_dB(1)/20) * (los + scattered);
    for ii = 2:3
        alpha(ii) = 10^(power_dB(ii)/20) * (randn + 1j*randn)/sqrt(2);
    end
    H = zeros(size(k));
    for ii = 1:3
        H = H + alpha(ii) .* exp(-1j*2*pi*k*delay_samples(ii)/N);
    end
end

function H = local_cfr_from_cir(h, k, N)
    h = h(:).';
    n = 0:length(h)-1;
    H = zeros(size(k));
    for ii = 1:length(k)
        H(ii) = sum(h .* exp(-1j*2*pi*k(ii)*n/N));
    end
end

function [H_LS, noiseVarianceLS] = local_noisy_ls(H, X, SNR_dB)
    H = H(:);
    X = X(:);
    Y_clean = H .* X;
    signalPower = mean(abs(Y_clean).^2);
    noisePower = signalPower / 10^(SNR_dB/10);
    noise = sqrt(noisePower/2) .* (randn(size(Y_clean)) + 1j*randn(size(Y_clean)));
    Y = Y_clean + noise;
    H_LS = Y ./ X;
    noiseVarianceLS = mean(abs(noise ./ X).^2);
end

function g = local_fig2_transform(H_LS, M, windowType)
    H_LS = H_LS(:);
    % 101-point centered vector
    p101 = zeros(101,1);
    p101(1:50) = H_LS(1:50);
    p101(51) = 0.5 * (H_LS(50) + H_LS(51)); % DC interpolation
    p101(52:101) = H_LS(51:100);

    n = (-50:50).';
    switch lower(windowType)
        case 'rectangular'
            w = ones(101,1);
        case 'hann'
            w = 0.5 + 0.5*cos(2*pi*n/M);
        otherwise
            error('Unknown window: %s', windowType);
    end
    p_windowed = p101 .* w;

    m = (0:M-1).';
    q = (-50:50);
    G = (1/sqrt(M)) .* exp(1j*2*pi .* ((m - M/2) * q) / M);
    g = G * p_windowed;
    g = fftshift(g);
    g = g(:).';
end
