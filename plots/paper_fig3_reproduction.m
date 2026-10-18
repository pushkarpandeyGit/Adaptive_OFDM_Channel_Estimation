function paper_fig3_reproduction()
% PAPER_FIG3_REPRODUCTION
%
% Reproduction of Fig. 3 from:
% H.-R. Park et al., "A Low-Complexity Channel Estimation for OFDM Systems
% Based on CIR Length Adaptation," IEEE Access, 2022.
%
% Fig. 3: Distributions of the CIR length estimates obtained by
% SHT, SOT, and GAIC for SNR = 20 dB.
% (a) SUI-2 model
% (b) SUI-5 model

fprintf('\n============================================\n');
fprintf(' FIGURE 3 REPRODUCTION - SHT, SOT, GAIC\n');
fprintf('============================================\n');

%% 1. PROJECT SETUP
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

%% 2. PARAMETERS
SNRdB = 20;
numTrials = 1000;

N = 256;              % OFDM FFT size
M = 128;              % CIR transform size
G = 64;               % Maximum CIR length tested
K = 56;               % Noise-estimation taps (65:120)

T_SHT = 15;           % SHT threshold multiplier
T_SOT = 13.5;         % SOT theoretical multiplier at SNR = 20 dB
gamma_GAIC = 0.022;   % GAIC penalty parameter

kLTS = [-100:2:-2 2:2:100].';
Np = length(kLTS);    % 100

%% 3. KNOWN LTS
Pall = ieee80216_pall();
preamble = preamble_sequences(Pall);
[~, ltsTime] = preamble_to_ifft(preamble,p);
lts128 = get_known_lts(ltsTime);
X128 = fft(lts128,M);
ltsBins = mod(kLTS/2,M) + 1;
X = X128(ltsBins);
X = X(:);

%% 4. HANN WINDOW & CIR TRANSFORM MATRICES
q = (-50:50).';
w = 0.5 + 0.5*cos(2*pi*q/M);
w = w(:);

m = (0:M-1).';
Gmat = (1/sqrt(M)) * exp(1j*2*pi*((m-M/2)*q.')/M);

%% 5. PRECOMPUTE GAIC BASIS PROJECTIONS (High Speed)
A_full = exp(-1j*2*pi*(kLTS*(0:G-1))/N);
proj_matrices = cell(G,1);
for L = 1:G
    AL = A_full(:,1:L);
    proj_matrices{L} = AL * pinv(AL);
end

%% 6. PHYSICAL MULTIPATH PARAMETERS (CONTINUOUS DELAYS)
% SUI-2: [0, 0.4, 1.1] us -> [0, 1.6, 4.4] Ts
delays2 = [0, 0.4, 1.1] * 1e-6 / p.Ts;
powers2_lin = 10.^([0, -12, -15]/20);
K_lin = 10^(3/10); % 3 dB Ricean factor

% SUI-5: [0, 4.0, 10.0] us -> [0, 16.0, 40.0] Ts
delays5 = [0, 4.0, 10.0] * 1e-6 / p.Ts;
powers5_lin = 10.^([0, -5, -10]/20);

%% 7. STORAGE
L_SHT_SUI2  = zeros(numTrials,1);
L_SOT_SUI2  = zeros(numTrials,1);
L_GAIC_SUI2 = zeros(numTrials,1);

L_SHT_SUI5  = zeros(numTrials,1);
L_SOT_SUI5  = zeros(numTrials,1);
L_GAIC_SUI5 = zeros(numTrials,1);

rng(2026);

fprintf('SNR    = %d dB\n', SNRdB);
fprintf('Trials = %d\n\n', numTrials);

%% 8. MONTE-CARLO SIMULATION
for trial = 1:numTrials
    if mod(trial,250) == 0
        fprintf('  Trial %d / %d\n', trial, numTrials);
    end

    %% --- Continuous-Delay SUI-2 Channel ---
    los = sqrt(K_lin/(K_lin+1));
    nlos = sqrt(1/(K_lin+1)) * (randn + 1j*randn)/sqrt(2);
    a1 = powers2_lin(1) * (los + nlos);
    a2 = powers2_lin(2) * (randn + 1j*randn)/sqrt(2);
    a3 = powers2_lin(3) * (randn + 1j*randn)/sqrt(2);
    H2 = a1 * exp(-1j*2*pi*kLTS*delays2(1)/N) + ...
         a2 * exp(-1j*2*pi*kLTS*delays2(2)/N) + ...
         a3 * exp(-1j*2*pi*kLTS*delays2(3)/N);

    %% --- Continuous-Delay SUI-5 Channel ---
    b1 = powers5_lin(1) * (randn + 1j*randn)/sqrt(2);
    b2 = powers5_lin(2) * (randn + 1j*randn)/sqrt(2);
    b3 = powers5_lin(3) * (randn + 1j*randn)/sqrt(2);
    H5 = b1 * exp(-1j*2*pi*kLTS*delays5(1)/N) + ...
         b2 * exp(-1j*2*pi*kLTS*delays5(2)/N) + ...
         b3 * exp(-1j*2*pi*kLTS*delays5(3)/N);

    %% Received LTS & AWGN
    Y2_clean = X .* H2;
    Y5_clean = X .* H5;

    noisePower2 = mean(abs(Y2_clean).^2) / 10^(SNRdB/10);
    noisePower5 = mean(abs(Y5_clean).^2) / 10^(SNRdB/10);

    n2 = sqrt(noisePower2/2) .* (randn(size(Y2_clean)) + 1j*randn(size(Y2_clean)));
    n5 = sqrt(noisePower5/2) .* (randn(size(Y5_clean)) + 1j*randn(size(Y5_clean)));

    HLS2 = (Y2_clean + n2) ./ X;
    HLS5 = (Y5_clean + n5) ./ X;

    %% --- CIR Transforms ---
    % DC Interpolation
    H101_2 = [HLS2(1:50); 0.5*(HLS2(50)+HLS2(51)); HLS2(51:100)];
    H101_5 = [HLS5(1:50); 0.5*(HLS5(50)+HLS5(51)); HLS5(51:100)];

    % Windowed CIR for SHT
    cirSHT2 = fftshift(Gmat * (H101_2 .* w));
    cirSHT5 = fftshift(Gmat * (H101_5 .* w));

    % Unwindowed CIR for SOT
    cirSOT2 = fftshift(Gmat * H101_2);
    cirSOT5 = fftshift(Gmat * H101_5);

    %% --- SHT Evaluation ---
    sigma2_SHT_2 = mean(abs(cirSHT2(G+1:G+K)).^2);
    sigma2_SHT_5 = mean(abs(cirSHT5(G+1:G+K)).^2);

    L_SHT_SUI2(trial) = sht_find_length(cirSHT2, G, T_SHT * sigma2_SHT_2);
    L_SHT_SUI5(trial) = sht_find_length(cirSHT5, G, T_SHT * sigma2_SHT_5);

    %% --- SOT Evaluation ---
    sigma2_SOT_2 = mean(abs(cirSOT2(G+1:G+K)).^2);
    sigma2_SOT_5 = mean(abs(cirSOT5(G+1:G+K)).^2);

    L_SOT_SUI2(trial) = sot_find_length(cirSOT2, G, T_SOT * sigma2_SOT_2);
    L_SOT_SUI5(trial) = sot_find_length(cirSOT5, G, T_SOT * sigma2_SOT_5);

    %% --- GAIC Evaluation ---
    L_GAIC_SUI2(trial) = fast_gaic(HLS2, proj_matrices, G, gamma_GAIC);
    L_GAIC_SUI5(trial) = fast_gaic(HLS5, proj_matrices, G, gamma_GAIC);
end

%% 9. STATISTICS
fprintf('\n============================================\n');
fprintf('FIGURE 3 RESULTS\n');
fprintf('============================================\n');
print_stats('SUI-2 / SHT',  L_SHT_SUI2);
print_stats('SUI-2 / SOT',  L_SOT_SUI2);
print_stats('SUI-2 / GAIC', L_GAIC_SUI2);
fprintf('\n');
print_stats('SUI-5 / SHT',  L_SHT_SUI5);
print_stats('SUI-5 / SOT',  L_SOT_SUI5);
print_stats('SUI-5 / GAIC', L_GAIC_SUI5);

%% 10. PLOT FIGURE 3 (MATCHING IEEE ACCESS 2022 1:1)
fig3 = plot_fig3_master( ...
    L_SHT_SUI2, L_SOT_SUI2, L_GAIC_SUI2, ...
    L_SHT_SUI5, L_SOT_SUI5, L_GAIC_SUI5, ...
    projectRoot);

end

%% =================================================================
% LOCAL HELPER FUNCTIONS
%% =================================================================

function Lhat = sht_find_length(cir, G, threshold)
    idx = find(abs(cir(1:G)).^2 > threshold, 1, 'last');
    if isempty(idx)
        Lhat = 1;
    else
        Lhat = idx;
    end
end

function Lhat = sot_find_length(cir, G, threshold)
    idx = find(abs(cir(1:G)).^2 > threshold, 1, 'last');
    if isempty(idx)
        Lhat = 1;
    else
        Lhat = idx;
    end
end

function Lhat = fast_gaic(HLS, proj_matrices, G, gamma)
    HLS = HLS(:);
    crit = zeros(G,1);
    for L = 1:G
        res = HLS - proj_matrices{L} * HLS;
        mse = max(mean(abs(res).^2), 1e-12);
        crit(L) = log(mse) + gamma * L;
    end
    [~, Lhat] = min(crit);
end

function print_stats(name,x)
    fprintf('%-18s : mean = %7.3f   std = %7.3f\n', name, mean(x), std(x));
end

function fig = plot_fig3_master(L_SHT_SUI2, L_SOT_SUI2, L_GAIC_SUI2, ...
                                L_SHT_SUI5, L_SOT_SUI5, L_GAIC_SUI5, rootDir)

    maxL = 65;
    bins = 0:maxL;

    % Calculate PMFs
    p1_2 = histcounts(L_SHT_SUI2, [-0.5, bins+0.5], 'Normalization', 'probability');
    p2_2 = histcounts(L_SOT_SUI2, [-0.5, bins+0.5], 'Normalization', 'probability');
    p3_2 = histcounts(L_GAIC_SUI2, [-0.5, bins+0.5], 'Normalization', 'probability');

    p1_5 = histcounts(L_SHT_SUI5, [-0.5, bins+0.5], 'Normalization', 'probability');
    p2_5 = histcounts(L_SOT_SUI5, [-0.5, bins+0.5], 'Normalization', 'probability');
    p3_5 = histcounts(L_GAIC_SUI5, [-0.5, bins+0.5], 'Normalization', 'probability');

    fig = figure('Name', 'Figure 3 - CIR Length Distribution', 'Color', 'w', ...
                 'Units', 'pixels', 'Position', [100 50 680 840]);

    %% (a) SUI-2 model
    ax1 = subplot(2,1,1);
    plot(bins, p1_2, '-*k', 'LineWidth', 1.0, 'MarkerSize', 5); hold on;
    plot(bins, p2_2, '-sb', 'LineWidth', 1.0, 'MarkerSize', 5);
    plot(bins, p3_2, '-or', 'LineWidth', 1.0, 'MarkerSize', 5);
    grid on; box on;
    xlim([0 20]); ylim([0 1]);
    xticks(0:2:20); yticks(0:0.1:1);
    xlabel('CIR Length Estimates (T_s)', 'Color', 'k', 'FontSize', 10);
    ylabel('Probability', 'Color', 'k', 'FontSize', 10);
    legend('SHT with T=15', 'SOT', 'GAIC', 'Location', 'northeast', 'FontSize', 8, ...
           'TextColor', 'k', 'Color', 'w', 'EdgeColor', [0.7 0.7 0.7]);
    text(0.5, -0.20, '(a)', 'Units', 'normalized', ...
         'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'normal', 'Color', 'k');
    set(ax1, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
             'GridColor', [0.8 0.8 0.8], 'GridAlpha', 0.7, ...
             'FontName', 'Times New Roman', 'FontSize', 10);

    %% (b) SUI-5 model
    axMain = subplot(2,1,2);
    plot(bins, p1_5, '-*k', 'LineWidth', 1.0, 'MarkerSize', 5); hold on;
    plot(bins, p2_5, '-sb', 'LineWidth', 1.0, 'MarkerSize', 5);
    plot(bins, p3_5, '-or', 'LineWidth', 1.0, 'MarkerSize', 5);
    grid on; box on;
    xlim([0 65]); ylim([0 1]);
    xticks(0:10:60); yticks(0:0.1:1);
    xlabel('CIR Length Estimates (T_s)', 'Color', 'k', 'FontSize', 10);
    ylabel('Probability', 'Color', 'k', 'FontSize', 10);
    legend('SHT with T=15', 'SOT', 'GAIC', 'Location', 'northeast', 'FontSize', 8, ...
           'TextColor', 'k', 'Color', 'w', 'EdgeColor', [0.7 0.7 0.7]);
    text(0.5, -0.20, '(b)', 'Units', 'normalized', ...
         'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'normal', 'Color', 'k');
    set(axMain, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
                'GridColor', [0.8 0.8 0.8], 'GridAlpha', 0.7, ...
                'FontName', 'Times New Roman', 'FontSize', 10);

    %% Inset inside (b) zooming in on low probabilities
    pos = get(axMain, 'Position');
    axInset = axes('Position', [pos(1) + 0.08*pos(3), pos(2) + 0.38*pos(4), 0.36*pos(3), 0.46*pos(4)]);
    plot(bins, p1_5, '-*k', 'LineWidth', 0.8, 'MarkerSize', 3.5); hold on;
    plot(bins, p2_5, '-sb', 'LineWidth', 0.8, 'MarkerSize', 3.5);
    plot(bins, p3_5, '-or', 'LineWidth', 0.8, 'MarkerSize', 3.5);
    grid on; box on;
    xlim([0 65]); ylim([0 0.02]);
    xticks(0:20:60); yticks(0:0.005:0.02);
    set(axInset, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
                 'GridColor', [0.8 0.8 0.8], 'GridAlpha', 0.7, ...
                 'FontName', 'Times New Roman', 'FontSize', 8);

    %% Save output figures & results
    evalDir = fullfile(rootDir, 'evaluation');
    if ~exist(evalDir, 'dir'), mkdir(evalDir); end
    plotsDir = fullfile(rootDir, 'plots');
    if ~exist(plotsDir, 'dir'), mkdir(plotsDir); end

    exportgraphics(fig, fullfile(evalDir, 'paper_fig3_reproduction.png'), 'Resolution', 300);
    exportgraphics(fig, fullfile(plotsDir, 'paper_fig3_reproduction.png'), 'Resolution', 300);

    res_fig3.L_SHT_SUI2  = L_SHT_SUI2;
    res_fig3.L_SOT_SUI2  = L_SOT_SUI2;
    res_fig3.L_GAIC_SUI2 = L_GAIC_SUI2;
    res_fig3.L_SHT_SUI5  = L_SHT_SUI5;
    res_fig3.L_SOT_SUI5  = L_SOT_SUI5;
    res_fig3.L_GAIC_SUI5 = L_GAIC_SUI5;

    save(fullfile(evalDir, 'paper_fig3_results.mat'), '-struct', 'res_fig3');
    save(fullfile(plotsDir, 'paper_fig3_results.mat'), '-struct', 'res_fig3');

    fprintf('\nFigure 3 saved successfully to:\n');
    fprintf('  %s\n', fullfile(evalDir, 'paper_fig3_reproduction.png'));
end