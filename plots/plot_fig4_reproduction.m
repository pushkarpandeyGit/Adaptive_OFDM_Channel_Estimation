%% ============================================================
% plot_fig4_reproduction.m
%
% FIGURE 4 REPRODUCTION - SHT THRESHOLD MULTIPLIER COMPARISON
% Paper: H.-R. Park et al., "A Low-Complexity Channel Estimation for 
%        OFDM Systems Based on CIR Length Adaptation," IEEE Access, 2022.
%
% (a) SUI-2 model
% (b) SUI-5 model
%
% Distributions of CIR-length estimates obtained by SHT for:
%   T = 5, 15, 25 at SNR = 20 dB.
% =============================================================

function plot_fig4_reproduction()

fprintf('\n====================================================\n');
fprintf('       FIGURE 4 REPRODUCTION - SHT THRESHOLDS\n');
fprintf('====================================================\n');

%% 1. PROJECT PATH
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

%% 2. PARAMETERS
p = parameters();
rng(p.randomSeed);

SNR = 20;
T_values = [5 15 25];
numTrials = 2000;

G = p.G;              % 64
K = p.K;              % 56
noiseStart = G + 1;   % 65
noiseEnd   = G + K;   % 120

fprintf('SNR        = %d dB\n', SNR);
fprintf('Trials     = %d\n', numTrials);
fprintf('G          = %d\n', G);
fprintf('K          = %d (noise taps %d:%d)\n', K, noiseStart, noiseEnd);
fprintf('T values   = [%d %d %d]\n\n', T_values);

%% 3. BUILD IEEE 802.16 PREAMBLE
Pall = ieee80216_pall();
preamble = preamble_sequences(Pall);
[stsTime, ltsTime] = preamble_to_ifft(preamble, p);
txStruct = build_80216_preamble(stsTime, ltsTime, p);
tx = txStruct.complete;

knownLTS = get_known_lts(ltsTime);
knownLTS_FFT = known_lts_fft(knownLTS);

kLTS = [-100:2:-2 2:2:100];
ltsBins = mod(kLTS/2,128) + 1;

txPower = mean(abs(tx).^2);
noisePower = txPower / 10^(SNR/10);

%% 4. MONTE-CARLO SIMULATION
Lhat_SUI2 = zeros(numTrials, length(T_values));
Lhat_SUI5 = zeros(numTrials, length(T_values));

fprintf('Running Monte-Carlo trials...\n');
for trial = 1:numTrials
    %% SUI-2
    h2 = sui2_channel(p);
    rxClean2 = apply_channel_to_preamble(tx,h2);
    noise2 = sqrt(noisePower/2) .* (randn(size(rxClean2)) + 1j*randn(size(rxClean2)));
    rx2 = rxClean2 + noise2;

    [~, lts2_rx2] = extract_lts(rx2,p);
    receivedLTS2_FFT = lts_fft(lts2_rx2);
    [Y2_100,X2_100] = extract_lts_100(receivedLTS2_FFT, knownLTS_FFT);
    HLS2_100 = pilot_ls_estimate(Y2_100, X2_100);
    cir2 = local_paper_cir_transform(HLS2_100);
    sigma2_2 = mean(abs(cir2(noiseStart:noiseEnd)).^2);

    for t = 1:length(T_values)
        Lhat_SUI2(trial,t) = local_sht(cir2, sigma2_2, G, T_values(t));
    end

    %% SUI-5
    h5 = sui5_channel(p);
    rxClean5 = apply_channel_to_preamble(tx,h5);
    noise5 = sqrt(noisePower/2) .* (randn(size(rxClean5)) + 1j*randn(size(rxClean5)));
    rx5 = rxClean5 + noise5;

    [~, lts2_rx5] = extract_lts(rx5,p);
    receivedLTS5_FFT = lts_fft(lts2_rx5);
    [Y5_100,X5_100] = extract_lts_100(receivedLTS5_FFT, knownLTS_FFT);
    HLS5_100 = pilot_ls_estimate(Y5_100, X5_100);
    cir5 = local_paper_cir_transform(HLS5_100);
    sigma2_5 = mean(abs(cir5(noiseStart:noiseEnd)).^2);

    for t = 1:length(T_values)
        Lhat_SUI5(trial,t) = local_sht(cir5, sigma2_5, G, T_values(t));
    end

    if mod(trial, 500) == 0
        fprintf('  Completed %d / %d trials\n', trial, numTrials);
    end
end

%% 5. COMPUTE PROBABILITY DISTRIBUTIONS
L = 0:G;
P_SUI2 = zeros(length(L), length(T_values));
P_SUI5 = zeros(length(L), length(T_values));

for t = 1:length(T_values)
    for i = 1:length(L)
        P_SUI2(i,t) = mean(Lhat_SUI2(:,t) == L(i));
        P_SUI5(i,t) = mean(Lhat_SUI5(:,t) == L(i));
    end
end

fprintf('\nSummary SUI-2:\n');
for t = 1:length(T_values)
    fprintf('  T = %2d : mean = %.3f, std = %.3f\n', T_values(t), mean(Lhat_SUI2(:,t)), std(Lhat_SUI2(:,t)));
end
fprintf('Summary SUI-5:\n');
for t = 1:length(T_values)
    fprintf('  T = %2d : mean = %.3f, std = %.3f\n', T_values(t), mean(Lhat_SUI5(:,t)), std(Lhat_SUI5(:,t)));
end

%% 6. PLOT FIGURE 4 (EXACT IEEE ACCESS 2022 MATCH)
fig4 = figure('Name', 'Figure 4 - SHT Distributions for Different T', ...
              'Color', 'w', 'Units', 'pixels', 'Position', [100 50 650 820]);

% (a) SUI-2 model
axA = subplot(2,1,1);
plot(L, P_SUI2(:,1), '-*k', 'LineWidth', 1.0, 'MarkerSize', 5); hold on;
plot(L, P_SUI2(:,2), '-sb', 'LineWidth', 1.0, 'MarkerSize', 4);
plot(L, P_SUI2(:,3), '-or', 'LineWidth', 1.0, 'MarkerSize', 4);
grid on; box on;
xlim([0 20]); ylim([0 1]);
xticks(0:2:20); yticks(0:0.1:1);
xlabel('CIR Length Estimates (T_s)', 'FontSize', 10);
ylabel('Probability', 'FontSize', 10);
legend('T = 5', 'T = 15', 'T = 25', 'Location', 'northeast', 'FontSize', 8);
text(0.5, -0.20, '(a)', 'Units', 'normalized', ...
     'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'normal');
set(axA, 'FontName', 'Times New Roman', 'FontSize', 10);

% Inset in (a) showing tail at higher estimates
posA = get(axA, 'Position');
axAI = axes('Position', [posA(1) + 0.50*posA(3), posA(2) + 0.38*posA(4), 0.38*posA(3), 0.46*posA(4)]);
plot(L, P_SUI2(:,1), '-*k', 'LineWidth', 0.8, 'MarkerSize', 3.5); hold on;
plot(L, P_SUI2(:,2), '-sb', 'LineWidth', 0.8, 'MarkerSize', 3);
plot(L, P_SUI2(:,3), '-or', 'LineWidth', 0.8, 'MarkerSize', 3);
grid on; box on;
xlim([0 20]); ylim([0 0.06]);
xticks(0:5:20); yticks(0:0.02:0.06);
set(axAI, 'FontName', 'Times New Roman', 'FontSize', 8);

% (b) SUI-5 model
axB = subplot(2,1,2);
plot(L, P_SUI5(:,1), '-*k', 'LineWidth', 1.0, 'MarkerSize', 5); hold on;
plot(L, P_SUI5(:,2), '-sb', 'LineWidth', 1.0, 'MarkerSize', 4);
plot(L, P_SUI5(:,3), '-or', 'LineWidth', 1.0, 'MarkerSize', 4);
grid on; box on;
xlim([0 65]); ylim([0 1]);
xticks(0:10:60); yticks(0:0.1:1);
xlabel('CIR Length Estimates (T_s)', 'FontSize', 10);
ylabel('Probability', 'FontSize', 10);
legend('T = 5', 'T = 15', 'T = 25', 'Location', 'northeast', 'FontSize', 8);
text(0.5, -0.20, '(b)', 'Units', 'normalized', ...
     'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'normal');
set(axB, 'FontName', 'Times New Roman', 'FontSize', 10);

% Inset in (b) showing false-alarm tail
posB = get(axB, 'Position');
axBI = axes('Position', [posB(1) + 0.08*posB(3), posB(2) + 0.38*posB(4), 0.38*posB(3), 0.46*posB(4)]);
plot(L, P_SUI5(:,1), '-*k', 'LineWidth', 0.8, 'MarkerSize', 3.5); hold on;
plot(L, P_SUI5(:,2), '-sb', 'LineWidth', 0.8, 'MarkerSize', 3);
plot(L, P_SUI5(:,3), '-or', 'LineWidth', 0.8, 'MarkerSize', 3);
grid on; box on;
xlim([0 65]); ylim([0 0.02]);
xticks(0:20:60); yticks(0:0.005:0.02);
set(axBI, 'FontName', 'Times New Roman', 'FontSize', 8);

%% 7. EXPORT RESULTS
evalDir = fullfile(projectRoot, 'evaluation');
if ~exist(evalDir, 'dir'), mkdir(evalDir); end

exportgraphics(fig4, fullfile(evalDir, 'paper_fig4_reproduction.png'), 'Resolution', 300);
exportgraphics(fig4, fullfile(projectRoot, 'plots', 'paper_fig4_reproduction.png'), 'Resolution', 300);

save(fullfile(evalDir, 'paper_fig4_results.mat'), 'L', 'P_SUI2', 'P_SUI5', 'T_values', 'SNR');
save(fullfile(projectRoot, 'plots', 'paper_fig4_results.mat'), 'L', 'P_SUI2', 'P_SUI5', 'T_values', 'SNR');

fprintf('\nFigure 4 saved successfully to:\n');
fprintf('  %s\n', fullfile(evalDir, 'paper_fig4_reproduction.png'));

end

%% ============================================================
% LOCAL FUNCTIONS
% =============================================================

function cir = local_paper_cir_transform(HLS100)
    M = 128;
    Np = 100;
    HLS100 = HLS100(:);
    
    H101 = zeros(101,1);
    H101(1:50)   = HLS100(1:50);
    H101(51)     = 0.5 * (HLS100(50) + HLS100(51)); % DC interpolation
    H101(52:101) = HLS100(51:100);

    n = (-Np/2):(Np/2);
    w = 0.5 + 0.5*cos(2*pi*n/M);
    w = w(:);
    Hw = H101 .* w;

    m = (0:M-1).';
    q = (-Np/2):(Np/2);
    Gmat = (1/sqrt(M)) * exp(1j*2*pi/M * ((m-M/2)*q));
    gCentered = Gmat * Hw;
    cir = fftshift(gCentered);
end

function Lhat = local_sht(cir, sigma2, G, T)
    threshold = T * sigma2;
    tapEnergy = abs(cir(1:G)).^2;
    Lhat = 0;
    for n = 1:G
        if tapEnergy(n) > threshold
            Lhat = n;
        end
    end
    if Lhat == 0
        Lhat = 1;
    end
end