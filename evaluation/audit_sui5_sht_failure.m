clear;
clc;
close all;

fprintf('\n');
fprintf('====================================================\n');
fprintf(' SUI-5 SHT FAILURE DIAGNOSTIC\n');
fprintf('====================================================\n');

p = parameters();

SNR = 20;
numTrials = 1000;

G = 64;
K = 56;
T = 15;

noiseStart = 65;
noiseEnd   = 120;

% -----------------------------------------------------
% Storage
% -----------------------------------------------------

Lhat = zeros(1,numTrials);

sigmaHat = zeros(1,numTrials);
threshold = zeros(1,numTrials);

E16 = zeros(1,numTrials);
E17 = zeros(1,numTrials);
E18 = zeros(1,numTrials);

E40 = zeros(1,numTrials);
E41 = zeros(1,numTrials);
E42 = zeros(1,numTrials);

% -----------------------------------------------------
% Fixed transmitter
% -----------------------------------------------------

Pall = ieee80216_pall();

preSeq = preamble_sequences(Pall);

[stsTime,ltsTime] = ...
    preamble_to_ifft(preSeq,p);

pre = build_80216_preamble( ...
    stsTime,ltsTime,p);

knownLTS = get_known_lts(ltsTime);

X_LTS = known_lts_fft(knownLTS);

% -----------------------------------------------------
% Monte Carlo
% -----------------------------------------------------

for trial = 1:numTrials

    % NEW SUI-5 channel realization
    h = sui5_channel(p);

    % Channel
    rxClean = apply_channel( ...
        pre.complete,h);

    % AWGN
    rx = add_awgn( ...
        rxClean,SNR);

    % Extract LTS
    [~,lts2] = extract_lts(rx,p);

    % FFT
    Y_LTS = lts_fft(lts2);

    % Extract 100 training tones
    [Y_p,X_p,~,~] = ...
        extract_lts_100(Y_LTS,X_LTS);

    % LS channel estimate
    H_LS = pilot_ls_estimate(Y_p,X_p);

    % LS -> CIR
    [cir,~,~] = ...
        ls_to_cir(H_LS,128,'hann');

    % -------------------------------------------------
    % Noise variance
    % -------------------------------------------------

    sigmaHat(trial) = mean( ...
        abs(cir(noiseStart:noiseEnd)).^2);

    % -------------------------------------------------
    % Threshold
    % -------------------------------------------------

    threshold(trial) = ...
        T * sigmaHat(trial);

    % -------------------------------------------------
    % Important tap energies
    % -------------------------------------------------

    tapEnergy = abs(cir).^2;

    E16(trial) = tapEnergy(16);
    E17(trial) = tapEnergy(17);
    E18(trial) = tapEnergy(18);

    E40(trial) = tapEnergy(40);
    E41(trial) = tapEnergy(41);
    E42(trial) = tapEnergy(42);

    % -------------------------------------------------
    % SHT
    % -------------------------------------------------

    [Lhat(trial),~,~,~] = ...
        sht_cir_length( ...
            cir, ...
            sigmaHat(trial), ...
            T, ...
            G);

end

% =====================================================
% BASIC STATISTICS
% =====================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' OVERALL STATISTICS\n');
fprintf('====================================================\n');

fprintf('Mean Lhat       = %.4f\n',mean(Lhat));
fprintf('Std Lhat        = %.4f\n',std(Lhat));

fprintf('Mean sigma^2    = %.6e\n',mean(sigmaHat));
fprintf('Mean threshold  = %.6e\n',mean(threshold));

fprintf('\n');

% =====================================================
% FAILURE COUNTS
% =====================================================

fail17 = sum(Lhat == 17);
fail18 = sum(Lhat == 18);
fail40 = sum(Lhat == 40);
fail41 = sum(Lhat == 41);
fail42 = sum(Lhat == 42);

fprintf('Lhat = 17 : %d trials\n',fail17);
fprintf('Lhat = 18 : %d trials\n',fail18);

fprintf('Lhat = 40 : %d trials\n',fail40);
fprintf('Lhat = 41 : %d trials\n',fail41);
fprintf('Lhat = 42 : %d trials\n',fail42);

% =====================================================
% TAP ENERGY STATISTICS
% =====================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' TAP ENERGY STATISTICS\n');
fprintf('====================================================\n');

fprintf('\n');
fprintf('Tap     Mean Energy       Median Energy\n');

fprintf('16      %.6e       %.6e\n', ...
    mean(E16),median(E16));

fprintf('17      %.6e       %.6e\n', ...
    mean(E17),median(E17));

fprintf('18      %.6e       %.6e\n', ...
    mean(E18),median(E18));

fprintf('40      %.6e       %.6e\n', ...
    mean(E40),median(E40));

fprintf('41      %.6e       %.6e\n', ...
    mean(E41),median(E41));

fprintf('42      %.6e       %.6e\n', ...
    mean(E42),median(E42));

fprintf('\n');
fprintf('Mean threshold = %.6e\n', ...
    mean(threshold));

% =====================================================
% FAILURE CASE ANALYSIS
% =====================================================

failureMask = (Lhat == 17 | Lhat == 18);

fprintf('\n');
fprintf('====================================================\n');
fprintf(' FAILURE-CASE STATISTICS\n');
fprintf('====================================================\n');

fprintf('Failure trials = %d / %d (%.2f%%)\n', ...
    sum(failureMask), ...
    numTrials, ...
    100*mean(failureMask));

if any(failureMask)

    fprintf('\n');

    fprintf('Failure mean sigma^2 = %.6e\n', ...
        mean(sigmaHat(failureMask)));

    fprintf('Failure mean threshold = %.6e\n', ...
        mean(threshold(failureMask)));

    fprintf('\n');

    fprintf('Failure mean E17 = %.6e\n', ...
        mean(E17(failureMask)));

    fprintf('Failure mean E18 = %.6e\n', ...
        mean(E18(failureMask)));

    fprintf('Failure mean E40 = %.6e\n', ...
        mean(E40(failureMask)));

    fprintf('Failure mean E41 = %.6e\n', ...
        mean(E41(failureMask)));

    fprintf('Failure mean E42 = %.6e\n', ...
        mean(E42(failureMask)));

end

% =====================================================
% THRESHOLD COMPARISON
% =====================================================

figure;

plot(threshold,'LineWidth',1);

hold on;

plot(E40,'LineWidth',1);

plot(E41,'LineWidth',1);

plot(E42,'LineWidth',1);

xlabel('Trial');

ylabel('Energy');

title('SUI-5: SHT Threshold vs Final-Path Tap Energies');

legend( ...
    'SHT Threshold', ...
    'Tap 40 Energy', ...
    'Tap 41 Energy', ...
    'Tap 42 Energy');

grid on;

% =====================================================
% FAILURE CASE PLOT
% =====================================================

figure;

failureIndices = find(failureMask);

if ~isempty(failureIndices)

    plot(failureIndices,E40(failureMask),'.');

    hold on;

    plot(failureIndices, ...
         threshold(failureMask),'.');

    xlabel('Failure Trial');

    ylabel('Energy');

    title('SUI-5 Failure Trials: Tap 40 vs SHT Threshold');

    legend( ...
        'Tap 40 Energy', ...
        'SHT Threshold');

    grid on;

else

    text(0.2,0.5,'No Lhat = 17/18 failures detected');

    axis off;

end

fprintf('\n');
fprintf('Diagnostic complete.\n');