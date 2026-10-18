clear;
clc;
close all;

fprintf('\n====================================================\n');
fprintf(' SUI-5 SUCCESS / FAILURE SIDE-BY-SIDE AUDIT\n');
fprintf('====================================================\n\n');

p = parameters();

SNR = 20;
T = p.T;
G = p.G;

numTrials = 1000;

% ----------------------------------------------------
% Storage
% ----------------------------------------------------
foundSuccess = false;
foundFailure = false;

successData = struct();
failureData = struct();

% ----------------------------------------------------
% Monte Carlo
% ----------------------------------------------------
for trial = 1:numTrials

    fprintf('Trial %d / %d\n', trial, numTrials);

    % IMPORTANT:
    % Generate a NEW SUI-5 channel every trial.
    h = sui5_channel(p);

    % ------------------------------------------------
    % Transmitter
    % ------------------------------------------------
    Pall = ieee80216_pall();
    preSeq = preamble_sequences(Pall);

    [stsTime, ltsTime] = preamble_to_ifft(preSeq, p);

    pre = build_80216_preamble(stsTime, ltsTime, p);

    % ------------------------------------------------
    % Channel + AWGN
    % ------------------------------------------------
    rxClean = apply_channel(pre.complete, h);
    rx = add_awgn(rxClean, SNR);

    % ------------------------------------------------
    % Receiver
    % ------------------------------------------------
    [~, lts2] = extract_lts(rx, p);

    Y_LTS = lts_fft(lts2);

    knownLTS = get_known_lts(ltsTime);
    X_LTS = known_lts_fft(knownLTS);

    [Y_p, X_p, ~, ~] = extract_lts_100(Y_LTS, X_LTS);

    % ------------------------------------------------
    % LS frequency-domain channel estimate
    % ------------------------------------------------
    H_LS = pilot_ls_estimate(Y_p, X_p);

    % ------------------------------------------------
    % Convert to CIR
    % ------------------------------------------------
    [cir, ~, ~] = ls_to_cir(H_LS, 128, 'hann');

    % ------------------------------------------------
    % Noise variance and SHT
    % ------------------------------------------------
    noiseVariance = estimate_noise_variance( ...
        cir, ...
        p.noiseTapStart, ...
        p.noiseTapEnd);

    [L_hat, threshold, tapEnergy, decision] = ...
        sht_cir_length(cir, noiseVariance, T, G);

    % ------------------------------------------------
    % IMPORTANT INDEXING
    %
    % SUI-5 delays:
    %   0 samples  -> MATLAB tap 1
    %   16 samples -> MATLAB tap 17
    %   40 samples -> MATLAB tap 41
    %
    % ------------------------------------------------

    % Failure = final 40-sample-delay path not detected
    isFailure = L_hat < 41;

    % ------------------------------------------------
    % Save first failure
    % ------------------------------------------------
    if isFailure && ~foundFailure

        foundFailure = true;

        failureData.trial = trial;
        failureData.h = h;
        failureData.cir = cir;
        failureData.L_hat = L_hat;
        failureData.threshold = threshold;
        failureData.noiseVariance = noiseVariance;
        failureData.tapEnergy = tapEnergy;
        failureData.decision = decision;

        fprintf('\n*** FIRST FAILURE FOUND ***\n');
        fprintf('Trial = %d\n', trial);
        fprintf('L_hat = %d\n', L_hat);
        fprintf('Noise variance = %.6e\n', noiseVariance);
        fprintf('Threshold = %.6e\n', threshold);
        fprintf('E16 = %.6e\n', tapEnergy(16));
        fprintf('E17 = %.6e\n', tapEnergy(17));
        fprintf('E18 = %.6e\n', tapEnergy(18));
        fprintf('E40 = %.6e\n', tapEnergy(40));
        fprintf('E41 = %.6e\n', tapEnergy(41));
        fprintf('E42 = %.6e\n', tapEnergy(42));
    end

    % ------------------------------------------------
    % Save first success
    % ------------------------------------------------
    if ~isFailure && ~foundSuccess

        foundSuccess = true;

        successData.trial = trial;
        successData.h = h;
        successData.cir = cir;
        successData.L_hat = L_hat;
        successData.threshold = threshold;
        successData.noiseVariance = noiseVariance;
        successData.tapEnergy = tapEnergy;
        successData.decision = decision;

        fprintf('\n*** FIRST SUCCESS FOUND ***\n');
        fprintf('Trial = %d\n', trial);
        fprintf('L_hat = %d\n', L_hat);
        fprintf('Noise variance = %.6e\n', noiseVariance);
        fprintf('Threshold = %.6e\n', threshold);
        fprintf('E16 = %.6e\n', tapEnergy(16));
        fprintf('E17 = %.6e\n', tapEnergy(17));
        fprintf('E18 = %.6e\n', tapEnergy(18));
        fprintf('E40 = %.6e\n', tapEnergy(40));
        fprintf('E41 = %.6e\n', tapEnergy(41));
        fprintf('E42 = %.6e\n', tapEnergy(42));
    end

    % Stop once both are found
    if foundSuccess && foundFailure
        break;
    end
end

% ====================================================
% PRINT COMPARISON
% ====================================================

fprintf('\n====================================================\n');
fprintf(' SUCCESS VS FAILURE\n');
fprintf('====================================================\n\n');

if foundSuccess

    fprintf('SUCCESS TRIAL: %d\n', successData.trial);
    fprintf('L_hat         = %d\n', successData.L_hat);
    fprintf('sigma^2       = %.6e\n', successData.noiseVariance);
    fprintf('threshold     = %.6e\n', successData.threshold);

    fprintf('\nTap energies:\n');
    fprintf('E16 = %.6e\n', successData.tapEnergy(16));
    fprintf('E17 = %.6e\n', successData.tapEnergy(17));
    fprintf('E18 = %.6e\n', successData.tapEnergy(18));
    fprintf('E40 = %.6e\n', successData.tapEnergy(40));
    fprintf('E41 = %.6e\n', successData.tapEnergy(41));
    fprintf('E42 = %.6e\n', successData.tapEnergy(42));

else
    fprintf('No success found.\n');
end

fprintf('\n----------------------------------------------------\n');

if foundFailure

    fprintf('FAILURE TRIAL: %d\n', failureData.trial);
    fprintf('L_hat         = %d\n', failureData.L_hat);
    fprintf('sigma^2       = %.6e\n', failureData.noiseVariance);
    fprintf('threshold     = %.6e\n', failureData.threshold);

    fprintf('\nTap energies:\n');
    fprintf('E16 = %.6e\n', failureData.tapEnergy(16));
    fprintf('E17 = %.6e\n', failureData.tapEnergy(17));
    fprintf('E18 = %.6e\n', failureData.tapEnergy(18));
    fprintf('E40 = %.6e\n', failureData.tapEnergy(40));
    fprintf('E41 = %.6e\n', failureData.tapEnergy(41));
    fprintf('E42 = %.6e\n', failureData.tapEnergy(42));

else
    fprintf('No failure found.\n');
end

% ====================================================
% PLOT 1: ACTUAL SUI-5 CHANNEL
% ====================================================

figure;

stem(0:length(successData.h)-1, abs(successData.h), 'filled');

xlabel('Delay (samples)');
ylabel('|h[n]|');
title('Actual SUI-5 Channel - Successful Trial');
grid on;

% ====================================================
% PLOT 2: ESTIMATED CIR - SUCCESS
% ====================================================

figure;

stem(1:128, abs(successData.cir), 'filled');
hold on;

yline(sqrt(successData.threshold));

xline(17);
xline(41);

xlabel('MATLAB tap index');
ylabel('|CIR tap|');
title(sprintf('Estimated CIR - SUCCESS, Lhat = %d', ...
    successData.L_hat));

legend('Estimated CIR', 'sqrt(threshold)', ...
       'Tap 17', 'Tap 41');

grid on;
xlim([1 70]);

% ====================================================
% PLOT 3: ESTIMATED CIR - FAILURE
% ====================================================

figure;

stem(1:128, abs(failureData.cir), 'filled');
hold on;

yline(sqrt(failureData.threshold));

xline(17);
xline(41);

xlabel('MATLAB tap index');
ylabel('|CIR tap|');
title(sprintf('Estimated CIR - FAILURE, Lhat = %d', ...
    failureData.L_hat));

legend('Estimated CIR', 'sqrt(threshold)', ...
       'Tap 17', 'Tap 41');

grid on;
xlim([1 70]);

% ====================================================
% PLOT 4: SUCCESS VS FAILURE TAP ENERGY
% ====================================================

figure;

taps = [16 17 18 40 41 42];

successEnergy = successData.tapEnergy(taps);
failureEnergy = failureData.tapEnergy(taps);

bar(taps, [successEnergy(:), failureEnergy(:)]);

hold on;

yline(successData.threshold);

xlabel('CIR tap index');
ylabel('|g[n]|^2');
title('SUI-5 Tap Energy: Success vs Failure');

legend('Success', 'Failure', 'Threshold');

grid on;

% ====================================================
% ACTUAL CHANNEL TAP VALUES
% ====================================================

fprintf('\n====================================================\n');
fprintf(' ACTUAL SUI-5 CHANNEL COEFFICIENTS\n');
fprintf('====================================================\n');

fprintf('\nSUCCESS:\n');

for n = 1:length(successData.h)
    if abs(successData.h(n)) > 0
        fprintf('Delay %2d -> |h| = %.6e\n', ...
            n-1, abs(successData.h(n)));
    end
end

fprintf('\nFAILURE:\n');

for n = 1:length(failureData.h)
    if abs(failureData.h(n)) > 0
        fprintf('Delay %2d -> |h| = %.6e\n', ...
            n-1, abs(failureData.h(n)));
    end
end

fprintf('\n====================================================\n');
fprintf(' AUDIT COMPLETE\n');
fprintf('====================================================\n');