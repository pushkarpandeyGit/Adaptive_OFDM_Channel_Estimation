clear;
clc;
close all;

fprintf('\n====================================================\n');
fprintf(' SUI-5 NOISELESS CIR RECOVERY AUDIT\n');
fprintf('====================================================\n\n');

p = parameters();

numTrials = 100;

actual40 = zeros(1,numTrials);
estimated41 = zeros(1,numTrials);

relativeError = zeros(1,numTrials);

for trial = 1:numTrials

    fprintf('Trial %d / %d\n', trial, numTrials);

    % ------------------------------------------------
    % NEW stochastic SUI-5 channel
    % ------------------------------------------------
    h = sui5_channel(p);

    % Actual path at physical delay 40 samples
    actual40(trial) = abs(h(41))^2;

    % ------------------------------------------------
    % Generate IEEE 802.16 preamble
    % ------------------------------------------------
    Pall = ieee80216_pall();

    preSeq = preamble_sequences(Pall);

    [stsTime,ltsTime] = preamble_to_ifft(preSeq,p);

    pre = build_80216_preamble(stsTime,ltsTime,p);

    % ------------------------------------------------
    % Channel ONLY
    % NO AWGN
    % ------------------------------------------------
    rx = apply_channel(pre.complete,h);

    % ------------------------------------------------
    % Extract LTS
    % ------------------------------------------------
    [~,lts2] = extract_lts(rx,p);

    Y_LTS = lts_fft(lts2);

    knownLTS = get_known_lts(ltsTime);
    X_LTS = known_lts_fft(knownLTS);

    [Y_p,X_p,~,~] = extract_lts_100(Y_LTS,X_LTS);

    % ------------------------------------------------
    % LS channel estimate
    % ------------------------------------------------
    H_LS = pilot_ls_estimate(Y_p,X_p);

    % ------------------------------------------------
    % Convert to CIR
    % ------------------------------------------------
    [cir,~,~] = ls_to_cir(H_LS,128,'hann');

    % ------------------------------------------------
    % Estimated physical path
    %
    % Physical delay 40 samples
    % corresponds to MATLAB index 41
    % ------------------------------------------------
    estimated41(trial) = abs(cir(41))^2;

    % ------------------------------------------------
    % Relative recovery error
    % ------------------------------------------------
    if actual40(trial) > 0
        relativeError(trial) = ...
            abs(estimated41(trial)-actual40(trial)) ...
            / actual40(trial);
    end

end

% ====================================================
% STATISTICS
% ====================================================

fprintf('\n====================================================\n');
fprintf(' NOISELESS RECOVERY RESULTS\n');
fprintf('====================================================\n\n');

fprintf('Mean actual E40        = %.6e\n',mean(actual40));
fprintf('Mean estimated E41     = %.6e\n',mean(estimated41));
fprintf('Mean relative error    = %.6f\n',mean(relativeError));

fprintf('\n');

fprintf('Median actual E40      = %.6e\n',median(actual40));
fprintf('Median estimated E41   = %.6e\n',median(estimated41));
fprintf('Median relative error  = %.6f\n',median(relativeError));

% ====================================================
% CORRELATION
% ====================================================

R = corrcoef(actual40,estimated41);

fprintf('\nCorrelation coefficient = %.6f\n',R(1,2));

% ====================================================
% PLOT 1
% ====================================================

figure;

scatter(actual40,estimated41,'filled');

hold on;

maxValue = max([actual40 estimated41]);

plot([0 maxValue],[0 maxValue],'-');

xlabel('Actual |h(41)|^2');
ylabel('Estimated |\hat{g}(41)|^2');

title('SUI-5 Noiseless Path-Energy Recovery');

grid on;

% ====================================================
% PLOT 2
% ====================================================

figure;

plot(1:numTrials,actual40,'-');
hold on;

plot(1:numTrials,estimated41,'-');

xlabel('Trial');
ylabel('Path energy');

title('Actual vs Estimated SUI-5 40-Sample Path');

legend('Actual |h(41)|^2', ...
       'Estimated |\hat{g}(41)|^2');

grid on;

% ====================================================
% SHOW WORST CASE
% ====================================================

[~,idx] = max(relativeError);

fprintf('\n====================================================\n');
fprintf(' WORST RECOVERY TRIAL\n');
fprintf('====================================================\n');

fprintf('Trial                  = %d\n',idx);
fprintf('Actual E40             = %.6e\n',actual40(idx));
fprintf('Estimated E41          = %.6e\n',estimated41(idx));
fprintf('Relative error         = %.6f\n',relativeError(idx));

fprintf('\n====================================================\n');
fprintf(' AUDIT COMPLETE\n');
fprintf('====================================================\n');