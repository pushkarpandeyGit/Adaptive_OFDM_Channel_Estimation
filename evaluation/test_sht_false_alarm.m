function test_sht_false_alarm()

clc;
fprintf('\n');
fprintf('============================================================\n');
fprintf(' STEP 40B - SHT FALSE-ALARM PROBABILITY TEST\n');
fprintf('============================================================\n');

%% -----------------------------------------------------------
% Parameters
% ------------------------------------------------------------

K = 56;                 % Number of noise-only CIR samples
T = 15;                 % SHT scaling factor
G = 64;                 % Maximum channel tap region

numTrials = 1000000;   % 1 million Monte Carlo trials

rng(42);                % Reproducibility

%% -----------------------------------------------------------
% Theoretical false-alarm probability
%
% Paper:
%
%        Pf = (K/(T+K))^K
%
% ------------------------------------------------------------

theoreticalPf = (K/(T + K))^K;

%% -----------------------------------------------------------
% Monte Carlo test
% ------------------------------------------------------------

falseAlarms = 0;

for trial = 1:numTrials

    % -------------------------------------------------------
    % Generate complex Gaussian noise CIR
    %
    % E{|n|^2} = 1
    % -------------------------------------------------------

    noiseCir = (randn(1, G + K) + ...
                1j*randn(1, G + K)) / sqrt(2);

    % -------------------------------------------------------
    % Estimate noise variance from taps 65:120
    % -------------------------------------------------------

    noiseVariance = estimate_noise_variance( ...
        noiseCir, ...
        G + 1, ...
        G + K);

    % -------------------------------------------------------
    % Run Sequential Hypothesis Test
    %
    % Correct function argument order:
    %
    % sht_cir_length(cir, noiseVariance, T, G)
    % -------------------------------------------------------

    Lhat = sht_cir_length( ...
        noiseCir, ...
        noiseVariance, ...
        T, ...
        G);

    % -------------------------------------------------------
    % A false alarm occurs if SHT detects at least one
    % channel tap in the noise-only region.
    % -------------------------------------------------------

    if Lhat > 0
        falseAlarms = falseAlarms + 1;
    end

end

%% -----------------------------------------------------------
% Measured probability
% ------------------------------------------------------------

measuredPf = falseAlarms / numTrials;

%% -----------------------------------------------------------
% Expected number of false alarms
% ------------------------------------------------------------

expectedFalseAlarms = numTrials * theoreticalPf;

%% -----------------------------------------------------------
% Display results
% ------------------------------------------------------------

fprintf('\n');
fprintf('============================================================\n');
fprintf(' RESULT\n');
fprintf('============================================================\n');

fprintf('False alarms          = %d\n', falseAlarms);
fprintf('Total trials           = %d\n', numTrials);

fprintf('Measured Pf            = %.12e\n', measuredPf);
fprintf('Theoretical Pf         = %.12e\n', theoreticalPf);

fprintf('Expected false alarms  = %.4f\n', expectedFalseAlarms);

fprintf('\n');

%% -----------------------------------------------------------
% Basic interpretation
% ------------------------------------------------------------

if measuredPf == 0

    fprintf('Measured Pf is zero in this Monte Carlo run.\n');
    fprintf('This is possible because the theoretical Pf is very small.\n');

elseif measuredPf <= 5 * theoreticalPf

    fprintf('SHT false-alarm result is within a reasonable Monte Carlo range.\n');

else

    fprintf('WARNING: Measured false-alarm probability is larger\n');
    fprintf('than expected. Investigate SHT implementation.\n');

end

fprintf('============================================================\n');
fprintf('\n');

end