clc;
clear;

addpath(genpath('/MATLAB Drive/Adaptive_OFDM_Channel_Estimation'));

fprintf('\n');
fprintf('====================================================\n');
fprintf(' 100-TONE LTS EXTRACTION AUDIT\n');
fprintf('====================================================\n');

%% Load parameters

p = parameters;

%% Generate known LTS

Pall = ieee80216_pall;
preamble = preamble_sequences(Pall);

[stsTime,ltsTime] = preamble_to_ifft(preamble,p);

knownLTS = get_known_lts(ltsTime);

X_LTS = known_lts_fft(knownLTS);

%% Extract

[Y_p,X_p,ltsLocations,ltsBins] = ...
    extract_lts_100(X_LTS,X_LTS);

%% Tests

fprintf('\n[1] Number of LTS tones\n');

assert(length(ltsLocations) == 100);
fprintf('PASS: exactly 100 LTS tones.\n');

fprintf('\n[2] LTS locations\n');

expectedLocations = [-100:2:-2 2:2:100];

assert(isequal(ltsLocations,expectedLocations));
fprintf('PASS: LTS locations are correct.\n');

fprintf('\n[3] No DC subcarrier\n');

assert(~ismember(0,ltsLocations));
fprintf('PASS: DC is excluded.\n');

fprintf('\n[4] Bin count\n');

assert(length(ltsBins) == 100);
fprintf('PASS: exactly 100 FFT bins extracted.\n');

fprintf('\n[5] No duplicate FFT bins\n');

assert(length(unique(ltsBins)) == 100);
fprintf('PASS: all 100 FFT bins are unique.\n');

fprintf('\n[6] Bin range\n');

assert(all(ltsBins >= 1 & ltsBins <= 128));
fprintf('PASS: all bins lie inside 1...128.\n');

fprintf('\n[7] Known LTS extraction\n');

assert(max(abs(Y_p-X_p)) < 1e-12);
fprintf('PASS: extracted known LTS values match.\n');

fprintf('\n[8] Nonzero known values\n');

assert(all(abs(X_p) > 1e-12));
fprintf('PASS: all selected LTS training values are nonzero.\n');

fprintf('\n[9] Mapping table\n');

fprintf('\n');
fprintf(' Original k     MATLAB FFT bin\n');
fprintf(' -----------------------------\n');

for n = 1:length(ltsLocations)

    if n <= 5 || n > 95
        fprintf('%10d %18d\n', ...
            ltsLocations(n),ltsBins(n));
    end

end

fprintf('\n');
fprintf('====================================================\n');
fprintf(' ALL 100-TONE LTS TESTS PASSED\n');
fprintf('====================================================\n');