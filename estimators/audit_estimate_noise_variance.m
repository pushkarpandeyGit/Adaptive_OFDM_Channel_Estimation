clc;
clear;

addpath(genpath('/MATLAB Drive/Adaptive_OFDM_Channel_Estimation'));

fprintf('\n');
fprintf('====================================================\n');
fprintf(' NOISE VARIANCE ESTIMATION AUDIT\n');
fprintf('====================================================\n');

%% Test 1: Known noise variance

rng(1);

trueVariance = 0.25;

noise = sqrt(trueVariance/2) * ...
    (randn(1,56) + 1j*randn(1,56));

cir = zeros(1,128);
cir(65:120) = noise;

[estimatedVariance, noiseTaps] = ...
    estimate_noise_variance(cir,65,120);

fprintf('\n[1] Known noise variance test\n');
fprintf('True variance      : %.6f\n',trueVariance);
fprintf('Estimated variance : %.6f\n',estimatedVariance);

assert(length(noiseTaps) == 56);

% Finite-sample estimate does not need to be exactly equal.
assert(abs(estimatedVariance-trueVariance) < 0.1);

fprintf('PASS: noise variance estimated correctly.\n');

%% Test 2: Direct equation

fprintf('\n[2] Direct equation test\n');

expectedVariance = mean(abs(cir(65:120)).^2);

errorValue = abs(estimatedVariance-expectedVariance);

fprintf('Equation error : %.3e\n',errorValue);

assert(errorValue < 1e-12);

fprintf('PASS: variance equation implemented correctly.\n');

%% Test 3: Zero noise

fprintf('\n[3] Zero-noise test\n');

cirZero = zeros(1,128);

[varianceZero, tapsZero] = ...
    estimate_noise_variance(cirZero,65,120);

assert(varianceZero == 0);
assert(length(tapsZero) == 56);

fprintf('PASS: zero noise gives zero variance.\n');

%% Test 4: Isolation test

fprintf('\n[4] Noise-tap isolation test\n');

cirIsolation = zeros(1,128);

cirIsolation(65:120) = 2;

% Variance estimate is mean(|2|^2) = 4.
[varianceIsolation,~] = ...
    estimate_noise_variance(cirIsolation,65,120);

fprintf('Estimated variance : %.6f\n',varianceIsolation);

assert(abs(varianceIsolation-4) < 1e-12);

fprintf('PASS: correct taps are being used.\n');

%% Test 5: K = 56

fprintf('\n[5] K test\n');

[~,noiseTaps] = ...
    estimate_noise_variance(zeros(1,128),65,120);

K = length(noiseTaps);

fprintf('K = %d\n',K);

assert(K == 56);

fprintf('PASS: K = 56.\n');

%% Test 6: Row/column input

fprintf('\n[6] Row/column consistency test\n');

cirColumn = zeros(128,1);
cirColumn(65:120) = noise.';

[vRow,~] = estimate_noise_variance(cirColumn.',65,120);
[vColumn,~] = estimate_noise_variance(cirColumn,65,120);

assert(abs(vRow-vColumn) < 1e-12);

fprintf('PASS: row and column inputs give same result.\n');

%% Test 7: Invalid ranges

fprintf('\n[7] Invalid input tests\n');

caught = false;

try
    estimate_noise_variance(zeros(1,128),0,120);
catch
    caught = true;
end

assert(caught);
fprintf('PASS: invalid start tap rejected.\n');

caught = false;

try
    estimate_noise_variance(zeros(1,128),65,129);
catch
    caught = true;
end

assert(caught);
fprintf('PASS: invalid end tap rejected.\n');

caught = false;

try
    estimate_noise_variance(zeros(1,128),120,65);
catch
    caught = true;
end

assert(caught);
fprintf('PASS: reversed tap range rejected.\n');

%% Final

fprintf('\n');
fprintf('====================================================\n');
fprintf(' ALL NOISE VARIANCE TESTS PASSED\n');
fprintf('====================================================\n');