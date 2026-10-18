clc;
clear;

addpath(genpath('/MATLAB Drive/Adaptive_OFDM_Channel_Estimation'));

fprintf('\n');
fprintf('====================================================\n');
fprintf(' SHT CIR LENGTH ESTIMATION AUDIT\n');
fprintf('====================================================\n');

%% Parameters

G = 64;
K = 56;
T = 15;

%% ----------------------------------------------------
% Test 1: Threshold calculation
% -----------------------------------------------------

fprintf('\n[1] Threshold calculation\n');

noiseVariance = 0.1;

[cirLength, threshold, tapEnergy, decision] = ...
    sht_cir_length(zeros(1,128),noiseVariance,T,G);

expectedThreshold = T * noiseVariance;

fprintf('Expected threshold : %.6f\n',expectedThreshold);
fprintf('Actual threshold   : %.6f\n',threshold);

assert(abs(threshold-expectedThreshold) < 1e-12);

fprintf('PASS: threshold = T * noise variance.\n');

%% ----------------------------------------------------
% Test 2: No significant channel taps
% -----------------------------------------------------

fprintf('\n[2] No significant taps\n');

noiseVariance = 0.1;

cir = zeros(1,128);

[L_hat,threshold,tapEnergy,decision] = ...
    sht_cir_length(cir,noiseVariance,T,G);

fprintf('Estimated CIR length : %d\n',L_hat);

assert(L_hat == 0);
assert(all(decision == 0));

fprintf('PASS: no significant taps gives L_hat = 0.\n');

%% ----------------------------------------------------
% Test 3: Single significant tap
% -----------------------------------------------------

fprintf('\n[3] Single significant tap\n');

cir = zeros(1,128);

% Threshold = 1.5
% Tap energy = 4
% Therefore tap 5 must be detected.

noiseVariance = 0.1;
cir(5) = 2;

[L_hat,threshold,tapEnergy,decision] = ...
    sht_cir_length(cir,noiseVariance,T,G);

fprintf('Threshold          : %.6f\n',threshold);
fprintf('Tap-5 energy       : %.6f\n',tapEnergy(5));
fprintf('Estimated CIR length: %d\n',L_hat);

assert(L_hat == 5);
assert(decision(5) == 1);

fprintf('PASS: significant tap detected correctly.\n');

%% ----------------------------------------------------
% Test 4: Multiple significant taps
% -----------------------------------------------------

fprintf('\n[4] Multiple significant taps\n');

cir = zeros(1,128);

noiseVariance = 0.1;

cir(3)  = 2;
cir(7)  = 2;
cir(12) = 2;

[L_hat,threshold,tapEnergy,decision] = ...
    sht_cir_length(cir,noiseVariance,T,G);

fprintf('Detected taps      : %d\n',sum(decision));
fprintf('Estimated CIR length: %d\n',L_hat);

assert(decision(3) == 1);
assert(decision(7) == 1);
assert(decision(12) == 1);
assert(L_hat == 12);

fprintf('PASS: last significant tap determines L_hat.\n');

%% ----------------------------------------------------
% Test 5: Below-threshold tap
% -----------------------------------------------------

fprintf('\n[5] Below-threshold test\n');

cir = zeros(1,128);

noiseVariance = 0.1;

% Threshold = 1.5
% Energy = 1
% Therefore tap must NOT be detected.

cir(10) = 1;

[L_hat,threshold,tapEnergy,decision] = ...
    sht_cir_length(cir,noiseVariance,T,G);

fprintf('Threshold       : %.6f\n',threshold);
fprintf('Tap energy      : %.6f\n',tapEnergy(10));
fprintf('Decision        : %d\n',decision(10));

assert(decision(10) == 0);
assert(L_hat == 0);

fprintf('PASS: below-threshold tap rejected.\n');

%% ----------------------------------------------------
% Test 6: Above-threshold test
% -----------------------------------------------------

fprintf('\n[6] Above-threshold test\n');

cir = zeros(1,128);

% Threshold = 1.5
% Energy = 4

cir(10) = 2;

[L_hat,threshold,tapEnergy,decision] = ...
    sht_cir_length(cir,noiseVariance,T,G);

fprintf('Threshold       : %.6f\n',threshold);
fprintf('Tap energy      : %.6f\n',tapEnergy(10));
fprintf('Decision        : %d\n',decision(10));

assert(decision(10) == 1);
assert(L_hat == 10);

fprintf('PASS: above-threshold tap detected.\n');

%% ----------------------------------------------------
% Test 7: Boundary tap 64
% -----------------------------------------------------

fprintf('\n[7] G = 64 boundary test\n');

cir = zeros(1,128);

cir(64) = 2;

[L_hat,~,~,decision] = ...
    sht_cir_length(cir,noiseVariance,T,G);

assert(decision(64) == 1);
assert(L_hat == 64);

fprintf('PASS: tap 64 is tested.\n');

%% ----------------------------------------------------
% Test 8: Tap 65 must NOT be tested
% -----------------------------------------------------

fprintf('\n[8] Tap 65 exclusion test\n');

cir = zeros(1,128);

cir(65) = 100;

[L_hat,~,~,decision] = ...
    sht_cir_length(cir,noiseVariance,T,G);

assert(L_hat == 0);
assert(all(decision == 0));

fprintf('PASS: tap 65 is outside the SHT test region.\n');

%% ----------------------------------------------------
% Test 9: Complex-valued tap
% -----------------------------------------------------

fprintf('\n[9] Complex tap energy test\n');

cir = zeros(1,128);

% |3+4j|^2 = 25

cir(8) = 3 + 4j;

noiseVariance = 1;

% Threshold = 15
% Energy = 25 -> detected.

[L_hat,threshold,tapEnergy,decision] = ...
    sht_cir_length(cir,noiseVariance,T,G);

fprintf('Threshold       : %.6f\n',threshold);
fprintf('Tap energy      : %.6f\n',tapEnergy(8));

assert(abs(tapEnergy(8)-25) < 1e-12);
assert(decision(8) == 1);
assert(L_hat == 8);

fprintf('PASS: complex tap uses |g_n|^2 correctly.\n');

%% ----------------------------------------------------
% Test 10: Theoretical false-alarm probability
% -----------------------------------------------------

fprintf('\n[10] Theoretical false-alarm probability\n');

Pf = (K/(T+K))^K;

expectedPf = 1.690682214704e-6;

fprintf('K                 : %d\n',K);
fprintf('T                 : %d\n',T);
fprintf('Calculated P_f    : %.12e\n',Pf);
fprintf('Expected P_f      : %.12e\n',expectedPf);

assert(abs(Pf-expectedPf) < 1e-10);

fprintf('PASS: theoretical P_f matches paper value.\n');

%% ----------------------------------------------------
% Test 11: Invalid inputs
% -----------------------------------------------------

fprintf('\n[11] Invalid input tests\n');

caught = false;

try
    sht_cir_length([],0.1,15,64);
catch
    caught = true;
end

assert(caught);
fprintf('PASS: empty CIR rejected.\n');

caught = false;

try
    sht_cir_length(zeros(1,128),-1,15,64);
catch
    caught = true;
end

assert(caught);
fprintf('PASS: negative noise variance rejected.\n');

caught = false;

try
    sht_cir_length(zeros(1,128),0.1,0,64);
catch
    caught = true;
end

assert(caught);
fprintf('PASS: invalid T rejected.\n');

caught = false;

try
    sht_cir_length(zeros(1,64),0.1,15,65);
catch
    caught = true;
end

assert(caught);
fprintf('PASS: G larger than CIR rejected.\n');

%% ----------------------------------------------------
% Final
% -----------------------------------------------------

fprintf('\n');
fprintf('====================================================\n');
fprintf(' ALL SHT TESTS PASSED\n');
fprintf('====================================================\n');