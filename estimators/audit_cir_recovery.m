clc;
clear;

addpath(genpath('/MATLAB Drive/Adaptive_OFDM_Channel_Estimation'));

fprintf('\n');
fprintf('====================================================\n');
fprintf(' STEP 32B: CIR RECOVERY DIAGNOSTIC\n');
fprintf(' NO AWGN\n');
fprintf('====================================================\n');

%% 1. Parameters

p = parameters;

%% 2. Generate LTS

Pall = ieee80216_pall;

preamble = preamble_sequences(Pall);

[stsTime,ltsTime] = ...
    preamble_to_ifft(preamble,p);

knownLTS = get_known_lts(ltsTime);

X_LTS = known_lts_fft(knownLTS);

%% 3. Build complete preamble

preambleStruct = ...
    build_80216_preamble(stsTime,ltsTime,p);

tx = preambleStruct.complete;

%% 4. Create channel

hTrue = sui2_channel(p);

%% 5. Pass through channel

rx = apply_channel(tx,hTrue);

%% 6. Extract second LTS

[~,lts2_rx] = extract_lts(rx,p);

%% 7. FFT received LTS

Y_LTS = lts_fft(lts2_rx);

%% 8. Extract 100 LTS tones

[Y_p,X_p,ltsLocations,ltsBins] = ...
    extract_lts_100(Y_LTS,X_LTS);

%% 9. LS estimate

hLS_100 = pilot_ls_estimate(Y_p,X_p);

%% 10. Convert to CIR

[cir,windowedLS,window,Gmatrix] = ...
    ls100_to_cir128(hLS_100);

%% ====================================================
% Diagnostic information
% =====================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' TRUE CHANNEL\n');
fprintf('====================================================\n');

fprintf('Channel length = %d\n',length(hTrue));

for n = 1:length(hTrue)

    fprintf('h(%d) = %.6f %+.6fj   |h| = %.6f\n', ...
        n,hTrue(n),abs(hTrue(n)));

end

%% CIR energy

cirEnergy = abs(cir).^2;

[sortedEnergy,sortedIndex] = ...
    sort(cirEnergy,'descend');

fprintf('\n');
fprintf('====================================================\n');
fprintf(' STRONGEST ESTIMATED CIR SAMPLES\n');
fprintf('====================================================\n');

fprintf(' Rank     CIR index       Energy       Magnitude\n');
fprintf('----------------------------------------------------\n');

for n = 1:15

    idx = sortedIndex(n);

    fprintf('%4d %12d %14.6e %14.6e\n', ...
        n,idx,sortedEnergy(n),abs(cir(idx)));

end

%% First 20 CIR samples

fprintf('\n');
fprintf('====================================================\n');
fprintf(' FIRST 20 CIR SAMPLES\n');
fprintf('====================================================\n');

fprintf(' Index          Real             Imag          Magnitude\n');
fprintf('--------------------------------------------------------\n');

for n = 1:20

    fprintf('%5d %15.6f %15.6f %15.6f\n', ...
        n,real(cir(n)),imag(cir(n)),abs(cir(n)));

end

%% LS frequency-domain error

fprintf('\n');
fprintf('====================================================\n');
fprintf(' LS ESTIMATE MAGNITUDES\n');
fprintf('====================================================\n');

fprintf('Min |hLS| : %.6e\n',min(abs(hLS_100)));
fprintf('Max |hLS| : %.6e\n',max(abs(hLS_100)));
fprintf('Mean |hLS|: %.6e\n',mean(abs(hLS_100)));

%% Basic sanity

assert(length(hTrue) == 5);
assert(length(hLS_100) == 100);
assert(length(cir) == 128);

fprintf('\n');
fprintf('====================================================\n');
fprintf(' STRUCTURAL CHECKS PASSED\n');
fprintf('====================================================\n');

%% Plot CIR magnitude

figure;

stem(1:128,abs(cir),'filled');

grid on;

xlabel('CIR sample index');
ylabel('Magnitude');

title('Estimated CIR Magnitude — No AWGN');

xlim([1 128]);

%% Plot true channel

figure;

stem(1:length(hTrue),abs(hTrue),'filled');

grid on;

xlabel('Channel sample index');
ylabel('Magnitude');

title('True SUI-2 Channel');

xlim([1 10]);