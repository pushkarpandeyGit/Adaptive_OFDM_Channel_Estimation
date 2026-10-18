function audit_ls_frequency_response()

clc;

fprintf('\n');
fprintf('====================================================\n');
fprintf(' STEP 33: 100-TONE LS FREQUENCY RESPONSE AUDIT\n');
fprintf('====================================================\n');

%% 1. Load parameters
p = parameters();

fprintf('\nProject parameters loaded.\n');
fprintf('N  = %d\n', p.N);
fprintf('Ts = %.2f us\n', p.Ts*1e6);

%% 2. Build IEEE 802.16 preamble
Pall = ieee80216_pall();

preamble = preamble_sequences(Pall);

[stsTime, ltsTime] = preamble_to_ifft(preamble, p);

txPreamble = build_80216_preamble(stsTime, ltsTime, p);

%% 3. Generate SUI-2 channel
h = sui2_channel(p);

fprintf('\nTrue SUI-2 channel:\n');

for n = 1:length(h)

    fprintf('h(%d) = %.6f + j%.6f\n', ...
        n, real(h(n)), imag(h(n)));

end

%% 4. Pass preamble through channel
rx = apply_channel(txPreamble.complete, h);

%% 5. Extract second LTS
[~, lts2] = extract_lts(rx, p);

%% 6. FFT received LTS
Y_LTS = lts_fft(lts2);

%% 7. FFT known LTS
knownLTS = get_known_lts(ltsTime);

X_LTS = known_lts_fft(knownLTS);

%% 8. Extract the 100 LTS tones
%
% IMPORTANT:
% extract_lts_100 already knows the correct
% 100-tone LTS locations and 128-point mapping.
%
[Y_p, X_p, ltsLocations, ltsBins] = ...
    extract_lts_100(Y_LTS, X_LTS);

%% 9. Calculate 100-tone LS estimate
hLS_p = pilot_ls_estimate(Y_p, X_p);

%% 10. Calculate theoretical SUI-2 frequency response
H_true = zeros(size(ltsLocations));

for k = 1:length(ltsLocations)

    subcarrier = ltsLocations(k);

    H = 0;

    for n = 0:length(h)-1

        H = H + h(n+1) * ...
            exp(-1j*2*pi*subcarrier*n/p.N);

    end

    H_true(k) = H;

end

%% 11. Compare LS estimate with true response
errorVector = hLS_p(:).' - H_true;

MSE = mean(abs(errorVector).^2);

maxError = max(abs(errorVector));

relativeError = ...
    norm(errorVector) / norm(H_true);

fprintf('\n');
fprintf('====================================================\n');
fprintf(' FREQUENCY RESPONSE COMPARISON\n');
fprintf('====================================================\n');

fprintf('Number of LTS tones : %d\n', length(ltsLocations));

fprintf('LS MSE              : %.12e\n', MSE);

fprintf('Maximum absolute    : %.12e\n', maxError);

fprintf('Relative error      : %.12e\n', relativeError);

%% 12. Show selected tones
fprintf('\n');
fprintf('Sample tone comparison:\n');

fprintf('\n');
fprintf(' k       LS magnitude       TRUE magnitude       Error\n');
fprintf('---------------------------------------------------------\n');

sampleIndices = [1 2 3 10 25 50 75 99 100];

for i = sampleIndices

    fprintf('%4d    %14.8f    %14.8f    %.3e\n', ...
        ltsLocations(i), ...
        abs(hLS_p(i)), ...
        abs(H_true(i)), ...
        abs(errorVector(i)));

end

%% 13. Plot magnitude comparison
figure;

plot(ltsLocations, abs(H_true), 'o-');

hold on;

plot(ltsLocations, abs(hLS_p), 'x--');

grid on;

xlabel('Original LTS Subcarrier Index');

ylabel('Channel Magnitude');

title('True vs 100-Tone LS Channel Frequency Response');

legend('True SUI-2 Response', ...
       '100-Tone LS Estimate');

%% 14. Plot error
figure;

stem(ltsLocations, abs(errorVector), 'filled');

grid on;

xlabel('Original LTS Subcarrier Index');

ylabel('|LS - True|');

title('100-Tone LS Frequency Response Error');

%% 15. Final decision
fprintf('\n');
fprintf('====================================================\n');
fprintf(' STEP 33 DECISION\n');
fprintf('====================================================\n');

if relativeError < 1e-10

    fprintf('PASS: 100-tone LS matches the theoretical channel.\n');
    fprintf('LTS extraction and LS estimation are correct.\n');
    fprintf('The remaining issue is LS-to-CIR processing.\n');

elseif relativeError < 1e-6

    fprintf('PASS: LS response matches theoretical response\n');
    fprintf('within very small numerical error.\n');
    fprintf('The remaining issue is most likely LS-to-CIR processing.\n');

else

    fprintf('FAIL: 100-tone LS does NOT match theoretical response.\n');
    fprintf('Do NOT modify the CIR code yet.\n');
    fprintf('We must investigate the LTS frequency-bin mapping.\n');

end

fprintf('====================================================\n');

end