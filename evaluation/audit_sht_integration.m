clear
clc
close all

% =========================================================
% STEP 38
% SHT CIR-LENGTH ESTIMATION INTEGRATION
% =========================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' STEP 38 - SHT CIR LENGTH ESTIMATION\n');
fprintf('====================================================\n');

p = parameters();


% =========================================================
% 1. Generate IEEE 802.16 preamble
% =========================================================

Pall = ieee80216_pall();

preSeq = preamble_sequences(Pall);

[stsTime,ltsTime] = ...
    preamble_to_ifft(preSeq,p);

pre = build_80216_preamble( ...
    stsTime,ltsTime,p);


% =========================================================
% 2. Generate SUI-2 channel
% =========================================================

h = sui2_channel(p);


% =========================================================
% 3. Pass through channel
% =========================================================

rxClean = apply_channel( ...
    pre.complete,h);


% =========================================================
% 4. Add AWGN
%
% Paper CIR-length simulation:
% SNR = 20 dB
% =========================================================

SNRdB = 20;

rx = add_awgn( ...
    rxClean,SNRdB);


% =========================================================
% 5. Extract LTS
% =========================================================

[lts1,lts2] = ...
    extract_lts(rx,p);


% =========================================================
% 6. Received LTS FFT
% =========================================================

Y_LTS = lts_fft(lts2);


% =========================================================
% 7. Known LTS
% =========================================================

knownLTS = ...
    get_known_lts(ltsTime);

X_LTS = ...
    known_lts_fft(knownLTS);


% =========================================================
% 8. Extract 100 LTS training tones
% =========================================================

[Y_p,X_p,ltsLocations,ltsBins] = ...
    extract_lts_100(Y_LTS,X_LTS);


% =========================================================
% 9. LS channel estimate
% =========================================================

H_LS = ...
    pilot_ls_estimate(Y_p,X_p);


% =========================================================
% 10. Convert LS estimate to CIR
%
% Paper:
% M = 128
% Hann window
% =========================================================

[cir,windowedLS,window] = ...
    ls_to_cir(H_LS,128,'hann');


% =========================================================
% 11. Estimate noise variance
%
% G = 64
% K = 56
%
% samples:
% 65 ... 120
% =========================================================

noiseVariance = ...
    estimate_noise_variance( ...
        cir, ...
        p.noiseTapStart, ...
        p.noiseTapEnd);


% =========================================================
% 12. SHT
%
% threshold:
%
% eta0 = T * sigma_hat^2
% =========================================================

T = p.T;

[L_hat,threshold,testEnergy] = ...
    sht_cir_length( ...
        cir, ...
        noiseVariance, ...
        T, ...
        p.G);


% =========================================================
% 13. Display physical channel
% =========================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' PHYSICAL SUI-2 CHANNEL\n');
fprintf('====================================================\n');

for n = 1:length(h)

    if abs(h(n)) > 0

        fprintf('h(%d) = %.8f\n', ...
            n,h(n));

    end

end


% =========================================================
% 14. Expected CIR length
% =========================================================

pathDelays_samples = ...
    round(p.SUI2.pathDelay_us * ...
    1e-6 / p.Ts);

trueLastPath = ...
    max(pathDelays_samples) + 1;


fprintf('\n');
fprintf('Expected last physical path sample = %d\n', ...
    trueLastPath);


% =========================================================
% 15. Display SHT parameters
% =========================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' SHT PARAMETERS\n');
fprintf('====================================================\n');

fprintf('SNR                  = %d dB\n',SNRdB);

fprintf('G                    = %d\n',p.G);

fprintf('K                    = %d\n',p.K);

fprintf('Noise samples        = %d to %d\n', ...
    p.noiseTapStart, ...
    p.noiseTapEnd);

fprintf('T                    = %.2f\n',T);

fprintf('Estimated variance   = %.10e\n', ...
    noiseVariance);

fprintf('SHT threshold        = %.10e\n', ...
    threshold);


% =========================================================
% 16. Display CIR energy at each tap
% =========================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' CIR TAP ENERGY\n');
fprintf('====================================================\n');

fprintf('Tap       |g[n]|        |g[n]|^2       TEST\n');
fprintf('----------------------------------------------------\n');

for n = 1:p.G

    if testEnergy(n) > threshold

        result = 'PASS';

    else

        result = 'FAIL';

    end

    fprintf('%3d       %.8f    %.8e    %s\n', ...
        n, ...
        abs(cir(n)), ...
        testEnergy(n), ...
        result);

end


% =========================================================
% 17. Final SHT result
% =========================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' SHT RESULT\n');
fprintf('====================================================\n');

fprintf('Estimated CIR length = %d samples\n', ...
    L_hat);

fprintf('Expected last path   = %d samples\n', ...
    trueLastPath);


if L_hat >= trueLastPath

    fprintf(['PASS: SHT detected the physical channel ', ...
             'extent or beyond.\n']);

else

    fprintf(['INFO: SHT underestimated the physical ', ...
             'channel extent.\n']);

end


% =========================================================
% 18. Theoretical false alarm probability
% =========================================================

Pf = ...
    (p.K/(T+p.K))^p.K;

fprintf('\n');
fprintf('Theoretical false alarm probability:\n');

fprintf('Pf = %.10e\n',Pf);


% =========================================================
% 19. Plot CIR and threshold
% =========================================================

figure;

plot(1:p.G,abs(cir(1:p.G)).^2,'LineWidth',1.2);

hold on;

plot(1:p.G, ...
    threshold*ones(1,p.G), ...
    '--','LineWidth',1.2);

xlabel('CIR sample');
ylabel('|g[n]|^2');

title('SHT CIR-Length Detection');

legend( ...
    'CIR tap energy', ...
    'SHT threshold');

grid on;

xlim([1 p.G]);


% =========================================================
% 20. Mark detected CIR length
% =========================================================

xline(L_hat, ...
    '--', ...
    sprintf('L hat = %d',L_hat));

xline(trueLastPath, ...
    ':', ...
    sprintf('True last path = %d',trueLastPath));


% =========================================================
% END
% =========================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' STEP 38 COMPLETE\n');
fprintf('====================================================\n');