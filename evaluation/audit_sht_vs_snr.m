clear
clc
close all

% =========================================================
% STEP 39
% SHT CIR-LENGTH ESTIMATION VS SNR
% =========================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' STEP 39 - SHT VS SNR\n');
fprintf('====================================================\n');

p = parameters();

% =========================================================
% Test settings
% =========================================================

SNR_list = 0:5:40;

numTrials = 20;

estimatedL = zeros(length(SNR_list),numTrials);


% =========================================================
% True physical channel extent
% =========================================================

h = sui2_channel(p);

pathDelays_samples = ...
    round(p.SUI2.pathDelay_us * 1e-6 / p.Ts);

trueL = max(pathDelays_samples) + 1;


fprintf('\n');
fprintf('True physical CIR length = %d samples\n',trueL);

fprintf('Number of SNR points     = %d\n', ...
    length(SNR_list));

fprintf('Trials per SNR           = %d\n', ...
    numTrials);


% =========================================================
% Main SNR loop
% =========================================================

for s = 1:length(SNR_list)

    currentSNR = SNR_list(s);

    fprintf('\n');
    fprintf('--------------------------------------------\n');
    fprintf('SNR = %d dB\n',currentSNR);
    fprintf('--------------------------------------------\n');

    for trial = 1:numTrials

        % -------------------------------------------------
        % Generate IEEE 802.16 preamble
        % -------------------------------------------------

        Pall = ieee80216_pall();

        preSeq = preamble_sequences(Pall);

        [stsTime,ltsTime] = ...
            preamble_to_ifft(preSeq,p);

        pre = build_80216_preamble( ...
            stsTime,ltsTime,p);


        % -------------------------------------------------
        % Pass through SUI-2
        % -------------------------------------------------

        rxClean = apply_channel( ...
            pre.complete,h);


        % -------------------------------------------------
        % Add AWGN
        % -------------------------------------------------

        rx = add_awgn( ...
            rxClean,currentSNR);


        % -------------------------------------------------
        % Extract LTS
        % -------------------------------------------------

        [~,lts2] = ...
            extract_lts(rx,p);


        % -------------------------------------------------
        % FFT
        % -------------------------------------------------

        Y_LTS = ...
            lts_fft(lts2);


        % -------------------------------------------------
        % Known LTS
        % -------------------------------------------------

        knownLTS = ...
            get_known_lts(ltsTime);

        X_LTS = ...
            known_lts_fft(knownLTS);


        % -------------------------------------------------
        % Extract 100 training tones
        % -------------------------------------------------

        [Y_p,X_p,~,~] = ...
            extract_lts_100(Y_LTS,X_LTS);


        % -------------------------------------------------
        % LS channel estimate
        % -------------------------------------------------

        H_LS = ...
            pilot_ls_estimate(Y_p,X_p);


        % -------------------------------------------------
        % LS -> CIR
        % -------------------------------------------------

        [cir,~,~] = ...
            ls_to_cir(H_LS,128,'hann');


        % -------------------------------------------------
        % Noise variance
        % -------------------------------------------------

        noiseVariance = ...
            estimate_noise_variance( ...
                cir, ...
                p.noiseTapStart, ...
                p.noiseTapEnd);


        % -------------------------------------------------
        % SHT
        % -------------------------------------------------

        [L_hat,~,~] = ...
            sht_cir_length( ...
                cir, ...
                noiseVariance, ...
                p.T, ...
                p.G);


        estimatedL(s,trial) = L_hat;

    end

    fprintf('Mean L_hat = %.3f\n', ...
        mean(estimatedL(s,:)));

    fprintf('Min L_hat  = %d\n', ...
        min(estimatedL(s,:)));

    fprintf('Max L_hat  = %d\n', ...
        max(estimatedL(s,:)));

    fprintf('Std L_hat  = %.3f\n', ...
        std(estimatedL(s,:)));

end


% =========================================================
% Summary
% =========================================================

meanL = mean(estimatedL,2);

minL = min(estimatedL,[],2);

maxL = max(estimatedL,[],2);

stdL = std(estimatedL,0,2);


fprintf('\n');
fprintf('====================================================\n');
fprintf(' STEP 39 SUMMARY\n');
fprintf('====================================================\n');

fprintf('\n');
fprintf(' SNR    Mean L    Min L    Max L    Std L\n');
fprintf('--------------------------------------------\n');

for s = 1:length(SNR_list)

    fprintf('%3d     %.3f      %3d      %3d     %.3f\n', ...
        SNR_list(s), ...
        meanL(s), ...
        minL(s), ...
        maxL(s), ...
        stdL(s));

end


% =========================================================
% Detection accuracy
%
% We regard L_hat >= trueL as successful detection
% of the physical channel extent.
% =========================================================

accuracy = zeros(length(SNR_list),1);

for s = 1:length(SNR_list)

    successful = ...
        estimatedL(s,:) >= trueL;

    accuracy(s) = ...
        mean(successful)*100;

end


fprintf('\n');
fprintf('SHT DETECTION SUCCESS RATE\n');
fprintf('--------------------------------------------\n');

for s = 1:length(SNR_list)

    fprintf('SNR = %2d dB : %.1f %%\n', ...
        SNR_list(s), ...
        accuracy(s));

end


% =========================================================
% Plot 1 - Mean estimated CIR length
% =========================================================

figure;

plot(SNR_list,meanL,'o-','LineWidth',1.2);

hold on;

plot(SNR_list, ...
    trueL*ones(size(SNR_list)), ...
    '--','LineWidth',1.2);

xlabel('SNR (dB)');
ylabel('Estimated CIR length');

title('SHT Estimated CIR Length vs SNR');

legend( ...
    'Mean estimated L', ...
    'True physical length');

grid on;


% =========================================================
% Plot 2 - Standard deviation
% =========================================================

figure;

plot(SNR_list,stdL,'o-','LineWidth',1.2);

xlabel('SNR (dB)');
ylabel('Standard deviation of L hat');

title('SHT CIR-Length Variability vs SNR');

grid on;


% =========================================================
% Plot 3 - Detection success rate
% =========================================================

figure;

plot(SNR_list,accuracy,'o-','LineWidth',1.2);

xlabel('SNR (dB)');
ylabel('Detection success (%)');

title('SHT Channel-Extent Detection Success vs SNR');

grid on;

ylim([0 105]);


% =========================================================
% Final
% =========================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' STEP 39 COMPLETE\n');
fprintf('====================================================\n');