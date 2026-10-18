clear
clc
close all

% =========================================================
% STEP 40
% INDEPENDENT SHT VALIDATION USING KNOWN CIR
% =========================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' STEP 40 - INDEPENDENT SHT VALIDATION\n');
fprintf('====================================================\n');

p = parameters();

% =========================================================
% Test settings
% =========================================================

SNR_list = 0:5:40;

numTrials = 100;

trueLengths = [3 5 8];

results = zeros(length(trueLengths), ...
                length(SNR_list), ...
                numTrials);


% =========================================================
% Test each known CIR length
% =========================================================

for Lidx = 1:length(trueLengths)

    trueL = trueLengths(Lidx);

    fprintf('\n');
    fprintf('============================================\n');
    fprintf('TRUE CIR LENGTH = %d\n',trueL);
    fprintf('============================================\n');

    % -----------------------------------------------------
    % Create known deterministic CIR
    % -----------------------------------------------------

    cirClean = zeros(1,128);

    for n = 1:trueL
        cirClean(n) = 1/sqrt(n);
    end

    % Normalize total CIR energy
    cirClean = cirClean / norm(cirClean);


    % -----------------------------------------------------
    % SNR loop
    % -----------------------------------------------------

    for s = 1:length(SNR_list)

        SNRdB = SNR_list(s);

        Ltemp = zeros(1,numTrials);

        for trial = 1:numTrials

            % -------------------------------------------------
            % Generate complex AWGN
            % -------------------------------------------------

            signalPower = mean(abs(cirClean).^2);

            noisePower = ...
                signalPower / 10^(SNRdB/10);

            noise = sqrt(noisePower/2) * ...
                (randn(1,128) + ...
                 1j*randn(1,128));

            cirNoisy = cirClean + noise;


            % -------------------------------------------------
            % Estimate noise variance
            %
            % Use the same region as the project:
            % samples 65 to 120
            % -------------------------------------------------

            noiseVariance = ...
                estimate_noise_variance( ...
                    cirNoisy, ...
                    p.noiseTapStart, ...
                    p.noiseTapEnd);


            % -------------------------------------------------
            % Run SHT
            % -------------------------------------------------

            [L_hat,~,~] = ...
                sht_cir_length( ...
                    cirNoisy, ...
                    noiseVariance, ...
                    p.T, ...
                    p.G);


            Ltemp(trial) = L_hat;

        end

        results(Lidx,s,:) = Ltemp;


        % -------------------------------------------------
        % Display
        % -------------------------------------------------

        fprintf(['SNR = %2d dB : Mean L = %.2f, ', ...
                 'Min = %d, Max = %d\n'], ...
            SNRdB, ...
            mean(Ltemp), ...
            min(Ltemp), ...
            max(Ltemp));

    end

end


% =========================================================
% Summary
% =========================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' STEP 40 SUMMARY\n');
fprintf('====================================================\n');

for Lidx = 1:length(trueLengths)

    trueL = trueLengths(Lidx);

    fprintf('\n');
    fprintf('True L = %d\n',trueL);

    fprintf('SNR    Mean L    Std L\n');
    fprintf('-------------------------\n');

    for s = 1:length(SNR_list)

        Ltemp = squeeze(results(Lidx,s,:));

        fprintf('%3d    %.3f     %.3f\n', ...
            SNR_list(s), ...
            mean(Ltemp), ...
            std(Ltemp));

    end

end


% =========================================================
% Detection accuracy
%
% Exact detection:
%
% L_hat == trueL
% =========================================================

accuracy = zeros(length(trueLengths), ...
                 length(SNR_list));

for Lidx = 1:length(trueLengths)

    trueL = trueLengths(Lidx);

    for s = 1:length(SNR_list)

        Ltemp = squeeze(results(Lidx,s,:));

        accuracy(Lidx,s) = ...
            mean(Ltemp == trueL) * 100;

    end

end


% =========================================================
% Accuracy display
% =========================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' EXACT CIR-LENGTH DETECTION ACCURACY\n');
fprintf('====================================================\n');

for Lidx = 1:length(trueLengths)

    trueL = trueLengths(Lidx);

    fprintf('\nTrue L = %d\n',trueL);

    for s = 1:length(SNR_list)

        fprintf('SNR = %2d dB : %.1f %%\n', ...
            SNR_list(s), ...
            accuracy(Lidx,s));

    end

end


% =========================================================
% Plot accuracy
% =========================================================

figure;

for Lidx = 1:length(trueLengths)

    plot(SNR_list, ...
         accuracy(Lidx,:), ...
         'o-', ...
         'LineWidth',1.2);

    hold on;

end

xlabel('SNR (dB)');
ylabel('Exact CIR-length detection (%)');

title('Independent SHT CIR-Length Detection');

legend( ...
    'True L = 3', ...
    'True L = 5', ...
    'True L = 8');

grid on;

ylim([0 105]);


% =========================================================
% Final
% =========================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' STEP 40 COMPLETE\n');
fprintf('====================================================\n');