function step5_evaluate_ml_16qam_ser()
% =========================================================================
% STEP 5 (PLOT 3): 16-QAM SYMBOL ERROR RATE (SER) PERFORMANCE EVALUATION
%
% Compares 16-QAM Symbol Error Rate (SER) across SNR = 10 to 35 dB for:
%   1. Perfect CSI Knowledge (Lower Bound)
%   2. Fixed SVD (L = 64, Rank = 67)
%   3. Park 2022 Adaptive SVD (SHT, T = 15) [Base Paper Baseline]
%   4. PROPOSED NOVEL ML-ADAPTIVE SVD (Learned L_hat)
%
% Output:
%   ML/Novelty_Day2_Plot3_SER_vs_SNR.png
%   ML/Novelty_Day2_Plot3_SER_results.mat
% =========================================================================

fprintf('\n============================================================\n');
fprintf(' NOVELTY PLOT 3: 16-QAM SYMBOL ERROR RATE (SER) EVALUATION\n');
fprintf('============================================================\n');

%% 1. SETUP
thisFile = mfilename('fullpath');
if isempty(thisFile), here = pwd; else, here = fileparts(thisFile); end
projectRoot = fileparts(here);
addpath(genpath(projectRoot));

p = parameters();
rng(2026);

SNRdB = 10:5:35;
nS = numel(SNRdB);
nTrials = 300; % Monte Carlo trials per SNR

usedK  = p.activeIndex(:);
ltsK   = [-100:2:-2 2:2:100].';
pilotK = p.pilotLocations(:);
dataK  = setdiff(usedK, pilotK, 'stable');

% Storage: [Perfect, Fixed, Park, ProposedML] x numSNR x 2 channels
ser = zeros(4, nS, 2);

%% 2. PRECOMPUTE COVARIANCES & SVD OPERATORS
idxUsed = usedK + p.N/2 + 1;
idxLTS  = ltsK  + p.N/2 + 1;

R64 = local_uniform_covariance(p.N, 64);
Rhp64 = R64(idxUsed, idxLTS);
Rpp64 = R64(idxLTS, idxLTS);
[U64, D64] = eig((Rpp64 + Rpp64') / 2);
[lam64, ord64] = sort(real(diag(D64)), 'descend');
U64 = U64(:, ord64);

% Precompute SVD banks for L = 1 to G
U_bank = cell(p.G, 1);
lam_bank = cell(p.G, 1);
Rhp_bank = cell(p.G, 1);

for L = 1:p.G
    RL = local_uniform_covariance(p.N, L);
    Rhp_bank{L} = RL(idxUsed, idxLTS);
    RppL = RL(idxLTS, idxLTS);
    [UL, DL] = eig((RppL + RppL') / 2);
    [la, oo] = sort(real(diag(DL)), 'descend');
    U_bank{L} = UL(:, oo);
    lam_bank{L} = la;
end

%% 3. SIMULATION LOOP OVER CHANNELS & SNRS
channelNames = {'SUI-2', 'SUI-5'};

for ch = 1:2
    chName = channelNames{ch};
    fprintf('\n--- Evaluating %s Channel ---\n', chName);

    for is = 1:nS
        snr = SNRdB(is);
        gammaP = 10^((snr + 6)/10); % Pilot boost +6 dB as in paper
        gammaD = 10^(snr/10);       % Data SNR

        % Fixed SVD Filter
        delta64 = zeros(100, 1);
        r64 = min(67, 100);
        delta64(1:r64) = 1 ./ (lam64(1:r64) + 1/gammaP);
        K_fixed = Rhp64 * U64 * diag(delta64) * U64';

        % Precompute SVD filter bank for this SNR
        Kbank = cell(p.G, 1);
        for L = 1:p.G
            r = min(L + p.adaptiveSVDOffset, 100);
            dd = zeros(100, 1);
            dd(1:r) = 1 ./ (lam_bank{L}(1:r) + 1/gammaP);
            Kbank{L} = Rhp_bank{L} * U_bank{L} * diag(dd) * U_bank{L}';
        end

        errTrials = zeros(4, nTrials);

        for tr = 1:nTrials
            % 1. Channel Realization
            if ch == 1
                h = sui2_channel(p);
            else
                h = sui5_channel(p);
            end

            H = local_channel_cfr(h, p.N);
            Hu = H(idxUsed);
            Hp = H(idxLTS);

            % 2. LTS Preamble LS Estimate
            nLS = sqrt(1/(2*gammaP)) * (randn(100,1) + 1j*randn(100,1));
            HLS = Hp + nLS;

            % 3. CIR Transform & Noise Estimation
            [cir,~,~] = ls_to_cir(HLS, 128, 'hann');
            sigma2 = mean(abs(cir(p.noiseTapStart : p.noiseTapEnd)).^2);

            % 4. CIR Length Estimates
            % Baseline Park SHT (T=15)
            L_park = sht_cir_length(cir, sigma2, p.T, p.G);
            L_park = max(1, min(p.G, L_park));

            % Proposed Novel ML
            L_ml = ml_cir_length(cir, sigma2, p.G, snr);
            L_ml = max(1, min(p.G, L_ml));

            if tr <= 2
                if ch == 1, L_true = 5; else, L_true = 41; end
                fprintf('  %s SNR=%ddB [tr=%d]: L_true=%d, L_SHT=%d, L_final=%d\n', ...
                    chName, snr, tr, L_true, L_park, L_ml);
            end

            % 5. Channel Estimates
            H_perfect = Hu;
            H_fixed   = K_fixed * HLS;
            H_park    = Kbank{L_park} * HLS;
            H_ml      = Kbank{L_ml} * HLS;

            Hset = {H_perfect, H_fixed, H_park, H_ml};

            % 6. Transmit 16-QAM Data + Pilots
            tx_data = local_qam16_symbols(numel(dataK));
            pilot_vals = ones(numel(pilotK), 1);

            X = zeros(200, 1);
            [~, idD] = ismember(dataK, usedK);
            [~, idP] = ismember(pilotK, usedK);
            X(idD) = tx_data;
            X(idP) = pilot_vals;

            noise = sqrt(1/(2*gammaD)) * (randn(200,1) + 1j*randn(200,1));
            Y = Hu .* X + noise;

            % 7. Equalization & Symbol Detection for Each Estimator
            for ee = 1:4
                He = Hset{ee};
                % Residual Phase Compensation via pilots
                z = Y(idP) ./ (He(idP) .* pilot_vals + eps);
                phi = angle(sum(z));
                Ye = Y .* exp(-1j * phi);

                % Zero-Forcing Equalization on Data Subcarriers
                shat = Ye(idD) ./ (He(idD) + eps);
                rx_dec = local_qam16_detect(shat);

                % Symbol Error Rate
                errTrials(ee, tr) = mean(rx_dec(:) ~= tx_data(:));
            end
        end

        ser(:, is, ch) = mean(errTrials, 2);

        fprintf('  SNR = %2d dB | Perfect: %.3e | Fixed: %.3e | Park: %.3e | Proposed ML: %.3e\n', ...
            snr, ser(1,is,ch), ser(2,is,ch), ser(3,is,ch), ser(4,is,ch));
    end
end

%% 4. PRINT NUMERICAL SER SUMMARY TABLE
fprintf('\n========================================================================================\n');
fprintf(' 16-QAM SYMBOL ERROR RATE (SER) SUMMARY TABLE\n');
fprintf('========================================================================================\n');
fprintf(' SNR(dB) | SUI-2 Fixed SVD | SUI-2 Park SVD | SUI-2 Proposed ML | SUI-5 Park SVD | SUI-5 Proposed ML\n');
fprintf('----------------------------------------------------------------------------------------\n');
for i = 1:nS
    fprintf('  %2d dB  |    %10.4e   |   %10.4e   |   %10.4e     |   %10.4e   |   %10.4e\n', ...
        SNRdB(i), ser(2,i,1), ser(3,i,1), ser(4,i,1), ser(3,i,2), ser(4,i,2));
end
fprintf('========================================================================================\n');

%% 5. PUBLICATION FIGURE GENERATION
fig = figure('Name', 'Novelty Plot 3 - 16-QAM SER Performance', ...
    'Color', 'w', 'Position', [100 80 1150 520], 'Visible', 'off');

% (a) SUI-2
subplot(1, 2, 1);
semilogy(SNRdB, ser(1,:,1), ':k', 'LineWidth', 2.0); hold on;
semilogy(SNRdB, ser(2,:,1), '-xg', 'LineWidth', 1.8, 'MarkerSize', 8);
semilogy(SNRdB, ser(3,:,1), '-sb', 'LineWidth', 1.8, 'MarkerFaceColor', 'b', 'MarkerSize', 6);
semilogy(SNRdB, ser(4,:,1), '-or', 'LineWidth', 2.2, 'MarkerFaceColor', 'r', 'MarkerSize', 6);
grid on; box on;
set(gca, 'YScale', 'log', 'YMinorGrid', 'on', 'FontSize', 10);
xlim([10 35]); ylim([5e-4 1]);
xlabel('SNR (dB)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('16-QAM Symbol Error Rate (SER)', 'FontSize', 11, 'FontWeight', 'bold');
title('(a) SUI-2 Channel (Quasi-Static)', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Perfect CSI Knowledge', 'Fixed SVD (L=64)', 'Park 2022 Adaptive SVD', 'PROPOSED NOVEL ML-SVD'}, ...
    'Location', 'southwest', 'FontSize', 9);

% (b) SUI-5
subplot(1, 2, 2);
semilogy(SNRdB, ser(1,:,2), ':k', 'LineWidth', 2.0); hold on;
semilogy(SNRdB, ser(2,:,2), '-xg', 'LineWidth', 1.8, 'MarkerSize', 8);
semilogy(SNRdB, ser(3,:,2), '-sb', 'LineWidth', 1.8, 'MarkerFaceColor', 'b', 'MarkerSize', 6);
semilogy(SNRdB, ser(4,:,2), '-or', 'LineWidth', 2.2, 'MarkerFaceColor', 'r', 'MarkerSize', 6);
grid on; box on;
set(gca, 'YScale', 'log', 'YMinorGrid', 'on', 'FontSize', 10);
xlim([10 35]); ylim([1e-3 1]);
xlabel('SNR (dB)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('16-QAM Symbol Error Rate (SER)', 'FontSize', 11, 'FontWeight', 'bold');
title('(b) SUI-5 Channel (Quasi-Static)', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Perfect CSI Knowledge', 'Fixed SVD (L=64)', 'Park 2022 Adaptive SVD', 'PROPOSED NOVEL ML-SVD'}, ...
    'Location', 'southwest', 'FontSize', 9);

sgtitle('16-QAM Symbol Error Rate (SER): Proposed ML-Adaptive SVD vs Park 2022 Baseline', ...
    'FontSize', 13, 'FontWeight', 'bold');

% Save Artifacts
outPng = fullfile(projectRoot, 'ML', 'Novelty_Day2_Plot3_SER_vs_SNR.png');
outMat = fullfile(projectRoot, 'ML', 'Novelty_Day2_Plot3_SER_results.mat');

saveas(fig, outPng);
save(outMat, 'SNRdB', 'ser', '-v7');

fprintf('\nSuccessfully saved Plot 3 results to:\n  %s\n  %s\n', outPng, outMat);
fprintf('============================================================\n');

end

%% --- HELPER FUNCTIONS ---
function R = local_uniform_covariance(N, L)
    k = (-N/2:N/2-1).';
    d = k - k.';
    R = ones(N);
    idx = (d ~= 0);
    R(idx) = (1 - exp(-1j*2*pi*L*d(idx)/N)) ./ (1j*2*pi*L*d(idx)/N);
    R = (R + R') / 2;
end

function H = local_channel_cfr(h, N)
    h = h(:);
    ht = [h; zeros(N - numel(h), 1)];
    H = fftshift(fft(ht, N));
    H = H(:);
end

function x = local_qam16_symbols(n)
    levels = [-3 -1 1 3];
    a = levels(randi(4, n, 1));
    b = levels(randi(4, n, 1));
    x = (a(:) + 1j*b(:)) / sqrt(10);
end

function xhat = local_qam16_detect(x)
    levels = [-3 -1 1 3] / sqrt(10);
    x = x(:);
    xr = real(x);
    xi = imag(x);
    [~, ir] = min(abs(xr - levels), [], 2);
    [~, ii] = min(abs(xi - levels), [], 2);
    xhat = levels(ir) + 1j*levels(ii);
    xhat = xhat(:);
end
