function step9_evaluate_joint_system()
% =========================================================================
% STEP 9 (PLOT 7): JOINT SYSTEM EVALUATION (ML + ADAPTIVE WINDOWING)
%
% Evaluates the synergistic performance of the COMPLETE PROPOSED SYSTEM:
%   - Novelty 1: Data-Driven ML CIR Length Predictor
%   - Novelty 2: Channel-Aware Dynamic Windowing (Hamming/Hann/Blackman)
%
% Compares Channel Estimation MSE across SNR = 0 to 30 dB in SUI-2 and SUI-5:
%   1. Ideal LMMSE (Theoretical Lower Bound)
%   2. Fixed SVD (L = 64, Rank = 67)
%   3. Park 2022 Baseline (Fixed Hann + SHT T = 15)
%   4. Proposed Novelty 1 Alone (Fixed Hann + ML SVD)
%   5. PROPOSED COMPLETE DUAL-NOVELTY SYSTEM (Adaptive Window + ML SVD)
%
% Outputs:
%   novelty/Novelty_Day4_Plot7_Joint_System_MSE.png
%   novelty/Novelty_Day4_Plot7_Joint_System_results.mat
% =========================================================================

fprintf('\n============================================================\n');
fprintf(' NOVELTY PLOT 7: JOINT DUAL-NOVELTY SYSTEM (ML + ADAPTIVE WIN)\n');
fprintf('============================================================\n');

%% 1. SYSTEM SETUP
thisFile = mfilename('fullpath');
if isempty(thisFile), here = pwd; else, here = fileparts(thisFile); end
projectRoot = fileparts(here);
addpath(genpath(projectRoot));

p = parameters();
rng(2026);

N = p.N; % 256
G = p.G; % 64
M = 128;

SNRdB  = 0:5:30;
nS     = numel(SNRdB);
nTrials = 300; % Monte Carlo trials per SNR

usedK  = p.activeIndex(:);
ltsK   = [-100:2:-2 2:2:100].';
pilotK = p.pilotLocations(:);

Nu = numel(usedK); % 200
Np = numel(ltsK);  % 100

idxUsed = usedK + N/2 + 1;
idxLTS  = ltsK  + N/2 + 1;

% Storage: [Ideal, Fixed, Park, ML_Hann, Proposed_Joint] x numSNR x 2 channels
mse = zeros(5, nS, 2);

%% 2. PRECOMPUTE FIXED SVD & SVD FILTER BANKS
R64 = local_uniform_covariance(N, 64);
Rhp64 = R64(idxUsed, idxLTS);
Rpp64 = R64(idxLTS, idxLTS);
[U64, D64] = eig((Rpp64 + Rpp64') / 2);
[lam64, ord64] = sort(real(diag(D64)), 'descend');
U64 = U64(:, ord64);

U_bank   = cell(G, 1);
lam_bank = cell(G, 1);
Rhp_bank = cell(G, 1);

for L = 1:G
    RL = local_uniform_covariance(N, L);
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
    fprintf('\n--- Evaluating %s Channel (%d trials/SNR) ---\n', chName, nTrials);

    % True channel covariance for Ideal LMMSE
    if ch == 1
        % SUI-2 parameters
        delays = [0, 0.4, 1.1] * 1e-6 / p.Ts;
        powers = 10.^([0, -12, -15]/10); powers = powers / sum(powers);
    else
        % SUI-5 parameters
        delays = [0, 4.0, 10.0] * 1e-6 / p.Ts;
        powers = 10.^([0, -5, -10]/10); powers = powers / sum(powers);
    end

    d_used = usedK - usedK.';
    d_lts  = ltsK  - ltsK.';
    d_ul   = usedK - ltsK.';

    R_ideal_hp = zeros(Nu, Np);
    R_ideal_pp = zeros(Np, Np);
    for p_idx = 1:numel(delays)
        R_ideal_hp = R_ideal_hp + powers(p_idx) * exp(-1j * 2*pi * d_ul * delays(p_idx) / N);
        R_ideal_pp = R_ideal_pp + powers(p_idx) * exp(-1j * 2*pi * d_lts * delays(p_idx) / N);
    end

    for is = 1:nS
        snr = SNRdB(is);
        gammaP = 10^((snr + 6)/10); % Pilot boost +6 dB

        % Fixed SVD Filter
        delta64 = zeros(Np, 1);
        r64 = min(67, Np);
        delta64(1:r64) = 1 ./ (lam64(1:r64) + 1/gammaP);
        K_fixed = Rhp64 * U64 * diag(delta64) * U64';

        % SVD Filter Bank
        Kbank = cell(G, 1);
        for L = 1:G
            r = min(L + p.adaptiveSVDOffset, Np);
            dd = zeros(Np, 1);
            dd(1:r) = 1 ./ (lam_bank{L}(1:r) + 1/gammaP);
            Kbank{L} = Rhp_bank{L} * U_bank{L} * diag(dd) * U_bank{L}';
        end

        % Ideal LMMSE Filter
        K_ideal = R_ideal_hp / (R_ideal_pp + (1/gammaP) * eye(Np));

        errTrials = zeros(5, nTrials);

        for tr = 1:nTrials
            % 1. Channel Realization
            if ch == 1
                los = sqrt(10^(0.3)/(10^(0.3)+1));
                nlos = sqrt(1/(10^(0.3)+1)) * (randn + 1j*randn)/sqrt(2);
                a1 = sqrt(powers(1)) * (los + nlos);
                a2 = sqrt(powers(2)) * (randn + 1j*randn)/sqrt(2);
                a3 = sqrt(powers(3)) * (randn + 1j*randn)/sqrt(2);
                h_true = [a1; a2; a3];
            else
                b1 = sqrt(powers(1)) * (randn + 1j*randn)/sqrt(2);
                b2 = sqrt(powers(2)) * (randn + 1j*randn)/sqrt(2);
                b3 = sqrt(powers(3)) * (randn + 1j*randn)/sqrt(2);
                h_true = [b1; b2; b3];
            end

            H_full_u = zeros(Nu, 1);
            H_full_p = zeros(Np, 1);
            for p_idx = 1:numel(delays)
                H_full_u = H_full_u + h_true(p_idx) * exp(-1j * 2*pi * usedK * delays(p_idx) / N);
                H_full_p = H_full_p + h_true(p_idx) * exp(-1j * 2*pi * ltsK  * delays(p_idx) / N);
            end

            % 2. LTS Preamble LS Observation
            nLS = sqrt(1/(2*gammaP)) * (randn(Np,1) + 1j*randn(Np,1));
            HLS = H_full_p + nLS;

            % 3. Fixed Hann Window Transformation (Park 2022 & ML alone)
            [cir_hann,~,~] = ls_to_cir(HLS, M, 'hann');
            sigma2_hann = mean(abs(cir_hann(p.noiseTapStart : p.noiseTapEnd)).^2);

            % 4. Adaptive Window Transformation (Novelty 2)
            [winType_adapt, ~] = adaptive_window_selector(HLS, M, snr);
            [cir_adapt,~,~] = ls_to_cir(HLS, M, winType_adapt);
            sigma2_adapt = mean(abs(cir_adapt(p.noiseTapStart : p.noiseTapEnd)).^2);

            % 5. Length Estimations
            % Baseline: Park SHT (Fixed Hann)
            L_park = sht_cir_length(cir_hann, sigma2_hann, p.T, G);
            L_park = max(1, min(G, L_park));

            % Proposed Novelty 1 alone: ML (Fixed Hann)
            L_ml_hann = ml_cir_length(cir_hann, sigma2_hann, G, snr);
            L_ml_hann = max(1, min(G, L_ml_hann));

            % Proposed Complete System: ML with Adaptive Window CIR
            L_ml_joint = ml_cir_length(cir_adapt, sigma2_adapt, G, snr);
            L_ml_joint = max(1, min(G, L_ml_joint));

            if tr <= 2
                if ch == 1, L_true = 5; else, L_true = 41; end
                fprintf('  %s SNR=%ddB [tr=%d]: L_true=%d, L_SHT=%d, L_ML_hann=%d, L_ML_joint=%d, win=%s\n', ...
                    chName, snr, tr, L_true, L_park, L_ml_hann, L_ml_joint, winType_adapt);
            end

            % 6. Channel Estimates
            H_ideal    = K_ideal * HLS;
            H_fixed    = K_fixed * HLS;
            H_park     = Kbank{L_park} * HLS;
            H_ml_hann  = Kbank{L_ml_hann} * HLS;
            H_joint    = Kbank{L_ml_joint} * HLS;

            % 7. Squared Error
            errTrials(1, tr) = mean(abs(H_ideal   - H_full_u).^2);
            errTrials(2, tr) = mean(abs(H_fixed   - H_full_u).^2);
            errTrials(3, tr) = mean(abs(H_park    - H_full_u).^2);
            errTrials(4, tr) = mean(abs(H_ml_hann - H_full_u).^2);
            errTrials(5, tr) = mean(abs(H_joint   - H_full_u).^2);
        end

        mse(:, is, ch) = mean(errTrials, 2);

        fprintf('  SNR = %2d dB | Ideal: %.2e | Fixed: %.2e | Park: %.2e | ML-Hann: %.2e | Joint: %.2e\n', ...
            snr, mse(1,is,ch), mse(2,is,ch), mse(3,is,ch), mse(4,is,ch), mse(5,is,ch));
    end
end

%% 4. SUMMARY TABLE
fprintf('\n========================================================================================================\n');
fprintf(' JOINT SYSTEM (ML + ADAPTIVE WINDOWING) CHANNEL ESTIMATION MSE SUMMARY TABLE\n');
fprintf('========================================================================================================\n');
fprintf(' SNR(dB) | SUI-2 Fixed SVD | SUI-2 Park 2022 | SUI-2 JOINT NOVELTY | SUI-5 Fixed SVD | SUI-5 Park 2022 | SUI-5 JOINT NOVELTY\n');
fprintf('--------------------------------------------------------------------------------------------------------\n');
for i = 1:nS
    fprintf('  %2d dB  |    %10.3e   |   %10.3e    |     %10.3e     |    %10.3e   |   %10.3e    |     %10.3e\n', ...
        SNRdB(i), mse(2,i,1), mse(3,i,1), mse(5,i,1), mse(2,i,2), mse(3,i,2), mse(5,i,2));
end
fprintf('========================================================================================================\n');

%% 5. PUBLICATION FIGURE GENERATION
fig = figure('Name', 'Novelty Plot 7 - Joint System MSE Performance', ...
    'Color', 'w', 'Position', [100 80 1150 520], 'Visible', 'off');

% (a) SUI-2
subplot(1, 2, 1);
semilogy(SNRdB, mse(1,:,1), ':k', 'LineWidth', 2.0); hold on;
semilogy(SNRdB, mse(2,:,1), '-xg', 'LineWidth', 1.8, 'MarkerSize', 8);
semilogy(SNRdB, mse(3,:,1), '-sb', 'LineWidth', 1.8, 'MarkerFaceColor', 'b', 'MarkerSize', 6);
semilogy(SNRdB, mse(4,:,1), '--m', 'LineWidth', 1.8, 'MarkerFaceColor', 'm', 'MarkerSize', 5);
semilogy(SNRdB, mse(5,:,1), '-or', 'LineWidth', 2.2, 'MarkerFaceColor', 'r', 'MarkerSize', 7);
grid on; box on;
set(gca, 'YScale', 'log', 'YMinorGrid', 'on', 'FontSize', 10);
xlim([0 30]); ylim([1e-4 1]);
xlabel('Signal-to-Noise Ratio (dB)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Channel Estimation MSE', 'FontSize', 11, 'FontWeight', 'bold');
title('(a) SUI-2 Channel (Short Delay Spread)', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Ideal LMMSE Bound', 'Fixed SVD (L=64)', 'Park 2022 (Hann+SHT)', ...
        'Novelty 1 (Hann+ML)', 'PROPOSED JOINT (AdaptWin+ML)'}, ...
        'Location', 'southwest', 'FontSize', 9);

% (b) SUI-5
subplot(1, 2, 2);
semilogy(SNRdB, mse(1,:,2), ':k', 'LineWidth', 2.0); hold on;
semilogy(SNRdB, mse(2,:,2), '-xg', 'LineWidth', 1.8, 'MarkerSize', 8);
semilogy(SNRdB, mse(3,:,2), '-sb', 'LineWidth', 1.8, 'MarkerFaceColor', 'b', 'MarkerSize', 6);
semilogy(SNRdB, mse(4,:,2), '--m', 'LineWidth', 1.8, 'MarkerFaceColor', 'm', 'MarkerSize', 5);
semilogy(SNRdB, mse(5,:,2), '-or', 'LineWidth', 2.2, 'MarkerFaceColor', 'r', 'MarkerSize', 7);
grid on; box on;
set(gca, 'YScale', 'log', 'YMinorGrid', 'on', 'FontSize', 10);
xlim([0 30]); ylim([1e-4 1]);
xlabel('Signal-to-Noise Ratio (dB)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Channel Estimation MSE', 'FontSize', 11, 'FontWeight', 'bold');
title('(b) SUI-5 Channel (High Delay Spread)', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Ideal LMMSE Bound', 'Fixed SVD (L=64)', 'Park 2022 (Hann+SHT)', ...
        'Novelty 1 (Hann+ML)', 'PROPOSED JOINT (AdaptWin+ML)'}, ...
        'Location', 'southwest', 'FontSize', 9);

sgtitle('Channel Estimation MSE: Proposed Joint System (ML + Adaptive Windowing) vs Baselines', ...
    'FontSize', 13, 'FontWeight', 'bold');

% Save Artifacts
outPng = fullfile(projectRoot, 'novelty', 'Novelty_Day4_Plot7_Joint_System_MSE.png');
outMat = fullfile(projectRoot, 'novelty', 'Novelty_Day4_Plot7_Joint_System_results.mat');

saveas(fig, outPng);
save(outMat, 'SNRdB', 'mse', '-v7');

fprintf('\nSuccessfully saved Plot 7 results to:\n  %s\n  %s\n', outPng, outMat);
fprintf('============================================================\n');

end

%% --- HELPER FUNCTION ---
function R = local_uniform_covariance(N, L)
    k = (-N/2:N/2-1).';
    d = k - k.';
    R = ones(N);
    idx = (d ~= 0);
    R(idx) = (1 - exp(-1j*2*pi*L*d(idx)/N)) ./ (1j*2*pi*L*d(idx)/N);
    R = (R + R') / 2;
end
