function step6_evaluate_ml_time_varying()
% =========================================================================
% STEP 6 (PLOT 4): TIME-VARYING CHANNEL TRACKING PERFORMANCE
%
% Evaluates channel estimation and symbol detection under mobile time-varying
% conditions (ITU Vehicular-A channel, Doppler frequency fd = 30 Hz) with
% 2 Decision-Directed (DD) iterations across SNR = 10 to 40 dB.
%
% Compares:
%   1. Perfect CSI Knowledge (Lower Bound)
%   2. Fixed SVD (L = 64, Rank = 67)
%   3. Windowed DFT (L = 64) [Park 2022 Reference]
%   4. Park 2022 Adaptive SVD (SHT, T = 15) [Base Paper Baseline]
%   5. PROPOSED NOVEL ML-ADAPTIVE SVD (Learned L_hat)
%
% Outputs:
%   ML/Novelty_Day3_Plot4_TimeVarying_Tracking.png
%   ML/Novelty_Day3_Plot4_TimeVarying_results.mat
% =========================================================================

fprintf('\n============================================================\n');
fprintf(' NOVELTY PLOT 4: TIME-VARYING CHANNEL TRACKING (fd = 30 Hz)\n');
fprintf('============================================================\n');

%% 1. SYSTEM SETUP
thisFile = mfilename('fullpath');
if isempty(thisFile), here = pwd; else, here = fileparts(thisFile); end
projectRoot = fileparts(here);
addpath(genpath(projectRoot));

p = parameters();
rng(2026);

N       = p.N;                 % 256
G       = p.G;                 % 64
fd      = p.fd;                % 30 Hz

SNRdB   = 10:5:40;
numSNR  = numel(SNRdB);
NUM_PACKETS = 250;             % Monte Carlo packets per SNR
NUM_DD_ITER = 2;

usedK   = p.activeIndex(:);
ltsK    = [-100:2:-2 2:2:100].';
pilotK  = p.pilotLocations(:);
dataK   = setdiff(usedK, pilotK, 'stable');
dataK   = dataK(:);

Nu      = numel(usedK);        % 200
Np      = numel(ltsK);         % 100
Nd      = numel(dataK);        % 192

idxUsed = usedK + N/2 + 1;
idxLTS  = ltsK  + N/2 + 1;

[~, idD] = ismember(dataK, usedK);
[~, idP] = ismember(pilotK, usedK);

% Storage: [Park, Fixed, WDFT, ProposedML]
% MSE: 4 estimators x numSNR
% SER: [Perfect, Park, Fixed, WDFT, ProposedML] -> 5 x numSNR
mse = zeros(4, numSNR);
ser = zeros(5, numSNR);

%% 2. PRECOMPUTE STATIC OPERATORS
% Windowed DFT filter
K_wdft = local_build_wdft_matrix(N, 128, ltsK, usedK);

% Fixed SVD covariance (L = 64)
R64 = local_uniform_covariance(N, 64);
Rhp64 = R64(idxUsed, idxLTS);
Rpp64 = R64(idxLTS, idxLTS);
[U64, D64] = eig((Rpp64 + Rpp64') / 2);
[lam64, ord64] = sort(real(diag(D64)), 'descend');
U64 = U64(:, ord64);

% SVD bank for L = 1 to G
U_bank = cell(G, 1);
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

%% 3. SIMULATION LOOP OVER SNRS
fprintf('Evaluating %d SNR points across %d packets each...\n', numSNR, NUM_PACKETS);

for is = 1:numSNR
    snr = SNRdB(is);
    gammaP = 10^((snr + 6)/10); % Pilot boost +6 dB
    gammaD = 10^(snr/10);       % Data SNR

    % Fixed SVD Filter at this SNR
    delta64 = zeros(Np, 1);
    r64 = min(67, Np);
    delta64(1:r64) = 1 ./ (lam64(1:r64) + 1/gammaP);
    K_fixed = Rhp64 * U64 * diag(delta64) * U64';

    % SVD Filter bank at this SNR
    Kbank = cell(G, 1);
    for L = 1:G
        r = min(L + p.adaptiveSVDOffset, Np);
        dd = zeros(Np, 1);
        dd(1:r) = 1 ./ (lam_bank{L}(1:r) + 1/gammaP);
        Kbank{L} = Rhp_bank{L} * U_bank{L} * diag(dd) * U_bank{L}';
    end

    msePkt = zeros(4, NUM_PACKETS);
    serPkt = zeros(5, NUM_PACKETS);

    for pkt = 1:NUM_PACKETS
        % A. Time-Varying Channel Generation (Vehicular-A, fd = 30 Hz)
        h0 = local_vehicular_A(p);
        H0 = local_channel_cfr(h0, N);
        H0p = H0(idxLTS);

        % Evolve channel to data symbol
        h1 = local_time_evolve(h0, p, fd);
        H1 = local_channel_cfr(h1, N);
        H1u = H1(idxUsed);

        % B. LTS Preamble LS Observation
        noiseLS = sqrt(1/(2*gammaP)) * (randn(Np,1) + 1j*randn(Np,1));
        HLS = H0p + noiseLS;

        % C. CIR Transform & Length Estimation
        [cirLS,~,~] = ls_to_cir(HLS, 128, 'hann');
        noiseVar = mean(abs(cirLS(p.noiseTapStart : p.noiseTapEnd)).^2);

        % Park SHT
        L_park = sht_cir_length(cirLS, noiseVar, p.T, G);
        L_park = max(1, min(G, L_park));

        % Proposed ML
        L_ml = ml_cir_length(cirLS, noiseVar, G, snr);
        L_ml = max(1, min(G, L_ml));

        if pkt <= 2
            fprintf('  Doppler SNR=%ddB [pkt=%d]: L_SHT=%d, L_pred_ML=%d, L_final=%d\n', ...
                snr, pkt, L_park, L_ml, L_ml);
        end

        % D. Initial Channel Estimates at Preamble
        H_park0  = Kbank{L_park} * HLS;
        H_fixed0 = K_fixed * HLS;
        H_wdft0  = K_wdft * HLS;
        H_ml0    = Kbank{L_ml} * HLS;

        % E. Transmitted 16-QAM Data Symbol
        tx_data = local_qam16_symbols(Nd);
        pilot_vals = ones(numel(pilotK), 1);

        X = zeros(Nu, 1);
        X(idD) = tx_data;
        X(idP) = pilot_vals;

        noiseData = sqrt(1/(2*gammaD)) * (randn(Nu,1) + 1j*randn(Nu,1));
        Y = H1u .* X + noiseData;

        % Perfect CSI Lower Bound SER
        rxPerfect = local_qam16_detect(Y(idD) ./ (H1u(idD) + eps));
        serPkt(1, pkt) = mean(rxPerfect(:) ~= tx_data(:));

        % F. Decision-Directed Iterations for Each Estimator
        % Estimator ordering: 1: Park, 2: Fixed, 3: WDFT, 4: Proposed ML
        Hinit = {H_park0, H_fixed0, H_wdft0, H_ml0};
        Koper = {Kbank{L_park}, K_fixed, K_wdft, Kbank{L_ml}};

        for ee = 1:4
            He = Hinit{ee};
            K_filt = Koper{ee};

            for it = 1:NUM_DD_ITER
                % 1. Common Phase Error (CPE) correction using pilots
                zPilot = Y(idP) ./ (He(idP) .* pilot_vals + eps);
                phi = angle(sum(zPilot));
                Ycorr = Y .* exp(-1j * phi);

                % 2. Temporary Data Detection
                zData = Ycorr(idD) ./ (He(idD) + eps);
                detData = local_qam16_detect(zData);

                % 3. Pseudo Channel Observation
                Hobs = zeros(Nu, 1);
                Hobs(idP) = Ycorr(idP) ./ (pilot_vals + eps);
                Hobs(idD) = Ycorr(idD) ./ (detData + eps);

                % 4. Extract LTS Subcarriers & Re-filter
                Hpseudo = local_extract_lts(Hobs, usedK, ltsK);
                HeNew = K_filt * Hpseudo;

                % 5. Damped update to prevent decision-directed divergence
                He = 0.5 * He + 0.5 * HeNew;
            end

            % Final MSE on actual data channel H1u
            msePkt(ee, pkt) = mean(abs(He - H1u).^2);

            % Final SER with phase-corrected equalized symbols
            zPilot_final = Y(idP) ./ (He(idP) .* pilot_vals + eps);
            phi_final = angle(sum(zPilot_final));
            Y_final = Y .* exp(-1j * phi_final);
            shat = Y_final(idD) ./ (He(idD) + eps);
            rx_dec = local_qam16_detect(shat);
            serPkt(ee + 1, pkt) = mean(rx_dec(:) ~= tx_data(:));
        end
    end

    mse(:, is) = mean(msePkt, 2);
    ser(:, is) = mean(serPkt, 2);

    fprintf('  SNR = %2d dB | MSE: ML %.2e | Park %.2e | Fix %.2e | WDFT %.2e\n', ...
        snr, mse(4,is), mse(1,is), mse(2,is), mse(3,is));
    fprintf('               | SER: Perf %.2e | ML %.2e | Park %.2e | Fix %.2e | WDFT %.2e\n', ...
        ser(1,is), ser(5,is), ser(2,is), ser(3,is), ser(4,is));
end

%% 4. SUMMARY TABLE
fprintf('\n====================================================================================================\n');
fprintf(' TIME-VARYING (fd = 30 Hz) TRACKING PERFORMANCE SUMMARY TABLE\n');
fprintf('====================================================================================================\n');
fprintf(' SNR(dB) | MSE Proposed ML | MSE Park 2022 | MSE Fixed SVD | SER Proposed ML | SER Park 2022 | SER Fixed SVD\n');
fprintf('----------------------------------------------------------------------------------------------------\n');
for i = 1:numSNR
    fprintf('  %2d dB  |    %10.3e   |   %10.3e  |   %10.3e  |    %10.3e   |   %10.3e  |   %10.3e\n', ...
        SNRdB(i), mse(4,i), mse(1,i), mse(2,i), ser(5,i), ser(2,i), ser(3,i));
end
fprintf('====================================================================================================\n');

%% 5. PUBLICATION PLOT
fig = figure('Name', 'Novelty Plot 4 - Time-Varying Doppler Tracking', ...
    'Color', 'w', 'Position', [100 60 1100 520], 'Visible', 'off');

% (a) Channel MSE
subplot(1, 2, 1);
semilogy(SNRdB, mse(3,:), '-^c', 'LineWidth', 1.6, 'MarkerSize', 6); hold on;
semilogy(SNRdB, mse(2,:), '-xg', 'LineWidth', 1.8, 'MarkerSize', 8);
semilogy(SNRdB, mse(1,:), '-sb', 'LineWidth', 1.8, 'MarkerFaceColor', 'b', 'MarkerSize', 6);
semilogy(SNRdB, mse(4,:), '-or', 'LineWidth', 2.2, 'MarkerFaceColor', 'r', 'MarkerSize', 6);
grid on; box on;
set(gca, 'YScale', 'log', 'YMinorGrid', 'on', 'FontSize', 10);
xlim([10 40]); ylim([1e-4 1]);
xlabel('SNR (dB)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Channel Mean Square Error (MSE)', 'FontSize', 11, 'FontWeight', 'bold');
title('(a) Tracking MSE (fd = 30 Hz, 2 DD Iterations)', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Windowed DFT (L=64)', 'Fixed SVD (L=64)', 'Park 2022 Adaptive SVD', 'PROPOSED NOVEL ML-SVD'}, ...
    'Location', 'southwest', 'FontSize', 9);

% (b) 16-QAM SER
subplot(1, 2, 2);
semilogy(SNRdB, ser(1,:), ':k', 'LineWidth', 2.0); hold on;
semilogy(SNRdB, ser(4,:), '-^c', 'LineWidth', 1.6, 'MarkerSize', 6);
semilogy(SNRdB, ser(3,:), '-xg', 'LineWidth', 1.8, 'MarkerSize', 8);
semilogy(SNRdB, ser(2,:), '-sb', 'LineWidth', 1.8, 'MarkerFaceColor', 'b', 'MarkerSize', 6);
semilogy(SNRdB, ser(5,:), '-or', 'LineWidth', 2.2, 'MarkerFaceColor', 'r', 'MarkerSize', 6);
grid on; box on;
set(gca, 'YScale', 'log', 'YMinorGrid', 'on', 'FontSize', 10);
xlim([10 40]); ylim([5e-4 1]);
xlabel('SNR (dB)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('16-QAM Symbol Error Rate (SER)', 'FontSize', 11, 'FontWeight', 'bold');
title('(b) 16-QAM SER (fd = 30 Hz, 2 DD Iterations)', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Perfect CSI Knowledge', 'Windowed DFT (L=64)', 'Fixed SVD (L=64)', 'Park 2022 Adaptive SVD', 'PROPOSED NOVEL ML-SVD'}, ...
    'Location', 'southwest', 'FontSize', 9);

sgtitle('Time-Varying Rayleigh Channel Tracking: Proposed ML-Adaptive SVD vs Baselines', ...
    'FontSize', 13, 'FontWeight', 'bold');

% Save Artifacts
outPng = fullfile(projectRoot, 'ML', 'Novelty_Day3_Plot4_TimeVarying_Tracking.png');
outMat = fullfile(projectRoot, 'ML', 'Novelty_Day3_Plot4_TimeVarying_results.mat');

saveas(fig, outPng);
save(outMat, 'SNRdB', 'mse', 'ser', '-v7');

fprintf('\nSuccessfully saved Plot 4 results to:\n  %s\n  %s\n', outPng, outMat);
fprintf('============================================================\n');

end

%% --- HELPER FUNCTIONS ---
function HLS = local_extract_lts(Hused, usedK, ltsK)
    [~, loc] = ismember(ltsK, usedK);
    HLS = Hused(loc);
    HLS = HLS(:);
end

function R = local_uniform_covariance(N, L)
    k = (-N/2:N/2-1).';
    d = k - k.';
    R = ones(N);
    idx = (d ~= 0);
    R(idx) = (1 - exp(-1j*2*pi*L*d(idx)/N)) ./ (1j*2*pi*L*d(idx)/N);
    R = (R + R') / 2;
end

function K_WDFT = local_build_wdft_matrix(N, M, ltsLocations, usedLocations)
    Np = length(ltsLocations);
    L  = 64;
    evenK = ltsLocations(:);
    idx_m2 = find(evenK == -2);
    idx_p2 = find(evenK ==  2);

    D = zeros(Np+1, Np);
    D(1:50, 1:50)     = eye(50);
    D(51, idx_m2)     = 0.5;
    D(51, idx_p2)     = 0.5;
    D(52:101, 51:100) = eye(50);

    q = (-Np/2 : Np/2).';
    ell = 0 : L-1;
    Fp = exp(-1j * 2 * pi / M * (q * ell));
    A = exp(-1j * 2 * pi / N * (usedLocations(:) * ell));
    pinvFp = pinv(Fp);

    taperLen = 8;
    w_taper = ones(L, 1);
    n_taper = (0:taperLen-1).';
    w_taper(end-taperLen+1:end) = 0.5 * (1 + cos(pi * n_taper / taperLen));

    K_WDFT = A * diag(w_taper) * pinvFp * D;
end

function h = local_vehicular_A(p)
    delays = [0 0.31 0.71 1.09 1.73 2.01];
    pdB = [0 -1 -9 -10 -15 -20];
    delayS = round(delays * 1e-6 / p.Ts);
    amp = 10.^(pdB / 20);
    amp = amp / norm(amp);
    L = max(delayS) + 1;
    h = zeros(L, 1);
    for k = 1:numel(delayS)
        g = (randn + 1j*randn) / sqrt(2);
        h(delayS(k) + 1) = amp(k) * g;
    end
end

function h1 = local_time_evolve(h0, p, fd)
    h0 = h0(:);
    Tsym = (p.N + p.G) * p.Ts;
    rho = besselj(0, 2*pi*fd*Tsym);
    rho = max(-1, min(1, rho));
    h1 = zeros(size(h0));
    for k = 1:length(h0)
        q = (randn + 1j*randn) / sqrt(2);
        h1(k) = rho*h0(k) + sqrt(max(0, 1 - rho^2))*q;
    end
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
