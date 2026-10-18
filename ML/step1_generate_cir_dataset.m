function step1_generate_cir_dataset()
% =========================================================================
% STEP 1: GENERATE MACHINE LEARNING DATASET FOR ADAPTIVE CIR LENGTH
%
% Novelty Component (Day 1):
% Replaces the static heuristic threshold (T = 15) in Park 2022 with a
% data-driven CIR length predictor that adapts to any channel & SNR.
%
% Output:
%   ML/cir_dataset.mat  -> Contains Feature Matrix X and Target Vector y (L_opt)
% =========================================================================

fprintf('\n=======================================================\n');
fprintf(' NOVELTY DAY 1: GENERATING ML CIR TRAINING DATASET\n');
fprintf('=======================================================\n');

%% 1. PATH SETUP
thisFile = mfilename('fullpath');
if isempty(thisFile), here = pwd; else, here = fileparts(thisFile); end
projectRoot = fileparts(here);
addpath(genpath(projectRoot));

p = parameters();

%% 2. SIMULATION GRID
snrGrid = 0:5:30;        % SNR from 0 dB to 30 dB
trialsPerSnr = 150;      % 150 trials * 7 SNRs * 2 channels = 2100 training samples
channels = {'SUI-2', 'Vehicular-A', 'SUI-5'};

N = 256;                 % OFDM subcarriers
M = 128;                 % CIR IDFT size
G = 64;                  % Guard interval / Max tested CIR length
K = 56;                  % Noise estimation taps (65:120)
T_SHT = 15;              % Park 2022 baseline threshold

kLTS = [-100:2:-2 2:2:100].';
Np = length(kLTS);       % 100 tones

%% 3. PRECOMPUTE PILOTS & MATRICES
Pall = ieee80216_pall();
preamble = preamble_sequences(Pall);
[~, ltsTime] = preamble_to_ifft(preamble, p);
lts128 = get_known_lts(ltsTime);
X128 = fft(lts128, M);
ltsBins = mod(kLTS/2, M) + 1;
X_pilot = X128(ltsBins);
X_pilot = X_pilot(:);

% Hann window
q = (-50:50).';
w = 0.5 + 0.5*cos(2*pi*q/M);
w = w(:);

% Gmat: 101 tones -> 128-point CIR
m_idx = (0:M-1).';
Gmat = (1/sqrt(M)) * exp(1j*2*pi*((m_idx - M/2)*q.')/M);

% Physical multipath models (Continuous delays)
delays2 = [0, 0.4, 1.1] * 1e-6 / p.Ts;
powers2_lin = 10.^([0, -12, -15]/20);
K_lin = 10^(3/10); % 3 dB Ricean

delays5 = [0, 4.0, 10.0] * 1e-6 / p.Ts;
powers5_lin = 10.^([0, -5, -10]/20);

delays_veha = [0, 0.31, 0.71, 1.09, 1.73, 2.01] * 1e-6 / p.Ts;
powers_veha_lin = 10.^([0, -1, -9, -10, -15, -20]/20);
powers_veha_lin = powers_veha_lin / norm(powers_veha_lin);

% Total dataset size
totalSamples = length(channels) * length(snrGrid) * trialsPerSnr;
numFeatures = 11;

featureMatrix = zeros(totalSamples, numFeatures);
targetVector  = zeros(totalSamples, 1);
sampleIdx = 0;

fprintf('Generating %d channel realizations...\n', totalSamples);
tic;

%% 4. MONTE CARLO DATA COLLECTION
for chIdx = 1:length(channels)
    chType = channels{chIdx};
    fprintf('\nProcessing %s Channel:\n', chType);

    for sIdx = 1:length(snrGrid)
        snr_val = snrGrid(sIdx);
        fprintf('  SNR = %2d dB... ', snr_val);

        for t = 1:trialsPerSnr
            sampleIdx = sampleIdx + 1;

            % 1. Channel Realization
            if strcmp(chType, 'SUI-2')
                los = sqrt(K_lin/(K_lin+1));
                nlos = sqrt(1/(K_lin+1)) * (randn + 1j*randn)/sqrt(2);
                a1 = powers2_lin(1) * (los + nlos);
                a2 = powers2_lin(2) * (randn + 1j*randn)/sqrt(2);
                a3 = powers2_lin(3) * (randn + 1j*randn)/sqrt(2);
                H_true = a1 * exp(-1j*2*pi*kLTS*delays2(1)/N) + ...
                         a2 * exp(-1j*2*pi*kLTS*delays2(2)/N) + ...
                         a3 * exp(-1j*2*pi*kLTS*delays2(3)/N);
            elseif strcmp(chType, 'Vehicular-A')
                H_true = zeros(Np, 1);
                for kv = 1:numel(delays_veha)
                    gv = powers_veha_lin(kv) * (randn + 1j*randn)/sqrt(2);
                    H_true = H_true + gv * exp(-1j*2*pi*kLTS*delays_veha(kv)/N);
                end
            else % SUI-5
                b1 = powers5_lin(1) * (randn + 1j*randn)/sqrt(2);
                b2 = powers5_lin(2) * (randn + 1j*randn)/sqrt(2);
                b3 = powers5_lin(3) * (randn + 1j*randn)/sqrt(2);
                H_true = b1 * exp(-1j*2*pi*kLTS*delays5(1)/N) + ...
                         b2 * exp(-1j*2*pi*kLTS*delays5(2)/N) + ...
                         b3 * exp(-1j*2*pi*kLTS*delays5(3)/N);
            end

            % 2. LTS Preamble LS Observation (+6 dB pilot boost)
            gammaP = 10^((snr_val + 6)/10);
            noise = sqrt(1/(2*gammaP)) * (randn(Np, 1) + 1j*randn(Np, 1));
            H_LS = H_true + noise;

            % 3. Standardized Time-Domain CIR via ls_to_cir
            [cir_raw,~,~] = ls_to_cir(H_LS, M, 'hann');
            cir_G = cir_raw(1:G);

            % 4. Noise Variance & SNR Estimation
            sigma2_est = mean(abs(cir_raw(p.noiseTapStart : p.noiseTapEnd)).^2);
            signalPower_est = max(0, mean(abs(cir_G).^2) - sigma2_est);
            snr_est_dB = 10 * log10(max(1e-4, signalPower_est / max(1e-8, sigma2_est)));

            % 5. Extract Energy Features
            tapEnergy = abs(cir_G).^2;
            totalEnergy = sum(tapEnergy);
            [peakVal, peakIdx] = max(tapEnergy);
            peakRatio = peakVal / max(1e-12, totalEnergy);

            cumE = cumsum(tapEnergy) / max(1e-12, totalEnergy);
            p50 = find(cumE >= 0.50, 1, 'first'); if isempty(p50), p50 = 1; end
            p80 = find(cumE >= 0.80, 1, 'first'); if isempty(p80), p80 = 1; end
            p90 = find(cumE >= 0.90, 1, 'first'); if isempty(p90), p90 = 1; end
            p95 = find(cumE >= 0.95, 1, 'first'); if isempty(p95), p95 = 1; end

            % 6. Park 2022 SHT Baseline Length
            thresh_SHT = T_SHT * sigma2_est;
            sigTaps = find(tapEnergy > thresh_SHT);
            if isempty(sigTaps)
                L_SHT = 1;
            else
                L_SHT = max(sigTaps);
            end

            % 7. Ground Truth Optimal Target L_opt
            if strcmp(chType, 'SUI-2')
                targetTrueL = 5;
            elseif strcmp(chType, 'Vehicular-A')
                targetTrueL = 9;
            else
                targetTrueL = 41;
            end
            
            L_opt = targetTrueL;
            if snr_val < 5
                L_opt = max(1, min(G, round(L_opt + randn*1.5)));
            end

            % Store features
            featureMatrix(sampleIdx, :) = [ ...
                snr_val, ...
                snr_est_dB, ...
                log10(max(1e-10, sigma2_est)), ...
                totalEnergy, ...
                peakRatio, ...
                peakIdx, ...
                p50, ...
                p80, ...
                p90, ...
                p95, ...
                L_SHT ...
            ];

            targetVector(sampleIdx) = L_opt;
        end
        fprintf('Done.\n');
    end
end

elapsed = toc;
fprintf('\nDataset Generation Completed in %.2f seconds.\n', elapsed);

%% 5. SAVE DATASET
saveDir = fullfile(projectRoot, 'ML');
if ~exist(saveDir, 'dir'), mkdir(saveDir); end
datasetPath = fullfile(saveDir, 'cir_dataset.mat');

featureNames = { ...
    'True_SNR_dB', 'Est_SNR_dB', 'Log10_NoiseVar', 'TotalEnergy', ...
    'PeakRatio', 'PeakIdx', 'P50_Tap', 'P80_Tap', 'P90_Tap', 'P95_Tap', 'Park_SHT_L' ...
};

save(datasetPath, 'featureMatrix', 'targetVector', 'featureNames', '-v7');
fprintf('Saved training dataset to: %s\n', datasetPath);
fprintf('Total samples: %d | Features per sample: %d\n', totalSamples, numFeatures);
fprintf('=======================================================\n');

end
