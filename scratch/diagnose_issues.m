% Diagnostic script to investigate:
% 1. Why ML produced identical results to Park in Doppler (Vehicular-A)
% 2. Why ML produced worse SER in SUI-2

projectRoot = 'c:\Users\pushk\Downloads\Adaptive_OFDM_Channel_Estimation';
addpath(genpath(projectRoot));
p = parameters();
rng(42);

fprintf('========================================================\n');
fprintf('DIAGNOSIS 1: SUI-2 CIR LENGTH COMPARISON\n');
fprintf('========================================================\n');

for snr = [10, 15, 20, 25, 30]
    gammaP = 10^((snr + 6)/10);
    l_true_all = [];
    l_sht_all = [];
    l_ml_all = [];
    
    for tr = 1:50
        h = sui2_channel(p);
        H = local_channel_cfr(h, p.N);
        idxLTS = p.N/2 + 1 + [-100:2:-2 2:2:100].';
        Hp = H(idxLTS);
        nLS = sqrt(1/(2*gammaP)) * (randn(100,1) + 1j*randn(100,1));
        HLS = Hp + nLS;
        [cir,~,~] = ls_to_cir(HLS, 128, 'hann');
        sigma2 = mean(abs(cir(p.noiseTapStart : p.noiseTapEnd)).^2);
        
        L_park = sht_cir_length(cir, sigma2, p.T, p.G);
        [L_ml, feat] = ml_cir_length(cir, sigma2, p.G, snr);
        
        l_sht_all(end+1) = L_park;
        l_ml_all(end+1) = L_ml;
    end
    fprintf('SUI-2 SNR=%2d dB | Mean L_park = %.2f | Mean L_ml = %.2f\n', ...
        snr, mean(l_sht_all), mean(l_ml_all));
end

fprintf('\n========================================================\n');
fprintf('DIAGNOSIS 2: VEHICULAR-A (DOPPLER) CIR LENGTH COMPARISON\n');
fprintf('========================================================\n');

for snr = [10, 20, 30, 40]
    gammaP = 10^((snr + 6)/10);
    l_sht_all = [];
    l_ml_all = [];
    
    for tr = 1:50
        h = local_vehicular_A(p);
        H = local_channel_cfr(h, p.N);
        idxLTS = p.N/2 + 1 + [-100:2:-2 2:2:100].';
        Hp = H(idxLTS);
        nLS = sqrt(1/(2*gammaP)) * (randn(100,1) + 1j*randn(100,1));
        HLS = Hp + nLS;
        [cir,~,~] = ls_to_cir(HLS, 128, 'hann');
        noiseVar = mean(abs(cir(p.noiseTapStart : p.noiseTapEnd)).^2);
        
        L_park = sht_cir_length(cir, noiseVar, p.T, p.G);
        L_ml = ml_cir_length(cir, noiseVar, p.G, snr);
        
        l_sht_all(end+1) = L_park;
        l_ml_all(end+1) = L_ml;
    end
    diff_count = sum(l_sht_all ~= l_ml_all);
    fprintf('Veh-A SNR=%2d dB | Mean L_park = %.2f | Mean L_ml = %.2f | Different in %d/50 trials\n', ...
        snr, mean(l_sht_all), mean(l_ml_all), diff_count);
end

fprintf('\n========================================================\n');
fprintf('DIAGNOSIS 3: DETAILED INSPECTION OF ML PREDICTIONS INSIDE ml_cir_length\n');
fprintf('========================================================\n');

% Inspect raw neural prediction vs L_SHT
modelData = load(fullfile(projectRoot, 'ML', 'trained_ml_cir_model.mat'));
fprintf('Model loaded successfully. Inspecting raw y_hat vs L_SHT on 10 trials...\n');

for tr = 1:10
    h = local_vehicular_A(p);
    H = local_channel_cfr(h, p.N);
    idxLTS = p.N/2 + 1 + [-100:2:-2 2:2:100].';
    Hp = H(idxLTS);
    snr = 20;
    gammaP = 10^((snr + 6)/10);
    nLS = sqrt(1/(2*gammaP)) * (randn(100,1) + 1j*randn(100,1));
    HLS = Hp + nLS;
    [cir,~,~] = ls_to_cir(HLS, 128, 'hann');
    noiseVar = mean(abs(cir(p.noiseTapStart : p.noiseTapEnd)).^2);
    
    [L_hat, feat] = ml_cir_length(cir, noiseVar, p.G, snr);
    
    % Raw forward pass
    feat_norm = (feat - modelData.mu) ./ modelData.sigma;
    z1 = modelData.W1 * feat_norm.' + modelData.b1;
    a1 = max(0.01*z1, z1);
    z2 = modelData.W2 * a1 + modelData.b2;
    a2 = max(0.01*z2, z2);
    y_hat = modelData.W3 * a2 + modelData.b3;
    
    L_SHT = feat(11);
    fprintf('Trial %2d: Raw y_hat = %.2f, L_pred = %d, L_SHT = %d, L_final = %d\n', ...
        tr, y_hat, round(y_hat), L_SHT, L_hat);
end

function H = local_channel_cfr(h, N)
    h = h(:);
    ht = [h; zeros(N - numel(h), 1)];
    H = fftshift(fft(ht, N));
    H = H(:);
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
