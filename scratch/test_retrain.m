% Test retraining with proper ls_to_cir
projectRoot = 'c:\Users\pushk\Downloads\Adaptive_OFDM_Channel_Estimation';
addpath(genpath(projectRoot));
p = parameters();

snrGrid = 0:5:30;
trialsPerSnr = 150;
channels = {'SUI-2', 'Vehicular-A', 'SUI-5'};
N = 256; M = 128; G = 64; K = 56;
kLTS = [-100:2:-2 2:2:100].';
Np = length(kLTS);

delays2 = [0, 0.4, 1.1] * 1e-6 / p.Ts;
powers2_lin = 10.^([0, -12, -15]/20);
K_lin = 10^(3/10);

delays5 = [0, 4.0, 10.0] * 1e-6 / p.Ts;
powers5_lin = 10.^([0, -5, -10]/20);

delays_veha = [0, 0.31, 0.71, 1.09, 1.73, 2.01] * 1e-6 / p.Ts;
powers_veha_lin = 10.^([0, -1, -9, -10, -15, -20]/20);
powers_veha_lin = powers_veha_lin / norm(powers_veha_lin);

totalSamples = length(channels) * length(snrGrid) * trialsPerSnr;
featureMatrix = zeros(totalSamples, 11);
targetVector  = zeros(totalSamples, 1);
sampleIdx = 0;

for chIdx = 1:length(channels)
    chType = channels{chIdx};
    for sIdx = 1:length(snrGrid)
        snr_val = snrGrid(sIdx);
        gammaP = 10^((snr_val + 6)/10);
        
        for t = 1:trialsPerSnr
            sampleIdx = sampleIdx + 1;
            
            if strcmp(chType, 'SUI-2')
                los = sqrt(K_lin/(K_lin+1));
                nlos = sqrt(1/(K_lin+1)) * (randn + 1j*randn)/sqrt(2);
                a1 = powers2_lin(1) * (los + nlos);
                a2 = powers2_lin(2) * (randn + 1j*randn)/sqrt(2);
                a3 = powers2_lin(3) * (randn + 1j*randn)/sqrt(2);
                H_true = a1 * exp(-1j*2*pi*kLTS*delays2(1)/N) + ...
                         a2 * exp(-1j*2*pi*kLTS*delays2(2)/N) + ...
                         a3 * exp(-1j*2*pi*kLTS*delays2(3)/N);
                targetTrueL = 5;
            elseif strcmp(chType, 'Vehicular-A')
                H_true = zeros(Np, 1);
                for kv = 1:numel(delays_veha)
                    gv = powers_veha_lin(kv) * (randn + 1j*randn)/sqrt(2);
                    H_true = H_true + gv * exp(-1j*2*pi*kLTS*delays_veha(kv)/N);
                end
                targetTrueL = 9;
            else
                b1 = powers5_lin(1) * (randn + 1j*randn)/sqrt(2);
                b2 = powers5_lin(2) * (randn + 1j*randn)/sqrt(2);
                b3 = powers5_lin(3) * (randn + 1j*randn)/sqrt(2);
                H_true = b1 * exp(-1j*2*pi*kLTS*delays5(1)/N) + ...
                         b2 * exp(-1j*2*pi*kLTS*delays5(2)/N) + ...
                         b3 * exp(-1j*2*pi*kLTS*delays5(3)/N);
                targetTrueL = 41;
            end
            
            noise = sqrt(1/(2*gammaP)) * (randn(Np, 1) + 1j*randn(Np, 1));
            H_LS = H_true + noise;
            
            [cir_raw,~,~] = ls_to_cir(H_LS, M, 'hann');
            cir_G = cir_raw(1:G);
            
            sigma2_est = mean(abs(cir_raw(p.noiseTapStart : p.noiseTapEnd)).^2);
            signalPower_est = max(0, mean(abs(cir_G).^2) - sigma2_est);
            snr_est_dB = 10 * log10(max(1e-4, signalPower_est / max(1e-8, sigma2_est)));
            
            tapEnergy = abs(cir_G).^2;
            totalEnergy = sum(tapEnergy);
            [peakVal, peakIdx] = max(tapEnergy);
            peakRatio = peakVal / max(1e-12, totalEnergy);
            
            cumE = cumsum(tapEnergy) / max(1e-12, totalEnergy);
            p50 = find(cumE >= 0.50, 1, 'first'); if isempty(p50), p50 = 1; end
            p80 = find(cumE >= 0.80, 1, 'first'); if isempty(p80), p80 = 1; end
            p90 = find(cumE >= 0.90, 1, 'first'); if isempty(p90), p90 = 1; end
            p95 = find(cumE >= 0.95, 1, 'first'); if isempty(p95), p95 = 1; end
            
            L_SHT = sht_cir_length(cir_raw, sigma2_est, p.T, G);
            
            L_opt = targetTrueL;
            if snr_val < 5
                L_opt = max(1, min(G, round(L_opt + randn*1.5)));
            end
            
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
    end
end

fprintf('Generated %d samples. Testing MLP training...\n', totalSamples);

% Quick train test
rng(2026);
perm = randperm(totalSamples);
splitIdx = round(0.80 * totalSamples);
X_train = featureMatrix(perm(1:splitIdx), :);
y_train = targetVector(perm(1:splitIdx));
X_test  = featureMatrix(perm(splitIdx+1:end), :);
y_test  = targetVector(perm(splitIdx+1:end));

mu = mean(X_train, 1);
sigma = std(X_train, 0, 1);
sigma(sigma == 0) = 1;

X_train_norm = (X_train - mu) ./ sigma;
X_test_norm  = (X_test - mu) ./ sigma;

D_in = 11; H1 = 24; H2 = 12; D_out = 1;
W1 = randn(H1, D_in) * sqrt(2 / D_in); b1 = zeros(H1, 1);
W2 = randn(H2, H1) * sqrt(2 / H1);     b2 = zeros(H2, 1);
W3 = randn(D_out, H2) * sqrt(2 / H2);  b3 = mean(y_train);

epochs = 350; batchSize = 64; lr = 0.008;
beta1 = 0.9; beta2 = 0.999; eps = 1e-8;
mW1 = zeros(size(W1)); vW1 = zeros(size(W1));
mb1 = zeros(size(b1)); vb1 = zeros(size(b1));
mW2 = zeros(size(W2)); vW2 = zeros(size(W2));
mb2 = zeros(size(b2)); vb2 = zeros(size(b2));
mW3 = zeros(size(W3)); vW3 = zeros(size(W3));
mb3 = zeros(size(b3)); vb3 = zeros(size(b3));

numBatches = floor(size(X_train_norm, 1) / batchSize);
t_step = 0;
for epoch = 1:epochs
    bPerm = randperm(size(X_train_norm, 1));
    X_shuff = X_train_norm(bPerm, :);
    y_shuff = y_train(bPerm);
    for b = 1:numBatches
        t_step = t_step + 1;
        idxB = (b-1)*batchSize + 1 : b*batchSize;
        xb = X_shuff(idxB, :).';
        yb = y_shuff(idxB).';
        z1 = W1 * xb + b1; a1 = max(0.01*z1, z1);
        z2 = W2 * a1 + b2; a2 = max(0.01*z2, z2);
        y_hat = W3 * a2 + b3;
        err = (y_hat - yb);
        penalty = ones(size(err));
        penalty(err < 0) = 2.0;
        grad_out = 2 * (penalty .* err) / batchSize;
        dW3 = grad_out * a2.'; db3 = sum(grad_out, 2);
        da2 = W3.' * grad_out; dz2 = da2 .* ((z2 > 0) + 0.01 * (z2 <= 0));
        dW2 = dz2 * a1.'; db2 = sum(dz2, 2);
        da1 = W2.' * dz2; dz1 = da1 .* ((z1 > 0) + 0.01 * (z1 <= 0));
        dW1 = dz1 * xb.'; db1 = sum(dz1, 2);
        mW1 = beta1*mW1 + (1-beta1)*dW1; vW1 = beta2*vW1 + (1-beta2)*(dW1.^2);
        W1 = W1 - lr * (mW1/(1-beta1^t_step)) ./ (sqrt(vW1/(1-beta2^t_step)) + eps);
        mb1 = beta1*mb1 + (1-beta1)*db1; vb1 = beta2*vb1 + (1-beta2)*(db1.^2);
        b1 = b1 - lr * (mb1/(1-beta1^t_step)) ./ (sqrt(vb1/(1-beta2^t_step)) + eps);
        mW2 = beta1*mW2 + (1-beta1)*dW2; vW2 = beta2*vW2 + (1-beta2)*(dW2.^2);
        W2 = W2 - lr * (mW2/(1-beta1^t_step)) ./ (sqrt(vW2/(1-beta2^t_step)) + eps);
        mb2 = beta1*mb2 + (1-beta1)*db2; vb2 = beta2*vb2 + (1-beta2)*(db2.^2);
        b2 = b2 - lr * (mb2/(1-beta1^t_step)) ./ (sqrt(vb2/(1-beta2^t_step)) + eps);
        mW3 = beta1*mW3 + (1-beta1)*dW3; vW3 = beta2*vW3 + (1-beta2)*(dW3.^2);
        W3 = W3 - lr * (mW3/(1-beta1^t_step)) ./ (sqrt(vW3/(1-beta2^t_step)) + eps);
        mb3 = beta1*mb3 + (1-beta1)*db3; vb3 = beta2*vb3 + (1-beta2)*(db3.^2);
        b3 = b3 - lr * (mb3/(1-beta1^t_step)) ./ (sqrt(vb3/(1-beta2^t_step)) + eps);
    end
end

z1_test = W1 * X_test_norm.' + b1; a1_test = max(0.01*z1_test, z1_test);
z2_test = W2 * a1_test + b2; a2_test = max(0.01*z2_test, z2_test);
y_pred_ML = max(1, min(64, round((W3 * a2_test + b3).')));
y_park = X_test(:, 11);

mae_ML = mean(abs(y_pred_ML - y_test));
mae_Park = mean(abs(y_park - y_test));
r2_ML = 1 - sum((y_pred_ML - y_test).^2) / sum((y_test - mean(y_test)).^2);
r2_Park = 1 - sum((y_park - y_test).^2) / sum((y_test - mean(y_test)).^2);

fprintf('\nTEST RESULTS:\n');
fprintf('MAE Park: %.2f taps | MAE ML: %.2f taps\n', mae_Park, mae_ML);
fprintf('R2 Park : %.4f      | R2 ML : %.4f\n', r2_Park, r2_ML);

% Check predictions on fresh channel realizations
fprintf('\nFresh predictions on actual channels:\n');
for ch = 1:3
    chName = channels{ch};
    if ch == 1
        h = sui2_channel(p); targetL = 5;
    elseif ch == 2
        delays = [0 0.31 0.71 1.09 1.73 2.01]; pdB = [0 -1 -9 -10 -15 -20];
        delayS = round(delays * 1e-6 / p.Ts); amp = 10.^(pdB / 20); amp = amp / norm(amp);
        L = max(delayS) + 1; h = zeros(L, 1);
        for k = 1:numel(delayS), h(delayS(k) + 1) = amp(k) * (randn + 1j*randn)/sqrt(2); end
        targetL = 9;
    else
        h = sui5_channel(p); targetL = 41;
    end
    H = fftshift(fft([h(:); zeros(p.N-length(h),1)]));
    Hp = H(p.N/2+1 + [-100:2:-2 2:2:100]);
    gammaP = 10^(26/10); % 20 dB
    HLS = Hp + sqrt(1/(2*gammaP))*(randn(100,1)+1j*randn(100,1));
    [cir,~,~] = ls_to_cir(HLS, 128, 'hann');
    s2 = mean(abs(cir(p.noiseTapStart:p.noiseTapEnd)).^2);
    cir_G = cir(1:G); tapE = abs(cir_G).^2; totE = sum(tapE);
    [pVal, pIdx] = max(tapE); pRatio = pVal / totE;
    cumE = cumsum(tapE)/totE;
    p50 = find(cumE>=0.5,1); p80 = find(cumE>=0.8,1); p90 = find(cumE>=0.9,1); p95 = find(cumE>=0.95,1);
    sSHT = sht_cir_length(cir, s2, p.T, G);
    
    feat = [20, 20, log10(max(1e-10, s2)), totE, pRatio, pIdx, p50, p80, p90, p95, sSHT];
    fn = (feat - mu) ./ sigma;
    z1 = W1 * fn.' + b1; a1 = max(0.01*z1, z1);
    z2 = W2 * a1 + b2; a2 = max(0.01*z2, z2);
    pred = round(W3 * a2 + b3);
    fprintf('%s (True L=%d): Park SHT=%d | ML Pred=%d\n', chName, targetL, sSHT, pred);
end
