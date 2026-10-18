projectRoot = 'c:\Users\pushk\Downloads\Adaptive_OFDM_Channel_Estimation';
addpath(genpath(projectRoot));
p = parameters();
rng(2026);

SNRdB = 10:5:35;
nS = numel(SNRdB);
nTrials = 250;

usedK  = p.activeIndex(:);
ltsK   = [-100:2:-2 2:2:100].';
pilotK = p.pilotLocations(:);
dataK  = setdiff(usedK, pilotK, 'stable');
[~, idD] = ismember(dataK, usedK);
[~, idP] = ismember(pilotK, usedK);

idxUsed = usedK + p.N/2 + 1;
idxLTS  = ltsK  + p.N/2 + 1;

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

channels = {'SUI-2', 'SUI-5'};
for ch = 1:2
    chName = channels{ch};
    fprintf('\n================== %s ==================\n', chName);
    
    for is = 1:nS
        snr = SNRdB(is);
        gammaP = 10^((snr + 6)/10);
        gammaD = 10^(snr/10);
        
        Kbank = cell(p.G, 1);
        for L = 1:p.G
            r = min(L + 3, 100);
            dd = zeros(100, 1);
            dd(1:r) = 1 ./ (lam_bank{L}(1:r) + 1/gammaP);
            Kbank{L} = Rhp_bank{L} * U_bank{L} * diag(dd) * U_bank{L}';
        end
        
        err_park = 0;
        err_ml   = 0;
        
        for tr = 1:nTrials
            if ch == 1
                h = sui2_channel(p);
            else
                h = sui5_channel(p);
            end
            H = local_channel_cfr(h, p.N);
            Hu = H(idxUsed);
            Hp = H(idxLTS);
            
            nLS = sqrt(1/(2*gammaP)) * (randn(100,1) + 1j*randn(100,1));
            HLS = Hp + nLS;
            
            [cir,~,~] = ls_to_cir(HLS, 128, 'hann');
            sigma2 = mean(abs(cir(p.noiseTapStart : p.noiseTapEnd)).^2);
            
            L_park = sht_cir_length(cir, sigma2, p.T, p.G);
            L_park = max(1, min(p.G, L_park));
            
            % ML prediction
            [~, feat] = ml_cir_length(cir, sigma2, p.G, snr);
            modelData = load('ML/trained_ml_cir_model.mat');
            feat_norm = (feat - modelData.mu) ./ modelData.sigma;
            z1 = modelData.W1 * feat_norm.' + modelData.b1;
            a1 = max(0.01*z1, z1);
            z2 = modelData.W2 * a1 + modelData.b2;
            a2 = max(0.01*z2, z2);
            y_hat = modelData.W3 * a2 + modelData.b3;
            L_pred = round(y_hat);
            
            % Smart combination:
            % For short channels (L_pred <= 15): trust ML completely, avoiding Park's overestimation!
            % For long channels (L_pred > 15): use max(L_pred, L_SHT) to protect long multipath!
            if L_pred <= 15
                if L_pred >= L_park - 2
                    L_ml = L_pred;
                else
                    L_ml = ceil((L_pred + L_park) / 2);
                end
            else
                L_ml = max(L_pred, L_park);
            end
            L_ml = max(1, min(p.G, L_ml));
            
            H_park = Kbank{L_park} * HLS;
            H_ml   = Kbank{L_ml}   * HLS;
            
            tx_data = local_qam16_symbols(numel(dataK));
            pilot_vals = ones(numel(pilotK), 1);
            X = zeros(200, 1);
            X(idD) = tx_data;
            X(idP) = pilot_vals;
            noise = sqrt(1/(2*gammaD)) * (randn(200,1) + 1j*randn(200,1));
            Y = Hu .* X + noise;
            
            % Detect Park
            z = Y(idP) ./ (H_park(idP) .* pilot_vals + eps);
            Ye = Y .* exp(-1j * angle(sum(z)));
            rx = local_qam16_detect(Ye(idD) ./ (H_park(idD) + eps));
            err_park = err_park + mean(rx(:) ~= tx_data(:));
            
            % Detect ML
            z = Y(idP) ./ (H_ml(idP) .* pilot_vals + eps);
            Ye = Y .* exp(-1j * angle(sum(z)));
            rx = local_qam16_detect(Ye(idD) ./ (H_ml(idD) + eps));
            err_ml = err_ml + mean(rx(:) ~= tx_data(:));
        end
        
        ser_p = err_park / nTrials;
        ser_m = err_ml / nTrials;
        fprintf('SNR %2d dB: Park=%.4e | ML=%.4e | ML <= Park? %d\n', ...
            snr, ser_p, ser_m, ser_m <= ser_p);
    end
end

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
