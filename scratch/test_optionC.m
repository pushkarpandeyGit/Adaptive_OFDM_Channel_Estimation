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

ser_park = zeros(nS, 1);
ser_ml   = zeros(nS, 1);

for is = 1:nS
    snr = SNRdB(is);
    gammaP = 10^((snr + 6)/10);
    gammaD = 10^(snr/10);
    
    err_park = 0;
    err_ml   = 0;
    
    for tr = 1:nTrials
        h = sui2_channel(p);
        H = local_channel_cfr(h, p.N);
        Hu = H(idxUsed);
        Hp = H(idxLTS);
        
        nLS = sqrt(1/(2*gammaP)) * (randn(100,1) + 1j*randn(100,1));
        HLS = Hp + nLS;
        
        [cir,~,~] = ls_to_cir(HLS, 128, 'hann');
        sigma2 = mean(abs(cir(p.noiseTapStart : p.noiseTapEnd)).^2);
        
        L_park = sht_cir_length(cir, sigma2, p.T, p.G);
        L_park = max(1, min(p.G, L_park));
        
        L_ml = ml_cir_length(cir, sigma2, p.G, snr);
        L_ml = max(1, min(p.G, L_ml));
        
        % Park filter with rank = L_park + 3
        r_park = min(L_park + 3, 100);
        dd_park = zeros(100, 1);
        dd_park(1:r_park) = 1 ./ (lam_bank{L_park}(1:r_park) + 1/gammaP);
        K_park = Rhp_bank{L_park} * U_bank{L_park} * diag(dd_park) * U_bank{L_park}';
        
        % ML filter with adaptive margin: if L_ml <= 5, rank = L_ml + 2
        if L_ml <= 5
            r_ml = min(L_ml + 2, 100);
        else
            r_ml = min(L_ml + 3, 100);
        end
        dd_ml = zeros(100, 1);
        dd_ml(1:r_ml) = 1 ./ (lam_bank{L_ml}(1:r_ml) + 1/gammaP);
        K_ml = Rhp_bank{L_ml} * U_bank{L_ml} * diag(dd_ml) * U_bank{L_ml}';
        
        H_park = K_park * HLS;
        H_ml   = K_ml   * HLS;
        
        tx_data = local_qam16_symbols(numel(dataK));
        pilot_vals = ones(numel(pilotK), 1);
        X = zeros(200, 1);
        X(idD) = tx_data;
        X(idP) = pilot_vals;
        noise = sqrt(1/(2*gammaD)) * (randn(200,1) + 1j*randn(200,1));
        Y = Hu .* X + noise;
        
        % Detect Park
        z = Y(idP) ./ (H_park(idP) .* pilot_vals + eps);
        phi = angle(sum(z));
        Ye = Y .* exp(-1j * phi);
        shat = Ye(idD) ./ (H_park(idD) + eps);
        rx_park = local_qam16_detect(shat);
        err_park = err_park + mean(rx_park(:) ~= tx_data(:));
        
        % Detect ML
        z = Y(idP) ./ (H_ml(idP) .* pilot_vals + eps);
        phi = angle(sum(z));
        Ye = Y .* exp(-1j * phi);
        shat = Ye(idD) ./ (H_ml(idD) + eps);
        rx_ml = local_qam16_detect(shat);
        err_ml = err_ml + mean(rx_ml(:) ~= tx_data(:));
    end
    
    ser_park(is) = err_park / nTrials;
    ser_ml(is)   = err_ml / nTrials;
    
    fprintf('SUI-2 SNR %2d dB: Park SER = %.4e | ML SER = %.4e | ML <= Park? %d\n', ...
        snr, ser_park(is), ser_ml(is), ser_ml(is) <= ser_park(is));
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
