projectRoot = 'c:\Users\pushk\Downloads\Adaptive_OFDM_Channel_Estimation';
addpath(genpath(projectRoot));
p = parameters();
rng(2026);

SNRdB = 10:5:35;
nS = numel(SNRdB);
nTrials = 300;

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

fprintf('\nEvaluating SUI-2 SER for Option A with rank = L_ml+3 vs Park:\n');
for is = 1:nS
    snr = SNRdB(is);
    gammaP = 10^((snr + 6)/10);
    gammaD = 10^(snr/10);
    
    err_park = 0;
    err_ml_r3 = 0;
    err_ml_r2 = 0;
    err_ml_r1 = 0;
    
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
        
        % Park filter
        r_park = min(L_park + 3, 100);
        dd_park = zeros(100, 1);
        dd_park(1:r_park) = 1 ./ (lam_bank{L_park}(1:r_park) + 1/gammaP);
        K_park = Rhp_bank{L_park} * U_bank{L_park} * diag(dd_park) * U_bank{L_park}';
        
        % ML r+3
        r3 = min(L_ml + 3, 100);
        dd3 = zeros(100, 1);
        dd3(1:r3) = 1 ./ (lam_bank{L_ml}(1:r3) + 1/gammaP);
        K3 = Rhp_bank{L_ml} * U_bank{L_ml} * diag(dd3) * U_bank{L_ml}';
        
        % ML r+2
        r2 = min(L_ml + 2, 100);
        dd2 = zeros(100, 1);
        dd2(1:r2) = 1 ./ (lam_bank{L_ml}(1:r2) + 1/gammaP);
        K2 = Rhp_bank{L_ml} * U_bank{L_ml} * diag(dd2) * U_bank{L_ml}';
        
        % ML r+1
        r1 = min(L_ml + 1, 100);
        dd1 = zeros(100, 1);
        dd1(1:r1) = 1 ./ (lam_bank{L_ml}(1:r1) + 1/gammaP);
        K1 = Rhp_bank{L_ml} * U_bank{L_ml} * diag(dd1) * U_bank{L_ml}';
        
        H_park = K_park * HLS;
        H_m3   = K3 * HLS;
        H_m2   = K2 * HLS;
        H_m1   = K1 * HLS;
        
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
        
        % Detect ML r+3
        z = Y(idP) ./ (H_m3(idP) .* pilot_vals + eps);
        Ye = Y .* exp(-1j * angle(sum(z)));
        rx = local_qam16_detect(Ye(idD) ./ (H_m3(idD) + eps));
        err_ml_r3 = err_ml_r3 + mean(rx(:) ~= tx_data(:));
        
        % Detect ML r+2
        z = Y(idP) ./ (H_m2(idP) .* pilot_vals + eps);
        Ye = Y .* exp(-1j * angle(sum(z)));
        rx = local_qam16_detect(Ye(idD) ./ (H_m2(idD) + eps));
        err_ml_r2 = err_ml_r2 + mean(rx(:) ~= tx_data(:));
        
        % Detect ML r+1
        z = Y(idP) ./ (H_m1(idP) .* pilot_vals + eps);
        Ye = Y .* exp(-1j * angle(sum(z)));
        rx = local_qam16_detect(Ye(idD) ./ (H_m1(idD) + eps));
        err_ml_r1 = err_ml_r1 + mean(rx(:) ~= tx_data(:));
    end
    
    fprintf('SNR %2d dB: Park=%.4e | ML(r+3)=%.4e | ML(r+2)=%.4e | ML(r+1)=%.4e\n', ...
        snr, err_park/nTrials, err_ml_r3/nTrials, err_ml_r2/nTrials, err_ml_r1/nTrials);
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
