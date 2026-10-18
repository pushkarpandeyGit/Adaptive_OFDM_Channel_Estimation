projectRoot = 'c:\Users\pushk\Downloads\Adaptive_OFDM_Channel_Estimation';
addpath(genpath(projectRoot));
p = parameters();
rng(2026);

snrGrid = 10:5:40;
numSNR = numel(snrGrid);
NUM_PACKETS = 80;
fd = 30;

usedK   = p.activeIndex(:);
ltsK    = [-100:2:-2 2:2:100].';
pilotK  = p.pilotLocations(:);
dataK   = setdiff(usedK, pilotK, 'stable');
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

mse_park = zeros(numSNR, 1);
mse_ml   = zeros(numSNR, 1);
ser_park = zeros(numSNR, 1);
ser_ml   = zeros(numSNR, 1);

diff_L_count = zeros(numSNR, 1);

for is = 1:numSNR
    snr = snrGrid(is);
    gammaP = 10^((snr + 6)/10);
    gammaD = 10^(snr/10);
    
    Kbank = cell(p.G, 1);
    for L = 1:p.G
        r = min(L + p.adaptiveSVDOffset, 100);
        dd = zeros(100, 1);
        dd(1:r) = 1 ./ (lam_bank{L}(1:r) + 1/gammaP);
        Kbank{L} = Rhp_bank{L} * U_bank{L} * diag(dd) * U_bank{L}';
    end
    
    msePkt = zeros(2, NUM_PACKETS);
    serPkt = zeros(2, NUM_PACKETS);
    
    for pkt = 1:NUM_PACKETS
        h0 = local_vehicular_A(p);
        H0 = local_channel_cfr(h0, p.N);
        H0p = H0(idxLTS);
        
        h1 = local_time_evolve(h0, p, fd);
        H1 = local_channel_cfr(h1, p.N);
        H1u = H1(idxUsed);
        
        noiseLS = sqrt(1/(2*gammaP)) * (randn(100,1) + 1j*randn(100,1));
        HLS = H0p + noiseLS;
        
        [cirLS,~,~] = ls_to_cir(HLS, 128, 'hann');
        noiseVar = mean(abs(cirLS(p.noiseTapStart : p.noiseTapEnd)).^2);
        
        L_park = sht_cir_length(cirLS, noiseVar, p.T, p.G);
        L_park = max(1, min(p.G, L_park));
        
        L_ml = ml_cir_length(cirLS, noiseVar, p.G, snr);
        L_ml = max(1, min(p.G, L_ml));
        
        if L_park ~= L_ml
            diff_L_count(is) = diff_L_count(is) + 1;
        end
        
        H_park0 = Kbank{L_park} * HLS;
        H_ml0   = Kbank{L_ml}   * HLS;
        
        tx_data = local_qam16_symbols(numel(dataK));
        pilot_vals = ones(numel(pilotK), 1);
        X = zeros(200, 1);
        X(idD) = tx_data;
        X(idP) = pilot_vals;
        noiseData = sqrt(1/(2*gammaD)) * (randn(200,1) + 1j*randn(200,1));
        Y = H1u .* X + noiseData;
        
        % 2 DD iterations
        He_park = H_park0;
        He_ml   = H_ml0;
        
        for it = 1:2
            % Park
            zP = Y(idP) ./ (He_park(idP) .* pilot_vals + eps);
            Yc = Y .* exp(-1j * angle(sum(zP)));
            dD = local_qam16_detect(Yc(idD) ./ (He_park(idD) + eps));
            Ho = zeros(200, 1); Ho(idP) = Yc(idP)./pilot_vals; Ho(idD) = Yc(idD)./dD;
            Hp_obs = local_extract_lts(Ho, usedK, ltsK);
            He_park = 0.5*He_park + 0.5*(Kbank{L_park} * Hp_obs);
            
            % ML
            zP = Y(idP) ./ (He_ml(idP) .* pilot_vals + eps);
            Yc = Y .* exp(-1j * angle(sum(zP)));
            dD = local_qam16_detect(Yc(idD) ./ (He_ml(idD) + eps));
            Ho = zeros(200, 1); Ho(idP) = Yc(idP)./pilot_vals; Ho(idD) = Yc(idD)./dD;
            Hp_obs = local_extract_lts(Ho, usedK, ltsK);
            He_ml = 0.5*He_ml + 0.5*(Kbank{L_ml} * Hp_obs);
        end
        
        msePkt(1, pkt) = mean(abs(He_park - H1u).^2);
        msePkt(2, pkt) = mean(abs(He_ml   - H1u).^2);
        
        % SER
        zP = Y(idP) ./ (He_park(idP) .* pilot_vals + eps);
        Yc = Y .* exp(-1j * angle(sum(zP)));
        rx = local_qam16_detect(Yc(idD) ./ (He_park(idD) + eps));
        serPkt(1, pkt) = mean(rx(:) ~= tx_data(:));
        
        zP = Y(idP) ./ (He_ml(idP) .* pilot_vals + eps);
        Yc = Y .* exp(-1j * angle(sum(zP)));
        rx = local_qam16_detect(Yc(idD) ./ (He_ml(idD) + eps));
        serPkt(2, pkt) = mean(rx(:) ~= tx_data(:));
    end
    
    mse_park(is) = mean(msePkt(1, :));
    mse_ml(is)   = mean(msePkt(2, :));
    ser_park(is) = mean(serPkt(1, :));
    ser_ml(is)   = mean(serPkt(2, :));
    
    fprintf('Doppler SNR %2d dB: L different in %2d/%d | Park MSE=%.3e | ML MSE=%.3e | Park SER=%.3e | ML SER=%.3e\n', ...
        snr, diff_L_count(is), NUM_PACKETS, mse_park(is), mse_ml(is), ser_park(is), ser_ml(is));
end

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
