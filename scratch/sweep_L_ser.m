projectRoot = 'c:\Users\pushk\Downloads\Adaptive_OFDM_Channel_Estimation';
addpath(genpath(projectRoot));
p = parameters();
rng(2026);

snr = 20;
gammaP = 10^((snr + 6)/10);
gammaD = 10^(snr/10);

usedK  = p.activeIndex(:);
ltsK   = [-100:2:-2 2:2:100].';
pilotK = p.pilotLocations(:);
dataK  = setdiff(usedK, pilotK, 'stable');
[~, idD] = ismember(dataK, usedK);
[~, idP] = ismember(pilotK, usedK);

idxUsed = usedK + p.N/2 + 1;
idxLTS  = ltsK  + p.N/2 + 1;

fprintf('Testing SER in SUI-2 at 20 dB as a function of L_hat:\n');
for candidateL = [3, 4, 5, 6, 7, 8, 9, 10]
    RL = local_uniform_covariance(p.N, candidateL);
    RhpL = RL(idxUsed, idxLTS);
    RppL = RL(idxLTS, idxLTS);
    [UL, DL] = eig((RppL + RppL') / 2);
    [la, oo] = sort(real(diag(DL)), 'descend');
    UL = UL(:, oo);
    
    r = min(candidateL + 3, 100);
    dd = zeros(100, 1);
    dd(1:r) = 1 ./ (la(1:r) + 1/gammaP);
    KL = RhpL * UL * diag(dd) * UL';
    
    errList = zeros(150, 1);
    for tr = 1:150
        h = sui2_channel(p);
        H = local_channel_cfr(h, p.N);
        Hu = H(idxUsed);
        Hp = H(idxLTS);
        
        nLS = sqrt(1/(2*gammaP)) * (randn(100,1) + 1j*randn(100,1));
        HLS = Hp + nLS;
        
        He = KL * HLS;
        
        tx_data = local_qam16_symbols(numel(dataK));
        pilot_vals = ones(numel(pilotK), 1);
        X = zeros(200, 1);
        X(idD) = tx_data;
        X(idP) = pilot_vals;
        noise = sqrt(1/(2*gammaD)) * (randn(200,1) + 1j*randn(200,1));
        Y = Hu .* X + noise;
        
        z = Y(idP) ./ (He(idP) .* pilot_vals + eps);
        phi = angle(sum(z));
        Ye = Y .* exp(-1j * phi);
        shat = Ye(idD) ./ (He(idD) + eps);
        rx_dec = local_qam16_detect(shat);
        errList(tr) = mean(rx_dec(:) ~= tx_data(:));
    end
    fprintf('  L = %2d (Rank %2d) -> SER = %.4e\n', candidateL, r, mean(errList));
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
