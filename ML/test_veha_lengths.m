% Test script to check ML vs Park on Vehicular-A with max(L_pred, L_SHT)
p = parameters();
addpath(genpath(pwd));
delays = [0 0.31 0.71 1.09 1.73 2.01];
pdB = [0 -1 -9 -10 -15 -20];
delayS = round(delays * 1e-6 / p.Ts);
amp = 10.^(pdB / 20); amp = amp / norm(amp);
L_true = max(delayS) + 1;

N = p.N; G = p.G;
ltsK = [-100:2:-2 2:2:100].';
idxLTS = ltsK + N/2 + 1;

fprintf('Vehicular-A True CIR length: %d samples\n', L_true);
for snr = [10, 20, 30, 40]
    gammaP = 10^((snr+6)/10);
    L_p_list = zeros(50, 1);
    L_m_list = zeros(50, 1);
    for tr = 1:50
        h = zeros(L_true, 1);
        for k = 1:numel(delayS)
            h(delayS(k) + 1) = amp(k) * (randn + 1j*randn)/sqrt(2);
        end
        H = fftshift(fft([h; zeros(N-L_true,1)], N));
        H0p = H(idxLTS);
        noise = sqrt(1/(2*gammaP))*(randn(100,1)+1j*randn(100,1));
        HLS = H0p + noise;
        [cir,~,~] = ls_to_cir(HLS, 128, 'hann');
        sigma2 = mean(abs(cir(65:120)).^2);
        L_p = sht_cir_length(cir, sigma2, 15, G);
        [L_raw, ~] = ml_cir_length(cir, sigma2, G, snr);
        L_p_list(tr) = L_p;
        L_m_list(tr) = max(L_raw, L_p);
    end
    fprintf('SNR = %2d dB | True: %d | Park mean: %.2f (std: %.2f) | ML max mean: %.2f (std: %.2f)\n', ...
        snr, L_true, mean(L_p_list), std(L_p_list), mean(L_m_list), std(L_m_list));
end
