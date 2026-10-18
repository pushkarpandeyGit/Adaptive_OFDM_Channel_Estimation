projectRoot = 'c:\Users\pushk\Downloads\Adaptive_OFDM_Channel_Estimation';
addpath(genpath(projectRoot));
p = parameters();
rng(2026);

SNRs = [10, 15, 20, 25, 30, 35];
for s = 1:numel(SNRs)
    snr = SNRs(s);
    gammaP = 10^((snr + 6)/10);
    ltsK = [-100:2:-2 2:2:100].';
    idxLTS = ltsK + p.N/2 + 1;
    
    L_preds = zeros(50, 1);
    L_SHTs = zeros(50, 1);
    L_hats = zeros(50, 1);
    
    for tr = 1:50
        h = sui2_channel(p);
        H = fftshift(fft([h(:); zeros(p.N - numel(h), 1)], p.N));
        Hp = H(idxLTS);
        nLS = sqrt(1/(2*gammaP)) * (randn(100,1) + 1j*randn(100,1));
        HLS = Hp + nLS;
        [cir,~,~] = ls_to_cir(HLS, 128, 'hann');
        sigma2 = mean(abs(cir(p.noiseTapStart : p.noiseTapEnd)).^2);
        
        L_SHT = sht_cir_length(cir, sigma2, p.T, p.G);
        [L_hat, feat] = ml_cir_length(cir, sigma2, p.G, snr);
        
        L_preds(tr) = L_hat;
        L_SHTs(tr) = L_SHT;
        L_hats(tr) = L_hat;
    end
    fprintf('SNR %2d dB: L_SHT mean=%.2f, L_ML mean=%.2f\n', snr, mean(L_SHTs), mean(L_hats));
end
