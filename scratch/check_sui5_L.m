projectRoot = 'c:\Users\pushk\Downloads\Adaptive_OFDM_Channel_Estimation';
addpath(genpath(projectRoot));
p = parameters();
rng(2026);

for snr = [10, 20, 30]
    gammaP = 10^((snr + 6)/10);
    ltsK = [-100:2:-2 2:2:100].';
    idxLTS = ltsK + p.N/2 + 1;
    
    l_s = [];
    l_m = [];
    for tr = 1:50
        h = sui5_channel(p);
        H = fftshift(fft([h(:); zeros(p.N-numel(h),1)], p.N));
        Hp = H(idxLTS);
        nLS = sqrt(1/(2*gammaP))*(randn(100,1) + 1j*randn(100,1));
        [cir,~,~] = ls_to_cir(Hp+nLS, 128, 'hann');
        s2 = mean(abs(cir(p.noiseTapStart:p.noiseTapEnd)).^2);
        l_s(end+1) = sht_cir_length(cir, s2, p.T, p.G);
        l_m(end+1) = ml_cir_length(cir, s2, p.G, snr);
    end
    fprintf('SUI-5 SNR %ddB: SHT mean=%.2f (min=%d, max=%d) | ML mean=%.2f (min=%d, max=%d)\n', ...
        snr, mean(l_s), min(l_s), max(l_s), mean(l_m), min(l_m), max(l_m));
end
