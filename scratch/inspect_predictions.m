projectRoot = 'c:\Users\pushk\Downloads\Adaptive_OFDM_Channel_Estimation';
addpath(genpath(projectRoot));
p = parameters();
rng(2026);

for chName = {'SUI-2', 'Vehicular-A', 'SUI-5'}
    ch = chName{1};
    fprintf('\nChannel: %s\n', ch);
    for snr = [10, 20, 30]
        gammaP = 10^((snr + 6)/10);
        ltsK = [-100:2:-2 2:2:100].';
        idxLTS = ltsK + p.N/2 + 1;
        
        preds = [];
        shts = [];
        for tr = 1:20
            if strcmp(ch, 'SUI-2')
                h = sui2_channel(p);
            elseif strcmp(ch, 'Vehicular-A')
                h = vehic_a_channel(p);
            else
                h = sui5_channel(p);
            end
            H = fftshift(fft([h(:); zeros(p.N - numel(h), 1)], p.N));
            Hp = H(idxLTS);
            nLS = sqrt(1/(2*gammaP)) * (randn(100,1) + 1j*randn(100,1));
            HLS = Hp + nLS;
            [cir,~,~] = ls_to_cir(HLS, 128, 'hann');
            sigma2 = mean(abs(cir(p.noiseTapStart : p.noiseTapEnd)).^2);
            
            L_SHT = sht_cir_length(cir, sigma2, p.T, p.G);
            [L_hat, feat] = ml_cir_length(cir, sigma2, p.G, snr);
            preds(end+1) = L_hat;
            shts(end+1) = L_SHT;
        end
        fprintf('  SNR %ddB: L_ML mean=%.2f (min=%d, max=%d) | L_SHT mean=%.2f (min=%d, max=%d)\n', ...
            snr, mean(preds), min(preds), max(preds), mean(shts), min(shts), max(shts));
    end
end
