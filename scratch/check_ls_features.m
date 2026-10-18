projectRoot = 'c:\Users\pushk\Downloads\Adaptive_OFDM_Channel_Estimation';
addpath(genpath(projectRoot));
p = parameters();

snrGrid = [10, 20, 30];
channels = {'SUI-2', 'Vehicular-A', 'SUI-5'};

for ch = 1:3
    chName = channels{ch};
    fprintf('\n--- %s ---\n', chName);
    for s = 1:3
        snr = snrGrid(s);
        gammaP = 10^((snr + 6)/10);
        if ch == 1
            h = sui2_channel(p);
        elseif ch == 2
            delays = [0 0.31 0.71 1.09 1.73 2.01];
            pdB = [0 -1 -9 -10 -15 -20];
            delayS = round(delays * 1e-6 / p.Ts);
            amp = 10.^(pdB / 20); amp = amp / norm(amp);
            L = max(delayS) + 1; h = zeros(L, 1);
            for k = 1:numel(delayS)
                h(delayS(k) + 1) = amp(k) * (randn + 1j*randn)/sqrt(2);
            end
        else
            h = sui5_channel(p);
        end
        H = fftshift(fft([h(:); zeros(p.N-length(h),1)]));
        Hp = H(p.N/2+1 + [-100:2:-2 2:2:100]);
        HLS = Hp + sqrt(1/(2*gammaP))*(randn(100,1)+1j*randn(100,1));
        [cir,~,~] = ls_to_cir(HLS, 128, 'hann');
        s2 = mean(abs(cir(p.noiseTapStart:p.noiseTapEnd)).^2);
        
        cir_G = cir(1:p.G);
        tapE = abs(cir_G).^2;
        totE = sum(tapE);
        pRatio = max(tapE)/totE;
        cumE = cumsum(tapE)/totE;
        p50 = find(cumE>=0.5, 1);
        p80 = find(cumE>=0.8, 1);
        p90 = find(cumE>=0.9, 1);
        p95 = find(cumE>=0.95, 1);
        L_SHT = sht_cir_length(cir, s2, p.T, p.G);
        
        fprintf('SNR %2d dB: TotE=%.3f, LogS2=%.2f, P50=%d, P80=%d, P90=%d, P95=%d, SHT=%d\n', ...
            snr, totE, log10(s2), p50, p80, p90, p95, L_SHT);
    end
end
