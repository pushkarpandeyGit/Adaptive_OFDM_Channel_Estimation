projectRoot = 'c:\Users\pushk\Downloads\Adaptive_OFDM_Channel_Estimation';
addpath(genpath(projectRoot));
p = parameters();
rng(2026);

channels = {'SUI-2', 'Vehicular-A', 'SUI-5'};
snrs = [10, 20, 30];

for chIdx = 1:3
    chName = channels{chIdx};
    fprintf('\n========================================\n');
    fprintf('CHANNEL: %s\n', chName);
    fprintf('========================================\n');
    
    for snr = snrs
        gammaP = 10^((snr + 6)/10);
        l_pred_list = [];
        l_sht_list = [];
        l_optA_list = [];
        l_optB_list = [];
        
        for tr = 1:30
            if strcmp(chName, 'SUI-2')
                h = sui2_channel(p);
                L_true = 5;
            elseif strcmp(chName, 'Vehicular-A')
                h = local_vehicular_A(p);
                L_true = 9;
            else
                h = sui5_channel(p);
                L_true = 41;
            end
            
            H = local_channel_cfr(h, p.N);
            idxLTS = p.N/2 + 1 + [-100:2:-2 2:2:100].';
            Hp = H(idxLTS);
            nLS = sqrt(1/(2*gammaP)) * (randn(100,1) + 1j*randn(100,1));
            HLS = Hp + nLS;
            [cir,~,~] = ls_to_cir(HLS, 128, 'hann');
            sigma2 = mean(abs(cir(p.noiseTapStart : p.noiseTapEnd)).^2);
            
            L_SHT = sht_cir_length(cir, sigma2, p.T, p.G);
            
            % Get raw prediction
            [~, feat] = ml_cir_length(cir, sigma2, p.G, snr);
            modelData = load(fullfile(projectRoot, 'ML', 'trained_ml_cir_model.mat'));
            feat_norm = (feat - modelData.mu) ./ modelData.sigma;
            z1 = modelData.W1 * feat_norm.' + modelData.b1;
            a1 = max(0.01*z1, z1);
            z2 = modelData.W2 * a1 + modelData.b2;
            a2 = max(0.01*z2, z2);
            y_hat = modelData.W3 * a2 + modelData.b3;
            L_pred = round(y_hat);
            
            % Option A
            if L_pred >= L_SHT - 2
                L_optA = L_pred;
            else
                L_optA = ceil((L_pred + L_SHT) / 2);
            end
            L_optA = max(1, min(p.G, L_optA));
            
            % Option B
            L_optB = max(1, min(p.G, L_pred));
            
            l_pred_list(end+1) = L_pred;
            l_sht_list(end+1)  = L_SHT;
            l_optA_list(end+1) = L_optA;
            l_optB_list(end+1) = L_optB;
        end
        
        fprintf('SNR=%2d dB | True L=%2d | Mean SHT=%.2f | Mean Pred=%.2f | OptA=%.2f | OptB=%.2f\n', ...
            snr, L_true, mean(l_sht_list), mean(l_pred_list), mean(l_optA_list), mean(l_optB_list));
    end
end

function H = local_channel_cfr(h, N)
    h = h(:);
    ht = [h; zeros(N - numel(h), 1)];
    H = fftshift(fft(ht, N));
    H = H(:);
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
