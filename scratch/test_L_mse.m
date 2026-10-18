projectRoot = 'c:\Users\pushk\Downloads\Adaptive_OFDM_Channel_Estimation';
addpath(genpath(projectRoot));
p = parameters();
rng(2026);

SNRs = [10, 15, 20, 25, 30, 35];
usedK  = p.activeIndex(:);
ltsK   = [-100:2:-2 2:2:100].';
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

fprintf('\nTesting Channel Estimation MSE for fixed L on SUI-2:\n');
fprintf('SNR(dB) |   L=4       L=5       L=6       L=7\n');
fprintf('-----------------------------------------------\n');

for s = 1:numel(SNRs)
    snr = SNRs(s);
    gammaP = 10^((snr + 6)/10);
    
    mse_L = zeros(4, 1);
    Ls = [4, 5, 6, 7];
    
    for tr = 1:200
        h = sui2_channel(p);
        H = local_channel_cfr(h, p.N);
        Hu = H(idxUsed);
        Hp = H(idxLTS);
        nLS = sqrt(1/(2*gammaP)) * (randn(100,1) + 1j*randn(100,1));
        HLS = Hp + nLS;
        
        for il = 1:4
            L = Ls(il);
            r = min(L + 3, 100);
            dd = zeros(100, 1);
            dd(1:r) = 1 ./ (lam_bank{L}(1:r) + 1/gammaP);
            K = Rhp_bank{L} * U_bank{L} * diag(dd) * U_bank{L}';
            H_est = K * HLS;
            mse_L(il) = mse_L(il) + mean(abs(Hu - H_est).^2);
        end
    end
    mse_L = mse_L / 200;
    fprintf(' %2d dB  | %.4e  %.4e  %.4e  %.4e\n', snr, mse_L(1), mse_L(2), mse_L(3), mse_L(4));
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
