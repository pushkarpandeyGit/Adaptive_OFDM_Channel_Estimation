function paper_fig7_reproduction()
% PAPER_FIG7_REPRODUCTION Reproduction of Fig. 7 (Quasi-static 16-QAM SER)
% Reference: H.-R. Park, IEEE Access 2022
%   (a) SUI-2 channel model
%   (b) SUI-5 channel model

fprintf('\n============================================================\n');
fprintf(' FIGURE 7 REPRODUCTION - QUASI-STATIC 16-QAM SER\n');
fprintf('============================================================\n');

thisFile = mfilename('fullpath');
if isempty(thisFile), here = pwd; else, here = fileparts(thisFile); end

if exist(fullfile(here,'config','parameters.m'),'file')
    projectRoot = here;
elseif exist(fullfile(fileparts(here),'config','parameters.m'),'file')
    projectRoot = fileparts(here);
else
    projectRoot = pwd;
end

addpath(genpath(projectRoot));
p = parameters();

SNRdB = 10:5:40;
nS = numel(SNRdB);
ser = zeros(5,nS,2);

usedK = p.activeIndex(:);
ltsK = [-100:2:-2 2:2:100].';
dataK = setdiff(usedK,p.pilotLocations(:),'stable');
pilotK = p.pilotLocations(:);

nTrials = 200; % Default fast verification trials (increase to 1000 for final publish)
fprintf('Trials/SNR/channel = %d\n', nTrials);

for ch = 1:2
    if ch == 1
        name = 'SUI-2';
        Lmodel = 5;
    else
        name = 'SUI-5';
        Lmodel = 40;
    end

    fprintf('\n--- %s ---\n', name);

    for is = 1:nS
        snr = SNRdB(is);
        gammaP = 10^((snr+6)/10);
        gammaD = 10^(snr/10);

        % LMMSE
        Rm = local_uniform_covariance(p.N, Lmodel);
        idxUsed = usedK + p.N/2 + 1;
        idxLTS  = ltsK  + p.N/2 + 1;
        Rhp = Rm(idxUsed, idxLTS);
        Rpp = Rm(idxLTS, idxLTS);
        K_lmmse = Rhp / (Rpp + (1/gammaP)*eye(100));

        % Fixed SVD: L=64, rank=67
        R64 = local_uniform_covariance(p.N, 64);
        Rhp64 = R64(idxUsed, idxLTS);
        Rpp64 = R64(idxLTS, idxLTS);
        [U,D] = eig((Rpp64+Rpp64')/2);
        [lam,ord] = sort(real(diag(D)),'descend');
        U = U(:,ord);
        delta = zeros(100,1);
        delta(1:min(67,100)) = 1./(lam(1:min(67,100)) + 1/gammaP);
        K_fixed = Rhp64 * U * diag(delta) * U';

        % Adaptive SVD filter bank
        Kbank = cell(p.G,1);
        for L = 1:p.G
            RL = local_uniform_covariance(p.N, L);
            RhpL = RL(idxUsed, idxLTS);
            RppL = RL(idxLTS, idxLTS);
            [UL,DL] = eig((RppL+RppL')/2);
            [la,oo] = sort(real(diag(DL)),'descend');
            UL = UL(:,oo);
            r = min(L+p.adaptiveSVDOffset, 100);
            dd = zeros(100,1);
            dd(1:r) = 1./(la(1:r) + 1/gammaP);
            Kbank{L} = RhpL * UL * diag(dd) * UL';
        end

        % Calibrated WDFT estimator (unbiased pseudoinverse Fourier projection)
        K_wdft = local_build_wdft_matrix(p.N, 128, ltsK, usedK);

        err = zeros(5, nTrials);

        for tr = 1:nTrials
            if ch == 1
                h = sui2_channel(p);
            else
                h = sui5_channel(p);
            end

            H = local_channel_cfr(h, p.N);
            Hu = H(usedK + p.N/2 + 1);
            Hp = H(ltsK  + p.N/2 + 1);

            nLS = sqrt(1/(2*gammaP)) * (randn(100,1) + 1j*randn(100,1));
            HLS = Hp + nLS;

            % Paper SHT
            [cir,~,~] = ls_to_cir(HLS, 128, 'hann');
            sigma2 = mean(abs(cir(p.G+1 : p.G+p.K)).^2);
            Lhat = sht_cir_length(cir, sigma2, p.T, p.G);
            Lhat = max(1, min(p.G, Lhat));

            H1 = K_lmmse * HLS;
            H2 = Kbank{Lhat} * HLS;
            H3 = K_fixed * HLS;
            H4 = K_wdft * HLS;

            Hset = {H1(:), H2(:), H3(:), H4(:)};

            tx = local_qam16_symbols(numel(dataK));
            pilot = ones(numel(pilotK),1);

            X = zeros(200,1);
            [~,idD] = ismember(dataK, usedK);
            [~,idP] = ismember(pilotK, usedK);
            X(idD) = tx;
            X(idP) = pilot;

            noise = sqrt(1/(2*gammaD)) * (randn(200,1) + 1j*randn(200,1));
            Y = Hu .* X + noise;

            for ee = 1:4
                He = Hset{ee};
                z = Y(idP) ./ (He(idP).*pilot + eps);
                phi = angle(sum(z));
                Ye = Y .* exp(-1j*phi);
                shat = Ye(idD) ./ (He(idD) + eps);
                rx = local_qam16_detect(shat);
                err(ee, tr) = mean(rx(:) ~= tx(:));
            end

            % Perfect CSI
            z = Y(idP) ./ (Hu(idP).*pilot + eps);
            phi = angle(sum(z));
            Ye = Y .* exp(-1j*phi);
            shat = Ye(idD) ./ (Hu(idD) + eps);
            rx = local_qam16_detect(shat);
            err(5, tr) = mean(rx(:) ~= tx(:));
        end

        ser(:, is, ch) = mean(err, 2);

        fprintf('%s SNR=%2d dB: LMMSE %.3e | AD %.3e | FIX %.3e | WDFT %.3e | PERFECT %.3e\n', ...
            name, snr, ser(1,is,ch), ser(2,is,ch), ser(3,is,ch), ser(4,is,ch), ser(5,is,ch));
    end
end

%% Plot Figure 7
fig = figure('Name', 'Figure 7 - Quasi-Static 16-QAM SER', ...
             'Color', 'w', 'Units', 'pixels', 'Position', [100 50 680 850]);

% (a) SUI-2
axA = subplot(2,1,1);
semilogy(SNRdB, ser(1,:,1), '-*k', 'LineWidth', 1.0, 'MarkerSize', 5); hold on;
semilogy(SNRdB, ser(2,:,1), '-sb', 'LineWidth', 1.0, 'MarkerSize', 5);
semilogy(SNRdB, ser(3,:,1), '-or', 'LineWidth', 1.0, 'MarkerSize', 5);
semilogy(SNRdB, ser(4,:,1), '-vk', 'LineWidth', 1.0, 'MarkerSize', 5);
semilogy(SNRdB, ser(5,:,1), '-+m', 'LineWidth', 1.0, 'MarkerSize', 5);
grid on; box on;
xlim([10 40]); ylim([1e-4 1]);
xticks(10:5:40); yticks([1e-4 1e-3 1e-2 1e-1 1]);
xlabel('SNR (dB)', 'FontSize', 10);
ylabel('Symbol Error Rate', 'FontSize', 10);
legend('LMMSE with L=5', 'Adaptive SVD', 'Fixed SVD with L=64', ...
       'Windowed DFT with L=64', 'Perfect channel knowledge', ...
       'Location', 'northeast', 'FontSize', 8);
text(0.5, -0.20, '(a)', 'Units', 'normalized', ...
     'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'normal');
set(axA, 'FontName', 'Times New Roman', 'FontSize', 10);

% Inset (a) zoom around 27 - 30 dB
posA = get(axA, 'Position');
axInsetA = axes('Position', [posA(1) + 0.08*posA(3), posA(2) + 0.12*posA(4), 0.35*posA(3), 0.40*posA(4)]);
semilogy(SNRdB, ser(1,:,1), '-*k', 'LineWidth', 0.8, 'MarkerSize', 4); hold on;
semilogy(SNRdB, ser(2,:,1), '-sb', 'LineWidth', 0.8, 'MarkerSize', 3.5);
semilogy(SNRdB, ser(3,:,1), '-or', 'LineWidth', 0.8, 'MarkerSize', 3.5);
semilogy(SNRdB, ser(4,:,1), '-vk', 'LineWidth', 0.8, 'MarkerSize', 3.5);
semilogy(SNRdB, ser(5,:,1), '-+m', 'LineWidth', 0.8, 'MarkerSize', 4);
xlim([27 30]); xticks(27:30);
ylim([0.0008 0.003]); yticks([0.001 0.002 0.003]);
yticklabels({'0.001','0.002','0.003'});
grid on; box on; set(axInsetA, 'FontName', 'Times New Roman', 'FontSize', 8);

% (b) SUI-5
axB = subplot(2,1,2);
semilogy(SNRdB, ser(1,:,2), '-*k', 'LineWidth', 1.0, 'MarkerSize', 5); hold on;
semilogy(SNRdB, ser(2,:,2), '-sb', 'LineWidth', 1.0, 'MarkerSize', 5);
semilogy(SNRdB, ser(3,:,2), '-or', 'LineWidth', 1.0, 'MarkerSize', 5);
semilogy(SNRdB, ser(4,:,2), '-vk', 'LineWidth', 1.0, 'MarkerSize', 5);
semilogy(SNRdB, ser(5,:,2), '-+m', 'LineWidth', 1.0, 'MarkerSize', 5);
grid on; box on;
xlim([10 40]); ylim([1e-4 1]);
xticks(10:5:40); yticks([1e-4 1e-3 1e-2 1e-1 1]);
xlabel('SNR (dB)', 'FontSize', 10);
ylabel('Symbol Error Rate', 'FontSize', 10);
legend('LMMSE with L=40', 'Adaptive SVD', 'Fixed SVD with L=64', ...
       'Windowed DFT with L=64', 'Perfect channel knowledge', ...
       'Location', 'northeast', 'FontSize', 8);
text(0.5, -0.20, '(b)', 'Units', 'normalized', ...
     'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'normal');
set(axB, 'FontName', 'Times New Roman', 'FontSize', 10);

% Inset (b) zoom around 27 - 30 dB
posB = get(axB, 'Position');
axInsetB = axes('Position', [posB(1) + 0.08*posB(3), posB(2) + 0.12*posB(4), 0.35*posB(3), 0.40*posB(4)]);
semilogy(SNRdB, ser(1,:,2), '-*k', 'LineWidth', 0.8, 'MarkerSize', 4); hold on;
semilogy(SNRdB, ser(2,:,2), '-sb', 'LineWidth', 0.8, 'MarkerSize', 3.5);
semilogy(SNRdB, ser(3,:,2), '-or', 'LineWidth', 0.8, 'MarkerSize', 3.5);
semilogy(SNRdB, ser(4,:,2), '-vk', 'LineWidth', 0.8, 'MarkerSize', 3.5);
semilogy(SNRdB, ser(5,:,2), '-+m', 'LineWidth', 0.8, 'MarkerSize', 4);
xlim([27 30]); xticks(27:30);
ylim([0.008 0.022]); yticks([0.01 0.015 0.02]);
yticklabels({'0.01','0.015','0.02'});
grid on; box on; set(axInsetB, 'FontName', 'Times New Roman', 'FontSize', 8);

fprintf('\nFIGURE 7 COMPLETE\n');
end

%% Helper functions
function R = local_uniform_covariance(N, L)
    k = (-N/2:N/2-1).';
    d = k - k.';
    R = ones(N);
    idx = (d ~= 0);
    R(idx) = (1 - exp(-1j*2*pi*L*d(idx)/N)) ./ (1j*2*pi*L*d(idx)/N);
    R = (R + R') / 2;
end

function K_WDFT = local_build_wdft_matrix(N, M, ltsLocations, usedLocations)
    Np = length(ltsLocations);
    L  = 64;

    evenK = ltsLocations(:);
    idx_m2 = find(evenK == -2);
    idx_p2 = find(evenK ==  2);

    D = zeros(Np+1, Np);
    D(1:50, 1:50)     = eye(50);
    D(51, idx_m2)     = 0.5;
    D(51, idx_p2)     = 0.5;
    D(52:101, 51:100) = eye(50);

    q = (-Np/2 : Np/2).';
    ell = 0 : L-1;

    Fp = exp(-1j * 2 * pi / M * (q * ell));
    A = exp(-1j * 2 * pi / N * (usedLocations(:) * ell));

    pinvFp = (Fp' * Fp) \ (Fp');

    taperLen = 8;
    w_taper = ones(L, 1);
    n_taper = (0:taperLen-1).';
    w_taper(end-taperLen+1:end) = 0.5 * (1 + cos(pi * n_taper / taperLen));

    K_WDFT = A * diag(w_taper) * pinvFp * D;
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
