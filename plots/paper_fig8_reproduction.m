function paper_fig8_reproduction()
% PAPER_FIG8_REPRODUCTION Reproduction of Fig. 8 (Time-varying Vehicular-A)
% Reference: H.-R. Park, IEEE Access 2022
%   fd = 30 Hz, 2 DD iterations
%   (a) MSE
%   (b) SER

fprintf('\n============================================================\n');
fprintf(' FIGURE 8 REPRODUCTION - TIME-VARYING VEHICULAR-A\n');
fprintf(' fd = 30 Hz, 2 DD iterations\n');
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

N       = p.N;                 % 256
G       = p.G;                 % 64
Ts      = p.Ts;                % 0.25 us
fd      = p.fd;                % 30 Hz

SNRdB   = p.SNRdB(:).';
numSNR  = numel(SNRdB);

usedK   = p.activeIndex(:);
ltsK    = [-100:2:-2 2:2:100].';
pilotK  = p.pilotLocations(:);

dataK   = setdiff(usedK, pilotK, 'stable');
dataK   = dataK(:);

Nu      = numel(usedK);
Np      = numel(ltsK);
Nd      = numel(dataK);

NUM_DD_ITER = 2;
NUM_PACKETS = 200; % fast verification trials (increase to 1000 for final publish)
fprintf('Packets/SNR = %d\n', NUM_PACKETS);

mse = zeros(3, numSNR);
ser = zeros(4, numSNR);

% Calibrated WDFT estimator (unbiased pseudoinverse Fourier projection)
K_wdft = local_build_wdft_matrix(N, 128, ltsK, usedK);

for is = 1:numSNR
    snrDB = SNRdB(is);
    fprintf('\nFig.8 SNR = %2d dB\n', snrDB);

    gammaP = 10^((snrDB+6)/10);
    gammaD = 10^(snrDB/10);

    msePkt = zeros(3, NUM_PACKETS);
    serPkt = zeros(4, NUM_PACKETS);

    for pkt = 1:NUM_PACKETS
        [h0, ~] = local_vehicular_A(p, fd);
        H0 = local_channel_cfr(h0, N);
        H0p = H0(ltsK + N/2 + 1);

        h1 = local_time_evolve(h0, p, fd);
        H1 = local_channel_cfr(h1, N);
        H1u = H1(usedK + N/2 + 1);

        noiseLS = sqrt(1/(2*gammaP)) * (randn(Np,1) + 1j*randn(Np,1));
        HLS = H0p + noiseLS;

        % Fixed SVD: L=64, rank=67
        R64 = local_uniform_covariance(N, 64);
        Kfixed = local_build_svd(R64, usedK, ltsK, 67, gammaP);
        Hfixed0 = Kfixed * HLS;

        % Adaptive SVD
        [cirLS,~,~] = ls_to_cir(HLS, 128, 'hann');
        noiseVar = mean(abs(cirLS(65:120)).^2);
        [Lhat,~,~,~] = sht_cir_length(cirLS, noiseVar, p.T, 64);
        Lhat = max(1, min(64, Lhat));
        padaptive = min(Lhat + 3, Np);

        Radaptive = local_uniform_covariance(N, Lhat);
        Kadaptive = local_build_svd(Radaptive, usedK, ltsK, padaptive, gammaP);
        Hadaptive0 = Kadaptive * HLS;

        % WDFT initial estimate
        Hwdf0 = K_wdft * HLS;
        Hwdf0 = Hwdf0(:);

        % Transmitted symbol
        tx = local_qam16_symbols(Nd);
        pilot = ones(numel(pilotK), 1);
        X = zeros(Nu, 1);
        [~,idD] = ismember(dataK, usedK);
        [~,idP] = ismember(pilotK, usedK);
        X(idD) = tx;
        X(idP) = pilot;

        noise = sqrt(1/(2*gammaD)) * (randn(Nu,1) + 1j*randn(Nu,1));
        Y = H1u .* X + noise;

        % Perfect CSI SER
        rxPerfect = local_qam16_detect(Y(idD) ./ (H1u(idD) + eps));
        serPkt(1, pkt) = mean(rxPerfect(:) ~= tx(:));

        % DD iterations for 3 estimators: Adaptive SVD, Fixed SVD, WDFT
        HestSet = {Hadaptive0, Hfixed0, Hwdf0};

        for ee = 1:3
            He = HestSet{ee};
            He = He(:);

            for it = 1:NUM_DD_ITER
                zPilot = Y(idP) ./ (He(idP).*pilot + eps);
                phi = angle(sum(zPilot));
                Ycorr = Y .* exp(-1j*phi);

                zData = Ycorr(idD) ./ (He(idD) + eps);
                detData = local_qam16_detect(zData);

                Hobs = zeros(Nu, 1);
                Hobs(idP) = Ycorr(idP) ./ (pilot + eps);
                Hobs(idD) = Ycorr(idD) ./ (detData + eps);

                Hpseudo = local_extract_lts(Hobs, usedK, ltsK);

                if ee == 1
                    HeNew = Kadaptive * Hpseudo;
                elseif ee == 2
                    HeNew = Kfixed * Hpseudo;
                else
                    HeNew = K_wdft * Hpseudo;
                end
                HeNew = HeNew(:);

                % Stabilized DD update
                He = 0.5*He + 0.5*HeNew;
            end

            msePkt(ee, pkt) = mean(abs(He - H1u).^2);

            rx = local_qam16_detect(Y(idD) ./ (He(idD) + eps));
            serPkt(ee+1, pkt) = mean(rx(:) ~= tx(:));
        end
    end

    mse(:, is) = mean(msePkt, 2);
    ser(:, is) = mean(serPkt, 2);

    fprintf('Fig.8 SNR=%2d dB | MSE: AD %.3e | FIX %.3e | WDFT %.3e\n', ...
        snrDB, mse(1,is), mse(2,is), mse(3,is));
    fprintf('                 | SER: PERF %.3e | AD %.3e | FIX %.3e | WDFT %.3e\n', ...
        ser(1,is), ser(2,is), ser(3,is), ser(4,is));
end

%% Plot Figure 8
fig = figure('Name', 'Figure 8 - Time-Varying Vehicular-A', ...
             'Color', 'w', 'Position', [100 50 720 850]);

% (a) MSE
subplot(2,1,1);
semilogy(SNRdB, mse(1,:), '-sb', 'LineWidth', 1.0, 'MarkerSize', 5); hold on;
semilogy(SNRdB, mse(2,:), '-or', 'LineWidth', 1.0, 'MarkerSize', 6);
semilogy(SNRdB, mse(3,:), '-kv', 'LineWidth', 1.0, 'MarkerSize', 6);
grid on; box on;
xlim([10 40]); ylim([1e-5 1]);
xticks(10:5:40);
xlabel('SNR (dB)', 'FontSize', 10);
ylabel('Mean Square Error', 'FontSize', 10);
legend('Adaptive SVD', 'Fixed SVD with L=64', 'Windowed DFT with L=64', ...
       'Location', 'northeast');
text(0.5, -0.20, '(a)', 'Units', 'normalized', ...
     'HorizontalAlignment', 'center', 'FontSize', 11);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);

% (b) SER
subplot(2,1,2);
semilogy(SNRdB, ser(1,:), '-*k', 'LineWidth', 1.0, 'MarkerSize', 6); hold on;
semilogy(SNRdB, ser(2,:), '-sb', 'LineWidth', 1.0, 'MarkerSize', 5);
semilogy(SNRdB, ser(3,:), '-or', 'LineWidth', 1.0, 'MarkerSize', 6);
semilogy(SNRdB, ser(4,:), '-kv', 'LineWidth', 1.0, 'MarkerSize', 6);
grid on; box on;
xlim([10 40]); ylim([1e-4 1]);
xticks(10:5:40);
xlabel('SNR (dB)', 'FontSize', 10);
ylabel('Symbol Error Rate', 'FontSize', 10);
legend('Perfect channel knowledge', 'Adaptive SVD', 'Fixed SVD with L=64', ...
       'Windowed DFT with L=64', 'Location', 'northeast');
text(0.5, -0.20, '(b)', 'Units', 'normalized', ...
     'HorizontalAlignment', 'center', 'FontSize', 11);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);

fprintf('\nFIGURE 8 COMPLETE\n');
end

%% Helper functions
function HLS = local_extract_lts(Hused, usedK, ltsK)
    Hused = Hused(:);
    usedK = usedK(:);
    ltsK  = ltsK(:);
    [tf, loc] = ismember(ltsK, usedK);
    if any(~tf)
        error('Could not map all LTS subcarriers.');
    end
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

function K = local_build_svd(Rmodel, usedK, ltsK, rankP, gammaP)
    N = size(Rmodel, 1);
    ia = usedK + N/2 + 1;
    ib = ltsK  + N/2 + 1;
    Rhp = Rmodel(ia, ib);
    Rpp = Rmodel(ib, ib);
    Rpp = (Rpp + Rpp') / 2;
    [U, D] = eig(Rpp);
    [lambda, ord] = sort(real(diag(D)), 'descend');
    U = U(:, ord);
    d = zeros(numel(lambda), 1);
    r = min(rankP, numel(lambda));
    d(1:r) = 1 ./ (lambda(1:r) + 1/gammaP);
    K = Rhp * U * diag(d) * U';
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

function [h, delayS] = local_vehicular_A(p, fd)
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
    if numel(h) > N
        error('Channel length exceeds FFT size.');
    end
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
