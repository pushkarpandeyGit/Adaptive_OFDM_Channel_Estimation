%% FIGURE 6 - FULL DCS SIMULATION
% WDFT corrected version
clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
projectRoot = thisDir;
while ~exist(fullfile(projectRoot,'config','parameters.m'),'file')
    parentDir = fileparts(projectRoot);
    if strcmp(parentDir,projectRoot)
        error('Could not find config/parameters.m');
    end
    projectRoot = parentDir;
end
addpath(genpath(projectRoot));

p = parameters();
rng(p.randomSeed);

N = p.N; Nu = p.Nu; G = p.G;
M = 128; Np = 100; Lw = 64;
SNR_values = 0:5:40;
numSNR = numel(SNR_values);
numTrials = 1000;
T_SHT = 15;
noiseStart = p.noiseTapStart;
noiseEnd = p.noiseTapEnd;
fixedRank = G + p.fixedSVDOffset;

ltsLocations = [-100:2:-2 2:2:100];
usedLocations = p.activeIndex;

assert(numel(ltsLocations)==100,'LTS count must be 100');
assert(numel(usedLocations)==200,'Used-tone count must be 200');

Ltrue2 = 5;
Ltrue5 = 40;

MSE_LMMSE_2 = zeros(1,numSNR);
MSE_ADAPT_2 = zeros(1,numSNR);
MSE_FIXED_2 = zeros(1,numSNR);
MSE_WDFT_2  = zeros(1,numSNR);

MSE_LMMSE_5 = zeros(1,numSNR);
MSE_ADAPT_5 = zeros(1,numSNR);
MSE_FIXED_5 = zeros(1,numSNR);
MSE_WDFT_5  = zeros(1,numSNR);

R2 = local_uniform_covariance(N,Ltrue2);
R5 = local_uniform_covariance(N,Ltrue5);
Rf = local_uniform_covariance(N,G);

Rhh2 = local_subcov(R2,usedLocations,usedLocations,N);
Rhp2 = local_subcov(R2,usedLocations,ltsLocations,N);
Rpp2 = local_subcov(R2,ltsLocations,ltsLocations,N);

Rhh5 = local_subcov(R5,usedLocations,usedLocations,N);
Rhp5 = local_subcov(R5,usedLocations,ltsLocations,N);
Rpp5 = local_subcov(R5,ltsLocations,ltsLocations,N);

Rhpf = local_subcov(Rf,usedLocations,ltsLocations,N);
Rppf = local_subcov(Rf,ltsLocations,ltsLocations,N);

[Uf,Df] = eig(Rppf);
[lamf,ii] = sort(real(diag(Df)),'descend');
Uf = Uf(:,ii);
lamf = max(lamf,0);

% -------------------------------------------------------------------------
% Correct WDFT:
% 100 physical pilot tones -> 128-point compressed frequency grid
% -> IDFT -> zero-delay alignment -> 64 causal taps
% -> 256-point frequency reconstruction.
%
% IMPORTANT: WDFT is calibrated with a flat-channel test so that the
% MATLAB FFT normalization cannot create an artificial 1/M or M gain.
% -------------------------------------------------------------------------
K_WDFT = local_wdft_matrix(N,M,Lw,ltsLocations,usedLocations);

% Structural tests
H0p = ones(Np,1);
H0u = K_WDFT*H0p;
fprintf('\nWDFT flat-channel MSE = %.6e\n',mean(abs(H0u-1).^2));

delay = 10;
Hdp = exp(-1j*2*pi/N*(ltsLocations(:)*delay));
Hdu = exp(-1j*2*pi/N*(usedLocations(:)*delay));
Hde = K_WDFT*Hdp;
fprintf('WDFT 10-tap delay MSE = %.6e\n\n',mean(abs(Hde-Hdu).^2));

for s = 1:numSNR

    SNRdB = SNR_values(s);
    gamma = 10^((SNRdB+6)/10);

    % ============================= SUI-2 ===============================
    KL = Rhp2/(Rpp2+(1/gamma)*eye(Np));
    MSE_LMMSE_2(s) = local_mse(Rhh2,Rhp2,Rpp2,KL,gamma,Nu);

    Dp = local_svd_filter(lamf,fixedRank,gamma);
    KF = Rhpf*Uf*Dp*Uf';
    MSE_FIXED_2(s) = local_mse(Rhh2,Rhp2,Rpp2,KF,gamma,Nu);

    MSE_WDFT_2(s) = local_mse(Rhh2,Rhp2,Rpp2,K_WDFT,gamma,Nu);

    vals = zeros(numTrials,1);
    Lvals = zeros(numTrials,1);

    for t = 1:numTrials
        h = sui2_channel(p);
        H = local_cfr(h,N);

        Hu = H(usedLocations+N/2+1);
        Hp = H(ltsLocations+N/2+1);

        noise = sqrt(1/(2*gamma))*(randn(Np,1)+1j*randn(Np,1));
        HLS = Hp(:)+noise;

        cir = local_ls_to_cir(HLS,M);
        sigma2 = mean(abs(cir(noiseStart:noiseEnd)).^2);
        Lhat = local_sht(cir,sigma2,T_SHT,G);
        Lhat = max(1,min(G,Lhat));
        Lvals(t) = Lhat;

        Ra = local_uniform_covariance(N,Lhat);
        Rha = local_subcov(Ra,usedLocations,ltsLocations,N);
        Rpa = local_subcov(Ra,ltsLocations,ltsLocations,N);

        [Ua,Da] = eig(Rpa);
        [la,jj] = sort(real(diag(Da)),'descend');
        Ua = Ua(:,jj);
        la = max(la,0);

        rankA = min(Lhat+p.adaptiveSVDOffset,Np);
        Dpa = local_svd_filter(la,rankA,gamma);

        Hhat = Rha*Ua*Dpa*Ua'*HLS;
        vals(t) = mean(abs(Hhat(:)-Hu(:)).^2);
    end

    MSE_ADAPT_2(s) = mean(vals);

    % ============================= SUI-5 ===============================
    KL = Rhp5/(Rpp5+(1/gamma)*eye(Np));
    MSE_LMMSE_5(s) = local_mse(Rhh5,Rhp5,Rpp5,KL,gamma,Nu);

    KF = Rhpf*Uf*Dp*Uf';
    MSE_FIXED_5(s) = local_mse(Rhh5,Rhp5,Rpp5,KF,gamma,Nu);

    MSE_WDFT_5(s) = local_mse(Rhh5,Rhp5,Rpp5,K_WDFT,gamma,Nu);

    vals = zeros(numTrials,1);
    Lvals5 = zeros(numTrials,1);

    for t = 1:numTrials
        h = sui5_channel(p);
        H = local_cfr(h,N);

        Hu = H(usedLocations+N/2+1);
        Hp = H(ltsLocations+N/2+1);

        noise = sqrt(1/(2*gamma))*(randn(Np,1)+1j*randn(Np,1));
        HLS = Hp(:)+noise;

        cir = local_ls_to_cir(HLS,M);
        sigma2 = mean(abs(cir(noiseStart:noiseEnd)).^2);
        Lhat = local_sht(cir,sigma2,T_SHT,G);
        Lhat = max(1,min(G,Lhat));
        Lvals5(t) = Lhat;

        Ra = local_uniform_covariance(N,Lhat);
        Rha = local_subcov(Ra,usedLocations,ltsLocations,N);
        Rpa = local_subcov(Ra,ltsLocations,ltsLocations,N);

        [Ua,Da] = eig(Rpa);
        [la,jj] = sort(real(diag(Da)),'descend');
        Ua = Ua(:,jj);
        la = max(la,0);

        rankA = min(Lhat+p.adaptiveSVDOffset,Np);
        Dpa = local_svd_filter(la,rankA,gamma);

        Hhat = Rha*Ua*Dpa*Ua'*HLS;
        vals(t) = mean(abs(Hhat(:)-Hu(:)).^2);
    end

    MSE_ADAPT_5(s) = mean(vals);

    fprintf('%2d dB | SUI2 WDFT %.4e | SUI5 WDFT %.4e | Lhat %.2f / %.2f\n', ...
        SNRdB,MSE_WDFT_2(s),MSE_WDFT_5(s),mean(Lvals),mean(Lvals5));
end

% =============================== PLOT ==================================
figure('Color','w','Position',[100 50 700 850]);

subplot(2,1,1);
semilogy(SNR_values,MSE_LMMSE_2,'-k*','LineWidth',1,'MarkerSize',7); hold on;
semilogy(SNR_values,MSE_ADAPT_2,'-bs','LineWidth',1,'MarkerSize',5);
semilogy(SNR_values,MSE_FIXED_2,'-ro','LineWidth',1,'MarkerSize',6);
semilogy(SNR_values,MSE_WDFT_2,'-kv','LineWidth',1,'MarkerSize',6);
grid on; box on; xlim([0 40]); xticks(0:5:40);
ylim([1e-7 1e1]);
ylabel('Mean Square Error');
legend('LMMSE with L=5','Adaptive SVD (sim)', ...
    'Fixed SVD with L=64','Windowed DFT with L=64','Location','northeast');
text(.5,-.20,'(a)','Units','normalized','HorizontalAlignment','center');
set(gca,'FontName','Times New Roman','FontSize',10);

subplot(2,1,2);
semilogy(SNR_values,MSE_LMMSE_5,'-k*','LineWidth',1,'MarkerSize',7); hold on;
semilogy(SNR_values,MSE_ADAPT_5,'-bs','LineWidth',1,'MarkerSize',5);
semilogy(SNR_values,MSE_FIXED_5,'-ro','LineWidth',1,'MarkerSize',6);
semilogy(SNR_values,MSE_WDFT_5,'-kv','LineWidth',1,'MarkerSize',6);
grid on; box on; xlim([0 40]); xticks(0:5:40);
ylim([1e-7 1e1]);
xlabel('SNR (dB)'); ylabel('Mean Square Error');
legend('LMMSE with L=40','Adaptive SVD (sim)', ...
    'Fixed SVD with L=64','Windowed DFT with L=64','Location','northeast');
text(.5,-.20,'(b)','Units','normalized','HorizontalAlignment','center');
set(gca,'FontName','Times New Roman','FontSize',10);

resultDir = fullfile(projectRoot,'evaluation');
if ~exist(resultDir,'dir'), mkdir(resultDir); end

save(fullfile(resultDir,'paper_fig6_results.mat'), ...
    'SNR_values','MSE_LMMSE_2','MSE_ADAPT_2','MSE_FIXED_2','MSE_WDFT_2', ...
    'MSE_LMMSE_5','MSE_ADAPT_5','MSE_FIXED_5','MSE_WDFT_5');

exportgraphics(gcf,fullfile(resultDir,'paper_fig6_WDFT_fixed.png'),'Resolution',300);

fprintf('\nFigure 6 saved in: %s\n',resultDir);

%% ========================================================================
function K = local_wdft_matrix(N,M,L,ltsLocations,usedLocations)

    Np = numel(ltsLocations);

    % ==============================================================
    % STEP 1: Map physical LTS locations onto compressed 128-point
    % frequency grid.
    %
    % Physical locations:
    %
    %   -100,-98,...,-2,+2,...,+98,+100
    %
    % Since the spacing is 2:
    %
    %   q = k/2
    %
    % gives:
    %
    %   -50,...,-1,+1,...,+50
    % ==============================================================

    q = ltsLocations/2;

    X = zeros(M,Np);

    for n = 1:Np

        row = mod(q(n)+M/2,M)+1;

        X(row,n) = 1;

    end

    % ==============================================================
    % STEP 2: Convert each LTS basis vector into CIR.
    %
    % Each column represents the CIR obtained from one LTS
    % frequency observation.
    % ==============================================================

    C = zeros(M,Np);

    for n = 1:Np

        Xn = X(:,n);

        C(:,n) = ...
            fftshift( ...
            ifft( ...
            ifftshift(Xn), ...
            M));

    end

    % ==============================================================
    % STEP 3: Zero delay is at M/2 + 1 after fftshift.
    %
    % Keep L = 64 causal taps:
    %
    %   h[0], h[1], ..., h[63]
    % ==============================================================

    zeroDelay = M/2 + 1;

    causalIndex = ...
        zeroDelay:(zeroDelay+L-1);

    Ccausal = C(causalIndex,:);

    % ==============================================================
    % STEP 4: Apply the window in the CIR domain.
    %
    % IMPORTANT:
    %
    % We do NOT window the sparse LTS observations here.
    % The purpose of WDFT is to suppress the unwanted CIR tail.
    % ==============================================================

    w = ones(L,1);

    taperLength = max(2,round(0.15*L));

    n = (0:taperLength-1).';

    taper = ...
        0.5*(1+cos(pi*n/(taperLength-1)));

    w(end-taperLength+1:end) = taper;

    Ccausal = Ccausal .* w;

    % ==============================================================
    % STEP 5: Reconstruct frequency response on the actual 256-point
    % OFDM subcarriers.
    %
    % H(k) =
    %
    %   sum_{l=0}^{L-1}
    %       h(l) exp(-j 2*pi*k*l/N)
    % ==============================================================

    ell = 0:(L-1);

    A = ...
        exp( ...
        -1j*2*pi/N * ...
        (usedLocations(:)*ell));

    % ==============================================================
    % STEP 6: COMPLETE WDFT ESTIMATOR
    %
    % Ccausal is already a 64 x 100 mapping from LTS observations
    % to CIR taps.
    %
    % Therefore:
    %
    %       K = A * Ccausal
    %
    % gives:
    %
    %       (200 x 64) * (64 x 100)
    %
    %       = 200 x 100
    %
    % DO NOT multiply by S' here.
    % ==============================================================

    K = A*Ccausal;

    % ==============================================================
    % STEP 7: Flat-channel normalization
    %
    % For:
    %
    %       HLS = ones(100,1)
    %
    % we want approximately:
    %
    %       Hhat = ones(200,1)
    % ==============================================================

    testResponse = K*ones(Np,1);

    gain = mean(testResponse);

    if abs(gain) < 1e-12

        error( ...
            'WDFT normalization failed: calibration gain is zero.');

    end

    K = K/gain;

    % ==============================================================
    % STEP 8: FINAL DIMENSION CHECK
    % ==============================================================

    if size(K,1) ~= numel(usedLocations)

        error( ...
            'WDFT output dimension is incorrect.');

    end

    if size(K,2) ~= Np

        error( ...
            'WDFT input dimension is incorrect.');

    end

    if any(~isfinite(K(:)))

        error( ...
            'WDFT matrix contains NaN or Inf.');

    end

    fprintf('\n');
    fprintf('========== WDFT CHECK ==========\n');
    fprintf('LTS observations   : %d\n',Np);
    fprintf('DFT length         : %d\n',M);
    fprintf('CIR length         : %d\n',L);
    fprintf('Output subcarriers : %d\n',numel(usedLocations));
    fprintf('K dimensions       : %d x %d\n',size(K,1),size(K,2));
    fprintf('Matrix rank        : %d\n',rank(K));
    fprintf('Matrix norm        : %.6e\n',norm(K,'fro'));
    fprintf('================================\n');

end

%% ========================================================================
function R = local_uniform_covariance(N,L)

    k = (0:N-1).';
    d = k-k.';
    R = zeros(N,N);

    same = abs(d)<eps;
    R(same)=1;

    idx = ~same;
    z = 1j*2*pi*L.*d(idx)/N;
    R(idx)=(1-exp(-z))./z;

    R=(R+R')/2;
end

%% ========================================================================
function Rsub = local_subcov(R,rows,cols,N)

    ri = rows+N/2+1;
    ci = cols+N/2+1;
    Rsub=R(ri,ci);
end

%% ========================================================================
function Dp = local_svd_filter(lambda,r,gamma)

    n=numel(lambda);
    r=max(1,min(r,n));
    Dp=zeros(n,n);

    for k=1:r
        Dp(k,k)=1/(lambda(k)+1/gamma);
    end
end

%% ========================================================================
function mse=local_mse(Rhh,Rhp,Rpp,K,gamma,Nu)

    Np=size(Rpp,1);
    C=Rpp+(1/gamma)*eye(Np);

    E=Rhh-Rhp*K'-K*Rhp'+K*C*K';
    mse=max(real(trace(E))/Nu,eps);
end

%% ========================================================================
function cir=local_ls_to_cir(HLS,M)

    HLS=HLS(:);

    X=zeros(M,1);
    q=[-50:-1 1:50];
    rows=mod(q+M/2,M)+1;
    X(rows)=HLS;

    cir=fftshift(ifft(ifftshift(X),M));
    cir=cir(:);
end

%% ========================================================================
function Lhat=local_sht(cir,sigma2,T,G)

    threshold=T*sigma2;
    E=abs(cir(1:G)).^2;

    Lhat=0;
    for n=1:G
        if E(n)>threshold
            Lhat=n;
        end
    end

    Lhat=max(1,min(G,Lhat));
end

%% ========================================================================
function H=local_cfr(h,N)

    h=h(:);
    hp=[h;zeros(N-numel(h),1)];
    H=fftshift(fft(hp,N));
    H=H(:);
end
