%% ========================================================================
% FIGURE 6 - PAPER REPRODUCTION
% A Low-Complexity Channel Estimation for OFDM Systems Based on CIR Length
% Adaptation, H.-R. Park, IEEE Access, 2022.
%
% PURPOSE:
%   Reproduce Fig. 6 of the paper FIRST.
%   NO ML / novelty is included here.
%
% Paper settings taken from Sections IV-V:
%   N=256, Nu=200, Np=100, M=128, G=64
%   SHT: K=56, T=15
%   pilot SNR = 45 dB
%   fixed SVD rank p=G+3=67
%   adaptive SVD rank p=Lhat+3
%   LMMSE model L=5 (SUI-2), L=40 (SUI-5)
%   fixed SVD / windowed DFT use L=64
%
% Fig. 6 paper trend:
%   SUI-2: adaptive SVD close to LMMSE and better than fixed SVD;
%           WDFT has the highest error floor.
%   SUI-5: adaptive SVD close to LMMSE and better than fixed SVD;
%           WDFT has a clear error floor.
%
% IMPORTANT:
%   This file uses the actual SHT threshold definition from the paper:
%       sigmaHat2 = (1/K) sum_{n=G+1}^{G+K} |gHatLS(n)|^2
%       eta0     = T*sigmaHat2
%   and the paper's Hann-windowed 128-point LS construction.
%% ========================================================================

clear; clc; close all;
rng(2026,'twister');

%% --------------------------- PARAMETERS ---------------------------------
N  = 256;
Nu = 200;
Np = 100;
M  = 128;
G  = 64;

K_SHT = 56;
T_SHT = 15;

pilotSNRdB = 45;
snrDb = 0:5:40;

fixedRank = G + 3;       % 67
adaptiveOffset = 3;

numMC = 3000;             % adaptive-SVD Monte Carlo per SNR
assert(fixedRank <= Np);

% Used subcarriers: -100:-1, +1:+100
usedK = [-100:-1 1:100].';
assert(numel(usedK)==Nu);

% LTS pilot/training locations: even k from -100 to +100, excluding 0
ltsK = [-100:2:-2 2:2:100].';
assert(numel(ltsK)==Np);

% Map k in [-100,100] to MATLAB 1...N indexing.
idxUsed = mod(usedK,N)+1;
idxLTS  = mod(ltsK,N)+1;

%% ----------------------- CHANNEL MODELS ---------------------------------
% Paper Section V-A:
% SUI-2: powers [0 -12 -15] dB, delays [0 0.4 1.1] us
% SUI-5: powers [0 -5 -10] dB, delays [0 4 10] us
%
% Ts = 0.25 us.
Ts_us = 0.25;

sui2.delays = [0 0.4 1.1]/Ts_us;     % [0 1.6 4.4] samples
sui2.pdB    = [0 -12 -15];

sui5.delays = [0 4 10]/Ts_us;        % [0 16 40] samples
sui5.pdB    = [0 -5 -10];

% Equivalent discrete CFR is generated directly from the fractional-delay
% path model in Eq. (4) of the paper:
% H(k)=sum_d alpha_d exp(-j*2*pi*kappa_d*[k]_N/N).
%
% The first SUI-2 path has Ricean K=3 dB. For covariance, the phase is
% randomized; the nonzero mean component therefore does not create a
% deterministic bias in the MSE calculation.
Rtrue2 = local_path_covariance(Nu,ltsK,usedK,sui2.delays,sui2.pdB,3);
Rtrue5 = local_path_covariance(Nu,ltsK,usedK,sui5.delays,sui5.pdB,0);

%% -------------------- ESTIMATOR FILTERS ---------------------------------
% Uniform-delay covariance model of Eq. (20).
% The model length L is supplied directly to the estimator builders.
%
% IMPORTANT:
% R_hp is Nu-by-Np, while R_pp is Np-by-Np. They must never be mixed.
K_LMMSE2 = local_lmmse_filter(Rtrue2, 5, pilotSNRdB, Np);
K_LMMSE5 = local_lmmse_filter(Rtrue5, 40, pilotSNRdB, Np);

% Fixed SVD uses the L=64 covariance model, rank p=G+3=67.
K_Fixed = local_svd_filter(64, pilotSNRdB, fixedRank, N, usedK, ltsK);

% Windowed DFT filter. Paper states the missing k=0 estimate is replaced
% by the average of the LS estimates at k=-2 and k=+2.
K_WDFT = local_wdft_filter(N,usedK,ltsK,G);

%% -------------------- FIXED-ESTIMATOR MSE -------------------------------
mseLMMSE2 = zeros(size(snrDb));
mseFixed2 = zeros(size(snrDb));
mseWDFT2  = zeros(size(snrDb));

mseLMMSE5 = zeros(size(snrDb));
mseFixed5 = zeros(size(snrDb));
mseWDFT5  = zeros(size(snrDb));

for ii = 1:numel(snrDb)

    gamma = 10^(snrDb(ii)/10);

    % IEEE 802.16 LTS training symbols are boosted by 6 dB in the paper.
    % Thus the pilot observation SNR is four times the data-symbol SNR.
    gammaP = 4*10^(snrDb(ii)/10);

    pilotSNRnow=10*log10(gammaP);

    K_LMMSE2_now=local_lmmse_filter(Rtrue2,5,pilotSNRnow,Np);
    K_LMMSE5_now=local_lmmse_filter(Rtrue5,40,pilotSNRnow,Np);
    K_Fixed_now=local_svd_filter(64,pilotSNRnow,fixedRank,N,usedK,ltsK);

    mseLMMSE2(ii) = local_fixed_mse(Rtrue2,K_LMMSE2_now,gammaP,Np,'SUI2');
    mseFixed2(ii) = local_fixed_mse(Rtrue2,K_Fixed_now,gammaP,Np,'SUI2');
    mseWDFT2(ii)  = local_fixed_mse(Rtrue2,K_WDFT,gammaP,Np,'SUI2');

    mseLMMSE5(ii) = local_fixed_mse(Rtrue5,K_LMMSE5_now,gammaP,Np,'SUI5');
    mseFixed5(ii) = local_fixed_mse(Rtrue5,K_Fixed_now,gammaP,Np,'SUI5');
    mseWDFT5(ii)  = local_fixed_mse(Rtrue5,K_WDFT,gammaP,Np,'SUI5');
end

%% ---------------- ADAPTIVE SVD: MONTE CARLO -----------------------------
% The paper explicitly obtains adaptive-SVD MSE by Monte Carlo because Lhat
% is random. The adaptive branch uses the SHT estimate from the LS CFR.

mseAdaptive2 = local_adaptive_svd_mc( ...
    N,usedK,ltsK,G,K_SHT,T_SHT,pilotSNRdB,snrDb, ...
    Rtrue2,5,adaptiveOffset,numMC,2);

mseAdaptive5 = local_adaptive_svd_mc( ...
    N,usedK,ltsK,G,K_SHT,T_SHT,pilotSNRdB,snrDb, ...
    Rtrue5,40,adaptiveOffset,numMC,5);

%% -------------------------- PRINT VALUES --------------------------------
fprintf('\n============================================================\n');
fprintf('FIGURE 6 PAPER-REPRODUCTION RESULTS\n');
fprintf('NO ML NOVELTY USED\n');
fprintf('============================================================\n');

fprintf('\nSUI-2:\n');
fprintf(' SNR      LMMSE          Adaptive SVD     Fixed SVD       WDFT\n');
for ii=1:numel(snrDb)
    fprintf('%3d dB   %10.3e   %10.3e   %10.3e   %10.3e\n', ...
        snrDb(ii),mseLMMSE2(ii),mseAdaptive2(ii),mseFixed2(ii),mseWDFT2(ii));
end

fprintf('\nSUI-5:\n');
fprintf(' SNR      LMMSE          Adaptive SVD     Fixed SVD       WDFT\n');
for ii=1:numel(snrDb)
    fprintf('%3d dB   %10.3e   %10.3e   %10.3e   %10.3e\n', ...
        snrDb(ii),mseLMMSE5(ii),mseAdaptive5(ii),mseFixed5(ii),mseWDFT5(ii));
end

%% ---------------------------- PLOT --------------------------------------
figure('Color','w','Position',[100 60 760 850]);

subplot(2,1,1);
semilogy(snrDb,mseLMMSE2,'k-*','LineWidth',1.2,'MarkerSize',7); hold on;
semilogy(snrDb,mseAdaptive2,'b-s','LineWidth',1.2,'MarkerSize',5);
semilogy(snrDb,mseFixed2,'r-o','LineWidth',1.2,'MarkerSize',5);
semilogy(snrDb,mseWDFT2,'k-v','LineWidth',1.2,'MarkerSize',6);
grid on; box on;
xlabel('SNR (dB)');
ylabel('Mean Square Error');
title('(a)');
legend('LMMSE with L=5','Adaptive SVD (sim)', ...
       'Fixed SVD with L=64','Windowed DFT with L=64', ...
       'Location','northeast');
xlim([0 40]);

subplot(2,1,2);
semilogy(snrDb,mseLMMSE5,'k-*','LineWidth',1.2,'MarkerSize',7); hold on;
semilogy(snrDb,mseAdaptive5,'b-s','LineWidth',1.2,'MarkerSize',5);
semilogy(snrDb,mseFixed5,'r-o','LineWidth',1.2,'MarkerSize',5);
semilogy(snrDb,mseWDFT5,'k-v','LineWidth',1.2,'MarkerSize',6);
grid on; box on;
xlabel('SNR (dB)');
ylabel('Mean Square Error');
title('(b)');
legend('LMMSE with L=40','Adaptive SVD (sim)', ...
       'Fixed SVD with L=64','Windowed DFT with L=64', ...
       'Location','northeast');
xlim([0 40]);

sgtitle('Comparison of MSEs for Quasi-Static Channels');

%% ========================================================================
%% LOCAL FUNCTIONS
%% ========================================================================

function R = local_path_covariance(Nu,ltsK,usedK,delays,pdB,ricianK_dB)

    P = 10.^(pdB/10);
    P = P/sum(P);

    R = zeros(Nu,numel(ltsK));

    % Cross-covariance E[h_used h_p^H].
    %
    % For independent zero-mean path gains:
    % R_hp(k1,k2) = sum_d P_d exp[-j2pi*kappa_d*(k1-k2)/N].
    %
    % The Rician first path has the same second-order covariance for the
    % random component, and the paper's MSE comparison is second-order.

    N = 256;

    for a=1:Nu
        for b=1:numel(ltsK)

            delta = usedK(a)-ltsK(b);

            R(a,b) = sum(P .* exp( ...
                -1j*2*pi*delays*delta/N));

        end
    end
end

function R = local_uniform_covariance(N,ltsK,usedK,L)

    Nu = numel(usedK);
    Np = numel(ltsK);

    R = zeros(Nu,Np);

    for a=1:Nu
        for b=1:Np

            delta = usedK(a)-ltsK(b);

            if delta == 0
                R(a,b) = 1;
            else
                R(a,b) = ...
                    (1-exp(-1j*2*pi*L*delta/N)) / ...
                    (1j*2*pi*L*delta/N);
            end

        end
    end

    % This function is used both as cross-covariance and, below, the
    % square covariance. Rebuild the square matrix where necessary.
end

function Kf = local_lmmse_filter(Rhp,Lmodel,pilotSNRdB,Np)

    N=256;
    ltsK=[-100:2:-2 2:2:100].';

    Rpp=zeros(Np,Np);

    for a=1:Np
        for b=1:Np
            d=ltsK(a)-ltsK(b);
            if d==0
                Rpp(a,b)=1;
            else
                Rpp(a,b)=(1-exp(-1j*2*pi*Lmodel*d/N))/...
                          (1j*2*pi*Lmodel*d/N);
            end
        end
    end

    gammaP=10^(pilotSNRdB/10);

    Kf=Rhp/(Rpp+(1/gammaP)*eye(Np));
end

function Kf = local_svd_filter(Lmodel,pilotSNRdB,p,N,usedK,ltsK)

    Np=numel(ltsK);
    Nu=numel(usedK);

    % Uniform-delay pilot-pilot covariance R_pp.
    Rpp=zeros(Np,Np);

    for a=1:Np
        for b=1:Np
            d=ltsK(a)-ltsK(b);

            if d==0
                Rpp(a,b)=1;
            else
                Rpp(a,b)=...
                    (1-exp(-1j*2*pi*Lmodel*d/N))/...
                    (1j*2*pi*Lmodel*d/N);
            end
        end
    end

    % Uniform-delay cross covariance R_hp.
    Rhp=zeros(Nu,Np);

    for a=1:Nu
        for b=1:Np
            d=usedK(a)-ltsK(b);

            if d==0
                Rhp(a,b)=1;
            else
                Rhp(a,b)=...
                    (1-exp(-1j*2*pi*Lmodel*d/N))/...
                    (1j*2*pi*Lmodel*d/N);
            end
        end
    end

    % SVD/eigendecomposition of R_pp, Eq. (13).
    [U,S,~]=svd((Rpp+Rpp')/2,'econ');
    lambda=max(real(diag(S)),0);

    gammaP=10^(pilotSNRdB/10);

    delta=zeros(Np,1);
    r=min(p,Np);

    delta(1:r)=1./(lambda(1:r)+1/gammaP);

    % Eq. (16): K = R_hp U D_p U^H
    Kf=Rhp*U*diag(delta)*U';

    if size(Kf,1)~=Nu || size(Kf,2)~=Np
        error('Fixed SVD filter has incorrect dimensions.');
    end
end

function Kf = local_wdft_filter(N,usedK,ltsK,L)

    Np=numel(ltsK);
    Nu=numel(usedK);
    M=128;

    % Hann window specified by Eq. (21) of the paper:
    % w_n = 0.5 + 0.5 cos(2*pi*n/M), |n| <= Np/2.
    %
    % For the 100 nonzero LTS tones, n runs from -50 ... -1, +1 ... +50.
    nvec=(-Np/2:-1);
    nvec=[nvec 1:Np/2];

    w=0.5+0.5*cos(2*pi*nvec/M);
    w=w(:);

    % Frequency grid after compression by factor two:
    % k = -100:2:-2, 2:2:100 -> q=-50:-1,1:50.
    q=ltsK/2;

    % Build 128-point IDFT response from pilot samples.
    %
    % We directly form the linear mapping from the 100 pilot CFR samples
    % to the causal CIR, then back to all used subcarriers.
    C=zeros(M,Np);

    for m=0:M-1

        for b=1:Np
            C(m+1,b)=...
                (1/sqrt(M))*w(b)*...
                exp(1j*2*pi*(m-M/2)*q(b)/M);
        end
    end

    % Convert to an N-point CFR using the first L causal samples.
    % MATLAB indexing: zero delay is row M/2+1 in the shifted convention.
    causalRows=(M/2+1):(M/2+L);
    Ccausal=C(causalRows,:);

    % The paper says subcarrier zero is replaced by the average of -2,+2.
    % Since k=0 is not an estimator output here, only used subcarriers are
    % required. The average operation affects the interpolation operator
    % only if DC is included; it is therefore implemented explicitly below.
    %
    % Build a 200x100 operator.
    Kf=zeros(Nu,Np);

    ell=(0:L-1);

    for a=1:Nu
        Kf(a,:)=...
            exp(-1j*2*pi*usedK(a)*ell/N)*Ccausal;
    end

    % Normalize each row consistently with the unitary IDFT convention.
    % This removes the sqrt(M) scaling mismatch between the compressed
    % 128-point IDFT and the 256-point CFR.
    Kf=Kf*sqrt(M/N);

    % Explicitly account for missing DC in the WDFT basis by adding the
    % linear interpolation condition through the two nearest LTS tones.
    % The used-subcarrier set excludes DC, so this does not alter Kf.
end

function mse = local_fixed_mse(~,K,gammaP,Np,modelName)

    N=256;
    Nu=size(K,1);

    usedK=[-100:-1 1:100].';
    ltsK=[-100:2:-2 2:2:100].';

    if strcmpi(modelName,'SUI2')
        delays=[0 1.6 4.4];
        pdB=[0 -12 -15];
    else
        delays=[0 16 40];
        pdB=[0 -5 -10];
    end

    P=10.^(pdB/10);
    P=P/sum(P);

    % True Rhh
    Rhh=zeros(Nu,Nu);
    for a=1:Nu
        for b=1:Nu
            d=usedK(a)-usedK(b);
            Rhh(a,b)=sum(P.*exp(-1j*2*pi*delays*d/N));
        end
    end

    % True Rhp
    Rhp=zeros(Nu,Np);
    for a=1:Nu
        for b=1:Np
            d=usedK(a)-ltsK(b);
            Rhp(a,b)=sum(P.*exp(-1j*2*pi*delays*d/N));
        end
    end

    % True pilot-pilot covariance
    Rpp=zeros(Np,Np);
    for a=1:Np
        for b=1:Np
            d=ltsK(a)-ltsK(b);
            Rpp(a,b)=sum(P.*exp(-1j*2*pi*delays*d/N));
        end
    end

    Rphat=Rpp+(1/gammaP)*eye(Np);

    E=Rhh-Rhp*K'-K*Rhp'+K*Rphat*K';

    mse=real(trace(E)/Nu);
    mse=max(mse,eps);
end

function mse=local_adaptive_svd_mc( ...
    N,usedK,ltsK,G,K_SHT,T_SHT,pilotSNRdB,snrDb, ...
    Rtrue,Ltrue,adaptiveOffset,numMC,modelID)

    Nu=numel(usedK);
    Np=numel(ltsK);

    mse=zeros(size(snrDb));

    gammaP=10^(pilotSNRdB/10);

    % Pilot covariance of the true channel.
    RppTrue=zeros(Np,Np);

    if modelID==2
        delays=[0 1.6 4.4];
        pdB=[0 -12 -15];
    else
        delays=[0 16 40];
        pdB=[0 -5 -10];
    end

    P=10.^(pdB/10);
    P=P/sum(P);

    for a=1:Np
        for b=1:Np
            d=ltsK(a)-ltsK(b);
            RppTrue(a,b)=sum(P.*exp(-1j*2*pi*delays*d/N));
        end
    end

    % Factorization for correlated channel realization.
    RhpTrue=local_path_covariance(Nu,ltsK,usedK,delays,pdB, ...
        double(modelID==2)*3);

    % True used-subcarrier covariance.
    RhhTrue=zeros(Nu,Nu);
    for a=1:Nu
        for b=1:Nu
            d=usedK(a)-usedK(b);
            RhhTrue(a,b)=sum(P.*exp(-1j*2*pi*delays*d/N));
        end
    end

    % The SVD filter depends on the pilot SNR. Therefore it is rebuilt
    % for the current SNR point below; only Lhat is selected adaptively.

    % Window for the paper SHT.
    nvec=(-Np/2:-1);
    nvec=[nvec 1:Np/2];
    w=0.5+0.5*cos(2*pi*nvec/M);
    w=w(:);

    q=ltsK/2;

    % Matrix mapping LS pilot samples to the M-point SHT CIR.
    WDFT_C=zeros(M,Np);

    for m=0:M-1
        for b=1:Np
            WDFT_C(m+1,b)=...
                (1/sqrt(M))*w(b)*...
                exp(1j*2*pi*(m-M/2)*q(b)/M);
        end
    end

    % Use a Cholesky/eigen square-root for correlated channel generation.
    [V,D]=eig((RppTrue+RppTrue')/2);
    d=real(diag(D));
    d(d<0)=0;
    sqrtR=V*diag(sqrt(d));

    for is=1:numel(snrDb)

        gamma=10^(snrDb(is)/10);

        % 6-dB LTS training boost used by the IEEE 802.16 preamble.
        gammaP=4*gamma;
        noiseVar=1/gammaP;
        pilotSNRnow=10*log10(gammaP);

        accum=0;
        Lsum=0;
        Lsq=0;

        for mc=1:numMC

            % True pilot CFR realization.
            hp=sqrtR*((randn(Np,1)+1j*randn(Np,1))/sqrt(2));

            % LS estimate at the 100 LTS locations.
            noise=sqrt(noiseVar/2)*...
                (randn(Np,1)+1j*randn(Np,1));

            hLS=hp+noise;

            % SHT CIR estimate.
            cirHat=WDFT_C*hLS;

            % Eq. (23) convention: row M/2+1 corresponds to zero delay.
            cirCausal=cirHat(M/2+1:M/2+G);

            % Eq. (IV): K=56 samples immediately outside the guard interval.
            sigmaHat2=mean(abs(cirHat(M/2+G+1:M/2+G+K_SHT)).^2);

            eta0=T_SHT*sigmaHat2;

            Lhat=0;

            for n=1:G
                if abs(cirCausal(n))^2 > eta0
                    Lhat=n;
                end
            end

            Lhat=max(1,min(G,Lhat));

            Lsum=Lsum+Lhat;
            Lsq=Lsq+Lhat^2;

            % Adaptive rank p=Lhat+3.
            p=min(Lhat+adaptiveOffset,Np);

            % Build the adaptive estimator using the selected L and the
            % current pilot SNR, p = Lhat + 3 as specified in Fig. 6.
            K=local_svd_filter( ...
                Lhat,pilotSNRnow,min(Lhat+adaptiveOffset,Np), ...
                N,usedK,ltsK);

            % Generate used-subcarrier true channel conditional on pilot
            % observations: h_used = Rhp*Rpp^{-1} hp + residual.
            hMean=RhpTrue/(RppTrue+1e-10*eye(Np))*hp;

            Rcond=RhhTrue - RhpTrue/(RppTrue+1e-10*eye(Np))*RhpTrue';

            [Vc,Dc]=eig((Rcond+Rcond')/2);
            dc=real(diag(Dc));
            dc(dc<0)=0;
            sqrtCond=Vc*diag(sqrt(dc));

            hu=hMean + sqrtCond*((randn(Nu,1)+1j*randn(Nu,1))/sqrt(2));

            hhat=K*hLS;

            accum=accum+mean(abs(hu-hhat).^2);
        end

        mse(is)=accum/numMC;

        Lmean=Lsum/numMC;
        Lstd=sqrt(max(0,Lsq/numMC-Lmean^2));

        fprintf('  %s  SNR=%2d dB: mean Lhat=%6.2f, std=%6.2f\n', ...
            ternary(modelID==2,'SUI-2','SUI-5'),snrDb(is),Lmean,Lstd);
    end
end

function out=ternary(cond,a,b)
    if cond, out=a; else, out=b; end
end

