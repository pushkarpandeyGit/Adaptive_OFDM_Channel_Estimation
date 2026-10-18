%% ============================================================
% STEP 24 - UNIFIED OFDM CHANNEL ESTIMATOR COMPARISON
%
% Estimators:
%   1. LS
%   2. LMMSE
%   3. Fixed-Rank SVD
%   4. Adaptive-Rank SVD
%   5. WDFT
%
% Primary metric:
%   CIR NMSE
%
% IMPORTANT:
%   LS CIR is obtained using the validated:
%       ls_to_cir(HLS,128)
%
%   SVD is performed directly in the physical CFR observation
%   domain, but with controlled truncation.
%
% ============================================================

clear;
clc;
close all;

%% ============================================================
% 1. PARAMETERS
% ============================================================

p = parameters();
rng(p.randomSeed);

N  = p.N;                  % 256
G  = p.G;                  % 64
M  = p.M;                  % QAM order
numFrames = p.numFrames;

SNRdB = p.SNRdB;

fprintf('\n============================================================\n');
fprintf(' STEP 24 - UNIFIED ESTIMATOR COMPARISON\n');
fprintf('============================================================\n');

fprintf('N              = %d\n', N);
fprintf('CIR comparison = %d taps\n', G);
fprintf('SNR range      = %d to %d dB\n', ...
    SNRdB(1), SNRdB(end));
fprintf('Monte Carlo    = %d frames/SNR\n', numFrames);

%% ============================================================
% 2. PHYSICAL LTS SUBCARRIER LOCATIONS
% ============================================================

% IEEE 802.16 LTS tones:
% k = -100,-98,...,-2,2,...,98,100
%
% DC is not transmitted.

kLTS = -100:2:100;
kLTS = kLTS(kLTS ~= 0);

assert(length(kLTS) == 100, ...
    'Expected exactly 100 LTS tones.');

%% ============================================================
% 3. PHYSICAL CFR -> CIR MODEL MATRIX FOR SVD
% ============================================================

% Model:
%
% H(k) = sum_{l=0}^{G-1} h(l) exp(-j2*pi*k*l/N)
%
% Therefore:
%
% F(r,q) = exp(-j2*pi*k_r*(q-1)/N)

F = zeros(100,G);

for r = 1:100

    k = kLTS(r);

    for q = 1:G

        tap = q - 1;

        F(r,q) = exp( ...
            -1j*2*pi*k*tap/N );

    end

end

%% ============================================================
% 4. SVD OF PHYSICAL OBSERVATION MATRIX
% ============================================================

fprintf('\nComputing SVD of physical CFR matrix...\n');

[U,Sigma,V] = svd(F,'econ');

singularValues = diag(Sigma);

fprintf('F size        : %d x %d\n',size(F,1),size(F,2));
fprintf('Rank(F)       : %d\n',rank(F));
fprintf('Max singular  : %.6e\n',max(singularValues));
fprintf('Min singular  : %.6e\n',min(singularValues));
fprintf('Condition no. : %.6e\n', ...
    max(singularValues)/min(singularValues));

%% ============================================================
% 5. FIXED SVD RANK
% ============================================================

% Do NOT use full rank 64.
%
% Full-rank inversion activates extremely small singular-value
% directions and causes severe noise amplification.
%
% Rank 54 is used as the initial fixed-rank benchmark because
% the previous SVD stability experiment showed that ranks above
% approximately 54 begin entering the unstable region.

fixedRank = 54;

fixedRank = max(1,min(G,fixedRank));

fprintf('\nFixed SVD rank = %d\n',fixedRank);

%% ============================================================
% 6. SUI-2 CHANNEL COVARIANCE FOR LMMSE
% ============================================================

fprintf('\nBuilding SUI-2 covariance matrix...\n');

% Generate many stochastic SUI-2 CIRs to estimate R_hh.

numCovTrials = 3000;

Htraining = zeros(G,numCovTrials);

for ii = 1:numCovTrials

    htmp = sui2_channel(p);

    htmp = htmp(:);

    if length(htmp) < G
        htmp = [htmp; zeros(G-length(htmp),1)];
    else
        htmp = htmp(1:G);
    end

    Htraining(:,ii) = htmp;

end

Rhh = (Htraining * Htraining') / numCovTrials;

% Numerical Hermitian enforcement
Rhh = (Rhh + Rhh')/2;

fprintf('Rhh size = %d x %d\n',size(Rhh));

%% ============================================================
% WDFT WEIGHTS
% Use the exact weighting from the validated Step 23-J audit.
% ============================================================

wdftFile = 'step23J_actual_wdft_audit.mat';

if ~isfile(wdftFile)
    error('Cannot find step23J_actual_wdft_audit.mat');
end

S23 = load(wdftFile);

if ~isfield(S23,'step23J')
    error('Variable step23J not found.');
end

if ~isfield(S23.step23J,'theta')
    error('Field step23J.theta not found.');
end

wdftWeights = S23.step23J.theta(:);

if length(wdftWeights) ~= G
    error('WDFT theta must contain exactly %d elements.',G);
end

if any(~isfinite(wdftWeights))
    error('WDFT theta contains NaN or Inf.');
end

fprintf('\nWDFT weights loaded from validated Step 23-J:\n');
fprintf('  min       = %.12e\n',min(wdftWeights));
fprintf('  max       = %.12e\n',max(wdftWeights));
fprintf('  mean      = %.12e\n',mean(wdftWeights));
fprintf('  nonzero   = %d / %d\n', ...
    nnz(wdftWeights),length(wdftWeights));

%% ============================================================
% 8. WDFT MATRIX FROM PAPER EQUATION (23)
% ============================================================

% Paper Eq. (23):
%
% [Gtilde]_(m,n) =
% 1/sqrt(M) *
% exp(j*2*pi*(m-M/2)*(n-Np/2)/M)
%
% Here:
%
% M     = 128
% Np+1  = 101
%
% m = 0,...,127
% n = 0,...,100

Mwdft = 128;

Gwdft = zeros(Mwdft,101);

for m = 0:Mwdft-1

    for n = 0:100

        Gwdft(m+1,n+1) = ...
            (1/sqrt(Mwdft)) * ...
            exp(1j*2*pi* ...
            (m-Mwdft/2) * ...
            (n-100/2) / Mwdft);

    end

end

%% ============================================================
% 9. HANN WINDOW
% ============================================================

nWin = -50:50;

hannWindow = 0.5 + ...
    0.5*cos(2*pi*nWin/Mwdft);

hannWindow = hannWindow(:);

assert(length(hannWindow)==101);

%% ============================================================
% 10. STORAGE
% ============================================================

numSNR = length(SNRdB);

nmseLS        = zeros(1,numSNR);
nmseLMMSE     = zeros(1,numSNR);
nmseFixedSVD  = zeros(1,numSNR);
nmseAdaptive  = zeros(1,numSNR);
nmseWDFT      = zeros(1,numSNR);

meanLhat      = zeros(1,numSNR);
meanAdaptiveRank = zeros(1,numSNR);

%% ============================================================
% 11. MAIN MONTE CARLO LOOP
% ============================================================

for snrIndex = 1:numSNR

    currentSNR = SNRdB(snrIndex);

    fprintf('\n');
    fprintf('------------------------------------------------------------\n');
    fprintf(' SNR = %d dB\n',currentSNR);
    fprintf('------------------------------------------------------------\n');

    errorLS       = zeros(numFrames,1);
    errorLMMSE    = zeros(numFrames,1);
    errorFixed    = zeros(numFrames,1);
    errorAdaptive = zeros(numFrames,1);
    errorWDFT     = zeros(numFrames,1);

    LhatStore = zeros(numFrames,1);
    rankStore = zeros(numFrames,1);

    for frame = 1:numFrames

        %% ====================================================
        % A. GENERATE TRUE SUI-2 CHANNEL
        % ====================================================

        hTrue = sui2_channel(p);

        hTrue = hTrue(:);

        if length(hTrue) < G

            hTrue = [ ...
                hTrue;
                zeros(G-length(hTrue),1)];

        else

            hTrue = hTrue(1:G);

        end

        %% ====================================================
        % B. TRUE 256-POINT CFR
        % ====================================================

        h256 = [ ...
            hTrue;
            zeros(N-G,1)];

        Htrue256 = fft(h256,N);

        %% ====================================================
        % C. EXTRACT 100 PHYSICAL LTS TONES
        % ====================================================

        idx256 = kLTS + N/2 + 1;

        Htrue100 = Htrue256(idx256);

        Htrue100 = Htrue100(:);

        %% ====================================================
        % D. GENERATE CFR-DOMAIN NOISE
        % ====================================================

        signalPower = mean(abs(Htrue100).^2);

        noisePower = ...
            signalPower / ...
            10^(currentSNR/10);

        noise = sqrt(noisePower/2) * ...
            (randn(100,1) + ...
            1j*randn(100,1));

        Hnoisy100 = Htrue100 + noise;

        %% ====================================================
        % E. LS CFR ESTIMATE
        %
        % For the ideal known LTS:
        %
        % H_LS = Y/X
        %
        % Here Hnoisy100 represents that LS CFR observation.
        % ====================================================

        HLS = Hnoisy100;

        HLS = HLS(:);

        %% ====================================================
        % F. ESTIMATE NOISE VARIANCE
        %
        % The validated LS->CIR transformation gives a
        % 128-sample CIR.
        %
        % Noise region:
        %       taps 65:120
        %
        % sigma_n^2 =
        %       mean(|g[n]|^2)
        % ====================================================

        cir128 = ls_to_cir(HLS,128);

        cir128 = cir128(:);

        if length(cir128) < 120

            error('CIR output must contain at least 120 samples.');

        end

        estimatedNoiseVariance = ...
            mean(abs(cir128(65:120)).^2);

        %% ====================================================
        % G. SHT
        %
        % threshold = T * noise variance
        %
        % T = 15
        % G = 64
        % ====================================================

        threshold = ...
            p.T * estimatedNoiseVariance;

        tapPower = abs(cir128(1:G)).^2;

        aboveThreshold = ...
            find(tapPower > threshold);

        if isempty(aboveThreshold)

            Lhat = 1;

        else

            Lhat = max(aboveThreshold);

        end

        Lhat = max(1,min(G,Lhat));

        %% ====================================================
        % H. ADAPTIVE SVD RANK
        % ====================================================

        adaptiveRank = ...
            Lhat + p.adaptiveSVDOffset;

        adaptiveRank = ...
            max(1,min(G,adaptiveRank));

        LhatStore(frame) = Lhat;
        rankStore(frame) = adaptiveRank;

        %% ====================================================
        % I. LS CIR ESTIMATE
        %
        % IMPORTANT:
        %
        % DO NOT use:
        %
        %       pinv(F)*HLS
        %
        % because F is extremely ill-conditioned.
        %
        % Use the validated transformation instead.
        % ====================================================

        hLS = cir128(1:G);

        %% ====================================================
        % J. LMMSE ESTIMATE
        %
        % h_LMMSE =
        %
        % Rhh F^H
        % ------------------------- HLS
        % F Rhh F^H + Rnn
        %
        % ====================================================

        % Noise variance in LS-CFR domain.
        %
        % HLS = Htrue + n
        %
        % Therefore use the actual observation noise variance.

        Rnn = noisePower * eye(100);

        A_lmmse = ...
            F * Rhh * F' + Rnn;

        % Avoid explicit matrix inverse.
        %
        % h = Rhh F^H A^{-1} HLS

        hLMMSE = ...
            Rhh * F' * ...
            (A_lmmse \ HLS);

        hLMMSE = hLMMSE(:);

        %% ============================================================
% K. FIXED-RANK SVD
% ============================================================

% Full physical CIR model:
%
% HLS = F*h + n
%
% F is 100 x 64.
%
% Fixed SVD uses a predetermined rank.

rFixed = fixedRank;

hFixed = zeros(G,1);

for ii = 1:rFixed

    sigma_i = singularValues(ii);

    coefficient = ...
        (U(:,ii)' * HLS) / sigma_i;

    hFixed = ...
        hFixed + ...
        coefficient * V(:,ii);

end


%% ============================================================
% L. ADAPTIVE CIR-LENGTH SVD
% ============================================================

% IMPORTANT:
%
% Lhat is an estimated PHYSICAL CIR LENGTH.
%
% Therefore we should NOT interpret:
%
%     Lhat + offset
%
% as the number of singular vectors of the original
% 100 x 64 matrix.
%
% Instead, construct a reduced channel model using Lhat
% physical CIR taps.

Lmodel = Lhat + p.adaptiveSVDOffset;

Lmodel = max(1,min(G,Lmodel));


% ------------------------------------------------------------
% Construct reduced physical CFR matrix
% ------------------------------------------------------------

Fadaptive = F(:,1:Lmodel);


% ------------------------------------------------------------
% SVD of reduced physical channel model
% ------------------------------------------------------------

[Ua,Sa,Va] = svd(Fadaptive,'econ');

singularAdaptive = diag(Sa);


% ------------------------------------------------------------
% Adaptive SVD rank
% ------------------------------------------------------------

% The reduced model has Lmodel physical dimensions.
%
% We don't automatically have to retain every dimension.
%
% For this first clean implementation, retain the complete
% reduced physical model.
%
% The important adaptive quantity is therefore Lmodel itself.

rAdaptive = min(Lmodel,length(singularAdaptive));


% ------------------------------------------------------------
% Reconstruct adaptive channel
% ------------------------------------------------------------

hAdaptiveReduced = zeros(Lmodel,1);

for ii = 1:rAdaptive

    sigma_i = singularAdaptive(ii);

    coefficient = ...
        (Ua(:,ii)' * HLS) / sigma_i;

    hAdaptiveReduced = ...
        hAdaptiveReduced + ...
        coefficient * Va(:,ii);

end


% ------------------------------------------------------------
% Put adaptive result back into 64-tap CIR vector
% ------------------------------------------------------------

hAdaptive = zeros(G,1);

hAdaptive(1:Lmodel) = hAdaptiveReduced;

        %% ====================================================
        % M. WDFT ESTIMATE
        %
        % 100 LTS tones
        %      |
        %      v
        % Insert DC
        %      |
        %      v
        % 101-point vector
        %      |
        %      v
        % Hann window
        %      |
        %      v
        % 14 front zeros + 13 back zeros
        %      |
        %      v
        % Eq. (23) transform
        %      |
        %      v
        % fftshift
        %      |
        %      v
        % first 64 CIR samples
        %      |
        %      v
        % Wiener weighting
        %
        % ====================================================

        % DC substitution:
        %
        % Use average of nearest physical tones:
        % k = -2 and +2

        dcEstimate = ...
            0.5 * (HLS(50) + HLS(51));

        H101 = [ ...
            HLS(1:50);
            dcEstimate;
            HLS(51:100)];

        H101 = H101(:);

        % Hann window
        Hwindowed = ...
            H101 .* hannWindow;

        % WDFT transform
        gCentered = ...
            Gwdft * Hwindowed;

        % Reorder so physical CIR begins at tap 1
        gFull = fftshift(gCentered);

        gWDFT = ...
            gFull(1:G);

        % Wiener weighting
        hWDFT = ...
            wdftWeights .* gWDFT;

        hWDFT = hWDFT(:);

        %% ====================================================
        % N. CALCULATE CIR NMSE
        %
        % NMSE =
        %
        % ||h_hat-h||^2
        % -----------------
        % ||h||^2
        %
        % ====================================================

        trueEnergy = ...
            sum(abs(hTrue).^2);

        if trueEnergy < eps

            error('True CIR energy is too small.');

        end

        errorLS(frame) = ...
            sum(abs(hLS-hTrue).^2) / trueEnergy;

        errorLMMSE(frame) = ...
            sum(abs(hLMMSE-hTrue).^2) / trueEnergy;

        errorFixed(frame) = ...
            sum(abs(hFixed-hTrue).^2) / trueEnergy;

        errorAdaptive(frame) = ...
            sum(abs(hAdaptive-hTrue).^2) / trueEnergy;

        errorWDFT(frame) = ...
            sum(abs(hWDFT-hTrue).^2) / trueEnergy;

        %% ====================================================
        % O. FIRST-FRAME SANITY CHECK
        % ====================================================

        if frame == 1

            fprintf('\n--- FIRST FRAME SANITY CHECK ---\n');

            fprintf('True CIR energy       : %.6e\n', ...
                sum(abs(hTrue).^2));

            fprintf('LS CIR energy         : %.6e\n', ...
                sum(abs(hLS).^2));

            fprintf('LMMSE CIR energy      : %.6e\n', ...
                sum(abs(hLMMSE).^2));

            fprintf('Fixed-SVD CIR energy  : %.6e\n', ...
                sum(abs(hFixed).^2));

            fprintf('Adaptive-SVD energy   : %.6e\n', ...
                sum(abs(hAdaptive).^2));

            fprintf('WDFT CIR energy       : %.6e\n', ...
                sum(abs(hWDFT).^2));

            fprintf('Estimated noise      : %.6e\n', ...
                estimatedNoiseVariance);

            fprintf('SHT threshold        : %.6e\n', ...
                threshold);

            fprintf('Lhat                 : %d\n', ...
                Lhat);

            fprintf('Adaptive rank        : %d\n', ...
                adaptiveRank);

            fprintf('LS NMSE              : %.6e\n', ...
                errorLS(frame));

            fprintf('LMMSE NMSE           : %.6e\n', ...
                errorLMMSE(frame));

            fprintf('Fixed-SVD NMSE       : %.6e\n', ...
                errorFixed(frame));

            fprintf('Adaptive-SVD NMSE    : %.6e\n', ...
                errorAdaptive(frame));

            fprintf('WDFT NMSE            : %.6e\n', ...
                errorWDFT(frame));

            fprintf('--------------------------------\n\n');

        end

    end

    %% ========================================================
    % P. AVERAGE MONTE CARLO NMSE
    % ========================================================

    nmseLS(snrIndex) = ...
        mean(errorLS);

    nmseLMMSE(snrIndex) = ...
        mean(errorLMMSE);

    nmseFixedSVD(snrIndex) = ...
        mean(errorFixed);

    nmseAdaptive(snrIndex) = ...
        mean(errorAdaptive);

    nmseWDFT(snrIndex) = ...
        mean(errorWDFT);

    meanLhat(snrIndex) = ...
        mean(LhatStore);

    meanAdaptiveRank(snrIndex) = ...
        mean(rankStore);

    %% ========================================================
    % Q. SNR SUMMARY
    % ========================================================

    fprintf('\nSNR %2d dB summary:\n',currentSNR);

    fprintf('LS             : %.6e\n', ...
        nmseLS(snrIndex));

    fprintf('LMMSE          : %.6e\n', ...
        nmseLMMSE(snrIndex));

    fprintf('Fixed-SVD      : %.6e\n', ...
        nmseFixedSVD(snrIndex));

    fprintf('Adaptive-SVD   : %.6e\n', ...
        nmseAdaptive(snrIndex));

    fprintf('WDFT           : %.6e\n', ...
        nmseWDFT(snrIndex));

    fprintf('Mean Lhat      : %.4f\n', ...
        meanLhat(snrIndex));

    fprintf('Mean rank      : %.4f\n', ...
        meanAdaptiveRank(snrIndex));

end

%% ============================================================
% 12. FINAL RESULTS TABLE
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' STEP 24 FINAL CIR NMSE RESULTS\n');
fprintf('============================================================\n\n');

fprintf('%8s %15s %15s %15s %15s %15s\n', ...
    'SNR','LS','LMMSE','Fixed-SVD', ...
    'Adaptive-SVD','WDFT');

fprintf('%s\n',repmat('-',1,95));

for ii = 1:numSNR

    fprintf('%8d %15.6e %15.6e %15.6e %15.6e %15.6e\n', ...
        SNRdB(ii), ...
        nmseLS(ii), ...
        nmseLMMSE(ii), ...
        nmseFixedSVD(ii), ...
        nmseAdaptive(ii), ...
        nmseWDFT(ii));

end

%% ============================================================
% 13. ADAPTIVE SVD TABLE
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' ADAPTIVE SVD PARAMETERS\n');
fprintf('============================================================\n\n');

fprintf('%8s %15s %15s\n', ...
    'SNR','Mean Lhat','Mean Rank');

fprintf('%s\n',repmat('-',1,45));

for ii = 1:numSNR

    fprintf('%8d %15.4f %15.4f\n', ...
        SNRdB(ii), ...
        meanLhat(ii), ...
        meanAdaptiveRank(ii));

end

%% ============================================================
% 14. VALIDITY CHECKS
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' STEP 24 VALIDITY CHECKS\n');
fprintf('============================================================\n');

% Check for catastrophic numerical values.

allResults = [ ...
    nmseLS(:);
    nmseLMMSE(:);
    nmseFixedSVD(:);
    nmseAdaptive(:);
    nmseWDFT(:)];

if any(~isfinite(allResults))

    error('FAIL: Non-finite NMSE detected.');

end

% Catastrophic explosion check.
%
% A correctly behaving estimator should not produce the previous
% 1e25-1e48 numerical explosion.

if max(allResults) > 1e12

    warning(['Some estimator still shows extremely large NMSE. ', ...
        'Inspect the first-frame sanity check and SVD ranks.']);

else

    fprintf('Numerical magnitude check : PASS\n');

end

% Lhat sanity check.
%
% SUI-2 physical CIR has approximately 5-sample support.

if mean(meanLhat) >= 2 && mean(meanLhat) <= 15

    fprintf('Adaptive Lhat sanity      : PASS\n');

else

    warning(['Adaptive Lhat is outside expected SUI-2 range. ', ...
        'Do not proceed to ML yet.']);

end

% Rank sanity
if all(meanAdaptiveRank >= 1 & meanAdaptiveRank <= G)

    fprintf('Adaptive rank bounds      : PASS\n');

else

    error('FAIL: Adaptive rank outside valid range.');

end

%% ============================================================
% 15. PLOT
% ============================================================

figure;

semilogy(SNRdB,nmseLS,'-o','LineWidth',1.5);
hold on;

semilogy(SNRdB,nmseLMMSE,'-s','LineWidth',1.5);

semilogy(SNRdB,nmseFixedSVD,'-^','LineWidth',1.5);

semilogy(SNRdB,nmseAdaptive,'-d','LineWidth',1.5);

semilogy(SNRdB,nmseWDFT,'-x','LineWidth',1.5);

grid on;

xlabel('SNR (dB)');
ylabel('CIR NMSE');

title('OFDM Channel Estimator Comparison');

legend( ...
    'LS', ...
    'LMMSE', ...
    'Fixed-SVD', ...
    'Adaptive-SVD', ...
    'WDFT', ...
    'Location','best');

%% ============================================================
% 16. ADAPTIVE RANK PLOT
% ============================================================

figure;

plot(SNRdB,meanLhat,'-o','LineWidth',1.5);
hold on;

plot(SNRdB,meanAdaptiveRank,'-s','LineWidth',1.5);

grid on;

xlabel('SNR (dB)');
ylabel('Average CIR length / SVD rank');

title('Adaptive CIR Length and SVD Rank');

legend( ...
    'Mean Lhat', ...
    'Mean Adaptive Rank', ...
    'Location','best');

%% ============================================================
% 17. SAVE RESULTS
% ============================================================

save( ...
    'step24_estimator_comparison.mat', ...
    'SNRdB', ...
    'nmseLS', ...
    'nmseLMMSE', ...
    'nmseFixedSVD', ...
    'nmseAdaptive', ...
    'nmseWDFT', ...
    'meanLhat', ...
    'meanAdaptiveRank', ...
    'fixedRank', ...
    'singularValues', ...
    'Rhh');

%% ============================================================
% 18. COMPLETE
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' STEP 24 COMPLETE\n');
fprintf('============================================================\n');

fprintf('Primary metric: CIR NMSE\n');

fprintf('Results saved to:\n');
fprintf('step24_estimator_comparison.mat\n');

fprintf('============================================================\n');