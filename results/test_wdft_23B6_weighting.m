clear;
clc;

p = parameters();

rng(p.randomSeed);

%% =============================================================
% WDFT 23B-6
% Time-domain scalar Wiener weighting diagnostic
% =============================================================

N  = p.N;
M  = 128;
Np = 100;
L  = 64;

numTrials = 2000;

SNRdB = 45;

fprintf('\n');
fprintf('============================================\n');
fprintf(' WDFT 23B-6 WIENER WEIGHTING\n');
fprintf('============================================\n');

fprintf('Trials       = %d\n',numTrials);
fprintf('SNR          = %.1f dB\n',SNRdB);
fprintf('M            = %d\n',M);
fprintf('L            = %d\n',L);
fprintf('LTS tones    = %d\n',Np);

%% -------------------------------------------------------------
% LTS subcarrier locations
% --------------------------------------------------------------

evenK = -100:2:100;
evenK(evenK == 0) = [];

k101 = -100:2:100;

%% -------------------------------------------------------------
% Build DC substitution matrix
% --------------------------------------------------------------

D = zeros(101,100);

idx_m2 = find(evenK == -2);
idx_p2 = find(evenK ==  2);

for ii = 1:length(k101)

    kk = k101(ii);

    if kk == 0

        D(ii,idx_m2) = 0.5;
        D(ii,idx_p2) = 0.5;

    else

        idx = find(evenK == kk,1);

        D(ii,idx) = 1;

    end
end

%% -------------------------------------------------------------
% Hann window
% --------------------------------------------------------------

n = (-Np/2):(Np/2);

w = 0.5 + 0.5*cos(2*pi*n/M);

w = w(:);

W = diag(w);

%% -------------------------------------------------------------
% Eq. (23) matrix
% --------------------------------------------------------------

G = zeros(M,Np+1);

for m = 0:M-1

    for nn = 0:Np

        G(m+1,nn+1) = ...
            (1/sqrt(M)) * ...
            exp(1j*2*pi*(m-M/2)*(nn-Np/2)/M);

    end
end

%% -------------------------------------------------------------
% Keep causal 64 taps
% --------------------------------------------------------------

G_L = G(M/2+1:M,:);

%% -------------------------------------------------------------
% Complete transform from 100 LS tones to 64 effective taps
% --------------------------------------------------------------

B = G_L * W * D;

fprintf('\nB matrix size = %d x %d\n', ...
    size(B,1),size(B,2));

%% -------------------------------------------------------------
% Noise variance
% --------------------------------------------------------------

% Normalize according to the signal power of the pilot
% observation.

noiseVar = 10^(-SNRdB/10);

%% -------------------------------------------------------------
% Accumulators for Wiener coefficients
% --------------------------------------------------------------

crossAccum = zeros(L,1);
powerAccum = zeros(L,1);

%% -------------------------------------------------------------
% Monte-Carlo ensemble
% --------------------------------------------------------------

for trial = 1:numTrials

    %% Generate independent SUI-2 channel

    h = sui2_channel(p);

    h = h(:);

    %% ---------------------------------------------------------
    % True 256-point CFR
    % ----------------------------------------------------------

    hN = [h; zeros(N-length(h),1)];

    H = fft(hN,N);

    %% ---------------------------------------------------------
    % True LTS CFR
    % ----------------------------------------------------------

    idxLTS = evenK + N/2 + 1;

    pTrue = H(idxLTS);

    pTrue = pTrue(:);

    %% ---------------------------------------------------------
    % Ideal effective CIR
    % ----------------------------------------------------------

    gTrue = B * pTrue;

    %% ---------------------------------------------------------
    % Add complex AWGN to LS observations
    % ----------------------------------------------------------

    noise = sqrt(noiseVar/2) * ...
        (randn(Np,1) + 1j*randn(Np,1));

    pLS = pTrue + noise;

    %% ---------------------------------------------------------
    % Noisy effective CIR
    % ----------------------------------------------------------

    gHat = B * pLS;

    %% ---------------------------------------------------------
    % Accumulate statistics
    % ----------------------------------------------------------

    crossAccum = crossAccum + ...
        gTrue .* conj(gHat);

    powerAccum = powerAccum + ...
        abs(gHat).^2;

end

%% -------------------------------------------------------------
% Ensemble averages
% --------------------------------------------------------------

Rcross = crossAccum / numTrials;

Rhat = powerAccum / numTrials;

%% -------------------------------------------------------------
% Scalar Wiener weights
% --------------------------------------------------------------

wWiener = real(Rcross ./ Rhat);

%% -------------------------------------------------------------
% Safety check
% --------------------------------------------------------------

wWiener(~isfinite(wWiener)) = 0;

%% -------------------------------------------------------------
% Display weights
% --------------------------------------------------------------

fprintf('\n============================================\n');
fprintf(' WIENER WEIGHTS\n');
fprintf('============================================\n');

fprintf('Tap       Weight\n');
fprintf('----------------------\n');

for l = 1:L

    fprintf('%3d    %.10f\n', ...
        l,wWiener(l));

end

%% -------------------------------------------------------------
% Statistics
% --------------------------------------------------------------

fprintf('\n============================================\n');
fprintf(' WEIGHT STATISTICS\n');
fprintf('============================================\n');

fprintf('Minimum weight = %.10f\n',min(wWiener));
fprintf('Maximum weight = %.10f\n',max(wWiener));
fprintf('Mean weight    = %.10f\n',mean(wWiener));

%% -------------------------------------------------------------
% Effective CIR energy
% --------------------------------------------------------------

fprintf('\n============================================\n');
fprintf(' EFFECTIVE CIR ENERGY\n');
fprintf('============================================\n');

fprintf('Tap 1-10 energy = %.12e\n', ...
    sum(mean(abs(gHat(1:min(10,L))).^2)));

fprintf('Total energy not directly available from last trial.\n');

%% -------------------------------------------------------------
% Save diagnostic result
% --------------------------------------------------------------
% % 
save('wdft_23B6_wiener_weights.mat', ...
    'wWiener', ...
    'Rcross', ...
    'Rhat', ...
    'B', ...
    'M', ...
    'L', ...
    'Np', ...
    'SNRdB', ...
    'numTrials');

fprintf('\nSaved:\n');
fprintf('wdft_23B6_wiener_weights.mat\n');

fprintf('\n============================================\n');
fprintf(' 23B-6 COMPLETE\n');
fprintf('============================================\n');