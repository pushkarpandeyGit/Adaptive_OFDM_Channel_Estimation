%% ============================================================
% STEP 23B-2
% WINDOWED DFT -> 200 USED-SUBCARRIER CHANNEL -> MSE
%
% Paper:
% A Low-Complexity Channel Estimation for OFDM Systems
% Based on CIR Length Adaptation
%
% Purpose:
%   Complete the WDFT baseline after Step 23B.
%
% Flow:
%
%   IEEE 802.16 preamble
%          |
%          v
%      SUI-2 channel
%          |
%          v
%         AWGN
%          |
%          v
%       received LTS
%          |
%          v
%     100 LTS tones
%          |
%          v
%       LS estimate
%          |
%          v
%       DC insertion
%          |
%          v
%       Hann window
%          |
%          v
%       WDFT Eq.(23)
%          |
%          v
%       fftshift
%          |
%          v
%      retain L = 64
%          |
%          v
%   reconstruct 200 used tones
%          |
%          v
%        WDFT MSE
%
% ============================================================

clear;
clc;

fprintf('\n');
fprintf('============================================================\n');
fprintf('STEP 23B-2 - WINDOWED DFT CHANNEL ESTIMATION MSE\n');
fprintf('============================================================\n\n');


%% ============================================================
% 1. PARAMETERS
% ============================================================

p = parameters();

rng(p.randomSeed);

SNRdB = 45;

M = 128;

% Number of LTS training tones
Np = 100;

% Paper Fig. 6 uses L = 64 for WDFT
L = 64;

windowType = 'hann';


fprintf('--- SIMULATION PARAMETERS ---\n');
fprintf('N                 = %d\n', p.N);
fprintf('Nu                = %d\n', p.Nu);
fprintf('Np                = %d\n', Np);
fprintf('M                 = %d\n', M);
fprintf('WDFT CIR length L = %d\n', L);
fprintf('SNR               = %.1f dB\n', SNRdB);
fprintf('Window            = %s\n\n', windowType);


%% ============================================================
% 2. IEEE 802.16 SUBCARRIER INFORMATION
% ============================================================

usedLocations = p.activeIndex;

usedLocations = usedLocations(:).';

fprintf('Used subcarrier range:\n');
fprintf('Minimum used k = %d\n', min(usedLocations));
fprintf('Maximum used k = %d\n', max(usedLocations));
fprintf('Number used    = %d\n\n', length(usedLocations));


%% ============================================================
% 3. GENERATE IEEE 802.16 PALL
% ============================================================

Pall = ieee80216_pall();


%% ============================================================
% 4. GENERATE PREAMBLE SEQUENCES
% ============================================================

preambleFD = preamble_sequences(Pall);


%% ============================================================
% 5. PREAMBLE TO TIME DOMAIN
% ============================================================

[stsTime,ltsTime] = preamble_to_ifft(preambleFD,p);


%% ============================================================
% 6. BUILD COMPLETE PREAMBLE
% ============================================================

preamble = build_80216_preamble(stsTime,ltsTime,p);

txPreamble = preamble.complete;

fprintf('Preamble length = %d samples\n\n', ...
    length(txPreamble));


%% ============================================================
% 7. GENERATE SUI-2 CHANNEL
% ============================================================

h = sui2_channel(p);

h = h(:);

fprintf('============================================================\n');
fprintf('SUI-2 CHANNEL\n');
fprintf('============================================================\n');

fprintf('Channel length = %d samples\n',length(h));

fprintf('Channel taps:\n');
disp(h.');

fprintf('\n');


%% ============================================================
% 8. TRUE CHANNEL FREQUENCY RESPONSE
% ============================================================
%
% We need the true H(k) to calculate MSE.
%
% For a discrete-time channel:
%
%       H(k) = sum_l h(l) exp(-j*2*pi*k*l/N)
%
% where:
%
%       k = OFDM subcarrier index
%       l = physical delay tap
%       N = 256
%
% MATLAB indexing:
%
%       l = 0,...,length(h)-1

delayIndex = (0:length(h)-1).';

H_true_used = zeros(length(usedLocations),1);

for kk = 1:length(usedLocations)

    k = usedLocations(kk);

    H_true_used(kk) = sum( ...
        h .* exp(-1j*2*pi*k*delayIndex/p.N) ...
        );

end


fprintf('True channel response generated.\n');
fprintf('H_true_used size = %d x %d\n\n', ...
    size(H_true_used,1), ...
    size(H_true_used,2));


%% ============================================================
% 9. APPLY CHANNEL
% ============================================================

rxChannel = apply_channel(txPreamble,h);


%% ============================================================
% 10. ADD AWGN
% ============================================================

rx = add_awgn(rxChannel,SNRdB);


%% ============================================================
% 11. EXTRACT TWO LTSs
% ============================================================

[rxLTS1,rxLTS2] = extract_lts(rx,p);


%% ============================================================
% 12. KNOWN LTS
% ============================================================

knownLTS = get_known_lts(ltsTime);


%% ============================================================
% 13. FFT RECEIVED LTSs
% ============================================================

Y_LTS1 = lts_fft(rxLTS1);
Y_LTS2 = lts_fft(rxLTS2);


%% ============================================================
% 14. FFT KNOWN LTS
% ============================================================

X_LTS = lts_fft(knownLTS);


%% ============================================================
% 15. COMBINE LTSs
% ============================================================

Y_LTS = combine_lts(Y_LTS1,Y_LTS2);


%% ============================================================
% 16. EXTRACT 100 LTS TRAINING TONES
% ============================================================

[Y_p,X_p,ltsLocations,ltsBins] = ...
    extract_lts_100(Y_LTS,X_LTS);


fprintf('============================================================\n');
fprintf('100-TONE LTS TRAINING DATA\n');
fprintf('============================================================\n');

fprintf('Number of tones = %d\n',length(Y_p));

fprintf('First location  = %d\n',ltsLocations(1));
fprintf('Last location   = %d\n\n',ltsLocations(end));


%% ============================================================
% 17. LS CHANNEL ESTIMATION
% ============================================================

H_LS = ls_channel_estimate(Y_p,X_p);

H_LS = H_LS(:);


fprintf('H_LS length = %d\n\n',length(H_LS));


%% ============================================================
% 18. CONSTRUCT 101-POINT FREQUENCY VECTOR
% ============================================================
%
% The 100 LTS tones are:
%
%   -100,-98,...,-2,+2,...,+98,+100
%
% There is no LTS at DC.
%
% The paper substitutes DC by:
%
%       H(0) = [H(-2)+H(+2)]/2
%
% Therefore:
%
%       101-point vector =
%
%       [-100 ... -2, 0, +2 ... +100]

H_101 = zeros(Np+1,1);

% Negative frequencies
H_101(1:50) = H_LS(1:50);

% DC
H_101(51) = ...
    0.5 * (H_LS(50) + H_LS(51));

% Positive frequencies
H_101(52:101) = H_LS(51:100);


fprintf('101-point frequency vector created.\n');
fprintf('DC value obtained by -2/+2 averaging.\n\n');


%% ============================================================
% 19. HANN WINDOW
% ============================================================

n = (-Np/2:Np/2).';

window = 0.5 + ...
         0.5*cos(2*pi*n/M);

fprintf('Window length = %d\n',length(window));


%% ============================================================
% 20. WINDOW THE LS ESTIMATE
% ============================================================

H_windowed = H_101 .* window;


%% ============================================================
% 21. PAPER Eq. (23)
% ============================================================
%
%       [Gtilde]_(m,n)
%
%       = 1/sqrt(M)
%         exp[
%           j*2*pi*(m-M/2)*(n-Np/2)/M
%         ]
%
% m = 0,...,M-1
% n = 0,...,Np

m = (0:M-1).';

nIndex = 0:Np;

G = (1/sqrt(M)) .* ...
    exp( ...
    1j * 2*pi/M .* ...
    ((m-M/2) * (nIndex-Np/2)) ...
    );


fprintf('\n');
fprintf('============================================================\n');
fprintf('PAPER Eq. (23)\n');
fprintf('============================================================\n');

fprintf('G dimensions = %d x %d\n', ...
    size(G,1),size(G,2));


%% ============================================================
% 22. WINDOWED DFT
% ============================================================

g_centered = G * H_windowed;


fprintf('Centered WDFT CIR = %d x %d\n', ...
    size(g_centered,1), ...
    size(g_centered,2));


%% ============================================================
% 23. REORDER CIR
% ============================================================

g_wdft = fftshift(g_centered);

g_wdft = g_wdft(:);


fprintf('Causal/circular WDFT CIR = %d x %d\n\n', ...
    size(g_wdft,1), ...
    size(g_wdft,2));


%% ============================================================
% 24. RETAIN FIRST L CIR TAPS
% ============================================================
%
% Paper Fig. 6 uses L = 64 for WDFT.
%
% Therefore:
%
%       g_L = g_wdft(1:L)
%
% The remaining taps are discarded.

g_L = g_wdft(1:L);


fprintf('============================================================\n');
fprintf('CIR TRUNCATION\n');
fprintf('============================================================\n');

fprintf('Original CIR length = %d\n',M);
fprintf('Retained CIR length = %d\n',L);

fprintf('Retained energy = %.12e\n', ...
    sum(abs(g_L).^2));

fprintf('Discarded energy = %.12e\n\n', ...
    sum(abs(g_wdft(L+1:end)).^2));


%% ============================================================
% 25. RECONSTRUCT 200 USED SUBCARRIERS
% ============================================================
%
% The unitary WDFT convention uses:
%
%       1/sqrt(M)
%
% in the transform.
%
% Therefore the corresponding synthesis is:
%
%       H_hat(k)
%
%       = 1/sqrt(M)
%         sum_{l=0}^{L-1}
%         g_L(l)
%         exp(-j*2*pi*k*l/M)
%
% IMPORTANT:
%
% Here M=128 is the WDFT interpolation size.
%
% The 200 OFDM subcarriers are the actual N=256
% OFDM frequency locations.
%
% The physical OFDM channel response uses N=256.
%
% Therefore we explicitly map the retained CIR into
% the 256-point OFDM frequency grid.

H_WDFT_used = zeros(length(usedLocations),1);

for kk = 1:length(usedLocations)

    k = usedLocations(kk);

    l = (0:L-1).';

    H_WDFT_used(kk) = ...
        (1/sqrt(M)) * ...
        sum( ...
        g_L .* ...
        exp(-1j*2*pi*k*l/p.N) ...
        );

end


fprintf('WDFT reconstructed channel length = %d\n\n', ...
    length(H_WDFT_used));


%% ============================================================
% 26. CHECK FINITE VALUES
% ============================================================

if any(~isfinite(H_WDFT_used))

    error('WDFT reconstructed channel contains NaN/Inf.');

end


%% ============================================================
% 27. CALCULATE WDFT MSE
% ============================================================
%
%       MSE =
%
%       1/Nu sum_k |H_true(k)-H_hat(k)|^2

error_WDFT = ...
    H_true_used - H_WDFT_used;

MSE_WDFT = ...
    mean(abs(error_WDFT).^2);

maxError_WDFT = ...
    max(abs(error_WDFT));


%% ============================================================
% 28. NORMALIZED MSE
% ============================================================

channelPower = ...
    mean(abs(H_true_used).^2);

NMSE_WDFT = ...
    MSE_WDFT / channelPower;


%% ============================================================
% 29. RESULTS
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('WINDOWED DFT CHANNEL ESTIMATION RESULTS\n');
fprintf('============================================================\n');

fprintf('SNR                  = %.1f dB\n',SNRdB);
fprintf('WDFT M                = %d\n',M);
fprintf('WDFT retained L       = %d\n',L);
fprintf('Window                = %s\n',windowType);
fprintf('Training tones        = %d\n',Np);
fprintf('Used subcarriers      = %d\n',length(usedLocations));

fprintf('\n');

fprintf('True channel power    = %.12e\n',channelPower);
fprintf('WDFT MSE              = %.12e\n',MSE_WDFT);
fprintf('WDFT NMSE             = %.12e\n',NMSE_WDFT);
fprintf('Maximum abs. error    = %.12e\n',maxError_WDFT);


%% ============================================================
% 30. COMPARE ESTIMATED AND TRUE CHANNEL POWER
% ============================================================

estimatedPower = ...
    mean(abs(H_WDFT_used).^2);

fprintf('\n');
fprintf('Estimated channel power = %.12e\n', ...
    estimatedPower);


%% ============================================================
% 31. BASIC VALIDATION
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('VALIDATION\n');
fprintf('============================================================\n');

pass = true;


% ------------------------------------------------------------
% Check 1: WDFT CIR
% ------------------------------------------------------------

if length(g_wdft) ~= M

    fprintf('FAIL: WDFT CIR length incorrect.\n');

    pass = false;

else

    fprintf('PASS: WDFT CIR length = %d.\n',M);

end


% ------------------------------------------------------------
% Check 2: reconstructed channel
% ------------------------------------------------------------

if length(H_WDFT_used) ~= p.Nu

    fprintf('FAIL: reconstructed channel length incorrect.\n');

    pass = false;

else

    fprintf('PASS: reconstructed channel has %d used tones.\n', ...
        p.Nu);

end


% ------------------------------------------------------------
% Check 3: finite
% ------------------------------------------------------------

if any(~isfinite(H_WDFT_used))

    fprintf('FAIL: WDFT estimate contains NaN/Inf.\n');

    pass = false;

else

    fprintf('PASS: WDFT estimate contains no NaN/Inf.\n');

end


% ------------------------------------------------------------
% Check 4: nonzero
% ------------------------------------------------------------

if MSE_WDFT <= 0

    fprintf('FAIL: MSE is zero/invalid.\n');

    pass = false;

else

    fprintf('PASS: WDFT MSE calculated.\n');

end


%% ============================================================
% 32. FINAL STATUS
% ============================================================

fprintf('\n');
fprintf('============================================================\n');

if pass

    fprintf('STEP 23B-2 STATUS: PASS\n');

else

    fprintf('STEP 23B-2 STATUS: NOT YET PASS\n');

end

fprintf('============================================================\n');