%% ============================================================
% STEP 23B-3
% EXACT WINDOWED-DFT RECONSTRUCTION VERIFICATION
%
% Purpose:
%   Verify the complete WDFT estimator:
%
%       LTS
%        |
%        v
%       LS
%        |
%        v
%   101-point frequency vector
%        |
%        v
%      Hann
%        |
%        v
%   Paper Eq. (23)
%        |
%        v
%   centered CIR
%        |
%        v
%     fftshift
%        |
%        v
%     L = 64
%        |
%        v
%   zero-pad to N = 256
%        |
%        v
%     256-point FFT
%        |
%        v
%   200 used subcarriers
%        |
%        v
%       MSE
%
% Two tests are performed:
%
%   TEST A: noiseless
%   TEST B: SNR = 45 dB
%
% ============================================================

clear;
clc;

fprintf('\n');
fprintf('============================================================\n');
fprintf('STEP 23B-3 - EXACT WDFT RECONSTRUCTION VERIFICATION\n');
fprintf('============================================================\n\n');


%% ============================================================
% 1. PARAMETERS
% ============================================================

p = parameters();

rng(p.randomSeed);

N = p.N;          % 256
Nu = p.Nu;        % 200
M = 128;          % WDFT size
Np = 100;         % LTS tones
L = 64;           % paper Fig. 6 WDFT CIR length

SNRdB = 45;

windowType = 'hann';


fprintf('--- PARAMETERS ---\n');
fprintf('N       = %d\n',N);
fprintf('Nu      = %d\n',Nu);
fprintf('Np      = %d\n',Np);
fprintf('M       = %d\n',M);
fprintf('L       = %d\n',L);
fprintf('SNR     = %.1f dB\n',SNRdB);
fprintf('Window  = %s\n\n',windowType);


%% ============================================================
% 2. USED SUBCARRIERS
% ============================================================

usedLocations = p.activeIndex;

usedLocations = usedLocations(:).';

fprintf('Used subcarrier count = %d\n', ...
    length(usedLocations));

fprintf('Used range = %d to %d excluding DC\n\n', ...
    min(usedLocations), ...
    max(usedLocations));


%% ============================================================
% 3. PALL
% ============================================================

Pall = ieee80216_pall();


%% ============================================================
% 4. PREAMBLE SEQUENCES
% ============================================================

preambleFD = preamble_sequences(Pall);


%% ============================================================
% 5. PREAMBLE IFFT
% ============================================================

[stsTime,ltsTime] = ...
    preamble_to_ifft(preambleFD,p);


%% ============================================================
% 6. BUILD COMPLETE PREAMBLE
% ============================================================

preamble = ...
    build_80216_preamble(stsTime,ltsTime,p);

txPreamble = preamble.complete;

fprintf('Preamble length = %d samples\n\n', ...
    length(txPreamble));


%% ============================================================
% 7. SUI-2 CHANNEL
% ============================================================

h = sui2_channel(p);

h = h(:);

fprintf('============================================================\n');
fprintf('SUI-2 CHANNEL\n');
fprintf('============================================================\n');

fprintf('Channel length = %d\n',length(h));

fprintf('Channel taps:\n');
disp(h.');

fprintf('\n');


%% ============================================================
% 8. TRUE 256-POINT CHANNEL RESPONSE
% ============================================================
%
% The actual OFDM CFR is:
%
% H(k) =
%   sum_l h(l) exp(-j 2*pi*k*l/N)
%
% where:
%
%   N = 256
%
% and k is the signed OFDM subcarrier index.

delayIndex = (0:length(h)-1).';

H_true_used = zeros(Nu,1);

for ii = 1:Nu

    k = usedLocations(ii);

    H_true_used(ii) = ...
        sum( ...
        h .* ...
        exp(-1j*2*pi*k*delayIndex/N) ...
        );

end


fprintf('True CFR generated.\n');
fprintf('H_true_used = %d x %d\n\n', ...
    size(H_true_used,1), ...
    size(H_true_used,2));


%% ============================================================
% 9. CHANNEL OUTPUT
% ============================================================

rxChannel = ...
    apply_channel(txPreamble,h);


%% ============================================================
% 10. RUN TWO CONDITIONS
% ============================================================
%
% Test A:
%   noiseless
%
% Test B:
%   SNR = 45 dB

numTests = 2;

testNames = {
    'NOISELESS'
    '45 dB SNR'
    };

MSE_results = zeros(numTests,1);
NMSE_results = zeros(numTests,1);

%% ============================================================
% TEST LOOP
% ============================================================

for testID = 1:numTests

    fprintf('\n');
    fprintf('############################################################\n');
    fprintf('# TEST %d: %s\n',testID,testNames{testID});
    fprintf('############################################################\n\n');


    %% --------------------------------------------------------
    % RECEIVE SIGNAL
    % ---------------------------------------------------------

    if testID == 1

        % Completely noiseless
        rx = rxChannel;

        fprintf('Noise: NONE\n');

    else

        rx = add_awgn(rxChannel,SNRdB);

        fprintf('Noise: AWGN\n');
        fprintf('SNR  : %.1f dB\n',SNRdB);

    end


    %% --------------------------------------------------------
    % LTS EXTRACTION
    % ---------------------------------------------------------

    [rxLTS1,rxLTS2] = ...
        extract_lts(rx,p);


    %% --------------------------------------------------------
    % KNOWN LTS
    % ---------------------------------------------------------

    knownLTS = ...
        get_known_lts(ltsTime);


    %% --------------------------------------------------------
    % FFT RECEIVED LTS
    % ---------------------------------------------------------

    Y_LTS1 = ...
        lts_fft(rxLTS1);

    Y_LTS2 = ...
        lts_fft(rxLTS2);


    %% --------------------------------------------------------
    % FFT KNOWN LTS
    % ---------------------------------------------------------

    X_LTS = ...
        lts_fft(knownLTS);


    %% --------------------------------------------------------
    % COMBINE LTSs
    % ---------------------------------------------------------

    Y_LTS = ...
        combine_lts(Y_LTS1,Y_LTS2);


    %% --------------------------------------------------------
    % 100 LTS TONES
    % ---------------------------------------------------------

    [Y_p,X_p,ltsLocations,ltsBins] = ...
        extract_lts_100(Y_LTS,X_LTS);


    %% --------------------------------------------------------
    % LS ESTIMATE
    % ---------------------------------------------------------

    H_LS = ...
        ls_channel_estimate(Y_p,X_p);

    H_LS = H_LS(:);


    fprintf('\n100-tone LS estimate obtained.\n');


    %% ========================================================
    % 11. CREATE 101-POINT VECTOR
    % ========================================================
    %
    %       -100 ... -2 | 0 | +2 ... +100
    %
    % DC is unavailable from the LTS.
    %
    % Paper:
    %
    %       H(0) = [H(-2)+H(+2)]/2

    H_101 = zeros(101,1);

    H_101(1:50) = ...
        H_LS(1:50);

    H_101(51) = ...
        0.5*(H_LS(50)+H_LS(51));

    H_101(52:101) = ...
        H_LS(51:100);


    %% ========================================================
    % 12. HANN WINDOW
    % ========================================================

    n = (-50:50).';

    window = ...
        0.5 + 0.5*cos(2*pi*n/M);


    %% ========================================================
    % 13. APPLY WINDOW
    % ========================================================

    H_windowed = ...
        H_101 .* window;


    %% ========================================================
    % 14. PAPER Eq. (23)
    % ========================================================
    %
    % G(m,n) =
    %
    % 1/sqrt(M)
    % exp[
    %   j*2*pi*(m-M/2)*(n-Np/2)/M
    % ]
    %
    % m = 0,...,M-1
    % n = 0,...,Np

    m = (0:M-1).';

    nIndex = 0:Np;

    G = ...
        (1/sqrt(M)) .* ...
        exp( ...
        1j*2*pi/M .* ...
        ((m-M/2)*(nIndex-Np/2)) ...
        );


    %% ========================================================
    % 15. CENTERED WDFT CIR
    % ========================================================

    g_centered = ...
        G * H_windowed;


    %% ========================================================
    % 16. REORDER TO CAUSAL/CIRCULAR ORDER
    % ========================================================

    g_wdft = ...
        fftshift(g_centered);

    g_wdft = g_wdft(:);


    %% ========================================================
    % 17. ENERGY BEFORE TRUNCATION
    % ========================================================

    energyFull = ...
        sum(abs(g_wdft).^2);

    energyRetained = ...
        sum(abs(g_wdft(1:L)).^2);

    energyDiscarded = ...
        sum(abs(g_wdft(L+1:end)).^2);


    fprintf('\n');
    fprintf('WDFT CIR energy:\n');
    fprintf('Full       = %.12e\n',energyFull);
    fprintf('Retained   = %.12e\n',energyRetained);
    fprintf('Discarded  = %.12e\n',energyDiscarded);


    %% ========================================================
    % 18. RETAIN L=64 TAPS
    % ========================================================

    g_L = ...
        g_wdft(1:L);


    %% ========================================================
    % 19. EMBED 64-TAP CIR INTO 256-SAMPLE OFDM GRID
    % ========================================================
    %
    % IMPORTANT:
    %
    % The WDFT uses M=128.
    % The actual OFDM system uses N=256.
    %
    % We therefore place the retained physical CIR taps into
    % the first L samples of an N-point time-domain vector.
    %
    % The remaining N-L samples are zero.

    g_N = zeros(N,1);

    g_N(1:L) = ...
        g_L;


    %% ========================================================
    % 20. 256-POINT FFT
    % ========================================================
    %
    % Eq. (23) uses a UNITARY transform:
    %
    %       1/sqrt(M)
    %
    % Therefore the corresponding CFR reconstruction is:
    %
    %       H_hat =
    %       FFT_N(g_N) / sqrt(M)
    %
    % not simply FFT_N(g_N).

    H_hat_256 = ...
        fft(g_N,N) / sqrt(M);


    %% ========================================================
    % 21. EXTRACT 200 USED SUBCARRIERS
    % ========================================================
    %
    % MATLAB FFT ordering:
    %
    %   k = 0,...,N-1
    %
    % Signed OFDM index:
    %
    %   k = -128,...,-1,0,...,+127
    %
    % MATLAB index:
    %
    %   k + N/2 + 1

    usedBins = ...
        usedLocations + N/2 + 1;

    usedBins = ...
        round(usedBins);


    H_WDFT_used = ...
        H_hat_256(usedBins);


    H_WDFT_used = ...
        H_WDFT_used(:);


    %% ========================================================
    % 22. DIMENSION CHECK
    % ========================================================

    if length(H_WDFT_used) ~= Nu

        error( ...
            'WDFT reconstructed channel does not have 200 tones.');

    end


    %% ========================================================
    % 23. ERROR
    % ========================================================

    channelError = ...
        H_true_used - H_WDFT_used;


    %% ========================================================
    % 24. MSE
    % ========================================================

    MSE = ...
        mean(abs(channelError).^2);


    %% ========================================================
    % 25. NMSE
    % ========================================================

    truePower = ...
        mean(abs(H_true_used).^2);

    NMSE = ...
        MSE / truePower;


    %% ========================================================
    % 26. MAXIMUM ERROR
    % ========================================================

    maxError = ...
        max(abs(channelError));


    %% ========================================================
    % 27. ESTIMATED POWER
    % ========================================================

    estimatedPower = ...
        mean(abs(H_WDFT_used).^2);


    %% ========================================================
    % 28. SAVE RESULTS
    % ========================================================

    MSE_results(testID) = ...
        MSE;

    NMSE_results(testID) = ...
        NMSE;


    %% ========================================================
    % 29. PRINT RESULTS
    % ========================================================

    fprintf('\n');
    fprintf('------------------------------------------------------------\n');
    fprintf('RESULTS: %s\n',testNames{testID});
    fprintf('------------------------------------------------------------\n');

    fprintf('True channel power      = %.12e\n',truePower);
    fprintf('Estimated channel power = %.12e\n',estimatedPower);

    fprintf('WDFT MSE                = %.12e\n',MSE);
    fprintf('WDFT NMSE               = %.12e\n',NMSE);

    fprintf('Maximum abs error       = %.12e\n',maxError);


    %% ========================================================
    % 30. ERROR CHECKS
    % ========================================================

    if any(~isfinite(H_WDFT_used))

        fprintf('FAIL: NaN/Inf detected.\n');

    else

        fprintf('PASS: No NaN/Inf.\n');

    end


    if length(H_WDFT_used) == 200

        fprintf('PASS: 200 reconstructed used tones.\n');

    else

        fprintf('FAIL: Wrong number of reconstructed tones.\n');

    end


end


%% ============================================================
% 31. COMPARE NOISELESS VS 45 dB
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('TEST COMPARISON\n');
fprintf('============================================================\n');

fprintf('\n');

fprintf('%-20s %18s %18s\n', ...
    'Condition','MSE','NMSE');

fprintf('%-20s %18.12e %18.12e\n', ...
    'Noiseless', ...
    MSE_results(1), ...
    NMSE_results(1));

fprintf('%-20s %18.12e %18.12e\n', ...
    '45 dB', ...
    MSE_results(2), ...
    NMSE_results(2));


%% ============================================================
% 32. NOISE IMPACT
% ============================================================

fprintf('\n');

if MSE_results(1) > 0

    noiseIncrease = ...
        (MSE_results(2)-MSE_results(1)) / ...
        MSE_results(1) * 100;

    fprintf('MSE increase from noise = %.4f %%\n', ...
        noiseIncrease);

end


%% ============================================================
% 33. FINAL STRUCTURAL CHECKS
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('FINAL VALIDATION\n');
fprintf('============================================================\n');

pass = true;


% Check M
if M ~= 128

    fprintf('FAIL: M is not 128.\n');

    pass = false;

else

    fprintf('PASS: M = 128.\n');

end


% Check Np
if Np ~= 100

    fprintf('FAIL: Np is not 100.\n');

    pass = false;

else

    fprintf('PASS: Np = 100.\n');

end


% Check L
if L ~= 64

    fprintf('FAIL: L is not 64.\n');

    pass = false;

else

    fprintf('PASS: WDFT retained CIR length L = 64.\n');

end


% Check reconstructed dimension
if length(H_WDFT_used) ~= 200

    fprintf('FAIL: reconstructed CFR length is not 200.\n');

    pass = false;

else

    fprintf('PASS: reconstructed CFR has 200 used tones.\n');

end


% Check finite MSE
if ~isfinite(MSE_results(1)) || ...
   ~isfinite(MSE_results(2))

    fprintf('FAIL: MSE contains invalid value.\n');

    pass = false;

else

    fprintf('PASS: MSE values are finite.\n');

end


%% ============================================================
% 34. FINAL STATUS
% ============================================================

fprintf('\n');
fprintf('============================================================\n');

if pass

    fprintf('STEP 23B-3 STRUCTURAL STATUS: PASS\n');

else

    fprintf('STEP 23B-3 STRUCTURAL STATUS: NOT YET PASS\n');

end

fprintf('============================================================\n');