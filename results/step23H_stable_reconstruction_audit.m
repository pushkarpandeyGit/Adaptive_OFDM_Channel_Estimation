function results = step23H_stable_reconstruction_audit()

% ============================================================
% STEP 23-H : STABLE WDFT RECONSTRUCTION AUDIT
%
% Purpose:
%   Study the numerical stability of reconstructing a 64-tap
%   CIR from the FULL 128-output WDFT operator.
%
% Forward model:
%
%       h64 --> B --> g128
%
% where
%
%       B = Qshift * Gmat * Wmat * Dmat * Smat
%
% Reconstruction:
%
%       g128 --> B^+_r --> h64
%
% using TRUNCATED SVD (TSVD).
%
% IMPORTANT:
%   This script is an AUDIT.
%   It does not change the final WDFT estimator.
%   It does not select a final rank.
%   It does not use Wiener weights.
%
% ============================================================

clc;

fprintf('\n');
fprintf('============================================================\n');
fprintf(' STEP 23-H : STABLE WDFT RECONSTRUCTION AUDIT\n');
fprintf('============================================================\n\n');


%% ============================================================
% 1. PARAMETERS
% ============================================================

p = parameters();
rng(p.randomSeed);

N  = p.N;
M  = 128;
G  = 64;
Np = 100;

fprintf('N  = %d\n',N);
fprintf('M  = %d\n',M);
fprintf('G  = %d\n',G);
fprintf('Np = %d\n\n',Np);


%% ============================================================
% 2. LTS TONE LOCATIONS
% ============================================================

evenK = -100:2:100;
evenK_noDC = evenK(evenK ~= 0);

assert(length(evenK_noDC) == 100);

idx_m2 = find(evenK_noDC == -2);
idx_p2 = find(evenK_noDC == 2);

assert(~isempty(idx_m2));
assert(~isempty(idx_p2));


%% ============================================================
% 3. HANN WINDOW
% ============================================================

win = zeros(101,1);

for nn = -50:50

    ii = nn + 51;

    win(ii) = ...
        0.5 + 0.5*cos(2*pi*nn/M);

end


%% ============================================================
% 4. EQ. (23) MATRIX
% ============================================================

n = 0:Np;
m = (0:M-1).';

Gmat = (1/sqrt(M)) .* ...
    exp(1j*2*pi*(m-M/2).*(n-Np/2)/M);


%% ============================================================
% 5. CIR -> 100-TONE CFR MATRIX
% ============================================================

Smat = zeros(Np,G);

for r = 1:Np

    k = evenK_noDC(r);

    for q = 1:G

        tap = q-1;

        Smat(r,q) = ...
            exp(-1j*2*pi*k*tap/N);

    end

end


%% ============================================================
% 6. 100 -> 101 DC INSERTION
% ============================================================

Dmat = zeros(101,100);

Dmat(1:50,1:50) = eye(50);

Dmat(51,idx_m2) = 0.5;
Dmat(51,idx_p2) = 0.5;

Dmat(52:101,51:100) = eye(50);


%% ============================================================
% 7. WINDOW MATRIX
% ============================================================

Wmat = diag(win);


%% ============================================================
% 8. FFTSHIFT MATRIX
% ============================================================

Qshift = zeros(M,M);

I128 = eye(M);

for c = 1:M
    Qshift(:,c) = fftshift(I128(:,c));
end


%% ============================================================
% 9. FULL WDFT OPERATOR
% ============================================================

B = ...
    Qshift * ...
    Gmat * ...
    Wmat * ...
    Dmat * ...
    Smat;

fprintf('Full operator size = %d x %d\n', ...
    size(B,1),size(B,2));


%% ============================================================
% 10. SVD
% ============================================================

[U,S,V] = svd(B,'econ');

sigma = diag(S);

fprintf('\n');
fprintf('------------------------------------------------------------\n');
fprintf(' SINGULAR VALUES\n');
fprintf('------------------------------------------------------------\n');

for k = 1:G

    fprintf('%2d : %.12e\n', ...
        k,sigma(k));

end


%% ============================================================
% 11. TEST RANKS
% ============================================================
%
% We deliberately test several ranks.
%
% This lets us see the trade-off:
%
%   lower rank  -> more stable, less information
%   higher rank -> more information, more noise amplification
%
% ============================================================

testRanks = [ ...
    10 ...
    15 ...
    20 ...
    25 ...
    30 ...
    35 ...
    40 ...
    45 ...
    50 ...
    52 ...
    54 ...
    56 ...
    58 ...
    60 ...
    62 ...
    63 ...
    64 ];

numRanks = length(testRanks);


%% ============================================================
% 12. NOISE-FREE RANDOM CIR TEST
% ============================================================

numTests = 200;

mseNoiseFree = zeros(numRanks,1);
nmseNoiseFree = zeros(numRanks,1);

for rr = 1:numRanks

    r = testRanks(rr);

    Ur = U(:,1:r);
    Sr = S(1:r,1:r);
    Vr = V(:,1:r);

    % TSVD inverse:
    %
    % B_r^+ = V_r * inv(S_r) * U_r^H

    Bpinv_r = ...
        Vr * diag(1./diag(Sr)) * Ur';

    errors = zeros(numTests,1);

    for trial = 1:numTests

        h = ...
            randn(G,1) + 1j*randn(G,1);

        g = B*h;

        hHat = Bpinv_r*g;

        errors(trial) = ...
            norm(hHat-h)^2 / ...
            norm(h)^2;

    end

    mseNoiseFree(rr) = mean(errors);

    nmseNoiseFree(rr) = mean(errors);

end


%% ============================================================
% 13. NOISE-FREE SUI-2 TEST
% ============================================================

numSUI2 = 200;

mseSUI2 = zeros(numRanks,1);
nmseSUI2 = zeros(numRanks,1);

for rr = 1:numRanks

    r = testRanks(rr);

    Ur = U(:,1:r);
    Sr = S(1:r,1:r);
    Vr = V(:,1:r);

    Bpinv_r = ...
        Vr * diag(1./diag(Sr)) * Ur';

    mseTrial = zeros(numSUI2,1);
    nmseTrial = zeros(numSUI2,1);

    for trial = 1:numSUI2

        h = sui2_channel(p);
        h = h(:);

        h64 = zeros(G,1);

        Lh = min(length(h),G);

        h64(1:Lh) = h(1:Lh);

        g = B*h64;

        hHat = Bpinv_r*g;

        err = hHat-h64;

        mseTrial(trial) = ...
            mean(abs(err).^2);

        denom = mean(abs(h64).^2);

        if denom > 0

            nmseTrial(trial) = ...
                mseTrial(trial)/denom;

        else

            nmseTrial(trial) = NaN;

        end

    end

    mseSUI2(rr) = mean(mseTrial);

    nmseSUI2(rr) = mean(nmseTrial);

end


%% ============================================================
% 14. AWGN RECONSTRUCTION TEST
% ============================================================
%
% We now test the important practical issue:
%
%       noise amplification
%
% Several SNR values are tested.
%
% ============================================================

SNRlist = [0 5 10 15 20 25 30 35 40];

numSNR = length(SNRlist);

mseAWGN = zeros(numRanks,numSNR);
nmseAWGN = zeros(numRanks,numSNR);


for ss = 1:numSNR

    SNRdB = SNRlist(ss);

    fprintf('\nTesting SNR = %d dB\n',SNRdB);

    for rr = 1:numRanks

        r = testRanks(rr);

        Ur = U(:,1:r);
        Sr = S(1:r,1:r);
        Vr = V(:,1:r);

        Bpinv_r = ...
            Vr * diag(1./diag(Sr)) * Ur';

        mseTrial = zeros(numTests,1);
        nmseTrial = zeros(numTests,1);

        for trial = 1:numTests

            h = ...
                randn(G,1) + 1j*randn(G,1);

            gClean = B*h;

            signalPower = ...
                mean(abs(gClean).^2);

            noisePower = ...
                signalPower / ...
                (10^(SNRdB/10));

            noise = ...
                sqrt(noisePower/2) .* ...
                (randn(M,1) + ...
                 1j*randn(M,1));

            gNoisy = gClean + noise;

            hHat = Bpinv_r*gNoisy;

            err = hHat-h;

            mseTrial(trial) = ...
                mean(abs(err).^2);

            denom = mean(abs(h).^2);

            nmseTrial(trial) = ...
                mseTrial(trial)/denom;

        end

        mseAWGN(rr,ss) = mean(mseTrial);

        nmseAWGN(rr,ss) = mean(nmseTrial);

    end

end


%% ============================================================
% 15. NOISE AMPLIFICATION FACTOR
% ============================================================

%
% For the TSVD inverse:
%
%       1/sigma_i
%
% is the amplification applied to noise in singular direction i.
%
% ============================================================

noiseAmplification = zeros(numRanks,1);

for rr = 1:numRanks

    r = testRanks(rr);

    noiseAmplification(rr) = ...
        1/(sigma(r)^2);

end


%% ============================================================
% 16. PRINT MAIN RESULTS
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' TSVD RANK COMPARISON\n');
fprintf('============================================================\n');

fprintf('\n');
fprintf('%6s %18s %18s %18s\n', ...
    'Rank', ...
    'Sigma(rank)', ...
    'NF NMSE', ...
    'SUI2 NMSE');

fprintf('%s\n', ...
    repmat('-',1,68));

for rr = 1:numRanks

    fprintf('%6d %18.6e %18.6e %18.6e\n', ...
        testRanks(rr), ...
        sigma(testRanks(rr)), ...
        nmseNoiseFree(rr), ...
        nmseSUI2(rr));

end


%% ============================================================
% 17. AWGN TABLE
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' AWGN NMSE BY RANK\n');
fprintf('============================================================\n');

fprintf('%6s', 'Rank');

for ss = 1:numSNR
    fprintf('%14s',sprintf('%ddB',SNRlist(ss)));
end

fprintf('\n');

fprintf('%s\n', ...
    repmat('-',1,6+14*numSNR));

for rr = 1:numRanks

    fprintf('%6d',testRanks(rr));

    for ss = 1:numSNR

        fprintf('%14.4e', ...
            nmseAWGN(rr,ss));

    end

    fprintf('\n');

end


%% ============================================================
% 18. CONDITIONING / AMPLIFICATION TABLE
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' SVD STABILITY\n');
fprintf('============================================================\n');

fprintf('%6s %18s %22s\n', ...
    'Rank', ...
    'Sigma(rank)', ...
    '1/Sigma(rank)^2');

fprintf('%s\n', ...
    repmat('-',1,52));

for rr = 1:numRanks

    fprintf('%6d %18.6e %22.6e\n', ...
        testRanks(rr), ...
        sigma(testRanks(rr)), ...
        noiseAmplification(rr));

end


%% ============================================================
% 19. SAVE RESULTS
% ============================================================

results.B = B;

results.U = U;
results.S = S;
results.V = V;

results.singularValues = sigma;

results.testRanks = testRanks;

results.noiseFreeNMSE = nmseNoiseFree;
results.noiseFreeMSE = mseNoiseFree;

results.sui2NMSE = nmseSUI2;
results.sui2MSE = mseSUI2;

results.SNRlist = SNRlist;

results.awgnNMSE = nmseAWGN;
results.awgnMSE = mseAWGN;

results.noiseAmplification = ...
    noiseAmplification;

results.conditionNumber = ...
    sigma(1)/sigma(end);

save( ...
    'step23H_stable_reconstruction_audit.mat', ...
    'results');


%% ============================================================
% 20. FINAL AUDIT MESSAGE
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' STEP 23-H COMPLETE\n');
fprintf('============================================================\n');

fprintf('\n');
fprintf(['This experiment does NOT select the final WDFT rank.\n' ...
         'It only measures the reconstruction/stability trade-off.\n']);

fprintf('\n');
fprintf(['Use the AWGN NMSE table to determine where increasing\n' ...
         'the SVD rank starts causing serious noise amplification.\n']);

fprintf('\n');
fprintf('Saved:\n');
fprintf('step23H_stable_reconstruction_audit.mat\n');

fprintf('============================================================\n');

end