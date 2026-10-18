function results = step23I_regularized_wdft_audit()

% ============================================================
% STEP 23-I : REGULARIZED WDFT RECONSTRUCTION AUDIT
%
% Uses the verified WDFT forward operator B from Step 23-G.
%
% IMPORTANT:
%   This script computes the SVD of B itself.
%   It does NOT assume U/V were saved by Step 23-G.
%
% Forward model:
%
%       g = B*h
%
% Reconstruction:
%
%       TSVD
%       Tikhonov
%
% This is an AUDIT only.
% ============================================================

clc;

fprintf('\n');
fprintf('============================================================\n');
fprintf(' STEP 23-I : REGULARIZED WDFT RECONSTRUCTION AUDIT\n');
fprintf('============================================================\n\n');


%% ============================================================
% 1. LOAD STEP 23-G
% ============================================================

auditFile = 'step23G_basis_audit_results.mat';

if ~exist(auditFile,'file')

    error(['Missing ',auditFile, ...
        '. Run Step 23-G first.']);

end

S = load(auditFile);

if ~isfield(S,'results')
    error('The Step 23-G MAT file does not contain "results".');
end

if ~isfield(S.results,'B')
    error(['The Step 23-G results do not contain matrix B. ', ...
           'Please rerun Step 23-G and save B.']);
end

B = S.results.B;

fprintf('Loaded Step 23-G forward operator.\n');
fprintf('B size = %d x %d\n',size(B,1),size(B,2));


%% ============================================================
% 2. PARAMETERS
% ============================================================

p = parameters();

rng(p.randomSeed);

G = p.G;

fprintf('CIR dimension G = %d\n',G);


%% ============================================================
% 3. COMPUTE SVD DIRECTLY
% ============================================================

fprintf('\nComputing SVD of B...\n');

[U,Sigma,V] = svd(B,'econ');

sigma = diag(Sigma);

fprintf('Number of singular values = %d\n',length(sigma));

fprintf('\nSingular values:\n');

disp(sigma.');

fprintf('\nSingular-value ratio:\n');

fprintf('sigma_max / sigma_min = %.6e\n', ...
    sigma(1)/sigma(end));


%% ============================================================
% 4. TEST RANKS
% ============================================================

testRanks = [ ...
    30 35 40 45 50 52 54 56 58 60 62 63 64 ];

numRanks = length(testRanks);


%% ============================================================
% 5. TIKHONOV PARAMETERS
% ============================================================

lambdaList = 10.^(-14:-1);

numLambda = length(lambdaList);


%% ============================================================
% 6. NOISE-FREE RANDOM CIR TEST
% ============================================================

numTests = 200;

tsvdNF = zeros(numRanks,1);
tikNF  = zeros(numLambda,1);


%% ============================================================
% 7. TSVD NOISE-FREE
% ============================================================

for rr = 1:numRanks

    r = testRanks(rr);

    Ur = U(:,1:r);
    Vr = V(:,1:r);

    invSigma = diag(1 ./ sigma(1:r));

    Bpinv = Vr * invSigma * Ur';

    errValues = zeros(numTests,1);

    for trial = 1:numTests

        h = randn(G,1) + 1j*randn(G,1);

        g = B*h;

        hHat = Bpinv*g;

        errValues(trial) = ...
            norm(hHat-h)^2 / norm(h)^2;

    end

    tsvdNF(rr) = mean(errValues);

end


%% ============================================================
% 8. TIKHONOV NOISE-FREE
% ============================================================

for ll = 1:numLambda

    lambda = lambdaList(ll);

    filterFactor = ...
        sigma ./ (sigma.^2 + lambda);

    Breg = ...
        V * diag(filterFactor) * U';

    errValues = zeros(numTests,1);

    for trial = 1:numTests

        h = randn(G,1) + 1j*randn(G,1);

        g = B*h;

        hHat = Breg*g;

        errValues(trial) = ...
            norm(hHat-h)^2 / norm(h)^2;

    end

    tikNF(ll) = mean(errValues);

end


%% ============================================================
% 9. SUI-2 NOISE-FREE TEST
% ============================================================

numSUI2 = 200;

tsvdSUI2 = zeros(numRanks,1);
tikSUI2  = zeros(numLambda,1);


%% TSVD SUI-2

for rr = 1:numRanks

    r = testRanks(rr);

    Ur = U(:,1:r);
    Vr = V(:,1:r);

    invSigma = diag(1 ./ sigma(1:r));

    Bpinv = Vr * invSigma * Ur';

    errValues = zeros(numSUI2,1);

    for trial = 1:numSUI2

        h = sui2_channel(p);

        h = h(:);

        h64 = zeros(G,1);

        Lh = min(length(h),G);

        h64(1:Lh) = h(1:Lh);

        g = B*h64;

        hHat = Bpinv*g;

        errValues(trial) = ...
            norm(hHat-h64)^2 / norm(h64)^2;

    end

    tsvdSUI2(rr) = mean(errValues);

end


%% TIKHONOV SUI-2

for ll = 1:numLambda

    lambda = lambdaList(ll);

    filterFactor = ...
        sigma ./ (sigma.^2 + lambda);

    Breg = ...
        V * diag(filterFactor) * U';

    errValues = zeros(numSUI2,1);

    for trial = 1:numSUI2

        h = sui2_channel(p);

        h = h(:);

        h64 = zeros(G,1);

        Lh = min(length(h),G);

        h64(1:Lh) = h(1:Lh);

        g = B*h64;

        hHat = Breg*g;

        errValues(trial) = ...
            norm(hHat-h64)^2 / norm(h64)^2;

    end

    tikSUI2(ll) = mean(errValues);

end


%% ============================================================
% 10. AWGN TEST
% ============================================================

SNRlist = [0 5 10 15 20 25 30 35 40];

numSNR = length(SNRlist);

tsvdAWGN = zeros(numRanks,numSNR);

tikAWGN = zeros(numLambda,numSNR);

numAWGN = 200;


for ss = 1:numSNR

    SNRdB = SNRlist(ss);

    fprintf('\nTesting AWGN SNR = %d dB\n',SNRdB);


    %% --------------------------------------------------------
    % TSVD
    % --------------------------------------------------------

    for rr = 1:numRanks

        r = testRanks(rr);

        Ur = U(:,1:r);
        Vr = V(:,1:r);

        invSigma = diag(1 ./ sigma(1:r));

        Bpinv = Vr * invSigma * Ur';

        errValues = zeros(numAWGN,1);

        for trial = 1:numAWGN

            h = randn(G,1) + 1j*randn(G,1);

            gClean = B*h;

            signalPower = mean(abs(gClean).^2);

            noisePower = ...
                signalPower / 10^(SNRdB/10);

            noise = ...
                sqrt(noisePower/2) .* ...
                (randn(size(gClean)) + ...
                 1j*randn(size(gClean)));

            gNoisy = gClean + noise;

            hHat = Bpinv*gNoisy;

            errValues(trial) = ...
                norm(hHat-h)^2 / norm(h)^2;

        end

        tsvdAWGN(rr,ss) = mean(errValues);

    end


    %% --------------------------------------------------------
    % TIKHONOV
    % --------------------------------------------------------

    for ll = 1:numLambda

        lambda = lambdaList(ll);

        filterFactor = ...
            sigma ./ (sigma.^2 + lambda);

        Breg = ...
            V * diag(filterFactor) * U';

        errValues = zeros(numAWGN,1);

        for trial = 1:numAWGN

            h = randn(G,1) + 1j*randn(G,1);

            gClean = B*h;

            signalPower = mean(abs(gClean).^2);

            noisePower = ...
                signalPower / 10^(SNRdB/10);

            noise = ...
                sqrt(noisePower/2) .* ...
                (randn(size(gClean)) + ...
                 1j*randn(size(gClean)));

            gNoisy = gClean + noise;

            hHat = Breg*gNoisy;

            errValues(trial) = ...
                norm(hHat-h)^2 / norm(h)^2;

        end

        tikAWGN(ll,ss) = mean(errValues);

    end

end


%% ============================================================
% 11. DISPLAY TSVD NOISE-FREE
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' TSVD NOISE-FREE RANDOM-CIR RESULTS\n');
fprintf('============================================================\n');

fprintf('%8s %20s %20s\n', ...
    'Rank','Random CIR NMSE','SUI-2 NMSE');

fprintf('%s\n',repmat('-',1,55));

for rr = 1:numRanks

    fprintf('%8d %20.6e %20.6e\n', ...
        testRanks(rr), ...
        tsvdNF(rr), ...
        tsvdSUI2(rr));

end


%% ============================================================
% 12. DISPLAY TIKHONOV NOISE-FREE
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' TIKHONOV NOISE-FREE RESULTS\n');
fprintf('============================================================\n');

fprintf('%15s %20s %20s\n', ...
    'Lambda','Random CIR NMSE','SUI-2 NMSE');

fprintf('%s\n',repmat('-',1,65));

for ll = 1:numLambda

    fprintf('%15.0e %20.6e %20.6e\n', ...
        lambdaList(ll), ...
        tikNF(ll), ...
        tikSUI2(ll));

end


%% ============================================================
% 13. DISPLAY TSVD AWGN
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' TSVD AWGN NMSE\n');
fprintf('============================================================\n');

fprintf('%8s', 'Rank');

for ss = 1:numSNR

    fprintf('%14s', ...
        sprintf('%ddB',SNRlist(ss)));

end

fprintf('\n');

for rr = 1:numRanks

    fprintf('%8d',testRanks(rr));

    for ss = 1:numSNR

        fprintf('%14.4e', ...
            tsvdAWGN(rr,ss));

    end

    fprintf('\n');

end


%% ============================================================
% 14. DISPLAY TIKHONOV AWGN
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' TIKHONOV AWGN NMSE\n');
fprintf('============================================================\n');

fprintf('%15s','Lambda');

for ss = 1:numSNR

    fprintf('%14s', ...
        sprintf('%ddB',SNRlist(ss)));

end

fprintf('\n');

for ll = 1:numLambda

    fprintf('%15.0e',lambdaList(ll));

    for ss = 1:numSNR

        fprintf('%14.4e', ...
            tikAWGN(ll,ss));

    end

    fprintf('\n');

end


%% ============================================================
% 15. BEST OBSERVED TSVD
% ============================================================

[minTSVD,minTSVDidx] = min(tsvdAWGN(:));

[bestRankIdx,bestSNRIdx] = ...
    ind2sub(size(tsvdAWGN),minTSVDidx);


%% ============================================================
% 16. BEST OBSERVED TIKHONOV
% ============================================================

[minTik,minTikIdx] = min(tikAWGN(:));

[bestLambdaIdx,bestTikSNRIdx] = ...
    ind2sub(size(tikAWGN),minTikIdx);


%% ============================================================
% 17. DISPLAY BEST
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' BEST OBSERVED VALUES\n');
fprintf('============================================================\n');

fprintf('TSVD best observed NMSE     = %.6e\n',minTSVD);

fprintf('TSVD rank                   = %d\n', ...
    testRanks(bestRankIdx));

fprintf('TSVD SNR                    = %d dB\n', ...
    SNRlist(bestSNRIdx));

fprintf('\n');

fprintf('Tikhonov best observed NMSE = %.6e\n',minTik);

fprintf('Tikhonov lambda             = %.0e\n', ...
    lambdaList(bestLambdaIdx));

fprintf('Tikhonov SNR                = %d dB\n', ...
    SNRlist(bestTikSNRIdx));


%% ============================================================
% 18. SAVE RESULTS
% ============================================================

results.B = B;

results.U = U;
results.Sigma = Sigma;
results.V = V;
results.singularValues = sigma;

results.testRanks = testRanks;

results.lambdaList = lambdaList;

results.tsvdNoiseFree = tsvdNF;
results.tsvdSUI2 = tsvdSUI2;

results.tikhonovNoiseFree = tikNF;
results.tikhonovSUI2 = tikSUI2;

results.SNRlist = SNRlist;

results.tsvdAWGN = tsvdAWGN;
results.tikhonovAWGN = tikAWGN;

results.bestTSVD = minTSVD;
results.bestTSVDRank = ...
    testRanks(bestRankIdx);
results.bestTSVDSNR = ...
    SNRlist(bestSNRIdx);

results.bestTikhonov = minTik;
results.bestTikhonovLambda = ...
    lambdaList(bestLambdaIdx);
results.bestTikhonovSNR = ...
    SNRlist(bestTikSNRIdx);

save( ...
    'step23I_regularized_wdft_audit.mat', ...
    'results');


%% ============================================================
% 19. COMPLETE
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' STEP 23-I COMPLETE\n');
fprintf('============================================================\n');

fprintf('SVD was computed directly from B.\n');

fprintf('Results saved to:\n');
fprintf('step23I_regularized_wdft_audit.mat\n');

fprintf('============================================================\n\n');

end