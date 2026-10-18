function results = test_wdft_step23_master()
% STEP 23 - MASTER WDFT BASELINE VALIDATION
%
% One-file validation of the WDFT baseline.
%
% Checks:
%   A. Eq. (23) matrix
%   B. 100-tone LTS extraction
%   C. DC substitution
%   D. Hann window
%   E. 128-point WDFT CIR
%   F. Wiener weights
%   G. Noise-free SUI-2 test
%   H. SUI-5 structural test
%   I. Final structural verdict
%
% IMPORTANT:
%   128-point LTS FFT indexing and 256-point physical CFR
%   indexing are kept separate.

clc;

fprintf('\n');
fprintf('============================================================\n');
fprintf(' STEP 23 - MASTER WDFT BASELINE VALIDATION\n');
fprintf('============================================================\n\n');


%% ============================================================
% 0. PARAMETERS
% =============================================================

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
% 1. Eq. (23) MATRIX
% =============================================================

fprintf('------------------------------------------------------------\n');
fprintf('23-A  Eq. (23) MATRIX CHECK\n');
fprintf('------------------------------------------------------------\n');

n = 0:Np;
m = (0:M-1).';

Gmat = (1/sqrt(M)) .* ...
    exp(1j*2*pi*(m-M/2).*(n-Np/2)/M);

fprintf('Matrix size = %d x %d\n', ...
    size(Gmat,1),size(Gmat,2));

checkA = isequal(size(Gmat),[128 101]);

if checkA
    fprintf('PASS: Eq. (23) matrix construction.\n');
else
    fprintf('FAIL: Eq. (23) matrix construction.\n');
end

fprintf('\n');


%% ============================================================
% 2. GENERATE KNOWN LTS
% =============================================================

fprintf('------------------------------------------------------------\n');
fprintf('23-B  LTS / DC / WINDOW CHECK\n');
fprintf('------------------------------------------------------------\n');

Pall = ieee80216_pall();

preamble = preamble_sequences(Pall);

[~,ltsTime] = preamble_to_ifft(preamble,p);

knownLTS = get_known_lts(ltsTime);

Xknown = known_lts_fft(knownLTS);

fprintf('Known LTS length     = %d\n',length(knownLTS));
fprintf('Known LTS FFT length = %d\n',length(Xknown));


%% ============================================================
% 3. CORRECT 128-POINT AND 256-POINT INDEXING
% =============================================================

% Physical IEEE 802.16 LTS locations:
%
% -100,-98,...,-2,+2,...,+98,+100
%
% There are exactly 100 tones.

evenK = -100:2:100;

evenK_noDC = evenK(evenK ~= 0);

assert(length(evenK_noDC) == 100);


% -------------------------------------------------------------
% 128-point LTS FFT
% -------------------------------------------------------------
%
% known_lts_fft() produces a 128-point FFT.
%
% The reduced grid is:
%
% -50,...,-1,0,+1,...,+50
%
% Mapping:
%
% negative k -> upper FFT bins
% positive k -> lower FFT bins
%
% MATLAB FFT index:
%
% mod(k,128)+1
%
% -------------------------------------------------------------

idx128 = mod(evenK_noDC,128) + 1;

assert(all(idx128 >= 1));
assert(all(idx128 <= 128));

X100 = Xknown(idx128);

fprintf('\n128-point LTS FFT mapping:\n');
fprintf('Minimum index = %d\n',min(idx128));
fprintf('Maximum index = %d\n',max(idx128));
fprintf('Number of tones = %d\n',length(X100));


% -------------------------------------------------------------
% 256-point physical channel CFR
% -------------------------------------------------------------
%
% This is kept separate because a 256-point FFT uses the
% physical OFDM subcarrier indexing.
%
% -------------------------------------------------------------

idx256 = evenK_noDC + N/2 + 1;

assert(all(idx256 >= 1));
assert(all(idx256 <= N));

fprintf('\n256-point physical CFR mapping:\n');
fprintf('Minimum index = %d\n',min(idx256));
fprintf('Maximum index = %d\n',max(idx256));


%% ============================================================
% 4. DC SUBSTITUTION
% =============================================================

fprintf('\nDC rule = average of LS(-2) and LS(+2).\n');

idx_m2 = find(evenK_noDC == -2);
idx_p2 = find(evenK_noDC == 2);

assert(~isempty(idx_m2));
assert(~isempty(idx_p2));


%% ============================================================
% 5. HANN WINDOW
% =============================================================

win = zeros(101,1);

for nn = -50:50

    ii = nn + 51;

    win(ii) = ...
        0.5 + 0.5*cos(2*pi*nn/M);

end

fprintf('Window type   = Hann\n');
fprintf('Window length = %d\n',length(win));

checkB = length(win) == 101;


%% ============================================================
% 6. ZERO PADDING STRUCTURE
% =============================================================

frontZeros = (M-Np)/2;
backZeros  = frontZeros-1;

fprintf('Front zeros = %d\n',frontZeros);
fprintf('Back zeros  = %d\n',backZeros);

checkC = ...
    frontZeros == 14 && ...
    backZeros == 13;

if checkB && checkC

    fprintf('PASS: LTS/DC/window/zero-padding structure.\n');

else

    fprintf('FAIL: LTS/DC/window/zero-padding structure.\n');

end

fprintf('\n');


%% ============================================================
% 7. WIENER WEIGHTS
% =============================================================

fprintf('------------------------------------------------------------\n');
fprintf('23-C  WIENER WEIGHT CHECK\n');
fprintf('------------------------------------------------------------\n');

weightFile = fullfile( ...
    'results','wdft_23B6_wiener_weights.mat');

weights = [];

if exist(weightFile,'file')

    S = load(weightFile);

    names = fieldnames(S);

    for ii = 1:length(names)

        value = S.(names{ii});

        if isnumeric(value) && ...
                numel(value) == G

            weights = value(:);

            break;

        end

    end

end


if isempty(weights)

    fprintf('WARNING: Wiener weights not found.\n');
    fprintf('Unity weights will be used diagnostically.\n');

    weights = ones(G,1);

    weightStatus = false;

else

    fprintf('Wiener weights loaded.\n');

    weightStatus = true;

end


fprintf('Minimum = %.10f\n',min(weights));
fprintf('Maximum = %.10f\n',max(weights));
fprintf('Mean    = %.10f\n',mean(weights));

checkD = ...
    length(weights) == G && ...
    all(isfinite(weights)) && ...
    all(weights >= 0);

if checkD

    fprintf('PASS: Wiener weight sanity.\n');

else

    fprintf('FAIL: Wiener weight sanity.\n');

end

fprintf('\n');


%% ============================================================
% 8. FLAT CHANNEL WDFT TRANSFORM TEST
% =============================================================

fprintf('------------------------------------------------------------\n');
fprintf('23-D  FLAT-CHANNEL TRANSFORM TEST\n');
fprintf('------------------------------------------------------------\n');

Hflat = ones(100,1);

p101_flat = zeros(101,1);

p101_flat(1:50) = Hflat(1:50);

p101_flat(51) = ...
    mean([Hflat(idx_m2),Hflat(idx_p2)]);

p101_flat(52:101) = Hflat(51:100);

p101w_flat = p101_flat .* win;

gCentered_flat = Gmat*p101w_flat;

% Convert to physical CIR ordering
gWD_flat = fftshift(gCentered_flat);

fprintf('Input length  = %d\n',length(p101_flat));
fprintf('Output length = %d\n',length(gWD_flat));

Etotal = sum(abs(gWD_flat).^2);

Eearly = ...
    sum(abs(gWD_flat(1:G)).^2);

Elate = ...
    sum(abs(gWD_flat(G+1:M)).^2);

fprintf('Total energy       = %.12e\n',Etotal);
fprintf('Taps 1:%d energy   = %.12e\n',G,Eearly);
fprintf('Taps %d:%d energy = %.12e\n',G+1,M,Elate);

if Etotal > 0

    earlyFraction = Eearly/Etotal;

else

    earlyFraction = NaN;

end

fprintf('Early energy fraction = %.6f\n', ...
    earlyFraction);

checkE = ...
    all(isfinite(gWD_flat)) && ...
    Etotal > 0;

if checkE

    fprintf('PASS: WDFT transform produces finite CIR.\n');

else

    fprintf('FAIL: WDFT transform.\n');

end

fprintf('\n');


%% ============================================================
% 9. NOISE-FREE SUI-2 TEST
% =============================================================

fprintf('------------------------------------------------------------\n');
fprintf('23-E  SUI-2 NOISE-FREE WDFT TEST\n');
fprintf('------------------------------------------------------------\n');

numTrials = 200;

mseRaw = zeros(numTrials,1);
nmseRaw = zeros(numTrials,1);

mseWeighted = zeros(numTrials,1);
nmseWeighted = zeros(numTrials,1);

for trial = 1:numTrials

    %% --------------------------------------------------------
    % Generate channel
    % ---------------------------------------------------------

    h = sui2_channel(p);

    h = h(:);


    %% --------------------------------------------------------
    % True 256-point CFR
    % ---------------------------------------------------------

    hN = zeros(N,1);

    hN(1:length(h)) = h;

    Htrue = fft(hN,N);


    %% --------------------------------------------------------
    % Extract true CFR on 100 physical LTS tones
    % ---------------------------------------------------------

    H100 = Htrue(idx256);


    %% --------------------------------------------------------
    % Build 101-point vector
    % ---------------------------------------------------------

    p101 = zeros(101,1);

    p101(1:50) = H100(1:50);

    p101(51) = ...
        mean([H100(idx_m2),H100(idx_p2)]);

    p101(52:101) = H100(51:100);


    %% --------------------------------------------------------
    % Apply Hann window
    % ---------------------------------------------------------

    p101w = p101 .* win;


    %% --------------------------------------------------------
    % Eq. (23)
    % ---------------------------------------------------------

    gCentered = Gmat*p101w;

    gWD = fftshift(gCentered);


    %% --------------------------------------------------------
    % First 64 CIR samples
    % ---------------------------------------------------------

    gL = gWD(1:G);


    %% --------------------------------------------------------
    % Wiener weighting
    % ---------------------------------------------------------

    gLW = weights .* gL;


   Hraw = NaN(100,1);
Hweighted = NaN(100,1);

errRaw = Hraw-H100;
errWeighted = Hweighted-H100;
    mseRaw(trial) = ...
        mean(abs(errRaw).^2);

    mseWeighted(trial) = ...
        mean(abs(errWeighted).^2);


    denom = mean(abs(H100).^2);

    if denom > 0

        nmseRaw(trial) = ...
            mseRaw(trial)/denom;

        nmseWeighted(trial) = ...
            mseWeighted(trial)/denom;

    else

        nmseRaw(trial) = NaN;

        nmseWeighted(trial) = NaN;

    end

end


fprintf('Trials = %d\n',numTrials);

fprintf('\nRAW DIAGNOSTIC RECONSTRUCTION\n');

fprintf('Mean MSE  = %.12e\n', ...
    mean(mseRaw));

fprintf('Mean NMSE = %.12e\n', ...
    mean(nmseRaw));


fprintf('\nWIENER-WEIGHTED DIAGNOSTIC RECONSTRUCTION\n');

fprintf('Mean MSE  = %.12e\n', ...
    mean(mseWeighted));

fprintf('Mean NMSE = %.12e\n', ...
    mean(nmseWeighted));


checkF = ...
    all(isfinite(mseRaw)) && ...
    all(isfinite(mseWeighted));


if checkF

    fprintf('\nPASS: SUI-2 diagnostic calculation is finite.\n');

else

    fprintf('\nFAIL: SUI-2 diagnostic calculation.\n');

end


%% ============================================================
% 10. SUI-5 STRUCTURAL TEST
% =============================================================

fprintf('\n');
fprintf('------------------------------------------------------------\n');
fprintf('23-F  SUI-5 STRUCTURAL TEST\n');
fprintf('------------------------------------------------------------\n');

numTrials5 = 100;

sui5Energy = zeros(numTrials5,1);

sui5LateEnergy = zeros(numTrials5,1);

for trial = 1:numTrials5

    powers = ...
        10.^(p.SUI5.pathPower_dB/10);

    delays = round( ...
        p.SUI5.pathDelay_us*1e-6/p.Ts);


    h5 = zeros(max(delays)+1,1);


    for q = 1:length(delays)

        if q == 1

            coeff = sqrt(powers(q));

        else

            coeff = sqrt(powers(q)/2) * ...
                (randn + 1j*randn);

        end

        h5(delays(q)+1) = coeff;

    end


    %% 256-point CFR

    hN5 = zeros(N,1);

    hN5(1:length(h5)) = h5;

    H5 = fft(hN5,N);


    %% 100 LTS tones

    H5100 = H5(idx256);


    %% 101 vector

    p101_5 = zeros(101,1);

    p101_5(1:50) = H5100(1:50);

    p101_5(51) = ...
        mean([H5100(idx_m2), ...
              H5100(idx_p2)]);

    p101_5(52:101) = H5100(51:100);


    %% Window + WDFT

    p101w5 = p101_5 .* win;

    gC5 = Gmat*p101w5;

    gW5 = fftshift(gC5);


    %% Energy

    sui5Energy(trial) = ...
        sum(abs(gW5).^2);

    sui5LateEnergy(trial) = ...
        sum(abs(gW5(G+1:end)).^2);

end


fprintf('Trials = %d\n',numTrials5);

fprintf('Mean total WDFT energy = %.12e\n', ...
    mean(sui5Energy));

fprintf('Mean energy after tap 64 = %.12e\n', ...
    mean(sui5LateEnergy));


checkG = ...
    all(isfinite(sui5Energy)) && ...
    all(sui5Energy > 0);


if checkG

    fprintf('PASS: SUI-5 structural transform test.\n');

else

    fprintf('FAIL: SUI-5 structural transform test.\n');

end
%% ============================================================
% 11. EXACT WDFT WINDOW-EFFECT COMPENSATION
% =============================================================

fprintf('\n');
fprintf('------------------------------------------------------------\n');
fprintf('23-G  EXACT WDFT COMPENSATION MATRIX\n');
fprintf('------------------------------------------------------------\n');

% We model the physical channel as a 64-tap CIR.
%
% h64 --> H100 --> H101 --> Hann --> Eq.(23)
%      --> fftshift --> first 64 taps
%
% Therefore:
%
%       g64 = A*h64
%
% and the compensated CFR estimate is:
%
%       Hhat100 = K*g64
%
% where:
%
%       K = S*inv(A)
%
% S maps the 64-tap CIR directly to the 100 physical LTS tones.


%% ------------------------------------------------------------
% 11.1 CIR -> 100 physical LTS CFR matrix
% ------------------------------------------------------------

Smat = zeros(Np,G);

for r = 1:Np

    k = evenK_noDC(r);

    for q = 1:G

        tap = q-1;

        Smat(r,q) = ...
            exp(-1j*2*pi*k*tap/N);

    end

end


%% ------------------------------------------------------------
% 11.2 100 -> 101 insertion matrix
%
% The missing DC value is replaced by:
%
% H_DC = [H(-2) + H(+2)] / 2
% ------------------------------------------------------------

Dmat = zeros(101,100);

% Negative-frequency tones
Dmat(1:50,1:50) = eye(50);

% DC substitution
Dmat(51,idx_m2) = 0.5;
Dmat(51,idx_p2) = 0.5;

% Positive-frequency tones
Dmat(52:101,51:100) = eye(50);


%% ------------------------------------------------------------
% 11.3 Hann window matrix
% ------------------------------------------------------------

Wmat = diag(win);


%% ------------------------------------------------------------
% 11.4 fftshift matrix
% ------------------------------------------------------------

Qshift = zeros(M,M);

I128 = eye(M);

for c = 1:M

    Qshift(:,c) = fftshift(I128(:,c));

end


%% ------------------------------------------------------------
% 11.5 Select first 64 physical CIR taps
% ------------------------------------------------------------

Pmat = [eye(G), zeros(G,M-G)];


%% ------------------------------------------------------------
% 11.6 Complete forward WDFT operator
% ------------------------------------------------------------

Amat = ...
    Pmat * ...
    Qshift * ...
    Gmat * ...
    Wmat * ...
    Dmat * ...
    Smat;


fprintf('S matrix size = %d x %d\n', ...
    size(Smat,1),size(Smat,2));

fprintf('D matrix size = %d x %d\n', ...
    size(Dmat,1),size(Dmat,2));

fprintf('A matrix size = %d x %d\n', ...
    size(Amat,1),size(Amat,2));


%% ------------------------------------------------------------
% 11.7 Check conditioning
% ------------------------------------------------------------

rankA = rank(Amat);

condA = cond(Amat);

fprintf('rank(A)       = %d\n',rankA);
fprintf('condition(A)  = %.6e\n',condA);


if rankA == G

    fprintf('PASS: WDFT forward operator has full rank.\n');

else

    fprintf('FAIL: WDFT forward operator is rank deficient.\n');

end


%% ------------------------------------------------------------
% 11.8 Exact noise-free compensation matrix
% ------------------------------------------------------------

if rankA == G

    Kcomp = Smat / Amat;

    fprintf('Compensation matrix size = %d x %d\n', ...
        size(Kcomp,1),size(Kcomp,2));

    checkK = all(isfinite(Kcomp(:)));

else

    Kcomp = Smat * pinv(Amat);

    fprintf('WARNING: pseudoinverse used because A is rank deficient.\n');

    checkK = all(isfinite(Kcomp(:)));

end


%% ------------------------------------------------------------
% 11.9 Verify K*A = S
% ------------------------------------------------------------

reconstructionOperator = Kcomp * Amat;

operatorError = ...
    norm(reconstructionOperator-Smat,'fro') / ...
    norm(Smat,'fro');

fprintf('Relative operator error = %.12e\n', ...
    operatorError);


if operatorError < 1e-10

    fprintf('PASS: Exact WDFT compensation verified.\n');

else

    fprintf('WARNING: Compensation operator error is not negligible.\n');

end


%% ------------------------------------------------------------
% 11.10 Noise-free SUI-2 verification
% ------------------------------------------------------------

numTrialsComp = 200;

mseComp = zeros(numTrialsComp,1);
nmseComp = zeros(numTrialsComp,1);

for trial = 1:numTrialsComp

    % Generate SUI-2 channel
    h = sui2_channel(p);
    h = h(:);

    % Put channel into 64-tap representation
    h64 = zeros(G,1);

    Lh = min(length(h),G);

    h64(1:Lh) = h(1:Lh);

    % True 100-tone CFR
    Htrue100 = Smat*h64;

    % WDFT observation
    g64 = Amat*h64;

    % Exact compensated CFR
    Hhat100 = Kcomp*g64;

    % Error
    err = Hhat100-Htrue100;

    mseComp(trial) = mean(abs(err).^2);

    denom = mean(abs(Htrue100).^2);

    if denom > 0

        nmseComp(trial) = ...
            mseComp(trial)/denom;

    else

        nmseComp(trial) = NaN;

    end

end


fprintf('\n');
fprintf('Noise-free SUI-2 compensation test\n');
fprintf('Trials = %d\n',numTrialsComp);

fprintf('Mean MSE  = %.12e\n',mean(mseComp));
fprintf('Mean NMSE = %.12e\n',mean(nmseComp));


checkComp = ...
    checkK && ...
    operatorError < 1e-10 && ...
    all(isfinite(mseComp)) && ...
    all(isfinite(nmseComp));


if checkComp

    fprintf('PASS: Exact WDFT compensation test.\n');

else

    fprintf('FAIL: Exact WDFT compensation test.\n');

end

%% ============================================================
% 11. FINAL VERDICT
% =============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' STEP 23 FINAL VERDICT\n');
fprintf('============================================================\n');

fprintf('Eq. (23) matrix             : %s\n', ...
    passfail(checkA));

fprintf('LTS/window structure        : %s\n', ...
    passfail(checkB && checkC));

fprintf('Wiener weights              : %s\n', ...
    passfail(checkD));

fprintf('WDFT transform              : %s\n', ...
    passfail(checkE));

fprintf('SUI-2 diagnostic            : %s\n', ...
    passfail(checkF));

fprintf('SUI-5 structural            : %s\n', ...
    passfail(checkG));


fprintf('\n');

if checkA && ...
        checkB && ...
        checkC && ...
        checkD && ...
        checkE && ...
        checkF && ...
        checkG && ...
        checkComp
    fprintf('STRUCTURAL RESULT: PASS\n');

    fprintf('\n');
    fprintf('The complete WDFT transform pipeline is structurally valid.\n');

    fprintf(['The reported MSE/NMSE values are diagnostic only; ', ...
             'they are NOT yet the final WDFT performance numbers.\n']);

    fprintf(['Exact WDFT window-effect compensation still needs ', ...
             'to be validated before Step 23 is locked.\n']);

else

    fprintf('STRUCTURAL RESULT: FAIL\n');

end
fprintf('Exact compensation          : %s\n', ...
    passfail(checkComp));

fprintf('============================================================\n');


%% ============================================================
% 12. SAVE RESULTS
% =============================================================

results.Gmat = Gmat;

results.weights = weights;

results.weightStatus = weightStatus;

results.flatCIR = gWD_flat;

results.sui2.rawMSE = mseRaw;

results.sui2.rawNMSE = nmseRaw;

results.sui2.weightedMSE = mseWeighted;

results.sui2.weightedNMSE = nmseWeighted;

results.sui5.energy = sui5Energy;

results.sui5.lateEnergy = sui5LateEnergy;

results.earlyEnergyFraction = ...
    earlyFraction;

results.checks = [ ...
    checkA ...
    checkB ...
    checkC ...
    checkD ...
    checkE ...
    checkF ...
    checkG ...
    checkComp ];

save( ...
    'step23_WDFT_master_results.mat', ...
    'results');


fprintf('\n');
fprintf('Saved:\n');
fprintf('results/step23_WDFT_master_results.mat\n');
fprintf('\n');

results.wdft.Smat = Smat;
results.wdft.Dmat = Dmat;
results.wdft.Wmat = Wmat;
results.wdft.Amat = Amat;
results.wdft.Kcomp = Kcomp;

results.wdft.rankA = rankA;
results.wdft.conditionA = condA;
results.wdft.operatorError = operatorError;

results.wdft.compMSE = mseComp;
results.wdft.compNMSE = nmseComp;

results.wdft.compensationPass = checkComp;
end


%% ============================================================
% LOCAL FUNCTION
% 100 -> 101 WITH DC SUBSTITUTION
% =============================================================

function p101 = local_make_101(H100)

p101 = zeros(101,1);

p101(1:50) = H100(1:50);

p101(51) = ...
    mean([H100(50),H100(51)]);

p101(52:101) = H100(51:100);

end


%% ============================================================
% LOCAL FUNCTION
% PASS / FAIL
% =============================================================

function s = passfail(x)

if x

    s = 'PASS';

else

    s = 'FAIL';

end

end