function results = step23G_basis_audit()

% ============================================================
% STEP 23-G : WDFT BASIS / SINGULAR-VALUE AUDIT
%
% Purpose:
%   Determine whether the complete 128-output WDFT transform
%   preserves all 64 CIR dimensions and how well-conditioned
%   those dimensions are.
%
% IMPORTANT:
%   This is an AUDIT ONLY.
%   It does NOT change the WDFT estimator.
%   It does NOT apply Wiener weighting.
%   It does NOT use a pseudoinverse for reconstruction.
%
% Forward model:
%
%   h64
%      |
%      v
%   H100
%      |
%      v
%   H101
%      |
%      v
%   Hann window
%      |
%      v
%   Eq. (23)
%      |
%      v
%   fftshift
%      |
%      v
%   g128
%
% Therefore:
%
%   g128 = B * h64
%
% where B is 128 x 64.
%
% ============================================================

clc;

fprintf('\n');
fprintf('============================================================\n');
fprintf(' STEP 23-G : WDFT BASIS / SINGULAR-VALUE AUDIT\n');
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
% 2. PHYSICAL LTS TONE LOCATIONS
% ============================================================

evenK = -100:2:100;

evenK_noDC = evenK(evenK ~= 0);

assert(length(evenK_noDC) == 100);

fprintf('Number of physical LTS tones = %d\n', ...
    length(evenK_noDC));


%% ============================================================
% 3. 100 -> 101 DC SUBSTITUTION INDEX
% ============================================================

idx_m2 = find(evenK_noDC == -2);
idx_p2 = find(evenK_noDC == 2);

assert(~isempty(idx_m2));
assert(~isempty(idx_p2));


%% ============================================================
% 4. HANN WINDOW
% ============================================================

win = zeros(Np+1,1);

for nn = -50:50

    ii = nn + 51;

    win(ii) = ...
        0.5 + 0.5*cos(2*pi*nn/M);

end

assert(length(win) == 101);


%% ============================================================
% 5. EQ. (23) MATRIX
% ============================================================

n = 0:Np;
m = (0:M-1).';

Gmat = (1/sqrt(M)) .* ...
    exp(1j*2*pi*(m-M/2).*(n-Np/2)/M);

assert(isequal(size(Gmat),[128 101]));


%% ============================================================
% 6. CIR -> 100 PHYSICAL CFR MATRIX
% ============================================================
%
% For a 64-tap CIR:
%
%   h[q]
%
% the physical CFR at subcarrier k is:
%
%   H(k) =
%       sum h[q] exp(-j*2*pi*k*q/N)
%
% q = 0,...,63
%
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

fprintf('\nS matrix size = %d x %d\n', ...
    size(Smat,1),size(Smat,2));


%% ============================================================
% 7. 100 -> 101 INSERTION MATRIX
% ============================================================
%
% Missing DC is estimated as:
%
%   H(0) = [H(-2) + H(+2)]/2
%
% ============================================================

Dmat = zeros(101,100);

% Negative-frequency tones
Dmat(1:50,1:50) = eye(50);

% DC substitution
Dmat(51,idx_m2) = 0.5;
Dmat(51,idx_p2) = 0.5;

% Positive-frequency tones
Dmat(52:101,51:100) = eye(50);

fprintf('D matrix size = %d x %d\n', ...
    size(Dmat,1),size(Dmat,2));


%% ============================================================
% 8. HANN WINDOW MATRIX
% ============================================================

Wmat = diag(win);

fprintf('W matrix size = %d x %d\n', ...
    size(Wmat,1),size(Wmat,2));


%% ============================================================
% 9. FFTSHIFT MATRIX
% ============================================================

Qshift = zeros(M,M);

I128 = eye(M);

for c = 1:M

    Qshift(:,c) = fftshift(I128(:,c));

end

fprintf('Qshift size = %d x %d\n', ...
    size(Qshift,1),size(Qshift,2));


%% ============================================================
% 10. FULL 128-OUTPUT WDFT OPERATOR
% ============================================================
%
% IMPORTANT:
%
% We intentionally DO NOT use:
%
%   Pmat = [eye(G), zeros(G,M-G)]
%
% because that discards the last 64 WDFT outputs.
%
% Instead:
%
%   B = Qshift * Gmat * Wmat * Dmat * Smat
%
% giving:
%
%   128 x 64
%
% ============================================================

B = ...
    Qshift * ...
    Gmat * ...
    Wmat * ...
    Dmat * ...
    Smat;

fprintf('\n');
fprintf('Full WDFT operator B size = %d x %d\n', ...
    size(B,1),size(B,2));


%% ============================================================
% 11. FINITE-VALUE CHECK
% ============================================================

checkFinite = all(isfinite(B(:)));

if checkFinite

    fprintf('PASS: B contains only finite values.\n');

else

    fprintf('FAIL: B contains NaN or Inf.\n');

end


%% ============================================================
% 12. SVD
% ============================================================

[U,S,V] = svd(B,'econ');

singularValues = diag(S);

fprintf('\n');
fprintf('------------------------------------------------------------\n');
fprintf(' SINGULAR VALUES\n');
fprintf('------------------------------------------------------------\n');

for k = 1:length(singularValues)

    fprintf('%2d : %.12e\n', ...
        k,singularValues(k));

end


%% ============================================================
% 13. NORMALIZED SINGULAR VALUES
% ============================================================

normalizedSingularValues = ...
    singularValues / singularValues(1);

fprintf('\n');
fprintf('------------------------------------------------------------\n');
fprintf(' NORMALIZED SINGULAR VALUES\n');
fprintf('------------------------------------------------------------\n');

for k = 1:length(normalizedSingularValues)

    fprintf('%2d : %.12e\n', ...
        k,normalizedSingularValues(k));

end


%% ============================================================
% 14. CONDITION NUMBER
% ============================================================

sigmaMax = singularValues(1);
sigmaMin = singularValues(end);

condB = sigmaMax/sigmaMin;

fprintf('\n');
fprintf('------------------------------------------------------------\n');
fprintf(' CONDITIONING\n');
fprintf('------------------------------------------------------------\n');

fprintf('Largest singular value  = %.12e\n',sigmaMax);
fprintf('Smallest singular value = %.12e\n',sigmaMin);
fprintf('Condition number        = %.12e\n',condB);


%% ============================================================
% 15. RANK AT MULTIPLE TOLERANCES
% ============================================================
%
% This is important because MATLAB's numerical rank depends
% on the tolerance.
%
% We inspect several tolerances rather than trusting one value.
%
% ============================================================

tolerances = [ ...
    1e-8 ...
    1e-9 ...
    1e-10 ...
    1e-11 ...
    1e-12 ...
    1e-13 ...
    1e-14 ];

rankAtTolerance = zeros(size(tolerances));

fprintf('\n');
fprintf('------------------------------------------------------------\n');
fprintf(' NUMERICAL RANK VS TOLERANCE\n');
fprintf('------------------------------------------------------------\n');

for ii = 1:length(tolerances)

    tol = tolerances(ii);

    rankAtTolerance(ii) = ...
        sum(singularValues > sigmaMax*tol);

    fprintf('tol = %.0e   rank = %d / %d\n', ...
        tol,rankAtTolerance(ii),G);

end


%% ============================================================
% 16. DEFAULT MATLAB RANK
% ============================================================

rankDefault = rank(B);

fprintf('\n');
fprintf('MATLAB default numerical rank = %d / %d\n', ...
    rankDefault,G);


%% ============================================================
% 17. BASIS-VECTOR TEST
% ============================================================
%
% Every column of B is the WDFT output generated by one
% unit-energy CIR basis vector:
%
%   e1, e2, ..., e64
%
% Therefore:
%
%   B(:,q)
%
% is the WDFT response of CIR tap q.
%
% We verify that all basis responses are finite and nonzero.
%
% ============================================================

basisNorms = zeros(G,1);

for q = 1:G

    e = zeros(G,1);

    e(q) = 1;

    g = B*e;

    basisNorms(q) = norm(g);

end

checkBasis = ...
    all(isfinite(basisNorms)) && ...
    all(basisNorms > 0);

fprintf('\n');
fprintf('------------------------------------------------------------\n');
fprintf(' BASIS VECTOR TEST\n');
fprintf('------------------------------------------------------------\n');

fprintf('Minimum basis-response norm = %.12e\n', ...
    min(basisNorms));

fprintf('Maximum basis-response norm = %.12e\n', ...
    max(basisNorms));

if checkBasis

    fprintf('PASS: All 64 CIR basis vectors produce nonzero outputs.\n');

else

    fprintf('FAIL: At least one CIR basis vector is lost.\n');

end


%% ============================================================
% 18. ORTHOGONALITY / GRAM MATRIX CHECK
% ============================================================
%
% B^H B describes how distinguishable the CIR dimensions are
% after the WDFT transform.
%
% If columns are strongly correlated, inversion becomes unstable.
%
% ============================================================

Gram = B' * B;

gramEigenvalues = eig(Gram);

gramEigenvalues = real(gramEigenvalues);

gramEigenvalues = sort(gramEigenvalues,'descend');

fprintf('\n');
fprintf('------------------------------------------------------------\n');
fprintf(' GRAM MATRIX\n');
fprintf('------------------------------------------------------------\n');

fprintf('Largest eigenvalue  = %.12e\n', ...
    gramEigenvalues(1));

fprintf('Smallest eigenvalue = %.12e\n', ...
    gramEigenvalues(end));


%% ============================================================
% 19. FULL-OPERATOR RECONSTRUCTION TEST
% ============================================================
%
% IMPORTANT:
% This is NOT used to choose a reconstruction method.
%
% It only checks whether B contains enough information to
% reproduce arbitrary 64-dimensional CIR vectors numerically.
%
% We use MATLAB backslash on the overdetermined system:
%
%       B * h ~= g
%
% ============================================================

numTests = 100;

reconstructionError = zeros(numTests,1);

for trial = 1:numTests

    hTest = ...
        randn(G,1) + 1j*randn(G,1);

    gTest = B*hTest;

    hRecovered = B \ gTest;

    reconstructionError(trial) = ...
        norm(hRecovered-hTest)^2 / ...
        norm(hTest)^2;

end

fprintf('\n');
fprintf('------------------------------------------------------------\n');
fprintf(' FULL-OPERATOR NUMERICAL RECONSTRUCTION\n');
fprintf('------------------------------------------------------------\n');

fprintf('Tests = %d\n',numTests);

fprintf('Mean relative squared error = %.12e\n', ...
    mean(reconstructionError));

fprintf('Maximum relative squared error = %.12e\n', ...
    max(reconstructionError));


%% ============================================================
% 20. VERDICT
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' STEP 23-G AUDIT VERDICT\n');
fprintf('============================================================\n');

fprintf('Finite full operator       : %s\n', ...
    passfail(checkFinite));

fprintf('All basis vectors visible  : %s\n', ...
    passfail(checkBasis));

fprintf('Default numerical rank     : %d / %d\n', ...
    rankDefault,G);

fprintf('Condition number           : %.6e\n', ...
    condB);

fprintf('Mean reconstruction error  : %.6e\n', ...
    mean(reconstructionError));

fprintf('Maximum reconstruction err : %.6e\n', ...
    max(reconstructionError));

fprintf('\n');

if checkFinite && checkBasis && rankDefault == G

    fprintf('STRUCTURAL RESULT: FULL OPERATOR IS FULL RANK\n');

else

    fprintf('STRUCTURAL RESULT: FULL OPERATOR HAS LOST DIMENSIONS\n');

end

fprintf('\n');

if condB > 1e8

    fprintf('STABILITY WARNING:\n');
    fprintf(['The full WDFT operator is extremely ill-conditioned.\n' ...
             'A direct inverse may amplify noise strongly.\n']);

elseif condB > 1e5

    fprintf('STABILITY WARNING:\n');
    fprintf(['The full WDFT operator is poorly conditioned.\n']);

else

    fprintf('STABILITY: Operator conditioning is acceptable.\n');

end

fprintf('\n');
fprintf('DO NOT use this script to modify the estimator yet.\n');
fprintf('This is an operator audit only.\n');

fprintf('============================================================\n');


%% ============================================================
% 21. SAVE AUDIT RESULTS
% ============================================================

results.B = B;

results.Smat = Smat;
results.Dmat = Dmat;
results.Wmat = Wmat;
results.Qshift = Qshift;
results.Gmat = Gmat;

results.singularValues = singularValues;
results.normalizedSingularValues = ...
    normalizedSingularValues;

results.conditionNumber = condB;

results.rankDefault = rankDefault;

results.tolerances = tolerances;
results.rankAtTolerance = rankAtTolerance;

results.basisNorms = basisNorms;

results.gramEigenvalues = ...
    gramEigenvalues;

results.reconstructionError = ...
    reconstructionError;

results.checkFinite = checkFinite;
results.checkBasis = checkBasis;

save('step23G_basis_audit_results.mat','results');

fprintf('\nSaved:\n');
fprintf('step23G_basis_audit_results.mat\n\n');

end


%% ============================================================
% LOCAL FUNCTION
% ============================================================

function s = passfail(x)

if x
    s = 'PASS';
else
    s = 'FAIL';
end

end