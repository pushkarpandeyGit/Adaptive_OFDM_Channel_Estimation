clear;
clc;

fprintf('\n');
fprintf('============================================================\n');
fprintf('          SVD CHANNEL ESTIMATOR - FULL AUDIT\n');
fprintf('============================================================\n');


%% ============================================================
% [1] Identity test
% =============================================================

fprintf('\n[1] Identity test\n');

Rhp = eye(4);
U = eye(4);
Dp = eye(4);

hLS_p = [1+1j; 2-1j; -1+2j; 0.5-0.5j];

H = svd_channel_estimate(Rhp,U,Dp,hLS_p);

expected = hLS_p;

identityError = max(abs(H-expected));

fprintf('Maximum identity-test error = %.6e\n',identityError);

assert(identityError < 1e-12);

fprintf('PASS: Identity case returns LS estimate.\n');


%% ============================================================
% [2] Dimension test
% =============================================================

fprintf('\n[2] Dimension test\n');

Rhp = randn(6,3) + 1j*randn(6,3);

[U,~] = qr(randn(3,3) + 1j*randn(3,3));

Dp = diag([1 0.5 0.2]);

hLS_p = randn(3,1) + 1j*randn(3,1);

H = svd_channel_estimate(Rhp,U,Dp,hLS_p);

fprintf('R_hp size = %d x %d\n',size(Rhp,1),size(Rhp,2));
fprintf('Output size = %d x %d\n',size(H,1),size(H,2));

assert(isequal(size(H),[6 1]));

fprintf('PASS: Output dimensions are correct.\n');


%% ============================================================
% [3] Direct equation test
% =============================================================

fprintf('\n[3] Direct mathematical equation test\n');

H_direct = Rhp * U * Dp * U' * hLS_p;

equationError = max(abs(H-H_direct));

fprintf('Maximum equation error = %.6e\n',equationError);

assert(equationError < 1e-12);

fprintf('PASS: Implementation matches mathematical equation.\n');


%% ============================================================
% [4] Zero input test
% =============================================================

fprintf('\n[4] Zero-input test\n');

hZero = zeros(3,1);

H_zero = svd_channel_estimate(Rhp,U,Dp,hZero);

zeroError = max(abs(H_zero));

fprintf('Maximum output magnitude = %.6e\n',zeroError);

assert(zeroError < 1e-12);

fprintf('PASS: Zero LS input produces zero channel estimate.\n');


%% ============================================================
% [5] Complex-valued test
% =============================================================

fprintf('\n[5] Complex-valued test\n');

Rhp = randn(8,4) + 1j*randn(8,4);

A = randn(4,4) + 1j*randn(4,4);
[U,~] = qr(A);

Dp = diag([1.2 0.8 0.4 0.1]);

hLS_p = randn(4,1) + 1j*randn(4,1);

H = svd_channel_estimate(Rhp,U,Dp,hLS_p);

H_direct = Rhp * U * Dp * U' * hLS_p;

complexError = max(abs(H-H_direct));

fprintf('Maximum complex-case error = %.6e\n',complexError);

assert(complexError < 1e-12);

fprintf('PASS: Complex-valued calculation is correct.\n');


%% ============================================================
% [6] Rank truncation behavior
% =============================================================

fprintf('\n[6] Rank truncation behavior\n');

Rhp = eye(4);
U = eye(4);

Dp = diag([1 1 0 0]);

hLS_p = [1;2;3;4];

H = svd_channel_estimate(Rhp,U,Dp,hLS_p);

expected = [1;2;0;0];

rankError = max(abs(H-expected));

fprintf('Maximum rank-truncation error = %.6e\n',rankError);

assert(rankError < 1e-12);

fprintf('PASS: Components beyond selected rank are removed.\n');


%% ============================================================
% [7] Project integration test
% =============================================================

fprintf('\n[7] Project integration test\n');

p = parameters;

mapping = subcarrier_mapping(p);

R = build_uniform_covariance(p.N,5);

Rpp = build_Rpp( ...
    R, ...
    mapping.pilotLocations, ...
    p.N);

Rhp = build_Rhp( ...
    R, ...
    mapping.usedIndex, ...
    mapping.pilotLocations, ...
    p.N);
[U,Lambda] = decompose_Rpp(Rpp);

eigenvalues = diag(Lambda);

pRank = 4;

Dp = build_Dp( ...
    eigenvalues, ...
    pRank, ...
    1, ...
    10);

hLS_p = randn(p.Np,1) + 1j*randn(p.Np,1);

H_project = svd_channel_estimate( ...
    Rhp, ...
    U, ...
    Dp, ...
    hLS_p);

fprintf('R_hp size = %d x %d\n',size(Rhp,1),size(Rhp,2));
fprintf('R_pp size = %d x %d\n',size(Rpp,1),size(Rpp,2));
fprintf('Output size = %d x %d\n',size(H_project,1),size(H_project,2));

assert(isequal(size(Rhp),[p.Nu p.Np]));
assert(isequal(size(Rpp),[p.Np p.Np]));
assert(isequal(size(H_project),[p.Nu 1]));

fprintf('PASS: Project SVD estimator integration works.\n');


%% ============================================================
% [8] Invalid LS dimension test
% ============================================================

fprintf('\n[8] Invalid LS dimension test\n');

try

    svd_channel_estimate( ...
        Rhp,U,Dp, ...
        randn(p.Np+1,1));

    error('Test failed: invalid hLS dimension was accepted.');

catch ME

    if contains(ME.message,'Test failed')
        rethrow(ME);
    end

end

fprintf('PASS: Invalid LS dimension rejected.\n');


%% ============================================================
% [9] Invalid D_p dimension test
% ============================================================

fprintf('\n[9] Invalid D_p dimension test\n');

try

    svd_channel_estimate( ...
        Rhp,U, ...
        eye(p.Np+1), ...
        hLS_p);

    error('Test failed: invalid D_p dimension was accepted.');

catch ME

    if contains(ME.message,'Test failed')
        rethrow(ME);
    end

end

fprintf('PASS: Invalid D_p dimension rejected.\n');


%% ============================================================
% FINAL RESULT
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('       SVD CHANNEL ESTIMATOR AUDIT COMPLETE\n');
fprintf('============================================================\n');

fprintf('ALL SVD CHANNEL ESTIMATOR TESTS PASSED.\n');

fprintf('============================================================\n');