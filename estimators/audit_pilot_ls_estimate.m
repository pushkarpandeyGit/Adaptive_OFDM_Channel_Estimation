clear;
clc;

fprintf('\n');
fprintf('============================================================\n');
fprintf('            PILOT LS ESTIMATOR - FULL AUDIT\n');
fprintf('============================================================\n');


%% ============================================================
% [1] Basic noiseless test
% =============================================================

fprintf('\n[1] Basic noiseless LS test\n');

X_p = [1 -1 1 1 -1 1 -1 1];

H0 = 2 + 3j;

Y_p = H0 .* X_p;

hLS_p = pilot_ls_estimate(Y_p,X_p);

expected = H0 * ones(size(X_p));

errorValue = max(abs(hLS_p - expected));

fprintf('Maximum LS error = %.6e\n',errorValue);

assert(errorValue < 1e-12);

fprintf('PASS: Noiseless LS estimation is exact.\n');


%% ============================================================
% [2] Mathematical equation test
% =============================================================

fprintf('\n[2] Direct LS equation test\n');

directEstimate = Y_p ./ X_p;

equationError = max(abs(hLS_p - directEstimate));

fprintf('Maximum equation error = %.6e\n',equationError);

assert(equationError < 1e-12);

fprintf('PASS: Implementation matches Y_p ./ X_p.\n');


%% ============================================================
% [3] Complex channel test
% =============================================================

fprintf('\n[3] Complex channel test\n');

X_p = [ ...
    1+0j ...
    -1+0j ...
    1+0j ...
    -1+0j ...
    1+0j ...
    -1+0j ...
    1+0j ...
    -1+0j ];

H_complex = ...
    [ ...
    0.8+0.4j ...
    1.1-0.7j ...
    -0.5+1.3j ...
    0.2-0.9j ...
    1.4+0.2j ...
    -1.0-0.6j ...
    0.3+1.1j ...
    -0.7+0.8j ];

Y_p = H_complex .* X_p;

hLS_p = pilot_ls_estimate(Y_p,X_p);

complexError = max(abs(hLS_p - H_complex));

fprintf('Maximum complex LS error = %.6e\n',complexError);

assert(complexError < 1e-12);

fprintf('PASS: Complex channel estimation is correct.\n');


%% ============================================================
% [4] Noisy LS test
% =============================================================

fprintf('\n[4] Noisy LS test\n');

rng(42);

X_p = ones(1,8);

H0 = 1.5 - 0.8j;

noise = 0.01 * ...
    (randn(1,8) + 1j*randn(1,8));

Y_p = H0 .* X_p + noise;

hLS_p = pilot_ls_estimate(Y_p,X_p);

mse = mean(abs(hLS_p-H0).^2);

fprintf('Noisy LS MSE = %.6e\n',mse);

assert(isfinite(mse));
assert(mse > 0);

fprintf('PASS: Noisy LS estimate behaves correctly.\n');


%% ============================================================
% [5] Dimension test
% =============================================================

fprintf('\n[5] Dimension test\n');

X_p = randn(1,8) + 1j*randn(1,8);
Y_p = randn(1,8) + 1j*randn(1,8);

hLS_p = pilot_ls_estimate(Y_p,X_p);

fprintf('Y_p size   = %d x %d\n', ...
    size(Y_p,1),size(Y_p,2));

fprintf('X_p size   = %d x %d\n', ...
    size(X_p,1),size(X_p,2));

fprintf('hLS_p size = %d x %d\n', ...
    size(hLS_p,1),size(hLS_p,2));

assert(isequal(size(hLS_p),[1 8]));

fprintf('PASS: Output dimension is correct.\n');


%% ============================================================
% [6] Column-input consistency
% =============================================================

fprintf('\n[6] Row/column input consistency test\n');

X_row = [1 2 3 4];
Y_row = [2+1j 4+2j 6+3j 8+4j];

H_row = pilot_ls_estimate(Y_row,X_row);

H_col = pilot_ls_estimate(Y_row.',X_row.');

consistencyError = max(abs(H_row-H_col));

fprintf('Row/column consistency error = %.6e\n', ...
    consistencyError);

assert(consistencyError < 1e-12);

fprintf('PASS: Row and column inputs give identical estimates.\n');


%% ============================================================
% [7] Zero-pilot rejection
% =============================================================

fprintf('\n[7] Zero known-pilot rejection test\n');

X_bad = [1 1 0 1];
Y_bad = [1 2 3 4];

try

    pilot_ls_estimate(Y_bad,X_bad);

    error('Test failed: zero pilot was accepted.');

catch ME

    if contains(ME.message,'Test failed')
        rethrow(ME);
    end

end

fprintf('PASS: Zero known pilot rejected.\n');


%% ============================================================
% [8] Mismatched length rejection
% =============================================================

fprintf('\n[8] Mismatched-length rejection test\n');

try

    pilot_ls_estimate( ...
        [1 2 3], ...
        [1 2]);

    error('Test failed: mismatched lengths were accepted.');

catch ME

    if contains(ME.message,'Test failed')
        rethrow(ME);
    end

end

fprintf('PASS: Mismatched pilot lengths rejected.\n');


%% ============================================================
% [9] Empty input rejection
% =============================================================

fprintf('\n[9] Empty-input rejection test\n');

try

    pilot_ls_estimate([],[]);

    error('Test failed: empty input was accepted.');

catch ME

    if contains(ME.message,'Test failed')
        rethrow(ME);
    end

end

fprintf('PASS: Empty input rejected.\n');


%% ============================================================
% [10] Project 8-pilot integration
% =============================================================

fprintf('\n[10] Project 8-pilot integration test\n');

p = parameters;

mapping = subcarrier_mapping(p);

% Create known pilot values.

X_project = ones(1,p.Np);

H_project = ...
    0.75 - 1.25j;

Y_project = H_project .* X_project;

hLS_project = ...
    pilot_ls_estimate(Y_project,X_project);

projectError = ...
    max(abs(hLS_project-H_project));

fprintf('Number of project pilots = %d\n',p.Np);
fprintf('Project LS error          = %.6e\n',projectError);

assert(length(hLS_project) == p.Np);
assert(projectError < 1e-12);

fprintf('PASS: Project 8-pilot LS estimation works.\n');


%% ============================================================
% FINAL RESULT
% =============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('          PILOT LS ESTIMATOR AUDIT COMPLETE\n');
fprintf('============================================================\n');

fprintf('ALL PILOT LS ESTIMATOR TESTS PASSED.\n');

fprintf('============================================================\n');