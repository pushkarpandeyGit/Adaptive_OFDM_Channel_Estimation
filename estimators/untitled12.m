clear;
clc;

fprintf('\n');
fprintf('============================================================\n');
fprintf('             LTS PILOT EXTRACTION - FULL AUDIT\n');
fprintf('============================================================\n');


%% ============================================================
% [1] Basic project parameters
% =============================================================

fprintf('\n[1] Loading project parameters\n');

p = parameters;

fprintf('N  = %d\n',p.N);
fprintf('Np = %d\n',p.Np);

assert(p.N == 256);
assert(p.Np == 8);

fprintf('PASS: Project parameters loaded correctly.\n');


%% ============================================================
% [2] Construct a known 256-bin spectrum
% =============================================================

fprintf('\n[2] Construct known frequency-domain spectrum\n');

N = p.N;

X = zeros(1,N);

for k = -N/2:N/2-1

    bin = k + N/2 + 1;

    X(bin) = k + 1000;

end

% Make sure pilot values are non-zero.

fprintf('Spectrum length = %d\n',length(X));

assert(length(X) == N);

fprintf('PASS: Test spectrum created.\n');


%% ============================================================
% [3] Create received spectrum with known channel
% =============================================================

fprintf('\n[3] Create artificial received spectrum\n');

% Use a known complex channel gain.

H0 = 2 + 3j;

Y = H0 * X;

fprintf('Known channel gain = %.2f + %.2fj\n', ...
    real(H0),imag(H0));

fprintf('PASS: Artificial received spectrum created.\n');


%% ============================================================
% [4] Pilot extraction
% =============================================================

fprintf('\n[4] Extract pilot subcarriers\n');

[Y_p,X_p,pilotBins] = extract_lts_pilots( ...
    Y, ...
    X, ...
    p.pilotLocations, ...
    p.N);


%% ============================================================
% [5] Dimension test
% =============================================================

fprintf('\n[5] Dimension test\n');

fprintf('Y_p size = %d x %d\n', ...
    size(Y_p,1),size(Y_p,2));

fprintf('X_p size = %d x %d\n', ...
    size(X_p,1),size(X_p,2));

fprintf('pilotBins size = %d x %d\n', ...
    size(pilotBins,1),size(pilotBins,2));

assert(length(Y_p) == p.Np);
assert(length(X_p) == p.Np);
assert(length(pilotBins) == p.Np);

fprintf('PASS: Eight pilot values extracted.\n');


%% ============================================================
% [6] Pilot-bin mapping test
% =============================================================

fprintf('\n[6] Pilot-bin mapping test\n');

expectedBins = ...
    p.pilotLocations + p.N/2 + 1;

mappingError = max(abs(pilotBins - expectedBins));

fprintf('Maximum mapping error = %.6e\n',mappingError);

assert(mappingError < 1e-12);

fprintf('PASS: Signed subcarrier locations map correctly.\n');


%% ============================================================
% [7] Direct extraction test
% =============================================================

fprintf('\n[7] Direct extraction test\n');

expectedY = Y(expectedBins);
expectedX = X(expectedBins);

Y_error = max(abs(Y_p - expectedY));
X_error = max(abs(X_p - expectedX));

fprintf('Y_p extraction error = %.6e\n',Y_error);
fprintf('X_p extraction error = %.6e\n',X_error);

assert(Y_error < 1e-12);
assert(X_error < 1e-12);

fprintf('PASS: Extracted values match direct indexing.\n');


%% ============================================================
% [8] Known channel consistency test
% =============================================================

fprintf('\n[8] Known channel consistency test\n');

H_p = Y_p ./ X_p;

channelError = max(abs(H_p - H0));

fprintf('Maximum channel-gain error = %.6e\n',channelError);

assert(channelError < 1e-12);

fprintf('PASS: Pilot values preserve the known channel relationship.\n');


%% ============================================================
% [9] Pilot ordering test
% =============================================================

fprintf('\n[9] Pilot ordering test\n');

orderingError = max(abs( ...
    pilotBins - expectedBins));

assert(orderingError < 1e-12);

fprintf('Pilot locations:\n');
disp(p.pilotLocations);

fprintf('Pilot MATLAB bins:\n');
disp(pilotBins);

fprintf('PASS: Pilot ordering is preserved.\n');


%% ============================================================
% [10] Zero-pilot rejection test
% =============================================================

fprintf('\n[10] Zero known-pilot rejection test\n');

X_bad = X;

X_bad(expectedBins(3)) = 0;

try

    extract_lts_pilots( ...
        Y, ...
        X_bad, ...
        p.pilotLocations, ...
        p.N);

    error('Test failed: zero pilot was accepted.');

catch ME

    if contains(ME.message,'Test failed')
        rethrow(ME);
    end

end

fprintf('PASS: Zero known pilot rejected.\n');


%% ============================================================
% [11] Invalid pilot location test
% =============================================================

fprintf('\n[11] Invalid pilot-location test\n');

badLocations = [ ...
    -88 -63 -38 -13 ...
     13  38  63  200];

try

    extract_lts_pilots( ...
        Y, ...
        X, ...
        badLocations, ...
        N);

    error('Test failed: invalid pilot location was accepted.');

catch ME

    if contains(ME.message,'Test failed')
        rethrow(ME);
    end

end

fprintf('PASS: Invalid pilot location rejected.\n');


%% ============================================================
% [12] Duplicate pilot test
% =============================================================

fprintf('\n[12] Duplicate-pilot test\n');

duplicateLocations = [ ...
    -88 -63 -38 -13 ...
     13  38  63  63];

try

    extract_lts_pilots( ...
        Y, ...
        X, ...
        duplicateLocations, ...
        N);

    error('Test failed: duplicate pilot was accepted.');

catch ME

    if contains(ME.message,'Test failed')
        rethrow(ME);
    end

end

fprintf('PASS: Duplicate pilot location rejected.\n');


%% ============================================================
% [13] Complex channel test
% =============================================================

fprintf('\n[13] Complex channel test\n');

H_complex = ...
    0.73 - 1.42j;

Y_complex = H_complex * X;

[Y_pc, X_pc, ~] = extract_lts_pilots( ...
    Y_complex, ...
    X, ...
    p.pilotLocations, ...
    p.N);

H_estimated = Y_pc ./ X_pc;

complexError = max(abs(H_estimated - H_complex));

fprintf('Complex-channel estimation error = %.6e\n', ...
    complexError);

assert(complexError < 1e-12);

fprintf('PASS: Complex channel preserved correctly.\n');


%% ============================================================
% FINAL RESULT
% =============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('          LTS PILOT EXTRACTION AUDIT COMPLETE\n');
fprintf('============================================================\n');

fprintf('ALL LTS PILOT EXTRACTION TESTS PASSED.\n');

fprintf('============================================================\n');