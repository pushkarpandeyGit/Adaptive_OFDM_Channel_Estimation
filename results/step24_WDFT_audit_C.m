%% ============================================================
% STEP 24 - WDFT AUDIT C (FIXED)
%
% Compare:
%   1. Exact Wiener weights theta from successful Step 23-J
%   2. Wiener-weight-related variables from B6 MAT file
%
% IMPORTANT:
%   This version NEVER assumes the first 64-element vector is
%   the Wiener weight vector.
% ============================================================

clear;
clc;
close all;

fprintf('\n=============================================\n');
fprintf('STEP 24 - WDFT AUDIT C (FIXED)\n');
fprintf('WDFT Wiener Weight Comparison\n');
fprintf('=============================================\n\n');

%% ============================================================
% 1. LOAD SUCCESSFUL STEP 23-J RESULT
% ============================================================

file23J = 'step23J_actual_wdft_audit.mat';

if ~isfile(file23J)

    error(['Cannot find ',file23J]);

end

S23 = load(file23J);

if ~isfield(S23,'step23J')

    error('Variable "step23J" not found.');

end

step23J = S23.step23J;

if ~isfield(step23J,'theta')

    error('Field "theta" not found inside step23J.');

end

% EXACT Wiener weights used by successful Step 23-J
theta = step23J.theta(:);

if length(theta) ~= 64

    error( ...
        'Step 23-J theta must contain exactly 64 elements.');

end

%% ============================================================
% 2. DISPLAY STEP 23-J THETA
% ============================================================

fprintf('---------------------------------------------\n');
fprintf('STEP 23-J EXACT THETA\n');
fprintf('---------------------------------------------\n');

fprintf('Size       : %d x %d\n', ...
    size(theta,1), ...
    size(theta,2));

fprintf('Min        : %.12e\n',min(theta));

fprintf('Max        : %.12e\n',max(theta));

fprintf('Mean       : %.12e\n',mean(theta));

fprintf('Norm       : %.12e\n',norm(theta));

fprintf('Nonzero    : %d / %d\n', ...
    nnz(theta), ...
    length(theta));

fprintf('\nAll 64 theta values:\n');

disp(theta);

%% ============================================================
% 3. LOAD B6 MAT FILE
% ============================================================

fileB6 = 'wdft_23B6_wiener_weights.mat';

if ~isfile(fileB6)

    error(['Cannot find ',fileB6]);

end

SB6 = load(fileB6);

namesB6 = fieldnames(SB6);

%% ============================================================
% 4. DISPLAY ALL B6 VARIABLES
% ============================================================

fprintf('\n=============================================\n');
fprintf('B6 MAT FILE CONTENTS\n');
fprintf('=============================================\n');

for ii = 1:length(namesB6)

    name = namesB6{ii};

    value = SB6.(name);

    fprintf('\nVariable: %s\n',name);
    fprintf('Class   : %s\n',class(value));

    fprintf('Size    : ');

    sz = size(value);

    fprintf('%d',sz(1));

    for jj = 2:length(sz)

        fprintf(' x %d',sz(jj));

    end

    fprintf('\n');

    if isnumeric(value) && ~isempty(value)

        fprintf('Min abs : %.12e\n', ...
            min(abs(value(:))));

        fprintf('Max abs : %.12e\n', ...
            max(abs(value(:))));

        fprintf('Mean abs: %.12e\n', ...
            mean(abs(value(:))));

        fprintf('Norm    : %.12e\n', ...
            norm(value(:)));

    end

end

%% ============================================================
% 5. EXPLICITLY IDENTIFY B6 WEIGHT VARIABLE
% ============================================================

fprintf('\n=============================================\n');
fprintf('SEARCHING FOR B6 WIENER WEIGHT\n');
fprintf('=============================================\n');

% Possible names that may represent the final Wiener weights.
candidateNames = { ...
    'theta', ...
    'weights', ...
    'wienerWeights', ...
    'wdftWeights', ...
    'w', ...
    'W'};

weightsB6 = [];

selectedName = '';

% First: look for an explicit known weight name.
for ii = 1:length(candidateNames)

    candidateName = candidateNames{ii};

    if isfield(SB6,candidateName)

        candidate = SB6.(candidateName);

        if isnumeric(candidate) && ...
                isvector(candidate) && ...
                length(candidate) == 64

            weightsB6 = candidate(:);

            selectedName = candidateName;

            break;

        end

    end

end

%% ============================================================
% 6. IF NO EXPLICIT WEIGHT FOUND, REPORT CANDIDATES
% ============================================================

if isempty(weightsB6)

    fprintf('\nNo explicitly named 64-element Wiener-weight vector\n');
    fprintf('was found in the B6 file.\n\n');

    fprintf('64-element numeric variables found:\n');

    foundCandidate = false;

    for ii = 1:length(namesB6)

        name = namesB6{ii};

        value = SB6.(name);

        if isnumeric(value) && ...
                isvector(value) && ...
                length(value) == 64

            foundCandidate = true;

            fprintf('\nCandidate: %s\n',name);

            fprintf('  Min abs : %.12e\n', ...
                min(abs(value(:))));

            fprintf('  Max abs : %.12e\n', ...
                max(abs(value(:))));

            fprintf('  Mean abs: %.12e\n', ...
                mean(abs(value(:))));

            fprintf('  Norm    : %.12e\n', ...
                norm(value(:)));

        end

    end

    if ~foundCandidate

        fprintf('  NONE\n');

    end

    fprintf('\nIMPORTANT:\n');
    fprintf(['The B6 MAT file does not contain an explicitly ', ...
        'named final Wiener-weight vector.\n']);

else

    %% ========================================================
    % 7. DISPLAY IDENTIFIED B6 WEIGHTS
    % ========================================================

    fprintf('\nB6 Wiener-weight variable found: %s\n', ...
        selectedName);

    weightsB6 = real(weightsB6(:));

    fprintf('\n---------------------------------------------\n');
    fprintf('B6 IDENTIFIED WEIGHTS\n');
    fprintf('---------------------------------------------\n');

    fprintf('Size       : %d x %d\n', ...
        size(weightsB6,1), ...
        size(weightsB6,2));

    fprintf('Min        : %.12e\n',min(weightsB6));

    fprintf('Max        : %.12e\n',max(weightsB6));

    fprintf('Mean       : %.12e\n',mean(weightsB6));

    fprintf('Norm       : %.12e\n',norm(weightsB6));

    fprintf('Nonzero    : %d / %d\n', ...
        nnz(weightsB6), ...
        length(weightsB6));

    fprintf('\nFirst 20 values:\n');

    disp(weightsB6(1:20));

    %% ========================================================
    % 8. DIRECT COMPARISON
    % ========================================================

    difference = theta - weightsB6;

    fprintf('\n=============================================\n');
    fprintf('DIRECT COMPARISON\n');
    fprintf('=============================================\n');

    fprintf('||theta - B6||       : %.12e\n', ...
        norm(difference));

    fprintf('Relative difference  : %.12e\n', ...
        norm(difference)/max(norm(theta),eps));

    %% ========================================================
    % 9. CORRELATION
    % ========================================================

    thetaReal = real(theta);
    b6Real = real(weightsB6);

    if norm(thetaReal) > eps && norm(b6Real) > eps

        correlation = ...
            abs(thetaReal' * b6Real) / ...
            (norm(thetaReal)*norm(b6Real));

        fprintf('Absolute correlation : %.12e\n', ...
            correlation);

    end

    %% ========================================================
    % 10. PLOT
    % ========================================================

    figure;

    stem(1:64,thetaReal,'o');
    hold on;

    stem(1:64,b6Real,'x');

    grid on;

    xlabel('CIR tap index');

    ylabel('Wiener weight');

    title('Step 23-J vs B6 WDFT Wiener Weights');

    legend( ...
        'Step 23-J \theta', ...
        ['B6 ',selectedName], ...
        'Location','best');

end

%% ============================================================
% 11. SAVE AUDIT
% ============================================================

if isempty(weightsB6)

    save( ...
        'step24_WDFT_audit_C_fixed.mat', ...
        'theta');

else

    save( ...
        'step24_WDFT_audit_C_fixed.mat', ...
        'theta', ...
        'weightsB6', ...
        'selectedName', ...
        'difference');

end

%% ============================================================
% 12. COMPLETE
% ============================================================

fprintf('\n=============================================\n');
fprintf('AUDIT C (FIXED) COMPLETE\n');
fprintf('=============================================\n');

fprintf('\nIMPORTANT RESULT:\n');

fprintf(['Step 23-J theta is loaded explicitly from ', ...
    'step23J.theta.\n']);

if isempty(weightsB6)

    fprintf(['No final Wiener-weight vector was automatically ', ...
        'identified in the B6 file.\n']);

    fprintf(['Therefore B6 variables must NOT be used as WDFT ', ...
        'weights.\n']);

else

    fprintf('B6 selected variable: %s\n',selectedName);

end

fprintf('\n');