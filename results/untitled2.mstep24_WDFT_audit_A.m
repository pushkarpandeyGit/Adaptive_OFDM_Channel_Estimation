%% STEP 24 - WDFT AUDIT A
% Inspect the variables saved by the successful Step 23-J audit.
% DO NOT modify estimator code yet.

clear;
clc;

fprintf('\n=============================================\n');
fprintf('STEP 24 - WDFT AUDIT A\n');
fprintf('Inspecting Step 23-J saved variables\n');
fprintf('=============================================\n\n');

fileName = 'step23J_actual_wdft_audit.mat';

if ~isfile(fileName)
    error(['Cannot find ', fileName, ...
        '. Make sure the Step 23-J audit MAT file is in the current MATLAB folder.']);
end

fprintf('Found file: %s\n\n', fileName);

% Show all variables stored in the successful audit
info = whos('-file', fileName);

fprintf('Variables stored in Step 23-J MAT file:\n');
fprintf('---------------------------------------------\n');

for n = 1:length(info)

    fprintf('%-30s  %-12s  ', ...
        info(n).name, info(n).class);

    fprintf('size = ');

    sz = info(n).size;

    fprintf('%d', sz(1));

    for k = 2:length(sz)
        fprintf(' x %d', sz(k));
    end

    fprintf('\n');
end

fprintf('---------------------------------------------\n\n');

% Load the file
S = load(fileName);

names = fieldnames(S);

fprintf('Numeric variables and basic statistics:\n');
fprintf('---------------------------------------------\n');

for n = 1:length(names)

    name = names{n};
    value = S.(name);

    if isnumeric(value)

        fprintf('\nVariable: %s\n', name);

        fprintf('  Size      : ');
        fprintf('%d', size(value,1));

        if length(size(value)) > 1
            for k = 2:length(size(value))
                fprintf(' x %d', size(value,k));
            end
        end
        fprintf('\n');

        fprintf('  Class     : %s\n', class(value));

        if ~isempty(value)

            fprintf('  Min       : %.12g\n', min(abs(value(:))));
            fprintf('  Max       : %.12g\n', max(abs(value(:))));
            fprintf('  Mean abs  : %.12g\n', mean(abs(value(:))));
            fprintf('  Norm      : %.12g\n', norm(value(:)));

            if isvector(value)

                fprintf('  First 5 values:\n');
                disp(value(1:min(5,end)));

            end
        end
    end
end

fprintf('\n=============================================\n');
fprintf('AUDIT A COMPLETE\n');
fprintf('=============================================\n');