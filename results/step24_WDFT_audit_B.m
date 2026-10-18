%% ============================================================
% STEP 24 - WDFT AUDIT B
% Inspect fields inside the Step 23-J result structure
% ============================================================

clear;
clc;

fprintf('\n=============================================\n');
fprintf('STEP 24 - WDFT AUDIT B\n');
fprintf('Inspecting fields inside step23J\n');
fprintf('=============================================\n\n');

fileName = 'step23J_actual_wdft_audit.mat';

if ~isfile(fileName)
    error('Cannot find step23J_actual_wdft_audit.mat');
end

S = load(fileName);

% ------------------------------------------------------------
% Check top-level variable
% ------------------------------------------------------------

if ~isfield(S,'step23J')
    error('Variable "step23J" was not found.');
end

step23J = S.step23J;

fprintf('Top-level variable:\n');
fprintf('  Name  : step23J\n');
fprintf('  Class : %s\n',class(step23J));
fprintf('  Size  : ');

disp(size(step23J));

% ------------------------------------------------------------
% List fields
% ------------------------------------------------------------

fields = fieldnames(step23J);

fprintf('\nFields inside step23J:\n');
fprintf('------------------------------------------------------------\n');

for ii = 1:length(fields)

    name = fields{ii};
    value = step23J.(name);

    fprintf('\n[%d] %s\n',ii,name);
    fprintf('    Class : %s\n',class(value));
    fprintf('    Size  : ');

    sz = size(value);

    fprintf('%d',sz(1));

    for jj = 2:length(sz)
        fprintf(' x %d',sz(jj));
    end

    fprintf('\n');

    % Numeric information
    if isnumeric(value) && ~isempty(value)

        fprintf('    Numeric: YES\n');
        fprintf('    Min abs : %.12e\n',min(abs(value(:))));
        fprintf('    Max abs : %.12e\n',max(abs(value(:))));
        fprintf('    Mean abs: %.12e\n',mean(abs(value(:))));
        fprintf('    Norm    : %.12e\n',norm(value(:)));

        if isvector(value)

            fprintf('    First values:\n');

            nShow = min(10,length(value));

            for kk = 1:nShow

                fprintf('      (%d) ',kk);

                if isreal(value(kk))
                    fprintf('%.12e\n',value(kk));
                else
                    fprintf('%.12e %+.12ej\n', ...
                        real(value(kk)), ...
                        imag(value(kk)));
                end

            end

        end

    elseif isstruct(value)

        fprintf('    Struct: YES\n');

        subfields = fieldnames(value);

        fprintf('    Subfields:\n');

        for jj = 1:length(subfields)
            fprintf('      - %s\n',subfields{jj});
        end

    else

        fprintf('    Numeric: NO\n');

        if ischar(value) || isstring(value)
            fprintf('    Value: %s\n',string(value));
        end

    end

end

fprintf('\n------------------------------------------------------------\n');

% ------------------------------------------------------------
% Also recursively inspect one level deeper for structs
% ------------------------------------------------------------

fprintf('\nSTRUCT SUBFIELD DETAILS\n');
fprintf('------------------------------------------------------------\n');

for ii = 1:length(fields)

    name = fields{ii};
    value = step23J.(name);

    if isstruct(value)

        subfields = fieldnames(value);

        fprintf('\nStruct field: %s\n',name);

        for jj = 1:length(subfields)

            subname = subfields{jj};
            subvalue = value.(subname);

            fprintf('  %s\n',subname);
            fprintf('    Class: %s\n',class(subvalue));
            fprintf('    Size : ');

            sz = size(subvalue);

            fprintf('%d',sz(1));

            for kk = 2:length(sz)
                fprintf(' x %d',sz(kk));
            end

            fprintf('\n');

            if isnumeric(subvalue) && ~isempty(subvalue)

                fprintf('    Min abs : %.12e\n', ...
                    min(abs(subvalue(:))));

                fprintf('    Max abs : %.12e\n', ...
                    max(abs(subvalue(:))));

                fprintf('    Mean abs: %.12e\n', ...
                    mean(abs(subvalue(:))));

            end

        end

    end

end

fprintf('\n=============================================\n');
fprintf('AUDIT B COMPLETE\n');
fprintf('=============================================\n');