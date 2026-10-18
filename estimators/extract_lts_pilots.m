function [Y_p, X_p, pilotBins] = extract_lts_pilots( ...
    Y_LTS, X_LTS, pilotLocations, N)

% =========================================================
% Extract IEEE 802.16 pilot subcarriers from LTS spectrum
%
% Inputs:
%   Y_LTS          : received LTS frequency-domain vector
%   X_LTS          : known LTS frequency-domain vector
%   pilotLocations : pilot locations in signed subcarrier
%                    notation
%   N              : OFDM FFT size
%
% Outputs:
%   Y_p            : received pilot values
%   X_p            : known pilot values
%   pilotBins      : MATLAB indices corresponding to pilots
%
% =========================================================

    % ---------- Force row vectors ----------

    Y_LTS = Y_LTS(:).';
    X_LTS = X_LTS(:).';
    pilotLocations = pilotLocations(:).';

    % ---------- Input checks ----------

    if length(Y_LTS) ~= N
        error('Y_LTS must contain N frequency bins.');
    end

    if length(X_LTS) ~= N
        error('X_LTS must contain N frequency bins.');
    end

    if isempty(pilotLocations)
        error('pilotLocations cannot be empty.');
    end

    if any(pilotLocations < -N/2) || ...
            any(pilotLocations > N/2-1)

        error('Pilot locations are outside the valid FFT range.');

    end

    if length(unique(pilotLocations)) ~= ...
            length(pilotLocations)

        error('Duplicate pilot locations detected.');

    end

    % ---------- Convert signed subcarrier locations
    %            to MATLAB FFT indices ----------
    %
    % Signed location:
    %
    %       k = -N/2,...,-1,0,1,...,N/2-1
    %
    % MATLAB index:
    %
    %       bin = k + N/2 + 1
    %
    % Examples for N = 256:
    %
    %       k = -128 -> bin 1
    %       k =    0 -> bin 129
    %       k =  127 -> bin 256

    pilotBins = pilotLocations + N/2 + 1;
    pilotBins = round(pilotBins);

    % ---------- Extract pilots ----------

    Y_p = Y_LTS(pilotBins);
    X_p = X_LTS(pilotBins);

    % ---------- Check that known pilots are non-zero ----------

    if any(abs(X_p) < 1e-12)
        error('At least one selected pilot has zero known LTS value.');
    end

    % ---------- Display ----------

    fprintf('\n');
    fprintf('============================================\n');
    fprintf(' LTS PILOT EXTRACTION\n');
    fprintf('============================================\n');

    fprintf('FFT size            : %d\n',N);
    fprintf('Number of pilots    : %d\n',length(pilotLocations));

    fprintf('\nPilot mapping:\n');

    for k = 1:length(pilotLocations)

        fprintf('Location %4d -> MATLAB bin %3d\n', ...
            pilotLocations(k), ...
            pilotBins(k));

    end

    fprintf('\nOutput dimensions:\n');
    fprintf('Y_p                 : %d x %d\n', ...
        size(Y_p,1),size(Y_p,2));

    fprintf('X_p                 : %d x %d\n', ...
        size(X_p,1),size(X_p,2));

    fprintf('============================================\n');

end