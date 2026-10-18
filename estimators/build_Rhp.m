function Rhp = build_Rhp(R, usedLocations, pilotLocations, N)

% =========================================================
% BUILD CHANNEL-TO-PILOT COVARIANCE MATRIX R_hp
%
% Rhp(i,j) = covariance between:
%
%   used subcarrier i
%
% and
%
%   pilot subcarrier j
%
% =========================================================


%% ---------------------------------------------------------
% Input validation
% ----------------------------------------------------------

if nargin < 4
    error('R, usedLocations, pilotLocations, and N are required.');
end

if ~ismatrix(R)
    error('R must be a matrix.');
end

if size(R,1) ~= N || size(R,2) ~= N
    error('R must be an N x N matrix.');
end

usedLocations = usedLocations(:).';
pilotLocations = pilotLocations(:).';

if isempty(usedLocations)
    error('usedLocations cannot be empty.');
end

if isempty(pilotLocations)
    error('pilotLocations cannot be empty.');
end

if any(usedLocations < -N/2) || ...
        any(usedLocations > N/2-1)

    error('Used subcarrier locations are outside the valid FFT range.');
end

if any(pilotLocations < -N/2) || ...
        any(pilotLocations > N/2-1)

    error('Pilot locations are outside the valid FFT range.');
end

if length(unique(usedLocations)) ~= length(usedLocations)
    error('Used subcarrier locations must be unique.');
end

if length(unique(pilotLocations)) ~= length(pilotLocations)
    error('Pilot locations must be unique.');
end


%% ---------------------------------------------------------
% Convert frequency locations to MATLAB indices
%
% k = -N/2,...,-1,0,1,...,N/2-1
%
% MATLAB index = k + N/2 + 1
% ----------------------------------------------------------

usedIndex = usedLocations + N/2 + 1;
pilotIndex = pilotLocations + N/2 + 1;

usedIndex = round(usedIndex);
pilotIndex = round(pilotIndex);


%% ---------------------------------------------------------
% Extract R_hp
% ----------------------------------------------------------

Rhp = R(usedIndex,pilotIndex);


%% ---------------------------------------------------------
% Information
% ----------------------------------------------------------

fprintf('\n');
fprintf('============================================\n');
fprintf(' Channel-Pilot Covariance Matrix R_hp\n');
fprintf('============================================\n');

fprintf('N                  : %d\n',N);

fprintf('Used subcarriers   : %d\n', ...
    length(usedLocations));

fprintf('Pilot subcarriers  : %d\n', ...
    length(pilotLocations));

fprintf('Matrix dimensions  : %d x %d\n', ...
    size(Rhp,1),size(Rhp,2));

fprintf('============================================\n');

end