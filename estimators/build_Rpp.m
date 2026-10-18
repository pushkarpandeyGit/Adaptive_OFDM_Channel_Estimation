function Rpp = build_Rpp(R, pilotLocations, N)

% =========================================================
% BUILD PILOT-PILOT CHANNEL COVARIANCE MATRIX
%
% Rpp(i,j) = channel covariance between pilot subcarriers
%            k_i and k_j
%
% R is the full N x N uniform channel covariance matrix.
%
% =========================================================


%% -----------------------------------------------------
% Input validation
% ------------------------------------------------------

if nargin < 3
    error('R, pilotLocations, and N are required.');
end

if ~ismatrix(R)
    error('R must be a matrix.');
end

if size(R,1) ~= N || size(R,2) ~= N
    error('R must be an N x N matrix.');
end

pilotLocations = pilotLocations(:).';

if isempty(pilotLocations)
    error('pilotLocations cannot be empty.');
end

if any(pilotLocations < -N/2) || ...
        any(pilotLocations > N/2-1)

    error('Pilot locations are outside the valid FFT range.');

end

if length(unique(pilotLocations)) ~= length(pilotLocations)
    error('Pilot locations must be unique.');
end


%% -----------------------------------------------------
% Convert subcarrier locations to MATLAB indices
%
% FFT bins:
%
% k = -N/2,...,-1,0,1,...,N/2-1
%
% MATLAB index:
%
% index = k + N/2 + 1
% ------------------------------------------------------

pilotIndex = pilotLocations + N/2 + 1;

pilotIndex = round(pilotIndex);


%% -----------------------------------------------------
% Extract pilot-pilot covariance
% ------------------------------------------------------

Rpp = R(pilotIndex,pilotIndex);


%% -----------------------------------------------------
% Information
% ------------------------------------------------------

fprintf('\n');
fprintf('============================================\n');
fprintf(' Pilot-Pilot Covariance Matrix R_pp\n');
fprintf('============================================\n');

fprintf('N                  : %d\n',N);
fprintf('Number of pilots   : %d\n',length(pilotLocations));

fprintf('Pilot locations    : ');

fprintf('%d ',pilotLocations);

fprintf('\n');

fprintf('Matrix dimensions  : %d x %d\n', ...
    size(Rpp,1),size(Rpp,2));

fprintf('============================================\n');

end