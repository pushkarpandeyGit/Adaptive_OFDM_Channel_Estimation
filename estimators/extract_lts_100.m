function [Y_p, X_p, ltsLocations, ltsBins] = ...
    extract_lts_100(Y_LTS, X_LTS)

% =========================================================
% Extract the 100 nonzero LTS training subcarriers
% for the paper's preamble-based LS channel estimation.
%
% The original LTS is defined on even subcarrier indices
% k = -100,...,100, excluding k = 0.
%
% The LTS is generated using a 256-point IFFT and produces
% a 128-sample periodic sequence.
%
% Therefore, after a 128-point FFT, the original even
% subcarriers map to 128-point FFT bins through:
%
% q = mod(k/2,128)
%
% MATLAB bin = q + 1
% =========================================================

Y_LTS = Y_LTS(:).';
X_LTS = X_LTS(:).';

if length(Y_LTS) ~= 128
    error('Received LTS must contain exactly 128 samples.');
end

if length(X_LTS) ~= 128
    error('Known LTS FFT must contain exactly 128 bins.');
end

% Original 100 nonzero LTS subcarrier locations
ltsLocations = [-100:2:-2 2:2:100];

% Convert original 256-grid signed subcarrier index
% into the corresponding 128-point FFT bin.

q = mod(ltsLocations/2,128);

ltsBins = q + 1;

% Extract received and known LTS values
Y_p = Y_LTS(ltsBins);
X_p = X_LTS(ltsBins);

% Safety check
if any(abs(X_p) < 1e-12)
    error('A selected LTS training value is zero.');
end

fprintf('\n');
fprintf('============================================\n');
fprintf(' 100-Tone LTS Extraction\n');
fprintf('============================================\n');

fprintf('Received LTS length : %d\n',length(Y_LTS));
fprintf('Known LTS length    : %d\n',length(X_LTS));
fprintf('Number of LTS tones : %d\n',length(ltsLocations));

fprintf('First LTS location  : %d\n',ltsLocations(1));
fprintf('Last LTS location   : %d\n',ltsLocations(end));

fprintf('First MATLAB bin    : %d\n',ltsBins(1));
fprintf('Last MATLAB bin     : %d\n',ltsBins(end));

fprintf('============================================\n');

end