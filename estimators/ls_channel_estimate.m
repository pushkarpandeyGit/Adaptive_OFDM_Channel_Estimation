function H_LS = ls_channel_estimate(Y, X)

% =========================================================
% Least-Squares (LS) Channel Estimation
% =========================================================
%
% Input:
%   Y -> received LTS in frequency domain
%   X -> known transmitted LTS in frequency domain
%
% Output:
%   H_LS -> LS channel estimate
%
% Mathematical equation:
%
%              Y[k]
%   H_LS[k] = -------
%              X[k]
%
% =========================================================

Y = Y(:).';
X = X(:).';

if length(Y) ~= length(X)
    error('Y and X must have the same length.');
end

% -----------------------------------------------------
% LS channel estimate
% -----------------------------------------------------

H_LS = Y ./ X;

% -----------------------------------------------------
% Display information
% -----------------------------------------------------

fprintf('\n');
fprintf('============================================\n');
fprintf(' LS Channel Estimation\n');
fprintf('============================================\n');

fprintf('Received LTS bins : %d\n', length(Y));
fprintf('Known LTS bins    : %d\n', length(X));
fprintf('Estimated H bins  : %d\n', length(H_LS));

fprintf('============================================\n');

end