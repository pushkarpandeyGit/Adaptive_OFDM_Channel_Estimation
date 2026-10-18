function hLS_p = pilot_ls_estimate(Y_p, X_p)

% =========================================================
% Pilot-based Least-Squares Channel Estimation
%
% hLS_p(k) = Y_p(k) / X_p(k)
%
% Inputs:
%   Y_p : received pilot values
%   X_p : known transmitted pilot values
%
% Output:
%   hLS_p : LS channel estimate at pilot locations
% =========================================================

Y_p = Y_p(:).';
X_p = X_p(:).';

% ---------- Input checks ----------

if isempty(Y_p) || isempty(X_p)
    error('Pilot vectors cannot be empty.');
end

if length(Y_p) ~= length(X_p)
    error('Y_p and X_p must have the same length.');
end

if any(~isfinite(Y_p)) || any(~isfinite(X_p))
    error('Pilot values must be finite.');
end

if any(abs(X_p) < 1e-12)
    error('Known pilot values cannot contain zeros.');
end

% ---------- LS estimation ----------

hLS_p = Y_p ./ X_p;

% ---------- Display ----------

fprintf('\n');
fprintf('============================================\n');
fprintf(' PILOT-BASED LS CHANNEL ESTIMATION\n');
fprintf('============================================\n');

fprintf('Number of pilots : %d\n',length(X_p));
fprintf('Y_p dimensions   : %d x %d\n', ...
    size(Y_p,1),size(Y_p,2));

fprintf('X_p dimensions   : %d x %d\n', ...
    size(X_p,1),size(X_p,2));

fprintf('hLS_p dimensions : %d x %d\n', ...
    size(hLS_p,1),size(hLS_p,2));

fprintf('============================================\n');

end