function H_SVD = svd_channel_estimate( ...
    Rhp, U, Dp, hLS_p)

% =========================================================
% Pilot-based SVD channel estimation
%
% H_SVD = R_hp * U * D_p * U^H * h_LS,p
% =========================================================

% ---------- Force column vectors where appropriate ----------

hLS_p = hLS_p(:);

% ---------- Input checks ----------

if ~ismatrix(Rhp)
    error('Rhp must be a matrix.');
end

[Nu,Np] = size(Rhp);

if ~ismatrix(U)
    error('U must be a matrix.');
end

if size(U,1) ~= Np || size(U,2) ~= Np
    error('U dimensions are inconsistent with Rhp.');
end

if ~ismatrix(Dp)
    error('Dp must be a matrix.');
end

if size(Dp,1) ~= Np || size(Dp,2) ~= Np
    error('Dp dimensions are inconsistent with Rhp.');
end

if length(hLS_p) ~= Np
    error('hLS_p length must equal the number of pilots.');
end

% ---------- SVD estimator ----------

H_SVD = Rhp * U * Dp * U' * hLS_p;

% ---------- Display ----------

fprintf('\n');
fprintf('============================================\n');
fprintf(' PILOT-BASED SVD CHANNEL ESTIMATION\n');
fprintf('============================================\n');

fprintf('R_hp dimensions       : %d x %d\n', ...
    size(Rhp,1),size(Rhp,2));

fprintf('U dimensions          : %d x %d\n', ...
    size(U,1),size(U,2));

fprintf('D_p dimensions        : %d x %d\n', ...
    size(Dp,1),size(Dp,2));

fprintf('LS pilot estimate     : %d x 1\n', ...
    length(hLS_p));

fprintf('Output channel        : %d x 1\n', ...
    length(H_SVD));

fprintf('============================================\n');

end