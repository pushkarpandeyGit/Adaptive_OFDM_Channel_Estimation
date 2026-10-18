function [cir, windowedLS, window, Gmatrix] = ...
    ls100_to_cir128(hLS_100)

% =========================================================
% Convert the 100-point LTS LS channel estimate into the
% 128-point CIR estimate used by the paper.
%
% Paper:
%   Np = 100 LTS training subcarriers
%   M  = 128 CIR/IDFT size
%
% Steps:
%   1. Construct Hann window
%   2. Apply window to 100-point LS estimate
%   3. Construct the centered transformation matrix G
%   4. Obtain 128-point CIR estimate
%
% Based on equations (21)-(23) of the paper.
% =========================================================

hLS_100 = hLS_100(:);

Np = 100;
M  = 128;

if length(hLS_100) ~= Np
    error('Input LS estimate must contain exactly 100 samples.');
end

% ---------------------------------------------------------
% Step 1: Hann window
% ---------------------------------------------------------

% 100 samples centered around zero.
n = (-Np/2):(Np/2-1);

window = 0.5 + 0.5*cos(2*pi*n/M);

window = window(:);

% ---------------------------------------------------------
% Step 2: Apply window
% ---------------------------------------------------------

windowedLS = window .* hLS_100;

% ---------------------------------------------------------
% Step 3: Construct G matrix from Eq. (23)
%
% [G]_(m,n) =
% 1/sqrt(M) *
% exp(j*2*pi*(m-M/2)*(n-Np/2)/M)
%
% MATLAB indices are converted to the paper's
% zero-based indices.
% ---------------------------------------------------------

m = (0:M-1).';

nG = 0:(Np-1);

Gmatrix = (1/sqrt(M)) * ...
    exp(1j*2*pi * ...
    ((m-M/2) * (nG-Np/2)) / M);

% ---------------------------------------------------------
% Step 4: CIR estimate
% ---------------------------------------------------------

cir = Gmatrix * windowedLS;

% ---------------------------------------------------------
% Display
% ---------------------------------------------------------

fprintf('\n');
fprintf('============================================\n');
fprintf(' 100-point LS -> 128-point CIR\n');
fprintf('============================================\n');

fprintf('LS estimate length : %d\n',length(hLS_100));
fprintf('Window length      : %d\n',length(window));
fprintf('CIR length         : %d\n',length(cir));
fprintf('G matrix size      : %d x %d\n', ...
    size(Gmatrix,1),size(Gmatrix,2));

fprintf('Window minimum     : %.6f\n',min(window));
fprintf('Window maximum     : %.6f\n',max(window));

fprintf('============================================\n');

end