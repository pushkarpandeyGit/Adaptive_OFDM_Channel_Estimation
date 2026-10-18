clear;
clc;

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' STEP 23A - PAPER WINDOWED DFT TRANSFORM\n');
fprintf('=============================================================\n');

% -------------------------------------------------------------
% Parameters
% -------------------------------------------------------------

M = 128;
L = 64;

fprintf('M                  : %d\n',M);
fprintf('L                  : %d\n',L);
fprintf('Window             : Hann\n');

% -------------------------------------------------------------
% Create a known physical CIR
% -------------------------------------------------------------

N = 256;

h_true = zeros(M,1);

h_true(1) = 1.0;

h_true(2) = 0.20 + 0.10j;

h_true(3) = 0.10 - 0.05j;

h_true(5) = 0.05 + 0.02j;

% -------------------------------------------------------------
% LTS frequency locations
% -------------------------------------------------------------

k = [-100:2:-2 2:2:100].';

% -------------------------------------------------------------
% Generate exact frequency response
% -------------------------------------------------------------

tap = (0:M-1).';

A = exp( ...
    -1j*2*pi * ...
    (k*tap.') / N );

H_true = A*h_true;

% -------------------------------------------------------------
% Run paper WDFT
% -------------------------------------------------------------

[H_est,cir_est,H_windowed,window] = ...
    windowed_dft_estimate(H_true,M,'hann',L);

% -------------------------------------------------------------
% MSE
% -------------------------------------------------------------

mse = mean(abs(H_est-H_true).^2);

% -------------------------------------------------------------
% Print
% -------------------------------------------------------------

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' RESULTS\n');
fprintf('=============================================================\n');

fprintf('Input H size       : %d x %d\n', ...
    size(H_true,1),size(H_true,2));

fprintf('Output H size      : %d x %d\n', ...
    size(H_est,1),size(H_est,2));

fprintf('CIR size            : %d x %d\n', ...
    size(cir_est,1),size(cir_est,2));

fprintf('Transform MSE       : %.12e\n',mse);

fprintf('Maximum CIR tail    : %.12e\n', ...
    max(abs(cir_est(L+1:end))));

fprintf('\n');
fprintf('First few true taps:\n');

for ii = 1:6

    fprintf('h_true(%d) = %.8f %+.8fj\n', ...
        ii,real(h_true(ii)),imag(h_true(ii)));

end

fprintf('\nFirst few estimated taps:\n');

for ii = 1:6

    fprintf('h_est(%d)  = %.8f %+.8fj\n', ...
        ii,real(cir_est(ii)),imag(cir_est(ii)));

end

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' VALIDATION\n');
fprintf('=============================================================\n');

pass = true;

if length(H_est) ~= 100
    fprintf('FAIL: output does not contain 100 tones.\n');
    pass = false;
else
    fprintf('PASS: output contains 100 tones.\n');
end

if length(cir_est) ~= 128
    fprintf('FAIL: CIR length is not 128.\n');
    pass = false;
else
    fprintf('PASS: CIR length is 128.\n');
end

if max(abs(cir_est(L+1:end))) > 1e-12
    fprintf('FAIL: CIR truncation incorrect.\n');
    pass = false;
else
    fprintf('PASS: CIR truncation correct.\n');
end

if mse < 1e-2
    fprintf('PASS: transform MSE is reasonable.\n');
else
    fprintf('FAIL: transform MSE is too large.\n');
    pass = false;
end

if pass
    fprintf('\nSTEP 23A: PASS\n');
else
    fprintf('\nSTEP 23A: FAIL\n');
end

fprintf('=============================================================\n');