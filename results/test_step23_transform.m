clear;
clc;

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' STEP 23A - WINDOWED DFT TRANSFORM SANITY TEST\n');
fprintf('=============================================================\n');

M = 128;
L = 64;

k = [-100:2:-2 2:2:100].';

% -------------------------------------------------------------
% Create a known CIR
% -------------------------------------------------------------

h_true = zeros(M,1);

h_true(1) = 1.0;

h_true(2) = 0.20 + 0.10j;

h_true(3) = 0.10 - 0.05j;

h_true(5) = 0.05 + 0.02j;

% -------------------------------------------------------------
% Generate frequency response at the 100 LTS tones
% -------------------------------------------------------------

H_true = zeros(100,1);

for ii = 1:100

    H_true(ii) = ...
        sum(h_true .* ...
        exp(-1j*2*pi*k(ii)*(0:M-1).'/256));

end

% -------------------------------------------------------------
% Add very small noise
% -------------------------------------------------------------

rng(42);

noisePower = 1e-10;

H_LS = H_true + ...
    sqrt(noisePower/2) * ...
    (randn(100,1) + 1j*randn(100,1));

% -------------------------------------------------------------
% Run windowed DFT
% -------------------------------------------------------------

[H_WDFT,cir_WDFT,~,window] = ...
    windowed_dft_estimate(H_LS,M,'hann',L);

% -------------------------------------------------------------
% MSE
% -------------------------------------------------------------

mse = mean(abs(H_WDFT(:)-H_true).^2);

fprintf('\n');
fprintf('M                 : %d\n',M);
fprintf('CIR length L      : %d\n',L);
fprintf('Number of tones   : %d\n',length(H_LS));
fprintf('Window             : Hann\n');

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' RESULTS\n');
fprintf('=============================================================\n');

fprintf('Input H size       : %d x %d\n', ...
    size(H_LS,1),size(H_LS,2));

fprintf('Output H size      : %d x %d\n', ...
    size(H_WDFT,1),size(H_WDFT,2));

fprintf('CIR size            : %d x %d\n', ...
    size(cir_WDFT,1),size(cir_WDFT,2));

fprintf('Transform MSE       : %.12e\n',mse);

fprintf('Maximum CIR tail    : %.12e\n', ...
    max(abs(cir_WDFT(L+1:end))));

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' VALIDATION\n');
fprintf('=============================================================\n');

pass = true;

if length(H_WDFT) ~= 100
    fprintf('FAIL: output does not contain 100 tones.\n');
    pass = false;
else
    fprintf('PASS: output contains 100 tones.\n');
end

if length(cir_WDFT) ~= 128
    fprintf('FAIL: CIR length is not 128.\n');
    pass = false;
else
    fprintf('PASS: CIR length is 128.\n');
end

if max(abs(cir_WDFT(L+1:end))) > 1e-12
    fprintf('FAIL: CIR truncation incorrect.\n');
    pass = false;
else
    fprintf('PASS: CIR truncation correct.\n');
end

if mse < 1e-2
    fprintf('PASS: transform produces reasonable channel estimate.\n');
else
    fprintf('FAIL: transform error is too large.\n');
    pass = false;
end

if pass
    fprintf('\nSTEP 23A: PASS\n');
else
    fprintf('\nSTEP 23A: FAIL\n');
end

fprintf('=============================================================\n');