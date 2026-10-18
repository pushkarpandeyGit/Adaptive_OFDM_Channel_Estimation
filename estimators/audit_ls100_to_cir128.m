clc;
clear;

addpath(genpath('/MATLAB Drive/Adaptive_OFDM_Channel_Estimation'));

fprintf('\n');
fprintf('====================================================\n');
fprintf(' 100-POINT LS -> 128-POINT CIR AUDIT\n');
fprintf('====================================================\n');

%% -----------------------------------------------------
% Test 1: Create a simple 100-point LS estimate
% ------------------------------------------------------

hLS = ones(100,1);

[cir,windowedLS,window,Gmatrix] = ...
    ls100_to_cir128(hLS);

%% -----------------------------------------------------
% Test 2: Dimensions
% ------------------------------------------------------

fprintf('\n[1] Dimension test\n');

assert(length(hLS) == 100);
assert(length(window) == 100);
assert(length(windowedLS) == 100);
assert(length(cir) == 128);

fprintf('PASS: LS = 100, CIR = 128.\n');

%% -----------------------------------------------------
% Test 3: Window multiplication
% ------------------------------------------------------

fprintf('\n[2] Window multiplication test\n');

expectedWindowed = window .* hLS;

windowError = max(abs(windowedLS - expectedWindowed));

fprintf('Window multiplication error : %.3e\n',windowError);

assert(windowError < 1e-12);

fprintf('PASS: Hann window applied correctly.\n');

%% -----------------------------------------------------
% Test 4: Hann window formula
% ------------------------------------------------------

fprintf('\n[3] Hann formula test\n');

Np = 100;
M = 128;

n = (-Np/2):(Np/2-1);

expectedWindow = ...
    0.5 + 0.5*cos(2*pi*n/M);

expectedWindow = expectedWindow(:);

windowError = ...
    max(abs(window - expectedWindow));

fprintf('Hann formula error : %.3e\n',windowError);

assert(windowError < 1e-12);

fprintf('PASS: Hann window matches paper equation.\n');

%% -----------------------------------------------------
% Test 5: G matrix dimensions
% ------------------------------------------------------

fprintf('\n[4] G matrix dimension test\n');

assert(isequal(size(Gmatrix),[128 100]));

fprintf('PASS: G matrix is 128 x 100.\n');

%% -----------------------------------------------------
% Test 6: G matrix equation
% ------------------------------------------------------

fprintf('\n[5] G matrix equation test\n');

m = (0:127).';
nG = 0:99;

expectedG = (1/sqrt(128)) * ...
    exp(1j*2*pi * ...
    ((m-64) * (nG-50)) / 128);

GError = max(max(abs(Gmatrix-expectedG)));

fprintf('G matrix equation error : %.3e\n',GError);

assert(GError < 1e-12);

fprintf('PASS: Eq. (23) implemented correctly.\n');

%% -----------------------------------------------------
% Test 7: Direct CIR equation
% ------------------------------------------------------

fprintf('\n[6] Direct CIR equation test\n');

expectedCIR = Gmatrix * windowedLS;

cirError = max(abs(cir-expectedCIR));

fprintf('CIR equation error : %.3e\n',cirError);

assert(cirError < 1e-12);

fprintf('PASS: CIR calculation is correct.\n');

%% -----------------------------------------------------
% Test 8: Complex LS estimate
% ------------------------------------------------------

fprintf('\n[7] Complex-valued LS test\n');

hComplex = ...
    (1 + 1j) * ones(100,1);

[cirComplex,~,~,~] = ...
    ls100_to_cir128(hComplex);

assert(length(cirComplex) == 128);
assert(all(isfinite(cirComplex)));

fprintf('PASS: complex LS estimate handled correctly.\n');

%% -----------------------------------------------------
% Test 9: Zero input
% ------------------------------------------------------

fprintf('\n[8] Zero-input test\n');

hZero = zeros(100,1);

[cirZero,~,~,~] = ...
    ls100_to_cir128(hZero);

zeroError = max(abs(cirZero));

fprintf('Zero-input CIR error : %.3e\n',zeroError);

assert(zeroError < 1e-12);

fprintf('PASS: zero LS input produces zero CIR.\n');

%% -----------------------------------------------------
% Test 10: Invalid length
% ------------------------------------------------------

fprintf('\n[9] Invalid input test\n');

caught = false;

try

    ls100_to_cir128(ones(99,1));

catch

    caught = true;

end

assert(caught);

fprintf('PASS: incorrect LS length rejected.\n');

%% -----------------------------------------------------
% Final
% ------------------------------------------------------

fprintf('\n');
fprintf('====================================================\n');
fprintf(' ALL LS -> CIR TESTS PASSED\n');
fprintf('====================================================\n');