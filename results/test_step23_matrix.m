clear;
clc;

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' STEP 23A - PAPER WDFT MATRIX TEST\n');
fprintf('=============================================================\n');

Np = 100;
M  = 128;

% -------------------------------------------------------------
% Construct 101-point frequency vector
% -------------------------------------------------------------

p = ones(Np+1,1);

% -------------------------------------------------------------
% Frequency index
% -------------------------------------------------------------

n = (-Np/2:Np/2).';

fprintf('\nNumber of frequency samples : %d\n',length(p));
fprintf('Number of n values          : %d\n',length(n));

% -------------------------------------------------------------
% Hann window
% -------------------------------------------------------------

w = 0.5 + 0.5*cos(2*pi*n/M);

fprintf('Window length               : %d\n',length(w));

% -------------------------------------------------------------
% Apply window
% -------------------------------------------------------------

p_windowed = p .* w;

% -------------------------------------------------------------
% Zero padding
% -------------------------------------------------------------

frontZeros = (M-Np)/2;
backZeros  = (M-Np)/2 - 1;

p_padded = [ ...
    zeros(frontZeros,1);
    p_windowed;
    zeros(backZeros,1)
];

fprintf('Front zeros                 : %d\n',frontZeros);
fprintf('Back zeros                  : %d\n',backZeros);
fprintf('Padded vector length        : %d\n',length(p_padded));

% -------------------------------------------------------------
% Paper transformation matrix
% -------------------------------------------------------------

m = (0:M-1).';
nn = (0:Np).';

G = (1/sqrt(M)) * ...
    exp( ...
    1j*2*pi * ...
    ((m-M/2)*(nn-Np/2).')/M ...
    );

fprintf('G dimensions               : %d x %d\n', ...
    size(G,1),size(G,2));

% -------------------------------------------------------------
% Apply transformation
% -------------------------------------------------------------

g = G * p_windowed;

fprintf('Output CIR dimensions       : %d x %d\n', ...
    size(g,1),size(g,2));

% -------------------------------------------------------------
% Checks
% -------------------------------------------------------------

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' VALIDATION\n');
fprintf('=============================================================\n');

pass = true;

if length(w) == 101
    fprintf('PASS: Hann window has 101 coefficients.\n');
else
    fprintf('FAIL: Window length incorrect.\n');
    pass = false;
end

if length(p_padded) == 128
    fprintf('PASS: zero-padded vector has 128 samples.\n');
else
    fprintf('FAIL: Padding length incorrect.\n');
    pass = false;
end

if size(G,1) == 128 && size(G,2) == 101
    fprintf('PASS: G is 128 x 101.\n');
else
    fprintf('FAIL: G dimensions incorrect.\n');
    pass = false;
end

if length(g) == 128
    fprintf('PASS: CIR output has 128 samples.\n');
else
    fprintf('FAIL: CIR output length incorrect.\n');
    pass = false;
end

if pass
    fprintf('\nSTEP 23A MATRIX TEST: PASS\n');
else
    fprintf('\nSTEP 23A MATRIX TEST: FAIL\n');
end

fprintf('=============================================================\n');