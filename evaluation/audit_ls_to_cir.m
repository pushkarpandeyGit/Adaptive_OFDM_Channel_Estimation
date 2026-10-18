function [cir, windowedLS, window] = ls_to_cir(H_LS, M, windowType)

% =========================================================
% 100-Tone LS Channel Estimate -> Time-Domain CIR
% =========================================================
%
% Input:
%   H_LS       -> 100-tone LS channel estimate
%                 ordered as:
%                 -100,-98,...,-2,+2,...,+98,+100
%
%   M          -> CIR IDFT size
%                 paper quasi-static case: M = 128
%
%   windowType -> 'hann', 'hamming', or 'blackman'
%
% Output:
%   cir        -> M-point time-domain CIR
%   windowedLS -> zero-padded, windowed M-point frequency vector
%   window    -> window applied to the 101 frequency positions
%
% =========================================================
% PAPER-STYLE PROCESS
%
% 100 nonzero LTS tones
%        |
%        v
% 101 positions including DC
%        |
%        v
% apply selected window
%        |
%        v
% 14 zeros + 101 positions + 13 zeros
%        |
%        v
% M-point IDFT
%
% =========================================================
% NOVELTY COMPATIBILITY
%
% windowType remains an input.
%
% Therefore:
%
%   'hann'     -> paper baseline
%   'hamming'  -> available for adaptive-window novelty
%   'blackman' -> available for adaptive-window novelty
%
% Later, the adaptive-window selector can decide which
% windowType to pass into this function.
%
% ML threshold is NOT placed here.
% ML threshold belongs to the SHT stage.
%
% =========================================================

%% ---------------------------------------------------------
% 1. INPUT FORMAT
% ---------------------------------------------------------

H_LS = H_LS(:).';

if nargin < 2 || isempty(M)

    M = 128;

end

if nargin < 3 || isempty(windowType)

    windowType = 'hann';

end

%% ---------------------------------------------------------
% 2. CHECK INPUT
% ---------------------------------------------------------

if length(H_LS) ~= 100

    error(['H_LS must contain exactly 100 LTS training ', ...
           'tones for the paper-based CIR estimator.']);

end

if M < 101 || mod(M,2) ~= 0

    error('M must be an even integer greater than or equal to 101.');

end

%% ---------------------------------------------------------
% 3. BUILD 101 FREQUENCY POSITIONS
% ---------------------------------------------------------
%
% The 100 measured LTS tones are:
%
% -100,-98,...,-2,+2,...,+98,+100
%
% The paper formulation includes the DC position:
%
% -100,-98,...,-2,0,+2,...,+98,+100
%
% Therefore:
%
% 100 measured tones + 1 DC position = 101 positions.
%
% DC is zero.

H101 = zeros(1,101);

H101(1:50) = H_LS(1:50);

H101(51) = 0;

H101(52:101) = H_LS(51:100);

%% ---------------------------------------------------------
% 4. BUILD SELECTED WINDOW
% ---------------------------------------------------------
%
% Window coordinate:
%
% n = -50,...,0,...,+50
%
% because Np = 100.
%
% Paper Hann:
%
% w[n] = 0.5 + 0.5*cos(2*pi*n/M)
%
% Our novelty-compatible implementation keeps the same
% coordinate system for Hann, Hamming, and Blackman.

n = -50:50;

switch lower(windowType)

    case 'hann'

        window = 0.5 + ...
                 0.5*cos(2*pi*n/M);

    case 'hamming'

        window = 0.54 + ...
                 0.46*cos(2*pi*n/M);

    case 'blackman'

        window = 0.42 ...
               + 0.50*cos(2*pi*n/M) ...
               + 0.08*cos(4*pi*n/M);

    otherwise

        error('Unsupported window type: %s', windowType);

end

%% ---------------------------------------------------------
% 5. APPLY WINDOW TO 101 POSITIONS
% ---------------------------------------------------------

H101_windowed = H101 .* window;

%% ---------------------------------------------------------
% 6. ZERO PADDING
% ---------------------------------------------------------
%
% Paper:
%
% front zeros = (M-Np)/2
%
% back zeros = (M-Np)/2 - 1
%
% Here Np = 100.
%
% For M = 128:
%
% front = 14
% back  = 13
%
% Therefore:
%
% 14 + 101 + 13 = 128

frontZeros = (M - 100)/2;

backZeros = frontZeros - 1;

if frontZeros ~= floor(frontZeros)

    error('M does not produce integer paper-style zero padding.');

end

windowedLS = [ ...
    zeros(1,frontZeros), ...
    H101_windowed, ...
    zeros(1,backZeros) ...
];

%% ---------------------------------------------------------
% 7. VERIFY FINAL VECTOR LENGTH
% ---------------------------------------------------------

if length(windowedLS) ~= M

    error('Internal error: zero-padded CIR vector has wrong length.');

end

%% ---------------------------------------------------------
% 8. M-POINT IDFT
% ---------------------------------------------------------

cir = ifft(windowedLS, M);

%% ---------------------------------------------------------
% 9. DISPLAY INFORMATION
% ---------------------------------------------------------

fprintf('\n');
fprintf('============================================\n');
fprintf(' 100-Tone LS Estimate -> CIR\n');
fprintf('============================================\n');

fprintf('Input LS tones       : %d\n', length(H_LS));
fprintf('Frequency positions  : %d\n', length(H101));
fprintf('DC value             : %.4f\n', H101(51));
fprintf('Window               : %s\n', windowType);
fprintf('Front zero padding   : %d\n', frontZeros);
fprintf('Back zero padding    : %d\n', backZeros);
fprintf('IDFT size             : %d\n', M);
fprintf('CIR length             : %d samples\n', length(cir));

fprintf('============================================\n');

end