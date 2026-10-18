function [cir, windowedLS, window] = ls_to_cir(H_LS, M, windowType)

% =========================================================
% LS FREQUENCY RESPONSE -> CIR
%
% IEEE 802.16 LTS:
%   100 known training tones
%
% Paper-style implementation:
%   101-point frequency grid including the DC position
%   k = -50 : 50
%
% M-point IDFT
%
% Supported windows:
%   hann
%   hamming
%   blackman
% =========================================================

H_LS = H_LS(:).';

if nargin < 2 || isempty(M)
    M = 128;
end

if nargin < 3 || isempty(windowType)
    windowType = 'hann';
end


% =========================================================
% Input validation
% =========================================================

if length(H_LS) ~= 100
    error('H_LS must contain exactly 100 LTS training tones.');
end

if M < 101 || mod(M,2) ~= 0
    error('M must be an even integer >= 102.');
end


% =========================================================
% Construct 101-position frequency grid
%
% H_LS contains the 100 non-DC training tones:
%
%   -100, -98, ..., -2
%    +2, +4, ..., +100
%
% We insert the DC position between them:
%
%   -50, ..., -1, 0, +1, ..., +50
%
% The DC value is set to zero because DC is not a
% training tone.
% =========================================================

H_101 = zeros(1,101);

% Negative-frequency half
H_101(1:50) = H_LS(1:50);

% DC interpolation (eliminates Gibbs phenomenon oscillations)
H_101(51)   = 0.5 * (H_LS(50) + H_LS(51));

% Positive-frequency half
H_101(52:101) = H_LS(51:100);


% =========================================================
% Frequency-grid index
%
% 101 positions:
%
%   n = -50 : +50
%
% =========================================================

n = -50:50;


% =========================================================
% Window
% =========================================================

switch lower(windowType)

    case 'hann'

        window = ...
            0.5 + 0.5*cos(2*pi*n/M);

    case 'hamming'

        window = ...
            0.54 + 0.46*cos(2*pi*n/M);

    case 'blackman'

        window = ...
            0.42 ...
            + 0.50*cos(2*pi*n/M) ...
            + 0.08*cos(4*pi*n/M);

    otherwise

        error('Unsupported window type: %s',windowType);

end


% =========================================================
% Apply window
% =========================================================

H_windowed_101 = H_101 .* window;


% =========================================================
% Zero padding
%
% For M = 128:
%
%   128 - 101 = 27
%
% Paper-style asymmetric padding:
%
%   14 zeros in front
%   13 zeros at back
%
% Total:
%
%   14 + 101 + 13 = 128
% =========================================================

frontZeros = (M - 101)/2 + 0.5;
backZeros  = frontZeros - 1;

frontZeros = floor(frontZeros);
backZeros  = floor(backZeros);


windowedLS = [ ...
    zeros(1,frontZeros), ...
    H_windowed_101, ...
    zeros(1,backZeros)];


% =========================================================
% Final length check
% =========================================================

if length(windowedLS) ~= M

    error(['Internal error: frequency vector length ', ...
           'does not equal M.']);

end


% =========================================================
% M-point IDFT
% =========================================================

cir = ifft(windowedLS,M);


% Verbose print disabled during simulations
% =========================================================

end