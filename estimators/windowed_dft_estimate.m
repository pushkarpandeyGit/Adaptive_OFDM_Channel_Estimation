function [H_WDFT, cir_WDFT, H_windowed, window] = ...
    windowed_dft_estimate(H_LS, M, windowType, L)

% ============================================================
% PAPER-BASED WINDOWED DFT CHANNEL ESTIMATOR
%
% Based on Park, IEEE Access 2022, Eqs. (21)-(23)
%
% H_LS:
%   100 LTS LS channel estimates
%
% M:
%   IDFT size = 128
%
% L:
%   assumed CIR length = 64 for paper baseline
%
% ============================================================

    % --------------------------------------------------------
    % Inputs
    % --------------------------------------------------------

    H_LS = H_LS(:);

    if nargin < 2 || isempty(M)
        M = 128;
    end

    if nargin < 3 || isempty(windowType)
        windowType = 'hann';
    end

    if nargin < 4 || isempty(L)
        L = 64;
    end

    Np = 100;

    % --------------------------------------------------------
    % Validation
    % --------------------------------------------------------

    if length(H_LS) ~= Np
        error('H_LS must contain exactly 100 LTS estimates.');
    end

    if M ~= 128
        error('For this paper reproduction M must be 128.');
    end

    if L < 1 || L > M
        error('Invalid CIR length L.');
    end

    % --------------------------------------------------------
    % Step 1:
    % Insert the missing DC estimate.
    %
    % Paper:
    %
    % H_LS(0) =
    % [H_LS(-2) + H_LS(+2)] / 2
    %
    % --------------------------------------------------------

    H_101 = zeros(101,1);

    % -100 ... -2
    H_101(1:50) = H_LS(1:50);

    % DC
    H_101(51) = ...
        0.5 * (H_LS(50) + H_LS(51));

    % +2 ... +100
    H_101(52:101) = H_LS(51:100);

    % --------------------------------------------------------
    % Step 2:
    % Construct window
    %
    % Eq. (21):
    %
    % w_n = 0.5 + 0.5 cos(2*pi*n/M)
    %
    % |n| <= Np/2
    %
    % --------------------------------------------------------

    n = (-Np/2:Np/2).';

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

            error('Unknown window type.');

    end

    % --------------------------------------------------------
    % Step 3:
    % Apply window
    % --------------------------------------------------------

    H_windowed = H_101 .* window;

    % --------------------------------------------------------
    % Step 4:
    % Zero padding
    %
    % Paper:
    %
    % (M-Np)/2 zeros at front
    %
    % (M-Np)/2 - 1 zeros at back
    %
    % --------------------------------------------------------

    frontZeros = (M-Np)/2;
    backZeros  = (M-Np)/2 - 1;

    p_windowed = [ ...
        zeros(frontZeros,1);
        H_windowed;
        zeros(backZeros,1)
    ];

    if length(p_windowed) ~= M
        error('Zero-padded vector does not have length M.');
    end

    % --------------------------------------------------------
    % Step 5:
    % Exact transformation matrix from Eq. (23)
    %
    % G(m,n) =
    % 1/sqrt(M) *
    % exp(j*2*pi*(m-M/2)*(n-Np/2)/M)
    %
    % MATLAB indices:
    % m = 0,...,M-1
    % n = 0,...,Np
    % --------------------------------------------------------

    m = (0:M-1).';
    nn = (0:Np).';

    G = (1/sqrt(M)) * ...
        exp( ...
        1j*2*pi * ...
        ((m-M/2) * (nn-Np/2).') / M );

    % --------------------------------------------------------
    % Step 6:
    % CIR estimate
    %
    % Eq. (22)
    % --------------------------------------------------------

    cir_WDFT = G * H_windowed;

    cir_WDFT = cir_WDFT(:);

    % --------------------------------------------------------
    % Step 7:
    % CIR truncation
    %
    % For the baseline:
    %
    % L = 64
    %
    % --------------------------------------------------------

    cir_WDFT(L+1:M) = 0;

    % --------------------------------------------------------
    % Step 8:
    % Reconstruct channel at the original 100 LTS tones
    %
    % This uses the same physical subcarrier locations.
    %
    % --------------------------------------------------------

    k100 = [-100:2:-2 2:2:100].';

    tap = (0:M-1).';

    A100 = exp( ...
        -1j*2*pi * ...
        (k100 * tap.') / 256 );

    H_WDFT = A100 * cir_WDFT;

    H_WDFT = H_WDFT(:);

end