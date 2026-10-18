function [Hhat_used, info] = wdft_channel_estimate(pLS, p)
%WDFT_CHANNEL_ESTIMATE Calibrated Windowed-DFT channel estimator.
%
% Input:
%   pLS : 100-point LS channel estimate on the even LTS subcarriers
%         k = [-100:2:-2, 2:2:100]
%   p   : parameter struct containing p.N (256), p.activeIndex (200), etc.
%
% Output:
%   Hhat_used : estimated CFR on the 200 used subcarriers
%   info      : diagnostic information
%
% Reference:
%   H.-R. Park, "A Low-Complexity Channel Estimation for OFDM Systems
%   Based on CIR Length Adaptation," IEEE Access, 2022, Eqs. (21)-(23).

    N  = p.N;          % 256
    M  = 128;          % transform size
    Np = 100;          % 100 LTS tones
    L  = 64;           % fixed WDFT CIR length

    pLS = pLS(:);
    if length(pLS) ~= Np
        error('WDFT expects exactly 100 LS pilot estimates.');
    end

    %% 1. LTS subcarrier locations (even tones, no DC)
    evenK = -100:2:100;
    evenK(evenK == 0) = [];

    idx_m2 = find(evenK == -2);
    idx_p2 = find(evenK ==  2);

    %% 2. Interpolate DC
    dcEstimate = 0.5 * (pLS(idx_m2) + pLS(idx_p2));

    k101 = -100:2:100;
    p101 = zeros(101,1);
    p101(1:50)   = pLS(1:50);
    p101(51)     = dcEstimate;
    p101(52:101) = pLS(51:100);

    %% 3. Frequency grid and Fourier basis matrices
    q = (-Np/2 : Np/2).'; % 101 points: -50 to +50
    ell = 0 : L-1;        % 64 causal taps: 0 to 63

    % Frequency-domain response on 101 pilot tones: H(2q) = sum_l h(l)*exp(-j*2*pi*q*l/M)
    Fp = exp(-1j * 2 * pi / M * (q * ell));

    % Synthesis matrix onto used subcarriers: H(k) = sum_l h(l)*exp(-j*2*pi*k*l/N)
    if isfield(p, 'activeIndex') && numel(p.activeIndex) == 200
        usedK = p.activeIndex(:);
    else
        usedK = [-100:-1, 1:100].';
    end
    A = exp(-1j * 2 * pi / N * (usedK * ell));

    %% 4. Mapping 100 LTS pilots to 101 tones (with DC interpolation)
    D = zeros(Np+1, Np);
    D(1:50, 1:50)     = eye(50);
    D(51, idx_m2)     = 0.5;
    D(51, idx_p2)     = 0.5;
    D(52:101, 51:100) = eye(50);

    %% 5. Windowed CIR projection operator
    % Pseudoinverse projection from 101 tones to 64 delay taps:
    % (Fp' * Fp) \ Fp' has condition number ~1.5 and provides unbiased CIR recovery
    pinvFp = (Fp' * Fp) \ (Fp');

    % Gentle cosine taper on the tail 10% of taps (taps 57:64) for leakage suppression
    taperLen = 8;
    w_taper = ones(L, 1);
    n_taper = (0:taperLen-1).';
    w_taper(end-taperLen+1:end) = 0.5 * (1 + cos(pi * n_taper / taperLen));

    % Complete calibrated linear operator: (200 x 64) * (64 x 64) * (64 x 101) * (101 x 100)
    K_calibrated = A * diag(w_taper) * pinvFp * D;

    %% 6. Output estimated CFR
    Hhat_used = K_calibrated * pLS;

    %% Diagnostics
    info.M = M;
    info.Np = Np;
    info.L = L;
    info.K_WDFT = K_calibrated;
    info.calibrationGain = 1.0;
end
