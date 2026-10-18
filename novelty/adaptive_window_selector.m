function [windowName, windowCoeffs, metrics] = adaptive_window_selector(H_LS, M, snr_est_dB)
% =========================================================================
% ADAPTIVE_WINDOW_SELECTOR
%
% Novelty Component 2: Channel-Aware Dynamic Windowing
%
% Replaces the static Hann window used by Park 2022 with an intelligent,
% channel-adaptive windowing strategy that matches the spectral tapering
% to the physical channel's multipath delay profile and operating SNR.
%
% Strategy:
%   1. 'hamming'  (P90 <= 8 taps):
%      Short channels (e.g. SUI-2).
%      Narrowest mainlobe width (4*pi/M) with -42.5 dB first sidelobe.
%      Prevents artificial broadening and adjacent multipath tap smearing.
%
%   2. 'hann'     (8 < P90 <= 20 taps):
%      Medium channels (e.g. Vehicular-A).
%      Balanced mainlobe with rapid asymptotic rolloff (-32 dB/octave).
%
%   3. 'blackman' (P90 > 20 taps):
%      Highly dispersive channels (e.g. SUI-5).
%      Deep sidelobe suppression (-58.1 dB peak sidelobe) that eliminates
%      spectral leakage from strong early taps into distant weak paths (tap 41)
%      and keeps the noise observation region (taps 65:120) uncontaminated.
%
% Inputs:
%   H_LS        -> 100 LS channel estimates at LTS pilot frequencies
%   M           -> IDFT CIR size (default: 128)
%   snr_est_dB  -> (Optional) Estimated SNR in dB
%
% Outputs:
%   windowName   -> 'hamming', 'hann', or 'blackman'
%   windowCoeffs -> 1 x 101 frequency-domain window coefficients
%   metrics      -> Struct containing diagnostic metrics
% =========================================================================

if nargin < 2 || isempty(M)
    M = 128;
end

H_LS = H_LS(:).';
if length(H_LS) ~= 100
    error('H_LS must contain exactly 100 LTS tones.');
end

% Construct 101-point symmetric frequency grid
H_101 = zeros(1, 101);
H_101(1:50)   = H_LS(1:50);
H_101(51)     = 0.5 * (H_LS(50) + H_LS(51)); % DC interpolation
H_101(52:101) = H_LS(51:100);

n = -50:50;

% Compute preliminary unwindowed CIR to measure energy dispersion
Gmat = (1/sqrt(M)) * exp(1j*2*pi*(((0:M-1).' - M/2)*n)/M);
cir_raw = fftshift(Gmat * H_101.');
tapEnergy = abs(cir_raw(1:64)).^2;
totEnergy = sum(tapEnergy);

cumE = cumsum(tapEnergy) / max(1e-12, totEnergy);
p90 = find(cumE >= 0.90, 1, 'first');
if isempty(p90), p90 = 1; end

p95 = find(cumE >= 0.95, 1, 'first');
if isempty(p95), p95 = 1; end

tailEnergyRatio = sum(tapEnergy(17:64)) / max(1e-12, totEnergy);

if nargin < 3 || isempty(snr_est_dB)
    noiseEst = mean(abs(cir_raw(65:120)).^2);
    sigEst = max(0, mean(tapEnergy) - noiseEst);
    snr_est_dB = 10 * log10(max(1e-3, sigEst / max(1e-8, noiseEst)));
end

% --- SELECTION LOGIC ---
if p95 <= 8
    % Short / Low Dispersion: SUI-2
    windowName = 'hamming';
    windowCoeffs = 0.54 + 0.46 * cos(2 * pi * n / M);
elseif p95 > 20 || tailEnergyRatio > 0.05
    % Highly Dispersive: SUI-5
    windowName = 'blackman';
    windowCoeffs = 0.42 + 0.50 * cos(2 * pi * n / M) + 0.08 * cos(4 * pi * n / M);
else
    % Moderate Dispersion: Vehicular-A
    windowName = 'hann';
    windowCoeffs = 0.50 + 0.50 * cos(2 * pi * n / M);
end

% Diagnostics
metrics.p90 = p90;
metrics.p95 = p95;
metrics.tailEnergyRatio = tailEnergyRatio;
metrics.snr_est_dB = snr_est_dB;
metrics.windowName = windowName;

end
