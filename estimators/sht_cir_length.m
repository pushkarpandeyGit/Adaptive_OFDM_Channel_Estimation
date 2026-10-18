function [L_hat, threshold, tapEnergy, decision] = ...
    sht_cir_length(cir, noiseVariance, T, G)

% =========================================================
% Sequential Hypothesis Test (SHT) for CIR Length
% =========================================================
%
% Input:
%   cir            -> estimated time-domain CIR
%   noiseVariance  -> estimated noise variance
%   T              -> SHT scaling factor
%   G              -> guard interval / maximum tested CIR length
%
% Output:
%   L_hat          -> estimated CIR length
%   threshold      -> SHT threshold eta_0
%   tapEnergy      -> |cir(n)|^2 for tested taps
%   decision       -> 1 if tap is significant, 0 otherwise
%
% Paper equations:
%
% eta_0 = T * sigma_hat^2
%
% sigma_hat^2 is estimated from the noise taps.
%
% The hypothesis test is:
%
% |g_hat(n)|^2 > eta_0
%
% If true, L_hat is replaced by n.
%
% =========================================================

cir = cir(:).';

% ---------------------------------------------------------
% Default parameters
% ---------------------------------------------------------

if nargin < 3 || isempty(T)
    T = 15;
end

if nargin < 4 || isempty(G)
    G = 64;
end

% ---------------------------------------------------------
% Input validation
% ---------------------------------------------------------

if isempty(cir)
    error('CIR cannot be empty.');
end

if ~isscalar(noiseVariance) || ~isfinite(noiseVariance)
    error('noiseVariance must be a finite scalar.');
end

if noiseVariance < 0
    error('noiseVariance cannot be negative.');
end

if ~isscalar(T) || ~isfinite(T) || T <= 0
    error('T must be a positive finite scalar.');
end

if ~isscalar(G) || ~isfinite(G) || G < 1 || G ~= floor(G)
    error('G must be a positive integer.');
end

if G > length(cir)
    error('G cannot exceed the CIR length.');
end

% ---------------------------------------------------------
% Calculate threshold
% ---------------------------------------------------------

threshold = T * noiseVariance;

% ---------------------------------------------------------
% Calculate tap energies
%
% MATLAB tap n corresponds to mathematical tap n.
% ---------------------------------------------------------

tapEnergy = abs(cir(1:G)).^2;

% ---------------------------------------------------------
% Sequential hypothesis testing
% ---------------------------------------------------------
%
% H0:
%   tap contains noise only
%
% H1:
%   tap contains a significant channel component
%
% The paper tests 1 <= n <= G and replaces L_hat
% with n whenever the tap energy exceeds threshold.
%
% Therefore the final L_hat is the largest significant
% tap found in the tested region.
% ---------------------------------------------------------

decision = zeros(1,G);

L_hat = 0;

for n = 1:G

    if tapEnergy(n) > threshold

        decision(n) = 1;

        L_hat = n;

    end

end

% ---------------------------------------------------------
% Display information
% ---------------------------------------------------------
% Display disabled during Monte Carlo simulations.
% The calling script will display the final results.

end