function [noiseVariance, noiseTaps] = estimate_noise_variance(cir, startTap, endTap)

% =========================================================
% CIR Noise Variance Estimation
% =========================================================
%
% Input:
%   cir      -> estimated channel impulse response
%   startTap -> first noise-estimation tap
%   endTap   -> last noise-estimation tap
%
% Output:
%   noiseVariance -> estimated noise variance
%   noiseTaps     -> CIR samples used for estimation
%
% Mathematical equation:
%
%                 1
% sigma_hat^2 = ----- sum |h_hat[n]|^2
%                 K
%
% where
%
% K = endTap - startTap + 1
%
% For the paper:
%
% startTap = 65
% endTap   = 120
% K        = 56
%
% =========================================================

cir = cir(:).';

if nargin < 2 || isempty(startTap)
    startTap = 65;
end

if nargin < 3 || isempty(endTap)
    endTap = 120;
end

% ---------------------------------------------------------
% Validate tap range
% ---------------------------------------------------------

if startTap < 1
    error('startTap must be at least 1.');
end

if endTap > length(cir)
    error('endTap exceeds the CIR length.');
end

if startTap > endTap
    error('startTap must be less than or equal to endTap.');
end

% ---------------------------------------------------------
% Extract noise-only taps
% ---------------------------------------------------------

noiseTaps = cir(startTap:endTap);

% Number of samples
K = length(noiseTaps);

% ---------------------------------------------------------
% Estimate noise variance
% ---------------------------------------------------------

noiseVariance = mean(abs(noiseTaps).^2);

% ---------------------------------------------------------
% Display information
% ---------------------------------------------------------

%fprintf('\n');
%fprintf('============================================\n');
%fprintf(' CIR Noise Variance Estimation\n');
%fprintf('============================================\n');

%fprintf('CIR length       : %d samples\n', length(cir));
%fprintf('Noise taps       : %d to %d\n', startTap, endTap);
%fprintf('Number of taps K : %d\n', K);
%fprintf('Noise variance   : %.6e\n', noiseVariance);

%fprintf('============================================\n');

end