function y = add_awgn(x, SNRdB)

% =========================================================
% AWGN Noise Addition
% =========================================================
%
% Input:
%   x      -> received signal before noise
%   SNRdB  -> desired signal-to-noise ratio in dB
%
% Output:
%   y      -> noisy received signal
%
% =========================================================

% Make sure signal is a row vector
x = x(:).';

% -----------------------------------------------------
% Calculate average signal power
% -----------------------------------------------------

signalPower = mean(abs(x).^2);

% -----------------------------------------------------
% Convert SNR from dB to linear
%
% SNR_linear = 10^(SNRdB/10)
% -----------------------------------------------------

SNRlinear = 10^(SNRdB/10);

% -----------------------------------------------------
% Calculate required noise power
%
% SNR = Signal Power / Noise Power
%
% Therefore:
%
% Noise Power = Signal Power / SNR
% -----------------------------------------------------

noisePower = signalPower / SNRlinear;

% -----------------------------------------------------
% Generate complex Gaussian noise
%
% Real part and imaginary part each have half
% of the total noise power.
% -----------------------------------------------------

noise = sqrt(noisePower/2) * ...
    (randn(size(x)) + 1j*randn(size(x)));

% -----------------------------------------------------
% Add noise to received signal
% -----------------------------------------------------

y = x + noise;

% -----------------------------------------------------
% Display information
% -----------------------------------------------------

fprintf('\n');
fprintf('============================================\n');
fprintf(' AWGN Noise Addition\n');
fprintf('============================================\n');

fprintf('Signal power : %.6f\n', signalPower);
fprintf('SNR           : %.2f dB\n', SNRdB);
fprintf('Noise power   : %.6f\n', noisePower);
fprintf('Samples       : %d\n', length(x));

fprintf('============================================\n');

end