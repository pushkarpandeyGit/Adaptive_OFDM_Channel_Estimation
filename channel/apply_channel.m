function y = apply_channel(x, h)

% =========================================================
% Apply Multipath Channel
% =========================================================
%
% Input:
%   x -> transmitted time-domain signal
%   h -> channel impulse response
%
% Output:
%   y -> received signal after multipath channel
%
% Mathematical operation:
%
%       y[n] = x[n] * h[n]
%
% where * represents linear convolution.
% =========================================================

% Make both signals row vectors
x = x(:).';
h = h(:).';

% Pass signal through channel
y = conv(x, h);

% Display information
fprintf('\n');
fprintf('============================================\n');
fprintf(' Channel Application\n');
fprintf('============================================\n');

fprintf('Input signal length    : %d samples\n', length(x));
fprintf('Channel length         : %d samples\n', length(h));
fprintf('Output signal length   : %d samples\n', length(y));

fprintf('============================================\n');

end