function [stsTime, ltsTime] = preamble_to_ifft(preamble, p)

% PREAMBLE_TO_IFFT Convert IEEE 802.16 frequency-domain
% STS/LTS sequences into 256-point time-domain signals.

% Create 256-point frequency-domain vectors
X_STS = zeros(1, p.N);
X_LTS = zeros(1, p.N);

% IEEE 802.16 subcarrier indices: -100 to +100
k = preamble.k;

% MATLAB FFT-bin mapping
% Negative frequencies -> bins near the end
% Positive frequencies -> bins near the beginning
for n = 1:length(k)

    if k(n) < 0
        bin = p.N + k(n) + 1;
    elseif k(n) > 0
        bin = k(n) + 1;
    else
        bin = 1;
    end

    X_STS(bin) = preamble.PSTS(n);
    X_LTS(bin) = preamble.PLTS(n);

end

% Convert frequency domain -> time domain
stsTime = ifft(X_STS, p.N);
ltsTime = ifft(X_LTS, p.N);

end