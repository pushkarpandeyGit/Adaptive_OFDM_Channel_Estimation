function X_LTS = known_lts_fft(knownLTS)

% =========================================================
% FFT of the known LTS reference
% =========================================================

    knownLTS = knownLTS(:).';

    if length(knownLTS) ~= 128
        error('Known LTS must contain exactly 128 samples.');
    end

    X_LTS = fft(knownLTS, 128);

    fprintf('\n');
    fprintf('============================================\n');
    fprintf(' Known LTS FFT\n');
    fprintf('============================================\n');

    fprintf('Input length  : %d samples\n', length(knownLTS));
    fprintf('FFT size      : %d\n', 128);
    fprintf('Output bins   : %d\n', length(X_LTS));

    fprintf('============================================\n');

end