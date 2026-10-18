function Y_LTS = lts_fft(lts2)

% =========================================================
% FFT of the received second LTS
% =========================================================

    lts2 = lts2(:).';

    if length(lts2) ~= 128
        error('Second LTS must contain exactly 128 samples.');
    end

    % The paper uses M = 128 for the quasi-static
    % second-LTS channel-estimation case.

    Y_LTS = fft(lts2, 128);

    fprintf('\n');
    fprintf('============================================\n');
    fprintf(' Received Second LTS FFT\n');
    fprintf('============================================\n');

    fprintf('Input LTS length : %d samples\n', length(lts2));
    fprintf('FFT size         : %d\n', 128);
    fprintf('Output length    : %d bins\n', length(Y_LTS));

    fprintf('============================================\n');

end