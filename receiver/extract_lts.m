function [lts1, lts2] = extract_lts(rx, p)

% =========================================================
% IEEE 802.16 LTS Extraction
% =========================================================

    rx = rx(:).';

    % Preamble structure:
    %
    % 64 CP
    % 256 STS
    % 64 CP
    % 256 LTS
    %
    % LTS section = 2 x 128 samples

    ltsStart = p.G + 4*64 + p.G + 1;

    % First LTS
    lts1Start = ltsStart;
    lts1End   = lts1Start + 127;

    % Second LTS
    lts2Start = lts1End + 1;
    lts2End   = lts2Start + 127;

    % Safety check
    if lts2End > length(rx)
        error('Received signal is too short to extract both LTS blocks.');
    end

    % Extract LTS blocks
    lts1 = rx(lts1Start:lts1End);
    lts2 = rx(lts2Start:lts2End);

    fprintf('\n');
    fprintf('============================================\n');
    fprintf(' IEEE 802.16 LTS Extraction\n');
    fprintf('============================================\n');

    fprintf('LTS1 samples : %d to %d\n', ...
        lts1Start, lts1End);

    fprintf('LTS2 samples : %d to %d\n', ...
        lts2Start, lts2End);

    fprintf('LTS1 length  : %d samples\n', length(lts1));
    fprintf('LTS2 length  : %d samples\n', length(lts2));

    fprintf('============================================\n');

end