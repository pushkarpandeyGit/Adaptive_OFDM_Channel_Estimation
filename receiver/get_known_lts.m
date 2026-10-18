function knownLTS = get_known_lts(ltsTime)

% =========================================================
% Obtain one known 128-sample LTS period
% =========================================================

    ltsTime = ltsTime(:).';

    if length(ltsTime) ~= 256
        error('ltsTime must contain exactly 256 samples.');
    end

    % One complete LTS period
    knownLTS = ltsTime(1:128);

    fprintf('\n');
    fprintf('============================================\n');
    fprintf(' Known LTS Reference\n');
    fprintf('============================================\n');

    fprintf('Original LTS length : %d samples\n', length(ltsTime));
    fprintf('Reference length    : %d samples\n', length(knownLTS));

    fprintf('============================================\n');

end