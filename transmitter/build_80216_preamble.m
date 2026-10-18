function preamble = build_80216_preamble(stsTime, ltsTime, p)

% =========================================================
% IEEE 802.16 OFDM PREAMBLE
% Based on the reference paper
% =========================================================

    stsTime = stsTime(:).';
    ltsTime = ltsTime(:).';

    % Extract one equivalent STS and one equivalent LTS
    sts64  = stsTime(1:64);
    lts128 = ltsTime(1:128);

    % IEEE 802.16 preamble structure:
    %
    % CP + 4 equivalent STSs + CP + 2 equivalent LTSs
    %
    % 64 + 256 + 64 + 256 = 640 samples

    sts = repmat(sts64, 1, 4);

    lts = repmat(lts128, 1, 2);

    % CP for the STS section
    cp_sts = sts64(end-p.G+1:end);

    % CP for the LTS section
    cp_lts = lts128(end-p.G+1:end);

    % Complete preamble
    complete = [cp_sts sts cp_lts lts];

    % Store everything
    preamble.sts64 = sts64;
    preamble.lts128 = lts128;

    preamble.sts = sts;
    preamble.lts = lts;

    preamble.cp_sts = cp_sts;
    preamble.cp_lts = cp_lts;

    preamble.complete = complete;

    fprintf('\n');
    fprintf('============================================\n');
    fprintf(' IEEE 802.16 Preamble Verification\n');
    fprintf('============================================\n');

    fprintf('STS section           : %d samples\n', length(sts));
    fprintf('Number of STS copies  : %d\n', length(sts)/64);

    fprintf('LTS section           : %d samples\n', length(lts));
    fprintf('Number of LTS copies  : %d\n', length(lts)/128);

    fprintf('STS CP                : %d samples\n', length(cp_sts));
    fprintf('LTS CP                : %d samples\n', length(cp_lts));

    fprintf('Complete preamble     : %d samples\n', length(complete));

    fprintf('============================================\n');

end