clc;
clear;

addpath(genpath('/MATLAB Drive/Adaptive_OFDM_Channel_Estimation'));

% Re-run the existing end-to-end diagnostic
p = parameters;

Pall = ieee80216_pall;
preamble = preamble_sequences(Pall);

[stsTime,ltsTime] = preamble_to_ifft(preamble,p);

knownLTS = get_known_lts(ltsTime);
X_LTS = known_lts_fft(knownLTS);

preambleStruct = build_80216_preamble(stsTime,ltsTime,p);
tx = preambleStruct.complete;

hTrue = sui2_channel(p);

rx = apply_channel(tx,hTrue);

[~,lts2_rx] = extract_lts(rx,p);

Y_LTS = lts_fft(lts2_rx);

[Y_p,X_p,~,~] = extract_lts_100(Y_LTS,X_LTS);

hLS_100 = pilot_ls_estimate(Y_p,X_p);

[cir,~,~,~] = ls100_to_cir128(hLS_100);

% Shift the centered CIR back to normal delay indexing
cirShifted = circshift(cir,-64);

energy = abs(cirShifted).^2;

[sortedEnergy,sortedIndex] = sort(energy,'descend');

fprintf('\n');
fprintf('====================================================\n');
fprintf(' CIR SHIFT DIAGNOSTIC\n');
fprintf('====================================================\n');

fprintf('\nStrongest shifted CIR samples:\n');
fprintf('Rank       Index        Energy       Magnitude\n');
fprintf('------------------------------------------------\n');

for n = 1:10
    idx = sortedIndex(n);

    fprintf('%4d %12d %14.6e %14.6e\n', ...
        n,idx,sortedEnergy(n),abs(cirShifted(idx)));
end

fprintf('\nFirst 10 shifted CIR samples:\n');
fprintf('Index          Magnitude\n');
fprintf('------------------------\n');

for n = 1:10
    fprintf('%5d %15.6f\n',n,abs(cirShifted(n)));
end

fprintf('\nTrue channel:\n');

for n = 1:length(hTrue)
    fprintf('h(%d) magnitude = %.6f\n',n,abs(hTrue(n)));
end

figure;

stem(1:128,abs(cirShifted),'filled');

grid on;

xlabel('CIR sample index');
ylabel('Magnitude');

title('CIR After 64-Sample Circular Shift');

xlim([1 128]);

fprintf('\n====================================================\n');
fprintf(' SHIFT DIAGNOSTIC COMPLETE\n');
fprintf('====================================================\n');