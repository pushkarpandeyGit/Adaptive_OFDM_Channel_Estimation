clear;
clc;
close all;

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' STEP 23 - WINDOWED DFT CHANNEL ESTIMATION\n');
fprintf('=============================================================\n');

% -------------------------------------------------------------
% Parameters
% -------------------------------------------------------------

p = parameters();
rng(p.randomSeed);

SNRdB = 45;
M = 128;
L = 64;

fprintf('SNR                 : %d dB\n', SNRdB);
fprintf('DFT/IFFT size       : %d\n', M);
fprintf('Fixed CIR length    : %d\n', L);
fprintf('Window              : Hann\n');

% -------------------------------------------------------------
% Generate OFDM preamble
% -------------------------------------------------------------

Pall = ieee80216_pall();

preamble = preamble_sequences(Pall);

[stsTime,ltsTime] = preamble_to_ifft(preamble,p);

preamble = build_80216_preamble(stsTime,ltsTime,p);

txPreamble = preamble.complete;

fprintf('Preamble length     : %d samples\n', length(txPreamble));

% -------------------------------------------------------------
% Generate stochastic SUI-2 channel
% -------------------------------------------------------------

h = sui2_channel(p);

fprintf('Channel length      : %d samples\n', length(h));

% -------------------------------------------------------------
% Pass preamble through channel
% -------------------------------------------------------------

rx = apply_channel(txPreamble,h);

% -------------------------------------------------------------
% Add AWGN
% -------------------------------------------------------------

rx = add_awgn(rx,SNRdB);

% -------------------------------------------------------------
% Extract LTS
% -------------------------------------------------------------

[rxLTS1,rxLTS2] = extract_lts(rx,p);

% -------------------------------------------------------------
% Combine the two LTS symbols
% -------------------------------------------------------------

rxLTS = combine_lts(rxLTS1,rxLTS2);

% -------------------------------------------------------------
% FFT
% -------------------------------------------------------------

Y_LTS = lts_fft(rxLTS);

% Known LTS
% -------------------------------------------------------------

knownLTS = get_known_lts(ltsTime);

X_LTS = known_lts_fft(knownLTS);

% -------------------------------------------------------------
% Extract 100 LTS training tones
% -------------------------------------------------------------

[Y_p,X_p,ltsLocations,ltsBins] = ...
    extract_lts_100(Y_LTS,X_LTS);

fprintf('Number of LTS tones  : %d\n',length(Y_p));

% -------------------------------------------------------------
% LS channel estimate
% -------------------------------------------------------------

H_LS = ls_channel_estimate(Y_p,X_p);

fprintf('LS estimate size     : %d x %d\n', ...
    size(H_LS,1),size(H_LS,2));

% -------------------------------------------------------------
% Windowed DFT
% -------------------------------------------------------------

[H_WDFT,cir_WDFT,windowedLS,window] = ...
    windowed_dft_estimate(H_LS,M,'hann',L);

fprintf('WDFT estimate size   : %d x %d\n', ...
    size(H_WDFT,1),size(H_WDFT,2));

fprintf('CIR size             : %d x %d\n', ...
    size(cir_WDFT,1),size(cir_WDFT,2));

% -------------------------------------------------------------
% True channel frequency response at LTS locations
% -------------------------------------------------------------

H_true_full = fft([h(:).' zeros(1,p.N-length(h))],p.N);

H_true_100 = zeros(1,100);

for ii = 1:100

    k = ltsLocations(ii);

    bin = mod(k,p.N) + 1;

    H_true_100(ii) = H_true_full(bin);

end

% -------------------------------------------------------------
% MSE
% -------------------------------------------------------------

mseLS = mean(abs(H_LS - H_true_100).^2);

mseWDFT = mean(abs(H_WDFT - H_true_100).^2);

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' MSE RESULTS\n');
fprintf('=============================================================\n');

fprintf('LS MSE              = %.12e\n',mseLS);
fprintf('Windowed DFT MSE    = %.12e\n',mseWDFT);

improvement = ...
    (mseLS-mseWDFT)/mseLS*100;

fprintf('MSE improvement     = %.4f %%\n',improvement);

% -------------------------------------------------------------
% CIR energy inspection
% -------------------------------------------------------------

cirEnergy = abs(cir_WDFT).^2;

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' CIR CHECK\n');
fprintf('=============================================================\n');

fprintf('Energy taps 1-%d     = %.12e\n', ...
    L,sum(cirEnergy(1:L)));

fprintf('Energy taps %d-%d  = %.12e\n', ...
    L+1,M,sum(cirEnergy(L+1:M)));

% -------------------------------------------------------------
% Check truncation
% -------------------------------------------------------------

tailMaximum = max(abs(cir_WDFT(L+1:end)));

fprintf('Maximum tail value   = %.12e\n',tailMaximum);

% -------------------------------------------------------------
% Basic validation
% -------------------------------------------------------------

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' VALIDATION\n');
fprintf('=============================================================\n');

pass = true;

if length(H_WDFT) ~= 100
    fprintf('FAIL: H_WDFT length is not 100.\n');
    pass = false;
else
    fprintf('PASS: H_WDFT has 100 tones.\n');
end

if length(cir_WDFT) ~= M
    fprintf('FAIL: CIR length is not %d.\n',M);
    pass = false;
else
    fprintf('PASS: CIR length is %d.\n',M);
end

if tailMaximum > 1e-12
    fprintf('FAIL: CIR tail was not properly truncated.\n');
    pass = false;
else
    fprintf('PASS: CIR tail is zero after tap %d.\n',L);
end

if mseWDFT < mseLS
    fprintf('PASS: Windowed DFT improves over LS for this trial.\n');
else
    fprintf('INFO: Windowed DFT did not improve LS in this trial.\n');
end

if pass
    fprintf('\nSTEP 23 BASIC VALIDATION: PASS\n');
else
    fprintf('\nSTEP 23 BASIC VALIDATION: FAIL\n');
end

fprintf('=============================================================\n');