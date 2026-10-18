
clear;
clc;

p = parameters();
rng(p.randomSeed);

%% -------------------------------------------------------------
% Generate one stochastic SUI-2 channel
% --------------------------------------------------------------

h = sui2_channel(p);

%% -------------------------------------------------------------
% Build noiseless LTS
% --------------------------------------------------------------

Pall = ieee80216_pall();

preamble = preamble_sequences(Pall);

[~, ltsTime] = preamble_to_ifft(preamble, p);

preambleWave = build_80216_preamble( ...
    preamble_sequences(Pall).PSTS, ...
    ltsTime, ...
    p);

%% -------------------------------------------------------------
% Use the known LTS frequency response
% --------------------------------------------------------------

Htrue = fftshift(fft(ltsTime(1:128),128));

% The actual channel CFR on the 256-point OFDM grid
Hchannel = fft( ...
    [h(:); zeros(p.N-length(h),1)], ...
    p.N);

%% -------------------------------------------------------------
% Generate ideal LS pilot observations
% --------------------------------------------------------------

kLTS = -100:2:100;
kLTS(kLTS == 0) = [];

used256 = p.subcarrierIndex;
used256 = used256(used256 ~= 0);

idx256 = used256 + p.N/2 + 1;

% Extract channel response at the 100 LTS tones.
HtrueLTS = Hchannel( ...
    kLTS + p.N/2 + 1);

% Ideal LS estimate.
pLS = HtrueLTS;

%% -------------------------------------------------------------
% Run WDFT
% --------------------------------------------------------------

[Hwdft, info] = wdft_channel_estimate(pLS,p);

%% -------------------------------------------------------------
% Reference CFR
% --------------------------------------------------------------

HtrueUsed = Hchannel(idx256);

%% -------------------------------------------------------------
% Metrics
% --------------------------------------------------------------

err = HtrueUsed - Hwdft;

MSE = mean(abs(err).^2);

NMSE = mean(abs(err).^2) / ...
       mean(abs(HtrueUsed).^2);

fprintf('\n============================================\n');
fprintf(' WDFT 23B-4 VALIDATION\n');
fprintf('============================================\n');

fprintf('N              = %d\n',p.N);
fprintf('Number LTS     = %d\n',length(pLS));
fprintf('M              = %d\n',info.M);
fprintf('Retained L     = %d\n',info.L);

fprintf('\nNoiseless WDFT MSE  = %.12e\n',MSE);
fprintf('Noiseless WDFT NMSE = %.12e\n',NMSE);

fprintf('\nTrue channel power = %.12e\n', ...
    mean(abs(HtrueUsed).^2));

fprintf('WDFT power         = %.12e\n', ...
    mean(abs(Hwdft).^2));

fprintf('Maximum error      = %.12e\n', ...
    max(abs(err)));

%% -------------------------------------------------------------
% Matrix consistency test
% --------------------------------------------------------------

K = info.K_WDFT;

H_matrix = K * pLS;

matrixError = max(abs(H_matrix - Hwdft));

fprintf('\nK matrix size = %d x %d\n', ...
    size(K,1),size(K,2));

fprintf('Matrix/direct WDFT difference = %.12e\n', ...
    matrixError);

if matrixError < 1e-12
    fprintf('PASS: K_WDFT reproduces direct estimator.\n');
else
    fprintf('FAIL: K_WDFT mismatch.\n');
end

%% -------------------------------------------------------------
% CIR diagnostic
% --------------------------------------------------------------

fprintf('\nCIR energy:\n');

fprintf('Total       = %.12e\n', ...
    sum(abs(info.gCentered).^2));

fprintf('Causal L=64 = %.12e\n', ...
    sum(abs(info.gCausal).^2));

fprintf('Discarded    = %.12e\n', ...
    sum(abs(info.gCentered).^2) - ...
    sum(abs(info.gCausal).^2));

%% -------------------------------------------------------------
% Show strongest causal taps
% --------------------------------------------------------------

[~,ord] = sort(abs(info.gCausal).^2,'descend');

fprintf('\nStrongest causal taps:\n');

for ii = 1:min(10,length(ord))

    tap = ord(ii);

    fprintf('tap %2d : |g|^2 = %.8e\n', ...
        tap,abs(info.gCausal(tap)).^2);

end
