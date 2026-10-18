clear;
clc;

p = parameters();
rng(p.randomSeed);

%% -------------------------------------------------------------
% Generate SUI-2 channel
% --------------------------------------------------------------

h = sui2_channel(p);

N = p.N;
M = 128;
Np = 100;
L = 64;

%% -------------------------------------------------------------
% True channel CFR
% --------------------------------------------------------------

hN = [h(:); zeros(N-length(h),1)];

Htrue = fft(hN,N);

usedK = p.subcarrierIndex;
usedK = usedK(usedK ~= 0);

usedIdx = usedK + N/2 + 1;

HtrueUsed = Htrue(usedIdx);

%% -------------------------------------------------------------
% LTS subcarriers
% --------------------------------------------------------------

kLTS = -100:2:100;
kLTS(kLTS == 0) = [];

HtrueLTS = Htrue(kLTS + N/2 + 1);

pLS = HtrueLTS;

%% -------------------------------------------------------------
% Run current WDFT
% --------------------------------------------------------------

[Hbase,info] = wdft_channel_estimate(pLS,p);

%% -------------------------------------------------------------
% Test different scalar normalizations
% --------------------------------------------------------------

scales = [0.5 1 1.5 2 sqrt(2) 2*sqrt(2)];

fprintf('\n');
fprintf('============================================\n');
fprintf(' WDFT 23B-5 SCALING DIAGNOSTIC\n');
fprintf('============================================\n');

fprintf('\nTrue channel power = %.12e\n', ...
    mean(abs(HtrueUsed).^2));

fprintf('\n');

for ii = 1:length(scales)

    a = scales(ii);

    Htest = a * Hbase;

    err = HtrueUsed - Htest;

    mse = mean(abs(err).^2);

    nmse = mse / mean(abs(HtrueUsed).^2);

    power = mean(abs(Htest).^2);

    fprintf('Scale = %8.5f   Power = %.8e   MSE = %.8e   NMSE = %.8e\n', ...
        a,power,mse,nmse);

end

%% -------------------------------------------------------------
% Least-squares optimal scalar
% --------------------------------------------------------------
%
% Find alpha that minimizes:
%
% || Htrue - alpha Hbase ||^2
%
% alpha = Re{Hbase^H Htrue}/||Hbase||^2
%
% --------------------------------------------------------------

alpha = real(Hbase' * HtrueUsed) / ...
        real(Hbase' * Hbase);

Hscaled = alpha * Hbase;

err = HtrueUsed - Hscaled;

mse = mean(abs(err).^2);

nmse = mse / mean(abs(HtrueUsed).^2);

fprintf('\n============================================\n');
fprintf(' OPTIMAL REAL SCALAR\n');
fprintf('============================================\n');

fprintf('Optimal alpha = %.12e\n',alpha);
fprintf('Scaled power  = %.12e\n', ...
    mean(abs(Hscaled).^2));

fprintf('Scaled MSE    = %.12e\n',mse);
fprintf('Scaled NMSE   = %.12e\n',nmse);