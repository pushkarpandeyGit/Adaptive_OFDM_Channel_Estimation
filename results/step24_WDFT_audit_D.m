%% ============================================================
% STEP 24 - WDFT AUDIT D
%
% Verify that Step 23-J theta is being used correctly.
%
% This does NOT run the complete Monte Carlo experiment.
% ============================================================

clear;
clc;

fprintf('\n=============================================\n');
fprintf('STEP 24 - WDFT AUDIT D\n');
fprintf('Correct Wiener Weight Verification\n');
fprintf('=============================================\n\n');

G = 64;
Mwdft = 128;

%% ------------------------------------------------------------
% Load Step 23-J result
% ------------------------------------------------------------

fileName = 'step23J_actual_wdft_audit.mat';

if ~isfile(fileName)
    error('Cannot find step23J_actual_wdft_audit.mat');
end

S = load(fileName);

step23J = S.step23J;

theta = step23J.theta(:);

fprintf('Loaded theta from Step 23-J.\n');

fprintf('Size       : %d x %d\n', ...
    size(theta,1),size(theta,2));

fprintf('Min        : %.12e\n',min(theta));
fprintf('Max        : %.12e\n',max(theta));
fprintf('Mean       : %.12e\n',mean(theta));
fprintf('Norm       : %.12e\n',norm(theta));
fprintf('Nonzero    : %d / %d\n',nnz(theta),length(theta));

%% ------------------------------------------------------------
% Reconstruct G matrix
% ------------------------------------------------------------

Gwdft = zeros(Mwdft,101);

for m = 0:Mwdft-1

    for n = 0:100

        Gwdft(m+1,n+1) = ...
            (1/sqrt(Mwdft)) * ...
            exp(1j*2*pi* ...
            (m-Mwdft/2) * ...
            (n-100/2) / Mwdft);

    end

end

%% ------------------------------------------------------------
% Hann window
% ------------------------------------------------------------

nWin = -50:50;

hannWindow = ...
    0.5 + 0.5*cos(2*pi*nWin/Mwdft);

hannWindow = hannWindow(:);

%% ------------------------------------------------------------
% Create a deterministic test CFR
% ------------------------------------------------------------

% Three-path test channel.
%
% This deliberately resembles the validated SUI-2 structure.

hTrue = zeros(G,1);

hTrue(1) = 1.0;

hTrue(3) = 0.25;

hTrue(5) = 0.15;

%% ------------------------------------------------------------
% Generate 256-point CFR
% ------------------------------------------------------------

N = 256;

h256 = [hTrue; zeros(N-G,1)];

Htrue256 = fft(h256,N);

%% ------------------------------------------------------------
% Extract 100 LTS tones
% ------------------------------------------------------------

kLTS = -100:2:100;
kLTS = kLTS(kLTS ~= 0);

idx256 = kLTS + N/2 + 1;

H100 = Htrue256(idx256);

H100 = H100(:);

%% ------------------------------------------------------------
% Insert DC
% ------------------------------------------------------------

Hdc = 0.5*(H100(50)+H100(51));

H101 = [ ...
    H100(1:50);
    Hdc;
    H100(51:100)];

H101 = H101(:);

%% ------------------------------------------------------------
% Apply Hann window
% ------------------------------------------------------------

Hwindowed = ...
    H101 .* hannWindow;

%% ------------------------------------------------------------
% WDFT transform
% ------------------------------------------------------------

gCentered = ...
    Gwdft * Hwindowed;

gFull = fftshift(gCentered);

gWDFT = gFull(1:G);

%% ------------------------------------------------------------
% Raw and Wiener-weighted output
% ------------------------------------------------------------

gRaw = gWDFT;

gWeighted = ...
    theta .* gWDFT;

%% ------------------------------------------------------------
% Display
% ------------------------------------------------------------

fprintf('\n=============================================\n');
fprintf('WDFT TEST RESULTS\n');
fprintf('=============================================\n\n');

fprintf('True CIR energy       : %.12e\n', ...
    sum(abs(hTrue).^2));

fprintf('Raw WDFT energy       : %.12e\n', ...
    sum(abs(gRaw).^2));

fprintf('Weighted WDFT energy  : %.12e\n', ...
    sum(abs(gWeighted).^2));

fprintf('\nFirst 10 raw WDFT taps:\n');
disp(gRaw(1:10));

fprintf('\nFirst 10 weighted WDFT taps:\n');
disp(gWeighted(1:10));

%% ------------------------------------------------------------
% NMSE
% ------------------------------------------------------------

trueEnergy = sum(abs(hTrue).^2);

rawNMSE = ...
    sum(abs(gRaw-hTrue).^2)/trueEnergy;

weightedNMSE = ...
    sum(abs(gWeighted-hTrue).^2)/trueEnergy;

fprintf('\nRaw WDFT NMSE      : %.12e\n',rawNMSE);

fprintf('Weighted WDFT NMSE : %.12e\n',weightedNMSE);

%% ------------------------------------------------------------
% Check theta
% ------------------------------------------------------------

if nnz(theta) == 3

    fprintf('\nTheta nonzero-count check : PASS\n');

else

    fprintf('\nTheta nonzero-count check : NOTE\n');
    fprintf('Nonzero count = %d\n',nnz(theta));

end

if all(isfinite(theta))

    fprintf('Theta finite check         : PASS\n');

else

    error('Theta contains NaN or Inf.');

end

fprintf('\n=============================================\n');
fprintf('AUDIT D COMPLETE\n');
fprintf('=============================================\n');