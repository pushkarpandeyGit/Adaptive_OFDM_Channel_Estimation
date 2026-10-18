%% ============================================================
% STEP 24 - WDFT AUDIT E
%
% Compare the ACTUAL B6 Wiener weights (wWiener) against
% the sparse theta stored by Step 23-J.
%
% Then test both on the same deterministic 3-path channel.
%
% IMPORTANT:
%   wWiener is explicitly identified from the B6 file.
%   Rcross and Rhat are NOT Wiener weights.
% ============================================================

clear;
clc;
close all;

fprintf('\n=============================================\n');
fprintf('STEP 24 - WDFT AUDIT E\n');
fprintf('Actual B6 Wiener Weight Verification\n');
fprintf('=============================================\n\n');

%% ============================================================
% 1. PARAMETERS
% ============================================================

N = 256;
G = 64;
Mwdft = 128;

%% ============================================================
% 2. LOAD STEP 23-J THETA
% ============================================================

file23J = 'step23J_actual_wdft_audit.mat';

if ~isfile(file23J)
    error('Cannot find step23J_actual_wdft_audit.mat');
end

S23 = load(file23J);

if ~isfield(S23,'step23J')
    error('step23J structure not found.');
end

if ~isfield(S23.step23J,'theta')
    error('step23J.theta not found.');
end

theta = real(S23.step23J.theta(:));

%% ============================================================
% 3. LOAD B6 wWiener EXPLICITLY
% ============================================================

fileB6 = 'wdft_23B6_wiener_weights.mat';

if ~isfile(fileB6)
    error('Cannot find wdft_23B6_wiener_weights.mat');
end

SB6 = load(fileB6);

if ~isfield(SB6,'wWiener')
    error('wWiener was not found in B6 MAT file.');
end

wWiener = real(SB6.wWiener(:));

if length(wWiener) ~= G
    error('wWiener must contain exactly 64 values.');
end

%% ============================================================
% 4. DISPLAY BOTH
% ============================================================

fprintf('---------------------------------------------\n');
fprintf('STEP 23-J THETA\n');
fprintf('---------------------------------------------\n');

fprintf('Min        : %.12e\n',min(theta));
fprintf('Max        : %.12e\n',max(theta));
fprintf('Mean       : %.12e\n',mean(theta));
fprintf('Norm       : %.12e\n',norm(theta));
fprintf('Nonzero    : %d / %d\n',nnz(theta),length(theta));

fprintf('\n---------------------------------------------\n');
fprintf('B6 wWiener\n');
fprintf('---------------------------------------------\n');

fprintf('Min        : %.12e\n',min(wWiener));
fprintf('Max        : %.12e\n',max(wWiener));
fprintf('Mean       : %.12e\n',mean(wWiener));
fprintf('Norm       : %.12e\n',norm(wWiener));
fprintf('Nonzero    : %d / %d\n',nnz(wWiener),length(wWiener));

fprintf('\nAll B6 wWiener values:\n');
disp(wWiener);

%% ============================================================
% 5. DIRECT COMPARISON
% ============================================================

difference = theta - wWiener;

fprintf('\n=============================================\n');
fprintf('WEIGHT COMPARISON\n');
fprintf('=============================================\n');

fprintf('||theta-wWiener|| : %.12e\n', ...
    norm(difference));

fprintf('Relative difference: %.12e\n', ...
    norm(difference)/max(norm(wWiener),eps));

fprintf('\nFirst 10 values:\n');

fprintf('       theta             wWiener\n');

for ii = 1:10

    fprintf('%2d   %.12e    %.12e\n', ...
        ii,theta(ii),wWiener(ii));

end

%% ============================================================
% 6. BUILD EXACT PAPER EQ. (23) MATRIX
% ============================================================

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

%% ============================================================
% 7. HANN WINDOW
% ============================================================

nWin = -50:50;

hannWindow = ...
    0.5 + ...
    0.5*cos(2*pi*nWin/Mwdft);

hannWindow = hannWindow(:);

%% ============================================================
% 8. SAME DETERMINISTIC 3-PATH CHANNEL
% ============================================================

hTrue = zeros(G,1);

hTrue(1) = 1.00;
hTrue(3) = 0.25;
hTrue(5) = 0.15;

%% ============================================================
% 9. TRUE CFR
% ============================================================

h256 = [hTrue; zeros(N-G,1)];

Htrue256 = fft(h256,N);

%% ============================================================
% 10. LTS LOCATIONS
% ============================================================

kLTS = -100:2:100;
kLTS = kLTS(kLTS ~= 0);

assert(length(kLTS)==100);

idx256 = kLTS + N/2 + 1;

H100 = Htrue256(idx256);

H100 = H100(:);

%% ============================================================
% 11. INSERT DC
% ============================================================

Hdc = 0.5*(H100(50)+H100(51));

H101 = [ ...
    H100(1:50);
    Hdc;
    H100(51:100)];

H101 = H101(:);

%% ============================================================
% 12. APPLY HANN WINDOW
% ============================================================

Hwindowed = ...
    H101 .* hannWindow;

%% ============================================================
% 13. WDFT
% ============================================================

gCentered = ...
    Gwdft * Hwindowed;

gFull = fftshift(gCentered);

gRaw = gFull(1:G);

%% ============================================================
% 14. APPLY BOTH WEIGHT VECTORS
% ============================================================

hTheta = ...
    theta .* gRaw;

hB6 = ...
    wWiener .* gRaw;

%% ============================================================
% 15. NMSE
% ============================================================

trueEnergy = ...
    sum(abs(hTrue).^2);

nmseRaw = ...
    sum(abs(gRaw-hTrue).^2) / trueEnergy;

nmseTheta = ...
    sum(abs(hTheta-hTrue).^2) / trueEnergy;

nmseB6 = ...
    sum(abs(hB6-hTrue).^2) / trueEnergy;

%% ============================================================
% 16. ENERGY
% ============================================================

energyRaw = ...
    sum(abs(gRaw).^2);

energyTheta = ...
    sum(abs(hTheta).^2);

energyB6 = ...
    sum(abs(hB6).^2);

%% ============================================================
% 17. RESULTS
% ============================================================

fprintf('\n=============================================\n');
fprintf('DETERMINISTIC WDFT TEST\n');
fprintf('=============================================\n');

fprintf('\nTrue CIR energy:\n');
fprintf('  %.12e\n',trueEnergy);

fprintf('\nRaw WDFT energy:\n');
fprintf('  %.12e\n',energyRaw);

fprintf('\nTheta-weighted energy:\n');
fprintf('  %.12e\n',energyTheta);

fprintf('\nB6 wWiener-weighted energy:\n');
fprintf('  %.12e\n',energyB6);

fprintf('\n---------------------------------------------\n');
fprintf('NMSE RESULTS\n');
fprintf('---------------------------------------------\n');

fprintf('Raw WDFT       : %.12e\n',nmseRaw);

fprintf('Theta weighted : %.12e\n',nmseTheta);

fprintf('B6 wWiener     : %.12e\n',nmseB6);

%% ============================================================
% 18. FIRST 10 CIR TAPS
% ============================================================

fprintf('\n=============================================\n');
fprintf('FIRST 10 CIR TAPS\n');
fprintf('=============================================\n');

fprintf('\nTap       True          Raw          Theta-weighted       B6-weighted\n');

for ii = 1:10

    fprintf('%2d   %10.6f   %10.6f   %16.6f   %16.6f\n', ...
        ii, ...
        real(hTrue(ii)), ...
        real(gRaw(ii)), ...
        real(hTheta(ii)), ...
        real(hB6(ii)));

end

%% ============================================================
% 19. PLOT WEIGHTS
% ============================================================

figure;

plot(1:G,theta,'-o','LineWidth',1.5);
hold on;

plot(1:G,wWiener,'-x','LineWidth',1.5);

grid on;

xlabel('CIR tap index');
ylabel('Wiener weight');

title('Step 23-J \theta vs B6 wWiener');

legend( ...
    'Step 23-J \theta', ...
    'B6 wWiener', ...
    'Location','best');

%% ============================================================
% 20. PLOT CIR
% ============================================================

figure;

stem(1:10,real(hTrue(1:10)),'o');
hold on;

stem(1:10,real(gRaw(1:10)),'x');

stem(1:10,real(hTheta(1:10)),'s');

stem(1:10,real(hB6(1:10)),'d');

grid on;

xlabel('CIR tap index');
ylabel('Real CIR amplitude');

title('WDFT Reconstruction Comparison');

legend( ...
    'True CIR', ...
    'Raw WDFT', ...
    'Theta weighted', ...
    'B6 wWiener weighted', ...
    'Location','best');

%% ============================================================
% 21. BASIC VALIDITY
% ============================================================

fprintf('\n=============================================\n');
fprintf('VALIDITY CHECKS\n');
fprintf('=============================================\n');

if all(isfinite(wWiener))

    fprintf('B6 wWiener finite : PASS\n');

else

    fprintf('B6 wWiener finite : FAIL\n');

end

if all(wWiener >= 0)

    fprintf('B6 wWiener nonnegative : PASS\n');

else

    fprintf('B6 wWiener nonnegative : FAIL\n');

end

if nmseB6 < nmseRaw

    fprintf('B6 weighting improves raw WDFT : PASS\n');

else

    fprintf('B6 weighting improves raw WDFT : FAIL\n');

end

fprintf('\n=============================================\n');
fprintf('AUDIT E COMPLETE\n');
fprintf('=============================================\n');