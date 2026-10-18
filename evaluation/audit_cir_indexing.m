clear
clc
close all

fprintf('\n');
fprintf('====================================================\n');
fprintf(' STEP 36 - CIR INDEXING / FREQUENCY GRID TEST\n');
fprintf('====================================================\n');

p = parameters();

% =========================================================
% 1. SUI-2 physical channel
% =========================================================

h = sui2_channel(p);

% =========================================================
% 2. 100 LTS training locations
% =========================================================

k = [-100:2:-2 2:2:100];

% =========================================================
% 3. Theoretical channel frequency response
% =========================================================

H = zeros(1,length(k));

for m = 1:length(k)

    for n = 0:length(h)-1

        H(m) = H(m) + ...
            h(n+1) * exp(-1j*2*pi*k(m)*n/p.N);

    end

end

fprintf('\nNumber of theoretical LTS tones = %d\n',length(H));


% =========================================================
% TEST A
% 100 tones only
%
% 14 zeros + 100 tones + 14 zeros = 128
% =========================================================

w100 = 0.5 + 0.5*cos(2*pi*(-49.5:49.5)/128);

Hw100 = H .* w100;

gridA = [zeros(1,14),Hw100,zeros(1,14)];

cirA = ifft(gridA,128);


% =========================================================
% TEST B
% 101 positions including DC
%
% 14 zeros + 101 positions + 13 zeros = 128
% =========================================================

H101 = zeros(1,101);

H101(1:50) = H(1:50);
H101(52:101) = H(51:100);

n101 = -50:50;

w101 = 0.5 + 0.5*cos(2*pi*n101/128);

Hw101 = H101 .* w101;

gridB = [zeros(1,14),Hw101,zeros(1,13)];

cirB = ifft(gridB,128);


% =========================================================
% TEST C
% 100 tones with integer reduced-frequency indexing
%
% k/2 = -50...-1,1...50
% =========================================================

r = k/2;

wC = 0.5 + 0.5*cos(2*pi*r/128);

HwC = H .* wC;

gridC = [zeros(1,14),HwC,zeros(1,14)];

cirC = ifft(gridC,128);


% =========================================================
% Function to display strongest samples
% =========================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' TEST A - 100 TONES, 14/14 ZERO PADDING\n');
fprintf('====================================================\n');

[pA,idxA] = sort(abs(cirA).^2,'descend');

for m = 1:10

    fprintf('Sample %3d : %.8f\n', ...
        idxA(m),abs(cirA(idxA(m))));

end


fprintf('\n');
fprintf('====================================================\n');
fprintf(' TEST B - 101 POSITIONS INCLUDING DC\n');
fprintf('====================================================\n');

[pB,idxB] = sort(abs(cirB).^2,'descend');

for m = 1:10

    fprintf('Sample %3d : %.8f\n', ...
        idxB(m),abs(cirB(idxB(m))));

end


fprintf('\n');
fprintf('====================================================\n');
fprintf(' TEST C - 100 TONES, REDUCED FREQUENCY INDEX\n');
fprintf('====================================================\n');

[pC,idxC] = sort(abs(cirC).^2,'descend');

for m = 1:10

    fprintf('Sample %3d : %.8f\n', ...
        idxC(m),abs(cirC(idxC(m))));

end


% =========================================================
% Energy comparison
% =========================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' ENERGY COMPARISON\n');
fprintf('====================================================\n');

fprintf('\nTest A:\n');
fprintf('1-64      = %.10e\n',sum(abs(cirA(1:64)).^2));
fprintf('65-120    = %.10e\n',sum(abs(cirA(65:120)).^2));
fprintf('121-128   = %.10e\n',sum(abs(cirA(121:128)).^2));
fprintf('Total     = %.10e\n',sum(abs(cirA).^2));

fprintf('\nTest B:\n');
fprintf('1-64      = %.10e\n',sum(abs(cirB(1:64)).^2));
fprintf('65-120    = %.10e\n',sum(abs(cirB(65:120)).^2));
fprintf('121-128   = %.10e\n',sum(abs(cirB(121:128)).^2));
fprintf('Total     = %.10e\n',sum(abs(cirB).^2));

fprintf('\nTest C:\n');
fprintf('1-64      = %.10e\n',sum(abs(cirC(1:64)).^2));
fprintf('65-120    = %.10e\n',sum(abs(cirC(65:120)).^2));
fprintf('121-128   = %.10e\n',sum(abs(cirC(121:128)).^2));
fprintf('Total     = %.10e\n',sum(abs(cirC).^2));


% =========================================================
% Plot all three
% =========================================================

figure;

plot(1:128,abs(cirA),'LineWidth',1.2);
hold on;

plot(1:128,abs(cirB),'LineWidth',1.2);

plot(1:128,abs(cirC),'LineWidth',1.2);

xlabel('CIR sample');
ylabel('|CIR|');

title('CIR Frequency-Grid Comparison');

legend( ...
    'A: 100 tones, 14/14', ...
    'B: 101 positions incl. DC', ...
    'C: reduced-frequency index');

grid on;

xlim([1 128]);


fprintf('\n');
fprintf('====================================================\n');
fprintf(' STEP 36 COMPLETE\n');
fprintf('====================================================\n');