clear
clc
close all

fprintf('\n');
fprintf('====================================================\n');
fprintf(' STEP 37 - CENTERED FREQUENCY / IFFTSHIFT TEST\n');
fprintf('====================================================\n');

p = parameters();

% =========================================================
% 1. Generate SUI-2 channel
% =========================================================

h = sui2_channel(p);

% =========================================================
% 2. 100 LTS training locations
% =========================================================

k = [-100:2:-2 2:2:100];

% Reduced-frequency index
r = k/2;

% =========================================================
% 3. Theoretical frequency response
% =========================================================

H = zeros(1,100);

for m = 1:100

    for n = 0:length(h)-1

        H(m) = H(m) + ...
            h(n+1) * exp(-1j*2*pi*k(m)*n/p.N);

    end

end

% =========================================================
% 4. Hann window
% =========================================================

window = ...
    0.5 + 0.5*cos(2*pi*r/128);

Hw = H .* window;

% =========================================================
% 5. Construct CENTERED 128-point frequency vector
%
% Frequency order:
%
% -64 ... -51
% -50 ... -1
%   1 ... 50
%  51 ... 63
%
% =========================================================

Hcentered = [ ...
    zeros(1,14), ...
    Hw, ...
    zeros(1,14)];

fprintf('\nCentered frequency vector length = %d\n', ...
    length(Hcentered));


% =========================================================
% TEST A
%
% WRONG for centered vector:
% direct IFFT
% =========================================================

cirDirect = ifft(Hcentered,128);


% =========================================================
% TEST B
%
% CORRECT handling of centered frequency vector:
%
% centered spectrum
%       ↓
% ifftshift
%       ↓
% IFFT
% =========================================================

cirShifted = ifft(ifftshift(Hcentered),128);


% =========================================================
% 6. Strongest samples - direct IFFT
% =========================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' DIRECT IFFT\n');
fprintf('====================================================\n');

[p1,idx1] = sort(abs(cirDirect).^2,'descend');

for m = 1:10

    fprintf('Sample %3d : %.8f\n', ...
        idx1(m), ...
        abs(cirDirect(idx1(m))));

end


% =========================================================
% 7. Strongest samples - IFFTSHIFT + IFFT
% =========================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' IFFTSHIFT + IFFT\n');
fprintf('====================================================\n');

[p2,idx2] = sort(abs(cirShifted).^2,'descend');

for m = 1:10

    fprintf('Sample %3d : %.8f\n', ...
        idx2(m), ...
        abs(cirShifted(idx2(m))));

end


% =========================================================
% 8. Expected physical paths
% =========================================================

pathDelays_samples = ...
    round(p.SUI2.pathDelay_us * 1e-6 / p.Ts);

fprintf('\n');
fprintf('====================================================\n');
fprintf(' EXPECTED PATH LOCATIONS\n');
fprintf('====================================================\n');

for m = 1:length(pathDelays_samples)

    n = pathDelays_samples(m) + 1;

    fprintf('Path %d -> sample %d\n',m,n);

end


% =========================================================
% 9. CIR values at physical paths
% =========================================================

fprintf('\n');
fprintf('DIRECT IFFT:\n');

for m = 1:length(pathDelays_samples)

    n = pathDelays_samples(m) + 1;

    fprintf('Sample %3d : %.8f\n', ...
        n,abs(cirDirect(n)));

end


fprintf('\n');
fprintf('IFFTSHIFT + IFFT:\n');

for m = 1:length(pathDelays_samples)

    n = pathDelays_samples(m) + 1;

    fprintf('Sample %3d : %.8f\n', ...
        n,abs(cirShifted(n)));

end


% =========================================================
% 10. Energy regions
% =========================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' ENERGY COMPARISON\n');
fprintf('====================================================\n');

fprintf('\nDIRECT IFFT:\n');

E1 = sum(abs(cirDirect(1:64)).^2);
E2 = sum(abs(cirDirect(65:120)).^2);
E3 = sum(abs(cirDirect(121:128)).^2);
ET = sum(abs(cirDirect).^2);

fprintf('1-64     = %.10e\n',E1);
fprintf('65-120   = %.10e\n',E2);
fprintf('121-128  = %.10e\n',E3);
fprintf('Total    = %.10e\n',ET);


fprintf('\nIFFTSHIFT + IFFT:\n');

E1 = sum(abs(cirShifted(1:64)).^2);
E2 = sum(abs(cirShifted(65:120)).^2);
E3 = sum(abs(cirShifted(121:128)).^2);
ET = sum(abs(cirShifted).^2);

fprintf('1-64     = %.10e\n',E1);
fprintf('65-120   = %.10e\n',E2);
fprintf('121-128  = %.10e\n',E3);
fprintf('Total    = %.10e\n',ET);


% =========================================================
% 11. Plot comparison
% =========================================================

figure;

plot(1:128,abs(cirDirect),'LineWidth',1.2);
hold on;

plot(1:128,abs(cirShifted),'LineWidth',1.2);

xlabel('CIR sample');
ylabel('|CIR|');

title('CIR: Direct IFFT vs IFFTSHIFT + IFFT');

legend( ...
    'Direct IFFT', ...
    'IFFTSHIFT + IFFT');

grid on;

xlim([1 128]);


% =========================================================
% END
% =========================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' STEP 37 COMPLETE\n');
fprintf('====================================================\n');