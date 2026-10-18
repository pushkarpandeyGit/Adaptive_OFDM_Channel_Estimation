clear
clc
close all

% =========================================================
% STEP 35
% THEORETICAL FREQUENCY RESPONSE -> CIR REFERENCE TEST
% =========================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' STEP 35 - CIR THEORETICAL REFERENCE TEST\n');
fprintf('====================================================\n');

p = parameters();


% =========================================================
% 1. Generate SUI-2 channel
% =========================================================

h = sui2_channel(p);

fprintf('\nPhysical channel:\n');

for n = 1:length(h)

    if abs(h(n)) > 0

        fprintf('h(%d) = %.8f\n', ...
            n, h(n));

    end

end


% =========================================================
% 2. Generate theoretical frequency response
%
% H(k) = sum h(n) exp(-j 2*pi*k*n/N)
%
% Here k corresponds to the original IEEE 802.16
% subcarrier locations.
% =========================================================

ltsLocations = [-100:2:-2 2:2:100];

H_theory = zeros(1,length(ltsLocations));

for m = 1:length(ltsLocations)

    k = ltsLocations(m);

    H = 0;

    for n = 0:length(h)-1

        H = H + h(n+1) * ...
            exp(-1j*2*pi*k*n/p.N);

    end

    H_theory(m) = H;

end


fprintf('\n');
fprintf('Theoretical H length = %d\n', ...
    length(H_theory));


% =========================================================
% 3. Convert theoretical frequency response to CIR
% =========================================================

[cirTheory,windowedTheory,window] = ...
    ls_to_cir(H_theory,128,'hann');


% =========================================================
% 4. Print strongest theoretical CIR samples
% =========================================================

[powerSorted,idxSorted] = ...
    sort(abs(cirTheory).^2,'descend');

fprintf('\n');
fprintf('STRONGEST THEORETICAL CIR SAMPLES\n');
fprintf('----------------------------------------------------\n');

for k = 1:10

    n = idxSorted(k);

    fprintf(['Sample %3d : magnitude = %.8f, ', ...
             'power = %.8f\n'], ...
             n, ...
             abs(cirTheory(n)), ...
             powerSorted(k));

end


% =========================================================
% 5. Expected physical path locations
% =========================================================

pathDelays_samples = ...
    round(p.SUI2.pathDelay_us * 1e-6 / p.Ts);

fprintf('\n');
fprintf('THEORETICAL CIR AT EXPECTED PATHS\n');
fprintf('----------------------------------------------------\n');

for k = 1:length(pathDelays_samples)

    n = pathDelays_samples(k) + 1;

    fprintf(['Path %d -> sample %d : magnitude = ', ...
             '%.8f\n'], ...
             k, ...
             n, ...
             abs(cirTheory(n)));

end


% =========================================================
% 6. Energy regions
% =========================================================

E1 = sum(abs(cirTheory(1:64)).^2);

E2 = sum(abs(cirTheory(65:120)).^2);

E3 = sum(abs(cirTheory(121:128)).^2);

ET = sum(abs(cirTheory).^2);


fprintf('\n');
fprintf('THEORETICAL CIR ENERGY\n');
fprintf('----------------------------------------------------\n');

fprintf('Samples 1-64     : %.10e\n',E1);
fprintf('Samples 65-120   : %.10e\n',E2);
fprintf('Samples 121-128  : %.10e\n',E3);
fprintf('Total            : %.10e\n',ET);
fprintf('Noise fraction   : %.6f\n',E2/ET);


% =========================================================
% 7. Plot theoretical CIR
% =========================================================

figure;

plot(1:128,abs(cirTheory),'LineWidth',1.2);

xlabel('CIR sample');
ylabel('|CIR|');

title('Theoretical SUI-2 CIR after Hann Window');

grid on;

xlim([1 128]);


% =========================================================
% 8. Compare with physical channel
% =========================================================

figure;

stem(1:length(h),abs(h),'filled');

xlabel('Physical channel sample');
ylabel('|h[n]|');

title('Original SUI-2 Physical Channel');

grid on;


fprintf('\n');
fprintf('====================================================\n');
fprintf(' STEP 35 COMPLETE\n');
fprintf('====================================================\n');