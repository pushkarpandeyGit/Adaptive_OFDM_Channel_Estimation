function h = sui2_channel(p)
% SUI2_CHANNEL
% Generate one stochastic SUI-2 channel realization.
%
% SUI-2:
%   Path delays  = [0, 0.4, 1.1] us
%   Path powers  = [0, -12, -15] dB
%   First path  = Ricean with K = 3 dB
%   Remaining paths = Rayleigh
%
% Output:
%   h = discrete-time channel impulse response

%% ---------------------------------------------------------
% 1. Read SUI-2 parameters
% ----------------------------------------------------------

pathDelays_us = p.SUI2.pathDelay_us;
pathPowers_dB = p.SUI2.pathPower_dB;
K_dB          = p.SUI2.RiceanK_dB;

%% ---------------------------------------------------------
% 2. Convert path delays from microseconds to samples
% ----------------------------------------------------------

pathDelays_samples = round( ...
    pathDelays_us * 1e-6 / p.Ts);

%% ---------------------------------------------------------
% 3. Convert Ricean K-factor from dB to linear
% ----------------------------------------------------------

K_linear = 10^(K_dB/10);

%% ---------------------------------------------------------
% 4. Determine required CIR length
% ----------------------------------------------------------

L = max(pathDelays_samples) + 1;

h = zeros(1,L);

%% ---------------------------------------------------------
% 5. Generate first path: Ricean fading
% ---------------------------------------------------------
%
% Ricean channel:
%
%   h = LOS component + scattered component
%
% The LOS power fraction is:
%
%   K/(K+1)
%
% The scattered power fraction is:
%
%   1/(K+1)
%
% Therefore:
%
%   E{|h|^2} = 1
%
% before applying the specified path power.

LOS_power_fraction = K_linear / (K_linear + 1);
NLOS_power_fraction = 1 / (K_linear + 1);

% Fixed LOS phase = 0
LOS_component = sqrt(LOS_power_fraction);

% Complex Gaussian scattered component
NLOS_component = sqrt(NLOS_power_fraction) * ...
    (randn + 1j*randn) / sqrt(2);

ricean_coeff = LOS_component + NLOS_component;

%% ---------------------------------------------------------
% 6. Apply specified average power to first path
% ----------------------------------------------------------

pathAmplitude1 = 10^(pathPowers_dB(1)/20);

h(pathDelays_samples(1)+1) = ...
    pathAmplitude1 * ricean_coeff;

%% ---------------------------------------------------------
% 7. Generate remaining paths: Rayleigh fading
% ----------------------------------------------------------

for k = 2:length(pathDelays_samples)

    % Unit-average-power complex Rayleigh fading
    rayleigh_coeff = ...
        (randn + 1j*randn) / sqrt(2);

    % Convert path power from dB to amplitude
    pathAmplitude = 10^(pathPowers_dB(k)/20);

    % Place path at its corresponding delay
    index = pathDelays_samples(k) + 1;

    h(index) = pathAmplitude * rayleigh_coeff;

end

%% ---------------------------------------------------------
% 8. Display channel information
% ----------------------------------------------------------

fprintf('\nSUI-2 stochastic channel generated:\n');

fprintf('Path delays (us): ');
fprintf('%.4f ', pathDelays_us);
fprintf('\n');

fprintf('Path delays (samples): ');
fprintf('%d ', pathDelays_samples);
fprintf('\n');

fprintf('Path powers (dB): ');
fprintf('%.2f ', pathPowers_dB);
fprintf('\n');

fprintf('Ricean K-factor: %.2f dB\n', K_dB);

fprintf('CIR length: %d samples\n', length(h));

end