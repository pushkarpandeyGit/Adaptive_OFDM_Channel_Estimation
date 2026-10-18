function h = sui5_channel(p)
% ==============================================================
% SUI5_CHANNEL
%
% Generate one stochastic SUI-5 channel realization.
%
% SUI-5 parameters used by the base-paper reproduction:
%   Path delays : [0, 4.0, 10.0] us
%   Path powers : [0, -5, -10] dB
%   Paths       : Rayleigh fading
%
% With Ts = 0.25 us:
%   delays = [0, 16, 40] samples
%
% Output:
%   h = discrete-time SUI-5 CIR
% ==============================================================

pathDelays_us = p.SUI5.pathDelay_us;
pathPowers_dB = p.SUI5.pathPower_dB;

pathDelays_samples = pathDelays_us * 1e-6 / p.Ts;

% SUI-5 delays are integer-valued for the project's Ts.
if any(abs(pathDelays_samples - round(pathDelays_samples)) > 1e-12)
    error('SUI-5 delays are not integer samples for the current Ts.');
end

pathDelays_samples = round(pathDelays_samples);

L = max(pathDelays_samples) + 1;
h = zeros(1,L);

for k = 1:length(pathDelays_samples)
    pathAmplitude = 10^(pathPowers_dB(k)/20);

    % Unit-average-power complex Rayleigh fading
    rayleigh_coeff = (randn + 1j*randn) / sqrt(2);

    index = pathDelays_samples(k) + 1;
    h(index) = pathAmplitude * rayleigh_coeff;
end

fprintf('\nSUI-5 stochastic channel generated:\n');

fprintf('Path delays (us): ');
fprintf('%.4f ', pathDelays_us);
fprintf('\n');

fprintf('Path delays (samples): ');
fprintf('%d ', pathDelays_samples);
fprintf('\n');

fprintf('Path powers (dB): ');
fprintf('%.2f ', pathPowers_dB);
fprintf('\n');

fprintf('CIR length: %d samples\n', length(h));

end
