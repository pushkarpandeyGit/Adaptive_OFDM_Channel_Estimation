projectRoot = 'c:\Users\pushk\Downloads\Adaptive_OFDM_Channel_Estimation';
addpath(genpath(projectRoot));
p = parameters();

H_LS = ones(100, 1);
[cir_real, winLS, w] = ls_to_cir(H_LS, 128, 'hann');

% From step1_generate_cir_dataset
M = 128;
q = (-50:50).';
w_step1 = 0.5 + 0.5*cos(2*pi*q/M);
m_idx = (0:M-1).';
Gmat = (1/sqrt(M)) * exp(1j*2*pi*((m_idx - M/2)*q.')/M);
H101 = [H_LS(1:50); 0.5*(H_LS(50)+H_LS(51)); H_LS(51:100)];
cir_step1 = fftshift(Gmat * (H101 .* w_step1));

fprintf('Scale check:\n');
fprintf('abs(cir_real(1:5)):  %s\n', mat2str(abs(cir_real(1:5)), 4));
fprintf('abs(cir_step1(1:5)): %s\n', mat2str(abs(cir_step1(1:5).'), 4));

fprintf('\nPeak position check:\n');
[~, p_real] = max(abs(cir_real));
[~, p_step1] = max(abs(cir_step1));
fprintf('Peak real: %d, Peak step1: %d\n', p_real, p_step1);
