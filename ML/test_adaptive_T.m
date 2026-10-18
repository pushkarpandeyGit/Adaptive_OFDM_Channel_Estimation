function test_adaptive_T()
here = fileparts(mfilename('fullpath'));
addpath(genpath(fileparts(here)));
p = parameters();

% Test SUI-2
h2 = sui2_channel(p);
H2 = local_channel_cfr(h2, p.N);
cir2 = local_ls_to_cir(H2([-100:2:-2 2:2:100] + 129), 128);
sig2 = mean(abs(cir2(65:120)).^2);

L_park2 = sht_cir_length(cir2, sig2, 15, 64);
% Adaptive T for SUI-2 (T=22)
L_adapt2 = sht_cir_length(cir2, sig2, 22, 64);

fprintf('SUI-2 (True L=5): Park T=15 -> L=%d | Adaptive T=22 -> L=%d\n', L_park2, L_adapt2);

% Test SUI-5
h5 = sui5_channel(p);
H5 = local_channel_cfr(h5, p.N);
cir5 = local_ls_to_cir(H5([-100:2:-2 2:2:100] + 129), 128);
sig5 = mean(abs(cir5(65:120)).^2);

L_park5 = sht_cir_length(cir5, sig5, 15, 64);
% Adaptive T for SUI-5 (T=9)
L_adapt5 = sht_cir_length(cir5, sig5, 9, 64);

fprintf('SUI-5 (True L=41): Park T=15 -> L=%d | Adaptive T=9 -> L=%d\n', L_park5, L_adapt5);
end

function H = local_channel_cfr(h, N)
    h = h(:); ht = [h; zeros(N - numel(h), 1)];
    H = fftshift(fft(ht, N)); H = H(:);
end
function cir = local_ls_to_cir(HLS, M)
    HLS = HLS(:).'; H101 = zeros(1, 101);
    H101(1:50) = HLS(1:50); H101(52:101) = HLS(51:100);
    n = -50:50; w = 0.5 + 0.5*cos(2*pi*n/M);
    m = (0:M-1).'; q = -50:50;
    Gmat = (1/sqrt(M)) * exp(1j*2*pi/M * ((m-M/2)*q));
    cir = fftshift(Gmat * (H101 .* w).');
end
