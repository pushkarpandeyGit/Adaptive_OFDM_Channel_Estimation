projectRoot = 'c:\Users\pushk\Downloads\Adaptive_OFDM_Channel_Estimation';
modelData = load(fullfile(projectRoot, 'ML', 'trained_ml_cir_model.mat'));
data = load(fullfile(projectRoot, 'ML', 'cir_dataset.mat'));

fprintf('mu:\n%s\n\n', mat2str(modelData.mu, 3));
fprintf('sigma:\n%s\n\n', mat2str(modelData.sigma, 3));

% Compare training sample 1 vs runtime feat for SUI-2
addpath(genpath(projectRoot));
p = parameters();
h = sui2_channel(p);
H = fftshift(fft([h(:); zeros(p.N-length(h),1)]));
Hp = H(p.N/2+1 + [-100:2:-2 2:2:100]);
snr = 20;
gammaP = 10^((snr+6)/10);
HLS = Hp + sqrt(1/(2*gammaP))*(randn(100,1)+1j*randn(100,1));
[cir,~,~] = ls_to_cir(HLS, 128, 'hann');
noiseVar = mean(abs(cir(p.noiseTapStart:p.noiseTapEnd)).^2);
[~, feat_runtime] = ml_cir_length(cir, noiseVar, p.G, snr);

fprintf('Training feature means:\n');
disp(modelData.mu);

fprintf('Runtime feature at 20dB SUI-2:\n');
disp(feat_runtime);
