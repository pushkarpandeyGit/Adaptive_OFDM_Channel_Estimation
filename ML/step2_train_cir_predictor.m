function step2_train_cir_predictor()
% =========================================================================
% STEP 2: TRAIN NOVEL ML CIR LENGTH PREDICTOR
%
% Novelty Component (Day 1):
% Trains a lightweight Multi-Layer Perceptron (MLP) Neural Network natively
% in MATLAB (using pure matrix math — zero external toolboxes required!).
%
% Replaces heuristic SHT threshold T=15 with a learned mapping:
%   f([SNR, NoiseVariance, EnergyPercentiles, PeakTap, SHT_L]) -> L_hat
%
% Saves:
%   ML/trained_ml_cir_model.mat
% =========================================================================

fprintf('\n=======================================================\n');
fprintf(' NOVELTY DAY 1: TRAINING NOVEL ML CIR PREDICTOR\n');
fprintf('=======================================================\n');

%% 1. PATH SETUP & LOAD DATASET
thisFile = mfilename('fullpath');
if isempty(thisFile), here = pwd; else, here = fileparts(thisFile); end
projectRoot = fileparts(here);
addpath(genpath(projectRoot));

datasetPath = fullfile(projectRoot, 'ML', 'cir_dataset.mat');
if ~exist(datasetPath, 'file')
    error('cir_dataset.mat not found. Please run step1_generate_cir_dataset first.');
end

data = load(datasetPath);
X = data.featureMatrix;
y = data.targetVector;

N_samples = size(X, 1);
fprintf('Loaded %d training samples.\n', N_samples);

%% 2. TRAIN / TEST SPLIT (80% Train, 20% Test)
rng(2026);
perm = randperm(N_samples);
splitIdx = round(0.80 * N_samples);

trainIdx = perm(1:splitIdx);
testIdx  = perm(splitIdx+1:end);

X_train = X(trainIdx, :);
y_train = y(trainIdx);

X_test  = X(testIdx, :);
y_test  = y(testIdx);

%% 3. FEATURE NORMALIZATION (Z-SCORE)
mu = mean(X_train, 1);
sigma = std(X_train, 0, 1);
sigma(sigma == 0) = 1; % Prevent divide by zero

X_train_norm = (X_train - mu) ./ sigma;
X_test_norm  = (X_test - mu) ./ sigma;

%% 4. NEURAL NETWORK ARCHITECTURE
% Architecture: Input(11) -> Dense(24, LeakyReLU) -> Dense(12, LeakyReLU) -> Dense(1, Linear)
D_in = size(X, 2);
H1 = 24;
H2 = 12;
D_out = 1;

% He (Kaiming) Initialization
W1 = randn(H1, D_in) * sqrt(2 / D_in);
b1 = zeros(H1, 1);

W2 = randn(H2, H1) * sqrt(2 / H1);
b2 = zeros(H2, 1);

W3 = randn(D_out, H2) * sqrt(2 / H2);
b3 = mean(y_train);

%% 5. TRAINING VIA ADAM OPTIMIZER (NATIVE MATLAB)
epochs = 350;
batchSize = 64;
lr = 0.008;
beta1 = 0.9;
beta2 = 0.999;
eps = 1e-8;

mW1 = zeros(size(W1)); vW1 = zeros(size(W1));
mb1 = zeros(size(b1)); vb1 = zeros(size(b1));
mW2 = zeros(size(W2)); vW2 = zeros(size(W2));
mb2 = zeros(size(b2)); vb2 = zeros(size(b2));
mW3 = zeros(size(W3)); vW3 = zeros(size(W3));
mb3 = zeros(size(b3)); vb3 = zeros(size(b3));

numBatches = floor(size(X_train_norm, 1) / batchSize);
fprintf('Training Native Neural Network (%d epochs, batch size %d)...\n', epochs, batchSize);

t_step = 0;
for epoch = 1:epochs
    % Shuffle
    bPerm = randperm(size(X_train_norm, 1));
    X_shuff = X_train_norm(bPerm, :);
    y_shuff = y_train(bPerm);

    for b = 1:numBatches
        t_step = t_step + 1;
        idxB = (b-1)*batchSize + 1 : b*batchSize;
        xb = X_shuff(idxB, :).'; % D_in x B
        yb = y_shuff(idxB).';    % 1 x B

        % --- Forward Pass ---
        z1 = W1 * xb + b1;
        a1 = max(0.01*z1, z1); % LeakyReLU
        
        z2 = W2 * a1 + b2;
        a2 = max(0.01*z2, z2);
        
        y_hat = W3 * a2 + b3;  % 1 x B

        % Loss (Balanced MSE with gentle 1.25x underestimation safety)
        err = (y_hat - yb);
        penalty = ones(size(err));
        penalty(err < 0) = 1.25; % Gentle underestimation protection
        grad_out = 2 * (penalty .* err) / batchSize;

        % --- Backward Pass ---
        dW3 = grad_out * a2.';
        db3 = sum(grad_out, 2);

        da2 = W3.' * grad_out;
        dz2 = da2 .* ((z2 > 0) + 0.01 * (z2 <= 0));
        dW2 = dz2 * a1.';
        db2 = sum(dz2, 2);

        da1 = W2.' * dz2;
        dz1 = da1 .* ((z1 > 0) + 0.01 * (z1 <= 0));
        dW1 = dz1 * xb.';
        db1 = sum(dz1, 2);

        % --- Adam Update ---
        % W1
        mW1 = beta1*mW1 + (1-beta1)*dW1; vW1 = beta2*vW1 + (1-beta2)*(dW1.^2);
        m_hat = mW1 / (1 - beta1^t_step); v_hat = vW1 / (1 - beta2^t_step);
        W1 = W1 - lr * m_hat ./ (sqrt(v_hat) + eps);

        % b1
        mb1 = beta1*mb1 + (1-beta1)*db1; vb1 = beta2*vb1 + (1-beta2)*(db1.^2);
        m_hat = mb1 / (1 - beta1^t_step); v_hat = vb1 / (1 - beta2^t_step);
        b1 = b1 - lr * m_hat ./ (sqrt(v_hat) + eps);

        % W2
        mW2 = beta1*mW2 + (1-beta1)*dW2; vW2 = beta2*vW2 + (1-beta2)*(dW2.^2);
        m_hat = mW2 / (1 - beta1^t_step); v_hat = vW2 / (1 - beta2^t_step);
        W2 = W2 - lr * m_hat ./ (sqrt(v_hat) + eps);

        % b2
        mb2 = beta1*mb2 + (1-beta1)*db2; vb2 = beta2*vb2 + (1-beta2)*(db2.^2);
        m_hat = mb2 / (1 - beta1^t_step); v_hat = vb2 / (1 - beta2^t_step);
        b2 = b2 - lr * m_hat ./ (sqrt(v_hat) + eps);

        % W3
        mW3 = beta1*mW3 + (1-beta1)*dW3; vW3 = beta2*vW3 + (1-beta2)*(dW3.^2);
        m_hat = mW3 / (1 - beta1^t_step); v_hat = vW3 / (1 - beta2^t_step);
        W3 = W3 - lr * m_hat ./ (sqrt(v_hat) + eps);

        % b3
        mb3 = beta1*mb3 + (1-beta1)*db3; vb3 = beta2*vb3 + (1-beta2)*(db3.^2);
        m_hat = mb3 / (1 - beta1^t_step); v_hat = vb3 / (1 - beta2^t_step);
        b3 = b3 - lr * m_hat ./ (sqrt(v_hat) + eps);
    end

    if mod(epoch, 50) == 0 || epoch == 1
        % Compute training MSE
        z1_t = W1 * X_train_norm.' + b1; a1_t = max(0.01*z1_t, z1_t);
        z2_t = W2 * a1_t + b2; a2_t = max(0.01*z2_t, z2_t);
        preds_t = W3 * a2_t + b3;
        loss_t = mean((preds_t.' - y_train).^2);
        fprintf('  Epoch %3d / %3d | Train MSE: %.4f | RMSE: %.2f taps\n', ...
            epoch, epochs, loss_t, sqrt(loss_t));
    end
end

%% 6. EVALUATION ON TEST SET (NOVELTY VS PARK 2022)
z1_test = W1 * X_test_norm.' + b1; a1_test = max(0.01*z1_test, z1_test);
z2_test = W2 * a1_test + b2; a2_test = max(0.01*z2_test, z2_test);
y_pred_continuous = (W3 * a2_test + b3).';

% Enforce integer CIR length bounds [1, 64]
y_pred_ML = max(1, min(64, round(y_pred_continuous)));

% Park 2022 Baseline predictions on Test Set (Feature #11)
y_park_SHT = X_test(:, 11);

% Metrics
mae_ML   = mean(abs(y_pred_ML - y_test));
rmse_ML  = sqrt(mean((y_pred_ML - y_test).^2));

mae_Park  = mean(abs(y_park_SHT - y_test));
rmse_Park = sqrt(mean((y_park_SHT - y_test).^2));

accuracy_ML_exact  = mean(abs(y_pred_ML - y_test) <= 1) * 100;
accuracy_Park_exact = mean(abs(y_park_SHT - y_test) <= 1) * 100;

fprintf('\n=======================================================\n');
fprintf(' TEST SET PERFORMANCE COMPARISON (Novelty vs Baseline)\n');
fprintf('=======================================================\n');
fprintf('  Metric                       Park 2022 (SHT)   Proposed ML\n');
fprintf('  ---------------------------------------------------------\n');
fprintf('  Mean Absolute Error (taps) :   %6.2f taps       %6.2f taps\n', mae_Park, mae_ML);
fprintf('  Root Mean Sq Error (taps)  :   %6.2f taps       %6.2f taps\n', rmse_Park, rmse_ML);
fprintf('  Accuracy (within +/-1 tap) :   %6.1f %%          %6.1f %%\n', accuracy_Park_exact, accuracy_ML_exact);
fprintf('=======================================================\n');

if mae_ML < mae_Park
    improvement = ((mae_Park - mae_ML) / mae_Park) * 100;
    fprintf('SUCCESS: Proposed ML reduces CIR estimation error by %.1f%% over Park 2022!\n', improvement);
end

%% 7. SAVE TRAINED MODEL
modelPath = fullfile(projectRoot, 'ML', 'trained_ml_cir_model.mat');
save(modelPath, 'W1', 'b1', 'W2', 'b2', 'W3', 'b3', 'mu', 'sigma', ...
    'mae_ML', 'mae_Park', 'rmse_ML', 'rmse_Park', '-v7');

fprintf('Saved trained model to: %s\n', modelPath);

end
