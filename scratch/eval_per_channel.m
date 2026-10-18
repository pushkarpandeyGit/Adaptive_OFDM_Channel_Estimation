projectRoot = 'c:\Users\pushk\Downloads\Adaptive_OFDM_Channel_Estimation';
load(fullfile(projectRoot, 'ML', 'cir_dataset.mat'));
model = load(fullfile(projectRoot, 'ML', 'trained_ml_cir_model.mat'));

X = featureMatrix;
y = targetVector;

% Split 80/20 with rng(2026) exactly as in step2
rng(2026);
perm = randperm(size(X, 1));
testIdx = perm(round(0.8*size(X, 1))+1 : end);

X_test = X(testIdx, :);
y_test = y(testIdx);

X_test_norm = (X_test - model.mu) ./ model.sigma;
z1 = model.W1 * X_test_norm.' + model.b1; a1 = max(0.01*z1, z1);
z2 = model.W2 * a1 + model.b2;            a2 = max(0.01*z2, z2);
y_pred = max(1, min(64, round((model.W3 * a2 + model.b3).')));

% Which channels are in test set?
is_sui2 = (testIdx <= 1050);
is_veha = (testIdx > 1050 & testIdx <= 2100);
is_sui5 = (testIdx > 2100);

fprintf('SUI-2 Test (N=%d): MAE SHT=%.2f, MAE ML=%.2f, Mean True=%.2f, Mean Pred=%.2f\n', ...
    sum(is_sui2), mean(abs(X_test(is_sui2, 11) - y_test(is_sui2))), mean(abs(y_pred(is_sui2) - y_test(is_sui2))), ...
    mean(y_test(is_sui2)), mean(y_pred(is_sui2)));

fprintf('Veh-A Test (N=%d): MAE SHT=%.2f, MAE ML=%.2f, Mean True=%.2f, Mean Pred=%.2f\n', ...
    sum(is_veha), mean(abs(X_test(is_veha, 11) - y_test(is_veha))), mean(abs(y_pred(is_veha) - y_test(is_veha))), ...
    mean(y_test(is_veha)), mean(y_pred(is_veha)));

fprintf('SUI-5 Test (N=%d): MAE SHT=%.2f, MAE ML=%.2f, Mean True=%.2f, Mean Pred=%.2f\n', ...
    sum(is_sui5), mean(abs(X_test(is_sui5, 11) - y_test(is_sui5))), mean(abs(y_pred(is_sui5) - y_test(is_sui5))), ...
    mean(y_test(is_sui5)), mean(y_pred(is_sui5)));
