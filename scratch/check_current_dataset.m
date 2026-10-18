% Check test errors per channel in test_retrain
projectRoot = 'c:\Users\pushk\Downloads\Adaptive_OFDM_Channel_Estimation';
addpath(genpath(projectRoot));

load(fullfile(projectRoot, 'ML', 'cir_dataset.mat'));
fprintf('Loaded cir_dataset with %d samples.\n', size(featureMatrix, 1));

% Check what channels are in the dataset
% Samples 1:1050 -> SUI-2
% Samples 1051:2100 -> Vehicular-A
% Samples 2101:3150 -> SUI-5

sui2_targets = targetVector(1:1050);
veha_targets = targetVector(1051:2100);
sui5_targets = targetVector(2101:3150);

fprintf('SUI-2 targets: min=%d, max=%d, mean=%.2f\n', min(sui2_targets), max(sui2_targets), mean(sui2_targets));
fprintf('Veh-A targets: min=%d, max=%d, mean=%.2f\n', min(veha_targets), max(veha_targets), mean(veha_targets));
fprintf('SUI-5 targets: min=%d, max=%d, mean=%.2f\n', min(sui5_targets), max(sui5_targets), mean(sui5_targets));

sui2_sht = featureMatrix(1:1050, 11);
veha_sht = featureMatrix(1051:2100, 11);
sui5_sht = featureMatrix(2101:3150, 11);

fprintf('SUI-2 SHT: mean=%.2f\n', mean(sui2_sht));
fprintf('Veh-A SHT: mean=%.2f\n', mean(veha_sht));
fprintf('SUI-5 SHT: mean=%.2f\n', mean(sui5_sht));
