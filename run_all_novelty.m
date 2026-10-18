function run_all_novelty()
% =========================================================================
% MASTER NOVELTY RUNNER: ADAPTIVE OFDM CHANNEL ESTIMATION
%
% Authors / Project: Adaptive OFDM Channel Estimation (Novelty Suite)
% Reference Base Paper: H.-R. Park, IEEE Access 2022
%
% Executes the complete novelty pipeline from start to finish:
%   [Day 1]
%     Step 1: ML Dataset Generation (SUI-2, Vehicular-A, SUI-5)
%     Step 2: Native MLP Neural Predictor Training
%     Step 3: Plot 1 - CIR Length Estimation Error Distribution
%   [Day 2]
%     Step 4: Plot 2 - Channel MSE vs SNR Evaluation
%     Step 5: Plot 3 - 16-QAM Symbol Error Rate (SER) Evaluation
%   [Day 3]
%     Step 6: Plot 4 - Time-Varying Doppler Channel Tracking (fd = 30 Hz)
%     Step 7: Plot 5 - ML Neural Diagnostics & Explainability Suite
%     Step 8: Plot 6 - Adaptive Windowing CIR Energy Compaction
%   [Day 4]
%     Step 9: Plot 7 - Joint System Evaluation (ML + Adaptive Windowing)
%     Step 10: Upgraded Table 2 - Computational Complexity & Timing
%
% Output: All PNG figures and .mat files in ML/ and novelty/
% =========================================================================

clc; close all;
fprintf('=========================================================================\n');
fprintf('  STARTING MASTER NOVELTY PIPELINE EXECUTION\n');
fprintf('=========================================================================\n');

thisFile = mfilename('fullpath');
if isempty(thisFile), here = pwd; else, here = fileparts(thisFile); end
addpath(genpath(here));

tStart = tic;

%% STEP 1: DATASET GENERATION
fprintf('\n>>> STEP 1: Generating Training CIR Dataset...\n');
step1_generate_cir_dataset();

%% STEP 2: NEURAL MODEL TRAINING
fprintf('\n>>> STEP 2: Training Native MATLAB Neural CIR Predictor...\n');
step2_train_cir_predictor();

%% STEP 3: PLOT 1 (CIR LENGTH DISTRIBUTION)
fprintf('\n>>> STEP 3: Plot 1 - Evaluating CIR Length Estimation Distribution...\n');
step3_plot_novelty_comparison();

%% STEP 4: PLOT 2 (CHANNEL MSE VS SNR)
fprintf('\n>>> STEP 4: Plot 2 - Evaluating Channel Estimation MSE vs SNR...\n');
step4_evaluate_ml_adaptive_svd();

%% STEP 5: PLOT 3 (16-QAM SER VS SNR)
fprintf('\n>>> STEP 5: Plot 3 - Evaluating 16-QAM Symbol Error Rate vs SNR...\n');
step5_evaluate_ml_16qam_ser();

%% STEP 6: PLOT 4 (TIME-VARYING TRACKING)
fprintf('\n>>> STEP 6: Plot 4 - Evaluating Time-Varying Doppler Tracking (fd=30Hz)...\n');
step6_evaluate_ml_time_varying();

%% STEP 7: PLOT 5 (ML NEURAL DIAGNOSTICS)
fprintf('\n>>> STEP 7: Plot 5 - Generating ML Neural Diagnostics Suite...\n');
step7_plot_ml_diagnostics();

%% STEP 8: PLOT 6 (ADAPTIVE WINDOWING LEAKAGE)
fprintf('\n>>> STEP 8: Plot 6 - Evaluating Adaptive Windowing on CIR Compaction...\n');
step8_plot_adaptive_window_effect();

%% STEP 9: PLOT 7 (JOINT SYSTEM MSE)
fprintf('\n>>> STEP 9: Plot 7 - Evaluating Joint System (ML + Adaptive Windowing)...\n');
step9_evaluate_joint_system();

%% STEP 10: UPGRADED TABLE 2 (COMPLEXITY)
fprintf('\n>>> STEP 10: Benchmarking Computational Complexity & Hardware Footprint...\n');
step10_complexity_analysis();

tElapsed = toc(tStart);

fprintf('\n=========================================================================\n');
fprintf('  MASTER NOVELTY PIPELINE COMPLETED SUCCESSFULLY IN %.1f SECONDS!\n', tElapsed);
fprintf('=========================================================================\n');
fprintf('Generated Figures:\n');
fprintf('  1. ML/Novelty_Day1_ML_vs_Park.png\n');
fprintf('  2. ML/Novelty_Day2_MSE_vs_SNR.png\n');
fprintf('  3. ML/Novelty_Day2_Plot3_SER_vs_SNR.png\n');
fprintf('  4. ML/Novelty_Day3_Plot4_TimeVarying_Tracking.png\n');
fprintf('  5. ML/Novelty_Day3_Plot5_ML_Diagnostics.png\n');
fprintf('  6. novelty/Novelty_Day3_Plot6_Adaptive_Windowing_Effect.png\n');
fprintf('  7. novelty/Novelty_Day4_Plot7_Joint_System_MSE.png\n');
fprintf('Generated Tables:\n');
fprintf('  1. novelty/Novelty_Table2_Complexity.csv\n');
fprintf('=========================================================================\n');

end
