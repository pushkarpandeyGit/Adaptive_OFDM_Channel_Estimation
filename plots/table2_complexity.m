%% ================================================================
% paper_table2_complexity.m
%
% TABLE 2
%
% Complexity of channel estimators per OFDM symbol.
%
% Values reproduced from the paper.
% No Monte-Carlo simulation is required.
% ================================================================

% clear;
clc;
close all;

%% ================================================================
% COMPLEXITY VALUES FROM PAPER TABLE 2
% ================================================================

estimators = { ...
    'LMMSE'; ...
    'Fixed SVD'; ...
    'Adaptive SVD'; ...
    'Windowed DFT' ...
    };

CE_via_preamble = [ ...
    19392; ...
    25920; ...
    3264; ...
    2000 ...
    ];

Iterative_CE = [ ...
    77184; ...
    51840; ...
    10368; ...
    5808 ...
    ];

%% ================================================================
% CREATE MATLAB TABLE
% ================================================================

T = table( ...
    estimators, ...
    CE_via_preamble, ...
    Iterative_CE, ...
    'VariableNames', { ...
        'Estimator', ...
        'CE_via_preamble', ...
        'Iterative_CE' ...
        });

%% ================================================================
% DISPLAY
% ================================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf('TABLE 2: COMPLEXITY OF CHANNEL ESTIMATORS\n');
fprintf('====================================================\n\n');

disp(T);

%% ================================================================
% PRINT PAPER-STYLE TABLE
% ================================================================

fprintf('\n');
fprintf('------------------------------------------------------------\n');
fprintf('%-20s %18s %18s\n', ...
    'Estimators', ...
    'CE via preamble', ...
    'Iterative CE');

fprintf('------------------------------------------------------------\n');

for k = 1:height(T)

    fprintf('%-20s %18d %18d\n', ...
        T.Estimator{k}, ...
        T.CE_via_preamble(k), ...
        T.Iterative_CE(k));

end

fprintf('------------------------------------------------------------\n');

%% ================================================================
% OPTIONAL: COMPLEXITY REDUCTION OF ADAPTIVE SVD
% COMPARED WITH FIXED SVD
% ================================================================

fixedPreamble = CE_via_preamble(2);
adaptivePreamble = CE_via_preamble(3);

fixedIterative = Iterative_CE(2);
adaptiveIterative = Iterative_CE(3);

reductionPreamble = ...
    100 * (fixedPreamble - adaptivePreamble) ...
    / fixedPreamble;

reductionIterative = ...
    100 * (fixedIterative - adaptiveIterative) ...
    / fixedIterative;

fprintf('\n');
fprintf('Adaptive SVD vs Fixed SVD:\n');

fprintf('CE via preamble reduction = %.2f %%\n', ...
    reductionPreamble);

fprintf('Iterative CE reduction    = %.2f %%\n', ...
    reductionIterative);

%% ================================================================
% OPTIONAL: SAVE TABLE
% ================================================================

save('paper_table2_complexity.mat', ...
    'T');

fprintf('\nTable 2 data saved as:\n');
fprintf('paper_table2_complexity.mat\n');

%% ================================================================
% OPTIONAL CSV EXPORT
% ================================================================

writetable(T, ...
    'paper_table2_complexity.csv');

fprintf('CSV saved as:\n');
fprintf('paper_table2_complexity.csv\n');

fprintf('\nTABLE 2 COMPLETE.\n');