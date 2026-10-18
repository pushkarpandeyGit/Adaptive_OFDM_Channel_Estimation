% paper_fig6_reproduction.m - Figure 6 Reproduction
% FIGURE 6 REPRODUCTION
%
% Comparison of MSEs for quasi-static channels:
%
%   (a) SUI-2
%   (b) SUI-5
%
% Estimators:
%
%   1. LMMSE
%   2. Adaptive SVD
%   3. Fixed SVD
%   4. Windowed DFT
%
% Paper settings:
%
%   N       = 256
%   Nu      = 200 used subcarriers
%   Np      = 100 LTS training tones
%   G       = 64
%   M       = 128
%   T_SHT   = 15
%
%   LMMSE covariance:
%       SUI-2 : L = 5
%       SUI-5 : L = 40
%
%   Fixed SVD:
%       covariance L = 64
%       rank p = G + 3 = 67
%
%   Adaptive SVD:
%       covariance L = L_hat
%       rank p = L_hat + 3
%
%   WDFT:
%       CIR length L = 64
%
% The fixed-estimator MSEs are calculated analytically using
% the paper's Appendix Eq. (31).
%
% The adaptive-SVD MSE is obtained by Monte Carlo simulation,
% as stated in the paper.
%
% ================================================================

% clear;
clc;
close all;

%% ================================================================
% 1. FIND PROJECT ROOT
% ================================================================

thisDir = fileparts(mfilename('fullpath'));

projectRoot = thisDir;

while ~exist(fullfile(projectRoot,'config','parameters.m'),'file')

    parentDir = fileparts(projectRoot);

    if strcmp(parentDir,projectRoot)

        error( ...
            'Could not find project root containing config/parameters.m.');

    end

    projectRoot = parentDir;

end

addpath(genpath(projectRoot));

fprintf('\n');
fprintf('============================================================\n');
fprintf(' FIGURE 6 REPRODUCTION\n');
fprintf('============================================================\n');

fprintf('Project root:\n%s\n\n',projectRoot);

%% ================================================================
% 2. PARAMETERS
% ================================================================

p = parameters();

rng(p.randomSeed);

%% OFDM parameters

N  = p.N;       % 256
Nu = p.Nu;      % 200
G  = p.G;       % 64

%% LTS parameters

M  = 128;
Np = 100;

%% SHT

T_SHT = 15;

noiseStart = p.noiseTapStart;     % 65
noiseEnd   = p.noiseTapEnd;       % 120

%% SVD ranks

fixedRank = G + p.fixedSVDOffset;

fprintf('N                = %d\n',N);
fprintf('Nu               = %d\n',Nu);
fprintf('Np               = %d\n',Np);
fprintf('G                = %d\n',G);
fprintf('M                = %d\n',M);
fprintf('SHT T            = %d\n',T_SHT);
fprintf('Fixed SVD rank   = %d\n',fixedRank);

%% ================================================================
% 3. SNR RANGE
% ================================================================

SNR_values = 0:5:40;

numSNR = length(SNR_values);

%% ================================================================
% 4. LTS SUBCARRIERS
% ================================================================

ltsLocations = ...
    [-100:2:-2 2:2:100];

if length(ltsLocations) ~= 100

    error('LTS locations must contain exactly 100 tones.');

end

%% ================================================================
% 5. USED SUBCARRIERS
% ================================================================

usedLocations = p.activeIndex;

if length(usedLocations) ~= 200

    error('Used subcarrier count must be 200.');

end

%% ================================================================
% 6. TRUE CHANNEL MODELS
% ================================================================

Ltrue_SUI2 = 5;

Ltrue_SUI5 = 40;

%% ================================================================
% 7. STORAGE
% ================================================================

%% SUI-2

MSE_LMMSE_SUI2 = zeros(1,numSNR);

MSE_FIXED_SUI2 = zeros(1,numSNR);

MSE_ADAPT_SUI2 = zeros(1,numSNR);

MSE_WDFT_SUI2 = zeros(1,numSNR);

%% SUI-5

MSE_LMMSE_SUI5 = zeros(1,numSNR);

MSE_FIXED_SUI5 = zeros(1,numSNR);

MSE_ADAPT_SUI5 = zeros(1,numSNR);

MSE_WDFT_SUI5 = zeros(1,numSNR);

%% ================================================================
% 8. BUILD TRUE COVARIANCE MATRICES
% ================================================================

fprintf('\n');
fprintf('Building covariance matrices...\n');

Rtrue_SUI2 = ...
    local_uniform_covariance(N,Ltrue_SUI2);

Rtrue_SUI5 = ...
    local_uniform_covariance(N,Ltrue_SUI5);

%% Fixed SVD uses L = 64

Rmodel_fixed = ...
    local_uniform_covariance(N,G);

%% ================================================================
% 9. BUILD COVARIANCE SUBMATRICES
% ================================================================

%% SUI-2 true covariance

Rhh_SUI2 = ...
    local_extract_covariance( ...
    Rtrue_SUI2, ...
    usedLocations, ...
    usedLocations, ...
    N);

Rhp_SUI2 = ...
    local_extract_covariance( ...
    Rtrue_SUI2, ...
    usedLocations, ...
    ltsLocations, ...
    N);

Rpp_SUI2 = ...
    local_extract_covariance( ...
    Rtrue_SUI2, ...
    ltsLocations, ...
    ltsLocations, ...
    N);

%% SUI-5 true covariance

Rhh_SUI5 = ...
    local_extract_covariance( ...
    Rtrue_SUI5, ...
    usedLocations, ...
    usedLocations, ...
    N);

Rhp_SUI5 = ...
    local_extract_covariance( ...
    Rtrue_SUI5, ...
    usedLocations, ...
    ltsLocations, ...
    N);

Rpp_SUI5 = ...
    local_extract_covariance( ...
    Rtrue_SUI5, ...
    ltsLocations, ...
    ltsLocations, ...
    N);

%% Fixed SVD covariance

Rhp_fixed = ...
    local_extract_covariance( ...
    Rmodel_fixed, ...
    usedLocations, ...
    ltsLocations, ...
    N);

Rpp_fixed = ...
    local_extract_covariance( ...
    Rmodel_fixed, ...
    ltsLocations, ...
    ltsLocations, ...
    N);

%% ================================================================
% 10. EIGENDECOMPOSITION FOR FIXED SVD
% ================================================================

fprintf('Decomposing fixed-SVD covariance...\n');

[U_fixed,D_fixed] = eig(Rpp_fixed);

[eigenvalues_fixed,sortIndex] = ...
    sort(real(diag(D_fixed)),'descend');

U_fixed = U_fixed(:,sortIndex);

%% Numerical cleanup

eigenvalues_fixed( ...
    eigenvalues_fixed < 0 & ...
    eigenvalues_fixed > -1e-12) = 0;

%% ================================================================
% 11. BUILD WDFT ESTIMATOR MATRIX
% ================================================================

fprintf('Building WDFT estimator matrix...\n');

K_WDFT = ...
    local_build_wdft_matrix( ...
    N, ...
    M, ...
    ltsLocations, ...
    usedLocations);

%% ================================================================
% 12. MONTE-CARLO SETTINGS
% ================================================================

numTrials = 1000;

fprintf('\n');
fprintf('Adaptive SVD Monte-Carlo trials = %d\n',numTrials);

%% ================================================================
% 13. MAIN SNR LOOP
% ================================================================

for iSNR = 1:numSNR

    SNRdB = SNR_values(iSNR);

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf(' SNR = %d dB\n',SNRdB);
    fprintf('============================================================\n');

    %% ------------------------------------------------------------
    % Pilot SNR
    %
    % IEEE 802.16 LTS training sequence is boosted by 6 dB.
    %
    % gamma_p = 10^((SNR + 6)/10)
    % -------------------------------------------------------------

    gamma_p = ...
        10^((SNRdB + 6)/10);

    %% ============================================================
    % 13A. FIXED ESTIMATORS — SUI-2
    % ============================================================

    %% LMMSE

    K_LMMSE_SUI2 = ...
        Rhp_SUI2 / ...
        (Rpp_SUI2 + ...
        (1/gamma_p)*eye(Np));

    MSE_LMMSE_SUI2(iSNR) = ...
        local_mse_eq31( ...
        Rhh_SUI2, ...
        Rhp_SUI2, ...
        Rpp_SUI2, ...
        K_LMMSE_SUI2, ...
        gamma_p, ...
        Nu);

    %% Fixed SVD

    Dp_fixed = ...
        local_build_Dp( ...
        eigenvalues_fixed, ...
        fixedRank, ...
        gamma_p);

    K_FIXED_SUI2 = ...
        Rhp_fixed * ...
        U_fixed * ...
        Dp_fixed * ...
        U_fixed';

    MSE_FIXED_SUI2(iSNR) = ...
        local_mse_eq31( ...
        Rhh_SUI2, ...
        Rhp_SUI2, ...
        Rpp_SUI2, ...
        K_FIXED_SUI2, ...
        gamma_p, ...
        Nu);

    %% WDFT

    MSE_WDFT_SUI2(iSNR) = ...
        local_mse_eq31( ...
        Rhh_SUI2, ...
        Rhp_SUI2, ...
        Rpp_SUI2, ...
        K_WDFT, ...
        gamma_p, ...
        Nu);

    %% ============================================================
    % 13B. FIXED ESTIMATORS — SUI-5
    % ============================================================

    %% LMMSE

    K_LMMSE_SUI5 = ...
        Rhp_SUI5 / ...
        (Rpp_SUI5 + ...
        (1/gamma_p)*eye(Np));

    MSE_LMMSE_SUI5(iSNR) = ...
        local_mse_eq31( ...
        Rhh_SUI5, ...
        Rhp_SUI5, ...
        Rpp_SUI5, ...
        K_LMMSE_SUI5, ...
        gamma_p, ...
        Nu);

    %% Fixed SVD

    K_FIXED_SUI5 = ...
        Rhp_fixed * ...
        U_fixed * ...
        Dp_fixed * ...
        U_fixed';

    MSE_FIXED_SUI5(iSNR) = ...
        local_mse_eq31( ...
        Rhh_SUI5, ...
        Rhp_SUI5, ...
        Rpp_SUI5, ...
        K_FIXED_SUI5, ...
        gamma_p, ...
        Nu);

    %% WDFT

    MSE_WDFT_SUI5(iSNR) = ...
        local_mse_eq31( ...
        Rhh_SUI5, ...
        Rhp_SUI5, ...
        Rpp_SUI5, ...
        K_WDFT, ...
        gamma_p, ...
        Nu);

    %% ============================================================
    % 13C. ADAPTIVE SVD — SUI-2
    % ============================================================

    mseAdaptiveTrials = zeros(numTrials,1);

    LhatTrials = zeros(numTrials,1);

    for trial = 1:numTrials

        %% --------------------------------------------------------
        % Generate SUI-2 channel
        % ---------------------------------------------------------

        h = sui2_channel(p);

        %% --------------------------------------------------------
        % True CFR
        % ---------------------------------------------------------

        Htrue_full = ...
            local_channel_frequency_response( ...
            h,N);

        Htrue_used = ...
            Htrue_full(usedLocations + N/2 + 1);

        Htrue_lts = ...
            Htrue_full(ltsLocations + N/2 + 1);

        Htrue_used = Htrue_used(:);

        Htrue_lts = Htrue_lts(:);

        %% --------------------------------------------------------
        % Generate LS noise directly at the pilot-estimate level
        %
        % Var(n_LS) = 1/gamma_p
        % ---------------------------------------------------------

        lsNoise = ...
            sqrt(1/(2*gamma_p)) * ...
            (randn(Np,1) + ...
             1j*randn(Np,1));

        HLS = ...
            Htrue_lts + lsNoise;

        %% --------------------------------------------------------
        % LS -> CIR
        % ---------------------------------------------------------

        cir = ...
            local_ls_to_cir( ...
            HLS, ...
            M);

        %% --------------------------------------------------------
        % Noise variance
        % ---------------------------------------------------------

        sigma2_hat = ...
            mean(abs(cir(noiseStart:noiseEnd)).^2);

        %% --------------------------------------------------------
        % SHT
        % ---------------------------------------------------------

        Lhat = ...
            local_sht( ...
            cir, ...
            sigma2_hat, ...
            T_SHT, ...
            G);

        %% Avoid zero-length covariance model

        Lmodel = max(1,Lhat);

        %% Adaptive rank

        adaptiveRank = ...
            min( ...
            Lhat + p.adaptiveSVDOffset, ...
            Np);

        adaptiveRank = max(1,adaptiveRank);

        LhatTrials(trial) = Lhat;

        %% --------------------------------------------------------
        % Build adaptive covariance
        % ---------------------------------------------------------

        Radaptive = ...
            local_uniform_covariance( ...
            N,Lmodel);

        Rhp_adapt = ...
            local_extract_covariance( ...
            Radaptive, ...
            usedLocations, ...
            ltsLocations, ...
            N);

        Rpp_adapt = ...
            local_extract_covariance( ...
            Radaptive, ...
            ltsLocations, ...
            ltsLocations, ...
            N);

        %% --------------------------------------------------------
        % Eigendecomposition
        % ---------------------------------------------------------

        [U_adapt,D_adapt] = ...
            eig(Rpp_adapt);

      [lambda_adapt,idx] = ...
            sort(real(diag(D_adapt)),'descend');

        U_adapt = ...
            U_adapt(:,idx);

        lambda_adapt( ...
            lambda_adapt < 0 & ...
            lambda_adapt > -1e-12) = 0;

        %% --------------------------------------------------------
        % Adaptive filtering matrix
        % ---------------------------------------------------------

        Dp_adapt = ...
            local_build_Dp( ...
            lambda_adapt, ...
            adaptiveRank, ...
            gamma_p);

        %% --------------------------------------------------------
        % Adaptive SVD estimate
        % ---------------------------------------------------------

        Hhat = ...
            Rhp_adapt * ...
            U_adapt * ...
            Dp_adapt * ...
            U_adapt' * ...
            HLS;

        %% --------------------------------------------------------
        % MSE
        % ---------------------------------------------------------

        mseAdaptiveTrials(trial) = ...
            mean(abs(Hhat - Htrue_used).^2);

    end

    MSE_ADAPT_SUI2(iSNR) = ...
        mean(mseAdaptiveTrials);

    %% ============================================================
    % 13D. ADAPTIVE SVD — SUI-5
    % ============================================================

    mseAdaptiveTrials = zeros(numTrials,1);

    LhatTrials5 = zeros(numTrials,1);

    for trial = 1:numTrials

        %% Generate SUI-5

        h = sui5_channel(p);

        %% True CFR

        Htrue_full = ...
            local_channel_frequency_response( ...
            h,N);

        Htrue_used = ...
            Htrue_full(usedLocations + N/2 + 1);

        Htrue_lts = ...
            Htrue_full(ltsLocations + N/2 + 1);

        Htrue_used = Htrue_used(:);

        Htrue_lts = Htrue_lts(:);

        %% LS noise

        lsNoise = ...
            sqrt(1/(2*gamma_p)) * ...
            (randn(Np,1) + ...
             1j*randn(Np,1));

        HLS = ...
            Htrue_lts + lsNoise;

        %% LS -> CIR

        cir = ...
            local_ls_to_cir( ...
            HLS, ...
            M);

        %% Noise variance

        sigma2_hat = ...
            mean(abs(cir(noiseStart:noiseEnd)).^2);

        %% SHT

        Lhat = ...
            local_sht( ...
            cir, ...
            sigma2_hat, ...
            T_SHT, ...
            G);

        Lmodel = max(1,Lhat);

        %% Adaptive rank

        adaptiveRank = ...
            min( ...
            Lhat + p.adaptiveSVDOffset, ...
            Np);

        adaptiveRank = max(1,adaptiveRank);

        LhatTrials5(trial) = Lhat;

        %% Adaptive covariance

        Radaptive = ...
            local_uniform_covariance( ...
            N,Lmodel);

        Rhp_adapt = ...
            local_extract_covariance( ...
            Radaptive, ...
            usedLocations, ...
            ltsLocations, ...
            N);

        Rpp_adapt = ...
            local_extract_covariance( ...
            Radaptive, ...
            ltsLocations, ...
            ltsLocations, ...
            N);

        %% Eigen decomposition

        [U_adapt,D_adapt] = ...
            eig(Rpp_adapt);

        [lambda_adapt,idx] = ...
            sort(real(diag(D_adapt)),'descend');

        U_adapt = ...
            U_adapt(:,idx);

        lambda_adapt( ...
            lambda_adapt < 0 & ...
            lambda_adapt > -1e-12) = 0;

        %% Filtering matrix

        Dp_adapt = ...
            local_build_Dp( ...
            lambda_adapt, ...
            adaptiveRank, ...
            gamma_p);

        %% Adaptive SVD estimate

        Hhat = ...
            Rhp_adapt * ...
            U_adapt * ...
            Dp_adapt * ...
            U_adapt' * ...
            HLS;

        %% MSE

        mseAdaptiveTrials(trial) = ...
            mean(abs(Hhat - Htrue_used).^2);

    end

    MSE_ADAPT_SUI5(iSNR) = ...
        mean(mseAdaptiveTrials);

    %% ============================================================
    % 13E. PRINT CURRENT RESULT
    % ============================================================

    fprintf('\nSUI-2:\n');

    fprintf('  LMMSE       = %.6e\n', ...
        MSE_LMMSE_SUI2(iSNR));

    fprintf('  Adaptive SVD= %.6e\n', ...
        MSE_ADAPT_SUI2(iSNR));

    fprintf('  Fixed SVD   = %.6e\n', ...
        MSE_FIXED_SUI2(iSNR));

    fprintf('  WDFT        = %.6e\n', ...
        MSE_WDFT_SUI2(iSNR));

    fprintf('\nSUI-5:\n');

    fprintf('  LMMSE       = %.6e\n', ...
        MSE_LMMSE_SUI5(iSNR));

    fprintf('  Adaptive SVD= %.6e\n', ...
        MSE_ADAPT_SUI5(iSNR));

    fprintf('  Fixed SVD   = %.6e\n', ...
        MSE_FIXED_SUI5(iSNR));

    fprintf('  WDFT        = %.6e\n', ...
        MSE_WDFT_SUI5(iSNR));

    fprintf('\n');

    fprintf('Average L_hat SUI-2 = %.3f\n', ...
        mean(LhatTrials));

    fprintf('Average L_hat SUI-5 = %.3f\n', ...
        mean(LhatTrials5));

end

%% ================================================================
% 14. DISPLAY RESULT TABLE
% ================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' FIGURE 6 NUMERICAL RESULTS\n');
fprintf('============================================================\n');

fprintf('\nSUI-2\n');

fprintf( ...
    'SNR     LMMSE          Adaptive       Fixed          WDFT\n');

fprintf( ...
    '------------------------------------------------------------\n');

for i = 1:numSNR

    fprintf( ...
        '%2d   %11.4e   %11.4e   %11.4e   %11.4e\n', ...
        SNR_values(i), ...
        MSE_LMMSE_SUI2(i), ...
        MSE_ADAPT_SUI2(i), ...
        MSE_FIXED_SUI2(i), ...
        MSE_WDFT_SUI2(i));

end

fprintf('\nSUI-5\n');

fprintf( ...
    'SNR     LMMSE          Adaptive       Fixed          WDFT\n');

fprintf( ...
    '------------------------------------------------------------\n');

for i = 1:numSNR

    fprintf( ...
        '%2d   %11.4e   %11.4e   %11.4e   %11.4e\n', ...
        SNR_values(i), ...
        MSE_LMMSE_SUI5(i), ...
        MSE_ADAPT_SUI5(i), ...
        MSE_FIXED_SUI5(i), ...
        MSE_WDFT_SUI5(i));

end

%% ================================================================
% 15. FIGURE 6(a) — SUI-2
% ================================================================

figure( ...
    'Name','Fig. 6(a) - SUI-2', ...
    'Color','w');

semilogy( ...
    SNR_values, ...
    MSE_LMMSE_SUI2, ...
    '-k*', ...
    'LineWidth',1.0, ...
    'MarkerSize',7);

hold on;

semilogy( ...
    SNR_values, ...
    MSE_ADAPT_SUI2, ...
    '-bs', ...
    'LineWidth',1.0, ...
    'MarkerSize',5);

semilogy( ...
    SNR_values, ...
    MSE_FIXED_SUI2, ...
    '-ro', ...
    'LineWidth',1.0, ...
    'MarkerSize',6);

semilogy( ...
    SNR_values, ...
    MSE_WDFT_SUI2, ...
    '-kv', ...
    'LineWidth',1.0, ...
    'MarkerSize',6);

grid on;
box on;

xlim([0 40]);

ylim([1e-6 3e-1]);

xticks(0:5:40);

xlabel('SNR (dB)');

ylabel('Mean Square Error');

legend( ...
    'LMMSE with L=5', ...
    'Adaptive SVD (sim)', ...
    'Fixed SVD with L=64', ...
    'Windowed DFT with L=64', ...
    'Location','northeast');

text( ...
    0.5, ...
    -0.22, ...
    '(a)', ...
    'Units','normalized', ...
    'HorizontalAlignment','center', ...
    'FontSize',11);

set(gca, ...
    'YScale','log', ...
    'FontName','Times New Roman', ...
    'FontSize',10);

%% ================================================================
% 16. FIGURE 6(b) — SUI-5
% ================================================================

figure( ...
    'Name','Fig. 6(b) - SUI-5', ...
    'Color','w');

semilogy( ...
    SNR_values, ...
    MSE_LMMSE_SUI5, ...
    '-k*', ...
    'LineWidth',1.0, ...
    'MarkerSize',7);

hold on;

semilogy( ...
    SNR_values, ...
    MSE_ADAPT_SUI5, ...
    '-bs', ...
    'LineWidth',1.0, ...
    'MarkerSize',5);

semilogy( ...
    SNR_values, ...
    MSE_FIXED_SUI5, ...
    '-ro', ...
    'LineWidth',1.0, ...
    'MarkerSize',6);

semilogy( ...
    SNR_values, ...
    MSE_WDFT_SUI5, ...
    '-kv', ...
    'LineWidth',1.0, ...
    'MarkerSize',6);

grid on;
box on;

xlim([0 40]);

ylim([1e-5 1]);

xticks(0:5:40);

xlabel('SNR (dB)');

ylabel('Mean Square Error');

legend( ...
    'LMMSE with L=40', ...
    'Adaptive SVD (sim)', ...
    'Fixed SVD with L=64', ...
    'Windowed DFT with L=64', ...
    'Location','northeast');

text( ...
    0.5, ...
    -0.22, ...
    '(b)', ...
    'Units','normalized', ...
    'HorizontalAlignment','center', ...
    'FontSize',11);

set(gca, ...
    'YScale','log', ...
    'FontName','Times New Roman', ...
    'FontSize',10);

%% ================================================================
% 17. COMBINED PAPER-STYLE FIGURE
% ================================================================

figure( ...
    'Name','Figure 6 - Paper Style', ...
    'Color','w', ...
    'Position',[100 50 700 850]);

%% ------------------------------------------------
% (a) SUI-2
% -------------------------------------------------

subplot(2,1,1);

semilogy( ...
    SNR_values, ...
    MSE_LMMSE_SUI2, ...
    '-k*', ...
    'LineWidth',1.0, ...
    'MarkerSize',7);

hold on;

semilogy( ...
    SNR_values, ...
    MSE_ADAPT_SUI2, ...
    '-bs', ...
    'LineWidth',1.0, ...
    'MarkerSize',5);

semilogy( ...
    SNR_values, ...
    MSE_FIXED_SUI2, ...
    '-ro', ...
    'LineWidth',1.0, ...
    'MarkerSize',6);

semilogy( ...
    SNR_values, ...
    MSE_WDFT_SUI2, ...
    '-kv', ...
    'LineWidth',1.0, ...
    'MarkerSize',6);

grid on;
box on;

xlim([0 40]);
ylim([1e-6 3e-1]);

xticks(0:5:40);

ylabel('Mean Square Error');

legend( ...
    'LMMSE with L=5', ...
    'Adaptive SVD (sim)', ...
    'Fixed SVD with L=64', ...
    'Windowed DFT with L=64', ...
    'Location','northeast');

text( ...
    0.5, ...
    -0.20, ...
    '(a)', ...
    'Units','normalized', ...
    'HorizontalAlignment','center', ...
    'FontSize',11);

%% ------------------------------------------------
% (b) SUI-5
% -------------------------------------------------

subplot(2,1,2);

semilogy( ...
    SNR_values, ...
    MSE_LMMSE_SUI5, ...
    '-k*', ...
    'LineWidth',1.0, ...
    'MarkerSize',7);

hold on;

semilogy( ...
    SNR_values, ...
    MSE_ADAPT_SUI5, ...
    '-bs', ...
    'LineWidth',1.0, ...
    'MarkerSize',5);

semilogy( ...
    SNR_values, ...
    MSE_FIXED_SUI5, ...
    '-ro', ...
    'LineWidth',1.0, ...
    'MarkerSize',6);

semilogy( ...
    SNR_values, ...
    MSE_WDFT_SUI5, ...
    '-kv', ...
    'LineWidth',1.0, ...
    'MarkerSize',6);

grid on;
box on;

xlim([0 40]);
ylim([1e-5 1]);

xticks(0:5:40);

xlabel('SNR (dB)');

ylabel('Mean Square Error');

legend( ...
    'LMMSE with L=40', ...
    'Adaptive SVD (sim)', ...
    'Fixed SVD with L=64', ...
    'Windowed DFT with L=64', ...
    'Location','northeast');

text( ...
    0.5, ...
    -0.20, ...
    '(b)', ...
    'Units','normalized', ...
    'HorizontalAlignment','center', ...
    'FontSize',11);

%% ================================================================
% 18. SAVE RESULTS
% ================================================================

resultFile = ...
    fullfile( ...
    projectRoot, ...
    'evaluation', ...
    'paper_fig6_results.mat');

save( ...
    resultFile, ...
    'SNR_values', ...
    'MSE_LMMSE_SUI2', ...
    'MSE_ADAPT_SUI2', ...
    'MSE_FIXED_SUI2', ...
    'MSE_WDFT_SUI2', ...
    'MSE_LMMSE_SUI5', ...
    'MSE_ADAPT_SUI5', ...
    'MSE_FIXED_SUI5', ...
    'MSE_WDFT_SUI5');

%% Save figure

figFile = ...
    fullfile( ...
    projectRoot, ...
    'evaluation', ...
    'paper_fig6_reproduction.png');

exportgraphics( ...
    gcf, ...
    figFile, ...
    'Resolution',300);

fprintf('\n');
fprintf('============================================================\n');
fprintf(' FIGURE 6 COMPLETE\n');
fprintf('============================================================\n');

fprintf('Results:\n%s\n',resultFile);

fprintf('\nFigure:\n%s\n',figFile);


%% =================================================================
% LOCAL FUNCTION 1
% Uniform covariance matrix — Eq. (20)
% =================================================================

function R = local_uniform_covariance(N,L)

    k = (0:N-1).';

    d = k-k.';

    R = zeros(N,N);

    same = abs(d) < eps;

    R(same) = 1;

    idx = ~same;

    z = ...
        1j * 2*pi*L .* d(idx) / N;

    R(idx) = ...
        (1-exp(-z))./z;

    %% Force Hermitian symmetry numerically

    R = ...
        (R + R')/2;

end


%% =================================================================
% LOCAL FUNCTION 2
% Extract covariance submatrix
% =================================================================

function Rsub = local_extract_covariance( ...
    R, ...
    rowLocations, ...
    colLocations, ...
    N)

    rowIndex = ...
        rowLocations + N/2 + 1;

    colIndex = ...
        colLocations + N/2 + 1;

    Rsub = ...
        R(rowIndex,colIndex);

end


%% =================================================================
% LOCAL FUNCTION 3
% Build SVD filtering matrix
% =================================================================

function Dp = local_build_Dp( ...
    eigenvalues, ...
    rankValue, ...
    gamma_p)

    Np = length(eigenvalues);

    Dp = zeros(Np,Np);

    beta = 1;

    rankValue = ...
        max(1,min(rankValue,Np));

    for k = 1:rankValue

        denominator = ...
            eigenvalues(k) + beta/gamma_p;

        Dp(k,k) = ...
            1/denominator;

    end

end


%% =================================================================
% LOCAL FUNCTION 4
% MSE from Appendix Eq. (31)
% =================================================================

function mse = local_mse_eq31( ...
    Rhh_true, ...
    Rhp_true, ...
    Rpp_true, ...
    K, ...
    gamma_p, ...
    Nu)

    %% Eq. (31):
    %
    % MSE = 1/Nu tr(
    %
    % Rhh
    % - Rhp K^H
    % - K Rhp^H
    % + K (Rpp + 1/gamma_p I) K^H
    %
    % )

    Np = size(Rpp_true,1);

    noiseCov = ...
        Rpp_true + ...
        (1/gamma_p)*eye(Np);

    E = ...
        Rhh_true ...
        - Rhp_true*K' ...
        - K*Rhp_true' ...
        + K*noiseCov*K';

    mse = ...
        real(trace(E))/Nu;

end


%% =================================================================
% LOCAL FUNCTION 5
% Channel frequency response
% =================================================================

function H = local_channel_frequency_response(h,N)

    h = h(:);

    hPadded = ...
        [h; zeros(N-length(h),1)];

    %% MATLAB FFT ordering

    Hfft = fft(hPadded,N);

    %% Convert to centered subcarrier ordering

    H = fftshift(Hfft);

    H = H(:);

end


%% =================================================================
% LOCAL FUNCTION 6
% LS estimate -> paper-style 128-point CIR
% =================================================================

function cir = local_ls_to_cir(HLS,M)

    HLS = HLS(:).';

    %% -------------------------------------------------------------
    % 100 LS tones
    %
    % Reconstruct:
    %
    %   -50,...,-1,0,+1,...,+50
    %
    % DC = 0
    % -------------------------------------------------------------

    H101 = zeros(1,101);

    H101(1:50) = HLS(1:50);

    H101(52:101) = HLS(51:100);

    %% -------------------------------------------------------------
    % Hann window
    % -------------------------------------------------------------

    n = -50:50;

    w = ...
        0.5 + ...
        0.5*cos(2*pi*n/M);

    Hwindowed = ...
        H101 .* w;

    %% -------------------------------------------------------------
    % Paper Eq. (23)
    %
    % centered transform
    % -------------------------------------------------------------

    m = (0:M-1).';

    q = -50:50;

    Gmat = ...
        (1/sqrt(M)) * ...
        exp( ...
        1j*2*pi/M * ...
        ((m-M/2)*q));

    gCentered = ...
        Gmat * Hwindowed(:);

    %% Convert to MATLAB-style ordering

    cir = fftshift(gCentered);

end


%% =================================================================
% LOCAL FUNCTION 7
% SHT
% =================================================================

function Lhat = local_sht( ...
    cir, ...
    sigma2, ...
    T, ...
    G)

    threshold = ...
        T*sigma2;

    tapEnergy = ...
        abs(cir(1:G)).^2;

    Lhat = 0;

    for n = 1:G

        if tapEnergy(n) > threshold

            Lhat = n;

        end

    end

end


%% =================================================================
% LOCAL FUNCTION 8
% Windowed DFT estimator matrix
%
% This implements the same structure used by the project's
% wdft_channel_estimate.m:
%
%   100 LS tones
%       ↓
%   insert DC estimate
%       ↓
%   Hann window
%       ↓
%   128-point centered transform
%       ↓
%   retain L = 64 CIR samples
%       ↓
%   reconstruct 200 used subcarriers
%
% =================================================================
function K_WDFT = local_build_wdft_matrix( ...
    N, ...
    M, ...
    ltsLocations, ...
    usedLocations)

    % =============================================================
    % WINDOWED DFT ESTIMATOR
    %
    % HLS on 100 LTS tones
    %       ↓
    % Hann window in the 100-tone pilot domain
    %       ↓
    % 64-tap CIR projection
    %       ↓
    % frequency-domain reconstruction on 200 used tones
    %
    % Hhat = K_WDFT * HLS
    % =============================================================

    Np = length(ltsLocations);
    L  = 64;

    % Interpolate DC
    evenK = ltsLocations(:);
    idx_m2 = find(evenK == -2);
    idx_p2 = find(evenK ==  2);

    D = zeros(Np+1, Np);
    D(1:50, 1:50)     = eye(50);
    D(51, idx_m2)     = 0.5;
    D(51, idx_p2)     = 0.5;
    D(52:101, 51:100) = eye(50);

    % Frequency grid and Fourier basis matrices
    q = (-Np/2 : Np/2).'; % 101 points: -50 to +50
    ell = 0 : L-1;        % 64 causal taps: 0 to 63

    % Frequency-domain response on 101 pilot tones: H(2q) = sum_l h(l)*exp(-j*2*pi*q*l/M)
    Fp = exp(-1j * 2 * pi / M * (q * ell));

    % Synthesis matrix onto used subcarriers: H(k) = sum_l h(l)*exp(-j*2*pi*k*l/N)
    A = exp(-1j * 2 * pi / N * (usedLocations(:) * ell));

    % Pseudoinverse projection from 101 tones to 64 delay taps
    pinvFp = (Fp' * Fp) \ (Fp');

    % Gentle cosine taper on the tail 10% of taps (taps 57:64)
    taperLen = 8;
    w_taper = ones(L, 1);
    n_taper = (0:taperLen-1).';
    w_taper(end-taperLen+1:end) = 0.5 * (1 + cos(pi * n_taper / taperLen));

    % Complete calibrated linear operator: (200 x 64) * (64 x 64) * (64 x 101) * (101 x 100)
    K_WDFT = A * diag(w_taper) * pinvFp * D;

    % -------------------------------------------------------------
    % Sanity checks
    % -------------------------------------------------------------

    if any(~isfinite(K_WDFT(:)))
        error('WDFT estimator contains NaN or Inf.');
    end

    if size(K_WDFT,1) ~= length(usedLocations)
        error('WDFT output dimension is incorrect.');
    end

    if size(K_WDFT,2) ~= length(ltsLocations)
        error('WDFT input dimension is incorrect.');
    end

    fprintf('\n');
    fprintf('========== WDFT CHECK ==========\n');
    fprintf('Input pilots       : %d\n',size(K_WDFT,2));
    fprintf('Output subcarriers : %d\n',size(K_WDFT,1));
    fprintf('CIR length         : %d\n',L);
    fprintf('Matrix rank        : %d\n',rank(K_WDFT));
    fprintf('Matrix norm        : %.6e\n',norm(K_WDFT,'fro'));
    fprintf('================================\n');

end