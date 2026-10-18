%% ========================================================================
% DCS_PAPER_REPRODUCTION_MASTER_FINAL_CORRECTED.m
%
% SINGLE ENTRY-POINT BASELINE REPRODUCTION
% Project: Adaptive OFDM Channel Estimation
% Paper : H.-R. Park, IEEE Access 2022
%
% Uses the user's existing project files for the baseline and keeps the
% ML/adaptive-window novelty completely OFF.
%
% FIG. 6 : quasi-static MSE, SUI-2 / SUI-5
% FIG. 7 : quasi-static 16-QAM SER, SUI-2 / SUI-5
% FIG. 8 : time-varying Vehicular-A, fd=30 Hz, 2 DD iterations
%
% IMPORTANT:
%   The original project Fig. 6 implementation is called inside a local
%   function. This prevents its CLEAR command from destroying the master
%   RUN_FIG7/RUN_FIG8 switches.
%% ========================================================================

% clear; 
clc; close all;

%% ----------------------------- SWITCHES ---------------------------------
if ~exist('RUN_FIG6','var'), RUN_FIG6 = true; end
if ~exist('RUN_FIG7','var'), RUN_FIG7 = true; end
if ~exist('RUN_FIG8','var'), RUN_FIG8 = true; end

if ~exist('NUM_TRIALS_FIG7','var'),  NUM_TRIALS_FIG7 = 1000; end
if ~exist('NUM_PACKETS_FIG8','var'), NUM_PACKETS_FIG8 = 500; end

rng(2026,'twister');

%% --------------------------- PROJECT ROOT -------------------------------
thisFile = mfilename('fullpath');
if isempty(thisFile)
    here = pwd;
else
    here = fileparts(thisFile);
end

if exist(fullfile(here,'config','parameters.m'),'file')
    projectRoot = here;
elseif exist(fullfile(fileparts(here),'config','parameters.m'),'file')
    projectRoot = fileparts(here);
else
    error('Place this master file in the project root or plots/ folder.');
end

addpath(genpath(projectRoot));
p = parameters();

fprintf('\n============================================================\n');
fprintf(' DCS PAPER REPRODUCTION MASTER - CORRECTED BASELINE\n');
fprintf('============================================================\n');
fprintf('Project root : %s\n',projectRoot);
fprintf('ML novelty   : OFF\n');
fprintf('Fig. 6       : %d\n',RUN_FIG6);
fprintf('Fig. 7       : %d\n',RUN_FIG7);
fprintf('Fig. 8       : %d\n',RUN_FIG8);
fprintf('============================================================\n');

%% =============================== FIG. 6 ================================
if RUN_FIG6
    fprintf('\n============================================================\n');
    fprintf(' FIGURE 6 - QUASI-STATIC MSE\n');
    fprintf('============================================================\n');

    % IMPORTANT: use the existing project implementation, not the broken
    % v4 script that passed a full Rhh matrix where Rhp was required.
    S6 = run_project_fig6(projectRoot);

    fprintf('\nFig. 6 loaded from the project reproduction:\n');
    fprintf('  SUI-2: LMMSE / Adaptive SVD / Fixed SVD / WDFT\n');
    fprintf('  SUI-5: LMMSE / Adaptive SVD / Fixed SVD / WDFT\n');

    % Save a stable copy for the master.
    save(fullfile(projectRoot,'evaluation', ...
        'DCS_MASTER_FIG6_RESULTS.mat'),'-struct','S6');
end

%% =============================== FIG. 7 ================================
if RUN_FIG7
    fprintf('\n============================================================\n');
    fprintf(' FIGURE 7 - QUASI-STATIC 16-QAM SER\n');
    fprintf('============================================================\n');

    [S7,fig7] = run_fig7(projectRoot,p,NUM_TRIALS_FIG7);

    save(fullfile(projectRoot,'evaluation', ...
        'DCS_MASTER_FIG7_RESULTS.mat'),'-struct','S7');

    exportgraphics(fig7,fullfile(projectRoot,'evaluation', ...
        'DCS_MASTER_FIG7.png'),'Resolution',300);
    exportgraphics(fig7,fullfile(projectRoot,'plots', ...
        'DCS_MASTER_FIG7.png'),'Resolution',300);
    exportgraphics(fig7,fullfile(projectRoot,'evaluation', ...
        'paper_fig7_reproduction.png'),'Resolution',300);
    exportgraphics(fig7,fullfile(projectRoot,'plots', ...
        'paper_fig7_reproduction.png'),'Resolution',300);
end

%% =============================== FIG. 8 ================================
if RUN_FIG8
    fprintf('\n============================================================\n');
    fprintf(' FIGURE 8 - TIME-VARYING VEHICULAR-A\n');
    fprintf('============================================================\n');

    [S8,fig8] = run_fig8(projectRoot,p,NUM_PACKETS_FIG8);

    save(fullfile(projectRoot,'evaluation', ...
        'DCS_MASTER_FIG8_RESULTS.mat'),'-struct','S8');

    exportgraphics(fig8,fullfile(projectRoot,'evaluation', ...
        'DCS_MASTER_FIG8.png'),'Resolution',300);
    exportgraphics(fig8,fullfile(projectRoot,'plots', ...
        'DCS_MASTER_FIG8.png'),'Resolution',300);
    exportgraphics(fig8,fullfile(projectRoot,'evaluation', ...
        'paper_fig8_reproduction.png'),'Resolution',300);
    exportgraphics(fig8,fullfile(projectRoot,'plots', ...
        'paper_fig8_reproduction.png'),'Resolution',300);
end

fprintf('\n============================================================\n');
fprintf(' MASTER RUN FINISHED\n');
fprintf('============================================================\n');

%% ========================================================================
%% FIG. 6 WRAPPER
%% ========================================================================
function S = run_project_fig6(projectRoot)

    f = fullfile(projectRoot,'plots','paper_fig6_reproduction.m');
    if ~isfile(f)
        error('Missing project Fig. 6 file: %s',f);
    end

    % This script is the original project reproduction. It contains its own
    % CLEAR command, so it MUST be executed in this local function rather
    % than in the master workspace.
    run(f);

    matFile = fullfile(projectRoot,'evaluation','paper_fig6_results.mat');
    if ~isfile(matFile)
        error('Fig. 6 script did not create %s',matFile);
    end

    S = load(matFile);

    required = {'SNR_values','MSE_LMMSE_SUI2','MSE_ADAPT_SUI2', ...
        'MSE_FIXED_SUI2','MSE_WDFT_SUI2', ...
        'MSE_LMMSE_SUI5','MSE_ADAPT_SUI5', ...
        'MSE_FIXED_SUI5','MSE_WDFT_SUI5'};

    for k=1:numel(required)
        if ~isfield(S,required{k})
            error('Fig. 6 result MAT is missing field %s.',required{k});
        end
    end

end

%% ========================================================================
%% FIG. 7
%% ========================================================================
function [S,fig] = run_fig7(projectRoot,p,nTrials)

    SNRdB = 10:5:40;
    nS = numel(SNRdB);
    ser = zeros(5,nS,2);

    usedK = p.activeIndex(:);
    ltsK = [-100:2:-2 2:2:100].';
    dataK = setdiff(usedK,p.pilotLocations(:),'stable');
    pilotK = p.pilotLocations(:);

    assert(numel(usedK)==200);
    assert(numel(ltsK)==100);
    assert(numel(dataK)==192);
    assert(numel(pilotK)==8);

    fprintf('Trials/SNR/channel = %d\n',nTrials);

    for ch=1:2
        if ch==1
            name='SUI-2';
            Lmodel=5;
        else
            name='SUI-5';
            Lmodel=40;
        end

        fprintf('\n%s\n',name);

        for is=1:nS
            snr = SNRdB(is);

            % IEEE 802.16 LTS training boost used by the project/paper
            % reproduction: +6 dB relative to the data SNR.
            gammaP = 10^((snr+6)/10);
            gammaD = 10^(snr/10);

            % ----------------------------------------------------------
            % LMMSE: covariance model matched to the channel.
            % ----------------------------------------------------------
            Rm = local_uniform_covariance(p.N,Lmodel);

idxUsed = usedK + p.N/2 + 1;
idxLTS  = ltsK  + p.N/2 + 1;

Rhp = Rm(idxUsed,idxLTS);
Rpp = Rm(idxLTS,idxLTS);
            K_lmmse = Rhp/(Rpp+(1/gammaP)*eye(100));

            % ----------------------------------------------------------
            % Fixed SVD: L=64, rank=G+3=67.
            % ----------------------------------------------------------
            R64 = local_uniform_covariance(p.N,64);
            idxUsed = usedK + p.N/2 + 1;
idxLTS  = ltsK  + p.N/2 + 1;

Rhp64 = R64(idxUsed,idxLTS);
            Rpp64 = R64(idxLTS,idxLTS);
            [U,D] = eig((Rpp64+Rpp64')/2);
            [lam,ord] = sort(real(diag(D)),'descend');
            U = U(:,ord);
            delta = zeros(100,1);
            delta(1:min(67,100)) = 1./(lam(1:min(67,100))+1/gammaP);
            K_fixed = Rhp64*U*diag(delta)*U';

            % ----------------------------------------------------------
            % Adaptive SVD filter bank.
            % ----------------------------------------------------------
            Kbank = cell(p.G,1);
            for L=1:p.G
                RL = local_uniform_covariance(p.N,L);
                RhpL = RL(idxUsed,idxLTS);
                RppL = RL(idxLTS,idxLTS);
                [UL,DL] = eig((RppL+RppL')/2);
                [la,oo] = sort(real(diag(DL)),'descend');
                UL = UL(:,oo);
                r = min(L+p.adaptiveSVDOffset,100);
                dd = zeros(100,1);
                dd(1:r) = 1./(la(1:r)+1/gammaP);
                Kbank{L} = RhpL*UL*diag(dd)*UL';
            end

            % WDFT is taken directly from the user's project estimator.
            pLocal = p;
            pLocal.activeIndex = usedK.';
            pLocal.subcarrierIndex = -100:100;
            [~, infoWDFT] = wdft_channel_estimate(ones(100,1), pLocal);
            K_wdft = infoWDFT.K_WDFT;

            err = zeros(5,nTrials);

            for tr=1:nTrials

                if ch==1
                    h = sui2_channel(pLocal);
                else
                    h = sui5_channel(pLocal);
                end

                H = local_channel_cfr(h,p.N);
                Hu = H(usedK+p.N/2+1);
                Hp = H(ltsK+p.N/2+1);

                % Second LTS LS estimate.
                nLS = sqrt(1/(2*gammaP))* ...
                    (randn(100,1)+1j*randn(100,1));
                HLS = Hp+nLS;

                % Paper SHT: Hann-windowed CIR, noise taps 65:120,
                % threshold T=15.
                [cir,~,~] = ls_to_cir(HLS,128,'hann');
                sigma2 = mean(abs(cir(p.G+1:p.G+p.K)).^2);
                Lhat = sht_cir_length(cir,sigma2,p.T,p.G);
                Lhat = max(1,min(p.G,Lhat));

                H1 = K_lmmse*HLS;
                H2 = Kbank{Lhat}*HLS;
                H3 = K_fixed*HLS;
                H4 = K_wdft*HLS;
                H4 = H4(:);

                Hset = {H1(:),H2(:),H3(:),H4};

                % ------------------------------------------------------
                % 16-QAM OFDM data symbol in frequency domain.
                % ------------------------------------------------------
                tx = local_qam16_symbols(numel(dataK));
                pilot = ones(numel(pilotK),1);

                X = zeros(200,1);
                [~,idD] = ismember(dataK,usedK);
                [~,idP] = ismember(pilotK,usedK);
                X(idD)=tx;
                X(idP)=pilot;

                noise = sqrt(1/(2*gammaD))* ...
                    (randn(200,1)+1j*randn(200,1));
                Y = Hu.*X+noise;

                % ------------------------------------------------------
                % Residual CPE tracking using the eight known pilots.
                % CFO itself is assumed corrected by the STS/LTS stage,
                % exactly as stated in the paper.
                % ------------------------------------------------------
                for ee=1:4
                    He = Hset{ee};
                    z = Y(idP)./(He(idP).*pilot+eps);
                    phi = angle(sum(z));
                    Ye = Y.*exp(-1j*phi);
                    shat = Ye(idD)./(He(idD)+eps);
                    rx = local_qam16_detect(shat);
                    err(ee,tr)=mean(rx(:) ~= tx(:));
                end

                % Perfect channel reference.
                z = Y(idP)./(Hu(idP).*pilot+eps);
                phi = angle(sum(z));
                Ye = Y.*exp(-1j*phi);
                shat = Ye(idD)./(Hu(idD)+eps);
                rx = local_qam16_detect(shat);
                err(5,tr)=mean(rx(:)~=tx(:));
            end

            ser(:,is,ch)=mean(err,2);

            fprintf('%s SNR=%2d dB: LMMSE %.3e | AD %.3e | FIX %.3e | WDFT %.3e | PERFECT %.3e\n', ...
                name,snr,ser(1,is,ch),ser(2,is,ch),ser(3,is,ch), ...
                ser(4,is,ch),ser(5,is,ch));
        end
    end

    % Guard against the exact failure seen in previous versions.
    if any(~isfinite(ser(:))) || all(ser(:)==1)
        error('Fig. 7 sanity check failed: SER output is invalid.');
    end

    fig = figure('Name', 'Figure 7 - Quasi-Static 16-QAM SER', ...
                 'Color', 'w', 'Units', 'pixels', 'Position', [100 50 680 850]);

    %% (a) SUI-2
    axA = subplot(2,1,1);
    semilogy(SNRdB, ser(1,:,1), '-*k', 'LineWidth', 1.0, 'MarkerSize', 5); hold on;
    semilogy(SNRdB, ser(2,:,1), '-sb', 'LineWidth', 1.0, 'MarkerSize', 5);
    semilogy(SNRdB, ser(3,:,1), '-or', 'LineWidth', 1.0, 'MarkerSize', 5);
    semilogy(SNRdB, ser(4,:,1), '-vk', 'LineWidth', 1.0, 'MarkerSize', 5);
    semilogy(SNRdB, ser(5,:,1), '-+m', 'LineWidth', 1.0, 'MarkerSize', 5);
    grid on; box on;
    xlim([10 40]); ylim([1e-4 1]);
    xticks(10:5:40); yticks([1e-4 1e-3 1e-2 1e-1 1]);
    xlabel('SNR (dB)', 'FontSize', 10);
    ylabel('Symbol Error Rate', 'FontSize', 10);
    legend('LMMSE with L=5', 'Adaptive SVD', 'Fixed SVD with L=64', ...
           'Windowed DFT with L=64', 'Perfect channel knowledge', ...
           'Location', 'northeast', 'FontSize', 8);
    text(0.5, -0.20, '(a)', 'Units', 'normalized', ...
         'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'normal');
    set(axA, 'FontName', 'Times New Roman', 'FontSize', 10);

    %% Inset (a) zoom around 27 - 30 dB
    posA = get(axA, 'Position');
    axInsetA = axes('Position', [posA(1) + 0.08*posA(3), posA(2) + 0.12*posA(4), 0.35*posA(3), 0.40*posA(4)]);
    semilogy(SNRdB, ser(1,:,1), '-*k', 'LineWidth', 0.8, 'MarkerSize', 4); hold on;
    semilogy(SNRdB, ser(2,:,1), '-sb', 'LineWidth', 0.8, 'MarkerSize', 3.5);
    semilogy(SNRdB, ser(3,:,1), '-or', 'LineWidth', 0.8, 'MarkerSize', 3.5);
    semilogy(SNRdB, ser(4,:,1), '-vk', 'LineWidth', 0.8, 'MarkerSize', 3.5);
    semilogy(SNRdB, ser(5,:,1), '-+m', 'LineWidth', 0.8, 'MarkerSize', 4);
    xlim([27 30]); xticks(27:30);
    ylim([3e-3 10e-3]); yticks([4e-3 6e-3 8e-3 10e-3]);
    yticklabels({'4','6','8','10'});
    text(0.02, 1.05, '\times 10^{-3}', 'Units', 'normalized', 'FontSize', 7, 'FontName', 'Times New Roman');
    grid on; box on; set(axInsetA, 'FontName', 'Times New Roman', 'FontSize', 8);

    %% (b) SUI-5
    axB = subplot(2,1,2);
    semilogy(SNRdB, ser(1,:,2), '-*k', 'LineWidth', 1.0, 'MarkerSize', 5); hold on;
    semilogy(SNRdB, ser(2,:,2), '-sb', 'LineWidth', 1.0, 'MarkerSize', 5);
    semilogy(SNRdB, ser(3,:,2), '-or', 'LineWidth', 1.0, 'MarkerSize', 5);
    semilogy(SNRdB, ser(4,:,2), '-vk', 'LineWidth', 1.0, 'MarkerSize', 5);
    semilogy(SNRdB, ser(5,:,2), '-+m', 'LineWidth', 1.0, 'MarkerSize', 5);
    grid on; box on;
    xlim([10 40]); ylim([1e-4 1]);
    xticks(10:5:40); yticks([1e-4 1e-3 1e-2 1e-1 1]);
    xlabel('SNR (dB)', 'FontSize', 10);
    ylabel('Symbol Error Rate', 'FontSize', 10);
    legend('LMMSE with L=40', 'Adaptive SVD', 'Fixed SVD with L=64', ...
           'Windowed DFT with L=64', 'Perfect channel knowledge', ...
           'Location', 'northeast', 'FontSize', 8);
    text(0.5, -0.20, '(b)', 'Units', 'normalized', ...
         'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'normal');
    set(axB, 'FontName', 'Times New Roman', 'FontSize', 10);

    %% Inset (b) zoom around 27 - 30 dB
    posB = get(axB, 'Position');
    axInsetB = axes('Position', [posB(1) + 0.08*posB(3), posB(2) + 0.12*posB(4), 0.35*posB(3), 0.40*posB(4)]);
    semilogy(SNRdB, ser(1,:,2), '-*k', 'LineWidth', 0.8, 'MarkerSize', 4); hold on;
    semilogy(SNRdB, ser(2,:,2), '-sb', 'LineWidth', 0.8, 'MarkerSize', 3.5);
    semilogy(SNRdB, ser(3,:,2), '-or', 'LineWidth', 0.8, 'MarkerSize', 3.5);
    semilogy(SNRdB, ser(4,:,2), '-vk', 'LineWidth', 0.8, 'MarkerSize', 3.5);
    semilogy(SNRdB, ser(5,:,2), '-+m', 'LineWidth', 0.8, 'MarkerSize', 4);
    xlim([27 30]); xticks(27:30);
    ylim([0.008 0.022]); yticks([0.01 0.015 0.02]);
    yticklabels({'0.01','0.015','0.02'});
    grid on; box on; set(axInsetB, 'FontName', 'Times New Roman', 'FontSize', 8);

    S.SNR7 = SNRdB;
    S.ser  = ser;
end

%% ========================================================================
%% FIG. 8 - TIME-VARYING RECONSTRUCTION
%% ========================================================================
function [S,fig] = run_fig8(projectRoot,p,NUM_PACKETS)
% ========================================================================
% FIGURE 8 - PAPER REPRODUCTION
%
% Time-varying Vehicular-A channel
% fd = 30 Hz
% 2 decision-directed iterations
%
% Paper Fig. 8:
%   (a) MSE:
%       Adaptive SVD
%       Fixed SVD
%       Windowed DFT
%
%   (b) SER:
%       Perfect channel knowledge
%       Adaptive SVD
%       Fixed SVD
%       Windowed DFT
%
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' FIGURE 8 : TIME-VARYING VEHICULAR-A CHANNEL\n');
fprintf(' fd = 30 Hz, 2 DD iterations\n');
fprintf('============================================================\n');

%% ------------------------------------------------------------------------
% Paths
% -------------------------------------------------------------------------

addpath(genpath(projectRoot));

%% ------------------------------------------------------------------------
% Parameters from paper/project
% -------------------------------------------------------------------------

N       = p.N;                 % 256
G       = p.G;                 % 64
Ts      = p.Ts;                % 0.25 us
fd      = p.fd;                % 30 Hz

SNRdB   = p.SNRdB(:).';
numSNR  = numel(SNRdB);

usedK   = p.activeIndex(:);
ltsK    = [-100:2:-2 2:2:100].';
pilotK  = p.pilotLocations(:);

dataK   = setdiff(usedK,pilotK,'stable');
dataK   = dataK(:);

Nu      = numel(usedK);
Np      = numel(ltsK);
Nd      = numel(dataK);

NUM_DD_ITER = 2;

%% ------------------------------------------------------------------------
% Storage
%
% Rows:
%   1 = Adaptive SVD
%   2 = Fixed SVD
%   3 = WDFT
%
% SER additionally stores perfect CSI separately.
% -------------------------------------------------------------------------

mse = zeros(3,numSNR);
ser = zeros(4,numSNR);

pLocal = p;
pLocal.N = N;
pLocal.subcarrierIndex = (-100:100);
[~, infoWDFT] = wdft_channel_estimate(ones(Np,1), pLocal);
K_wdft = infoWDFT.K_WDFT;

%% ------------------------------------------------------------------------
% Monte-Carlo simulation
% -------------------------------------------------------------------------

for is = 1:numSNR

    snrDB = SNRdB(is);

    fprintf('\nFig.8 SNR = %2d dB\n',snrDB);

    gammaP = 10^((snrDB+6)/10);
    gammaD = 10^(snrDB/10);

    msePkt = zeros(3,NUM_PACKETS);
    serPkt = zeros(4,NUM_PACKETS);

    for pkt = 1:NUM_PACKETS

        %% ================================================================
        % 1. Generate initial Vehicular-A channel
        % ================================================================

        [h0,~] = local_vehicular_A(p,fd);

        H0 = local_channel_cfr(h0,N);

        H0u = H0(usedK + N/2 + 1);
        H0p = H0(ltsK  + N/2 + 1);

        %% ================================================================
        % 2. Generate CHANNEL AT DATA SYMBOL
        %
        % The channel evolves according to J0 Doppler correlation.
        %
        % This is deliberately generated independently for every packet.
        % ================================================================

        h1 = local_time_evolve(h0,p,fd);

        H1 = local_channel_cfr(h1,N);

        H1u = H1(usedK + N/2 + 1);

        %% ================================================================
        % 3. LTS LS estimate
        % ================================================================

        noiseLS = sqrt(1/(2*gammaP)) * ...
            (randn(Np,1) + 1j*randn(Np,1));

        HLS = H0p + noiseLS;

        %% ================================================================
        % 4. FIXED SVD
        %
        % Paper:
        %       p = G + 3 = 67
        % ================================================================

        R64 = local_uniform_covariance(N,64);

        Kfixed = local_build_svd( ...
            R64,usedK,ltsK,67,gammaP);

        Hfixed0 = Kfixed * HLS;

        %% ================================================================
        % 5. ADAPTIVE SVD
        %
        % Estimate CIR length from the LTS LS estimate.
        %
        % p_hat = L_hat + 3
        % ================================================================

        [cirLS,~,~] = ls_to_cir(HLS,128,'hann');

        noiseVar = mean(abs(cirLS(65:120)).^2);

        [Lhat,~,~,~] = sht_cir_length( ...
            cirLS,noiseVar,p.T,64);

        Lhat = max(1,min(64,Lhat));

        padaptive = min(Lhat + 3,Np);

        Radaptive = local_uniform_covariance(N,Lhat);

        Kadaptive = local_build_svd( ...
            Radaptive,usedK,ltsK,padaptive,gammaP);

        Hadaptive0 = Kadaptive * HLS;

        %% ================================================================
        % 6. WDFT INITIAL ESTIMATE
        % ================================================================

        Hwdf0 = K_wdft * HLS;
        Hwdf0 = Hwdf0(:);

        %% ================================================================
        % 7. DATA SYMBOL
        % ================================================================

        tx = local_qam16_symbols(Nd);
        tx = tx(:);

        pilot = ones(numel(pilotK),1);

        X = zeros(Nu,1);

        [~,idD] = ismember(dataK,usedK);
        [~,idP] = ismember(pilotK,usedK);

        X(idD) = tx;
        X(idP) = pilot;

        %% ================================================================
        % 8. AWGN
        % ================================================================

        noise = sqrt(1/(2*gammaD)) * ...
            (randn(Nu,1) + 1j*randn(Nu,1));

        Y = H1u .* X + noise;

        %% ================================================================
        % 9. PERFECT CSI SER
        % ================================================================

        rxPerfect = local_qam16_detect( ...
            Y(idD)./(H1u(idD)+eps));

        rxPerfect = rxPerfect(:);

        serPkt(1,pkt) = mean(rxPerfect ~= tx);

        %% ================================================================
        % 10. PROCESS THE THREE ESTIMATORS
        % ================================================================

        HestSet = {
            Hadaptive0, ...
            Hfixed0, ...
            Hwdf0
        };

        for ee = 1:3

            He = HestSet{ee};
            He = He(:);

            %% ------------------------------------------------------------
            % Two decision-directed iterations
            % -------------------------------------------------------------

            for it = 1:NUM_DD_ITER

                %% --------------------------------------------------------
                % Pilot based phase correction
                % ---------------------------------------------------------

                zPilot = Y(idP) ./ ...
                    (He(idP).*pilot + eps);

                phi = angle(sum(zPilot));

                Ycorr = Y .* exp(-1j*phi);

                %% --------------------------------------------------------
                % Detect data
                % ---------------------------------------------------------

                zData = Ycorr(idD) ./ ...
                    (He(idD)+eps);

                detData = local_qam16_detect(zData);

                detData = detData(:);

                %% --------------------------------------------------------
                % Decision-directed observations
                %
                % Combine:
                %   8 known pilots
                %   192 detected data symbols
                %
                % This gives an observation over all 200 used tones.
                % ---------------------------------------------------------

                Hobs = zeros(Nu,1);

                Hobs(idP) = Ycorr(idP) ./ ...
                    (pilot + eps);

                Hobs(idD) = Ycorr(idD) ./ ...
                    (detData + eps);

                %% --------------------------------------------------------
                % Update estimator
                %
                % The update is intentionally estimator-specific.
                % ---------------------------------------------------------

                if ee == 1

                    % ====================================================
                    % ADAPTIVE SVD
                    %
                    % Build a 100-tone pseudo-LS vector from the
                    % decision-directed observations and reapply the
                    % adaptive SVD operator.
                    % ====================================================

                    Hpseudo = local_extract_lts( ...
                        Hobs,usedK,ltsK);

                    HeNew = Kadaptive * Hpseudo;

                elseif ee == 2

                    % ====================================================
                    % FIXED SVD
                    % ====================================================

                    Hpseudo = local_extract_lts( ...
                        Hobs,usedK,ltsK);

                    HeNew = Kfixed * Hpseudo;

                else

                    % ====================================================
                    % WINDOWED DFT
                    % ====================================================

                    Hpseudo = local_extract_lts( ...
                        Hobs,usedK,ltsK);

                    HeNew = K_wdft * Hpseudo;

                end

                HeNew = HeNew(:);

                %% --------------------------------------------------------
                % Stabilized DD update
                %
                % Avoid allowing one bad decision to completely replace
                % the previous estimate.
                % ---------------------------------------------------------

                He = 0.5*He + 0.5*HeNew;

            end

            %% ============================================================
            % 11. MSE
            % ============================================================

            msePkt(ee,pkt) = mean(abs(He-H1u).^2);

            %% ============================================================
            % 12. FINAL SER
            % ============================================================

            rx = local_qam16_detect( ...
                Y(idD)./(He(idD)+eps));

            rx = rx(:);
            tx = tx(:);

            if numel(rx) ~= numel(tx)

                error(['Fig.8 SER dimension mismatch: ', ...
                    'rx=%d, tx=%d'], ...
                    numel(rx),numel(tx));

            end

            serPkt(ee+1,pkt) = mean(rx ~= tx);

        end

    end

    %% --------------------------------------------------------------------
    % Average Monte-Carlo results
    % ---------------------------------------------------------------------

    mse(:,is) = mean(msePkt,2);

    ser(:,is) = mean(serPkt,2);

    fprintf(['  Adaptive SVD = %.4e | ', ...
             'Fixed SVD = %.4e | ', ...
             'WDFT = %.4e\n'], ...
             mse(1,is),mse(2,is),mse(3,is));

    fprintf(['  Perfect SER = %.4e | ', ...
             'Adaptive = %.4e | ', ...
             'Fixed = %.4e | ', ...
             'WDFT = %.4e\n'], ...
             ser(1,is),ser(2,is),ser(3,is),ser(4,is));

end

%% =========================================================================
% FIGURE
% =========================================================================

fig = figure( ...
    'Name','Figure 8 - Time-Varying Vehicular-A', ...
    'Color','w', ...
    'Position',[100 50 720 850]);

%% -------------------------------------------------------------------------
% (a) MSE
% --------------------------------------------------------------------------

subplot(2,1,1);

semilogy(SNRdB,mse(1,:),'-sb','LineWidth',1.0,'MarkerSize',5); hold on;
semilogy(SNRdB,mse(2,:),'-or','LineWidth',1.0,'MarkerSize',6);
semilogy(SNRdB,mse(3,:),'-kv','LineWidth',1.0,'MarkerSize',6);

grid on;
box on;

xlim([10 40]);
ylim([1e-5 1]);
xticks(10:5:40);

xlabel('SNR (dB)');
ylabel('Mean Square Error');

legend( ...
    'Adaptive SVD', ...
    'Fixed SVD with L=64', ...
    'Windowed DFT with L=64', ...
    'Location','northeast');

text(0.5,-0.20,'(a)', ...
    'Units','normalized', ...
    'HorizontalAlignment','center', ...
    'FontSize',11);

set(gca,'FontName','Times New Roman','FontSize',10);

%% -------------------------------------------------------------------------
% (b) SER
% --------------------------------------------------------------------------

subplot(2,1,2);

semilogy(SNRdB,ser(1,:),'-*k','LineWidth',1.0,'MarkerSize',6); hold on;
semilogy(SNRdB,ser(2,:),'-sb','LineWidth',1.0,'MarkerSize',5);
semilogy(SNRdB,ser(3,:),'-or','LineWidth',1.0,'MarkerSize',6);
semilogy(SNRdB,ser(4,:),'-kv','LineWidth',1.0,'MarkerSize',6);

grid on;
box on;

xlim([10 40]);
ylim([1e-4 1]);
xticks(10:5:40);

xlabel('SNR (dB)');
ylabel('Symbol Error Rate');

legend( ...
    'Perfect channel knowledge', ...
    'Adaptive SVD', ...
    'Fixed SVD with L=64', ...
    'Windowed DFT with L=64', ...
    'Location','northeast');

text(0.5,-0.20,'(b)', ...
    'Units','normalized', ...
    'HorizontalAlignment','center', ...
    'FontSize',11);

set(gca,'FontName','Times New Roman','FontSize',10);

%% -------------------------------------------------------------------------
% Save results
% --------------------------------------------------------------------------

S.SNR8 = SNRdB;
S.mse  = mse;
S.ser  = ser;

evalDir = fullfile(projectRoot, 'evaluation');
if ~exist(evalDir, 'dir'), mkdir(evalDir); end
save(fullfile(evalDir,'DCS_MASTER_FIG8_RESULTS.mat'),'S');
save(fullfile(projectRoot,'plots','DCS_MASTER_FIG8_RESULTS.mat'),'S');
save(fullfile(projectRoot,'DCS_MASTER_FIG8_RESULTS.mat'),'S');

fprintf('\n');
fprintf('============================================================\n');
fprintf(' FIGURE 8 FINISHED\n');
fprintf('============================================================\n');

end


%% =========================================================================
%% LOCAL: EXTRACT 100 LTS TONES FROM 200-TONE CFR
%% =========================================================================

function HLS = local_extract_lts(Hused,usedK,ltsK)

Hused = Hused(:);
usedK = usedK(:);
ltsK  = ltsK(:);

[tf,loc] = ismember(ltsK,usedK);

if any(~tf)
    error('Could not map all LTS subcarriers.');
end

HLS = Hused(loc);
HLS = HLS(:);

end


%% =========================================================================
%% LOCAL: UNIFORM COVARIANCE
%% =========================================================================

function R = local_uniform_covariance(N,L)

k = (-N/2:N/2-1).';

d = k-k.';

R = ones(N);

idx = (d~=0);

R(idx) = ...
    (1-exp(-1j*2*pi*L*d(idx)/N)) ./ ...
    (1j*2*pi*L*d(idx)/N);

R = (R+R')/2;

end


%% =========================================================================
%% LOCAL: SVD ESTIMATOR MATRIX
%% =========================================================================

function K = local_build_svd(Rmodel,usedK,ltsK,rankP,gammaP)

N = size(Rmodel,1);

ia = usedK + N/2 + 1;
ib = ltsK  + N/2 + 1;

Rhp = Rmodel(ia,ib);
Rpp = Rmodel(ib,ib);

Rpp = (Rpp+Rpp')/2;

[U,D] = eig(Rpp);

[lambda,ord] = sort(real(diag(D)),'descend');

U = U(:,ord);

d = zeros(numel(lambda),1);

r = min(rankP,numel(lambda));

d(1:r) = ...
    1 ./ (lambda(1:r)+1/gammaP);

K = Rhp * U * diag(d) * U';

end


%% =========================================================================
%% LOCAL: QAM16 SYMBOL GENERATOR
%% =========================================================================

function x = local_qam16_symbols(n)

levels = [-3 -1 1 3];

a = levels(randi(4,n,1));
b = levels(randi(4,n,1));

x = (a(:)+1j*b(:))/sqrt(10);

end


%% =========================================================================
%% LOCAL: QAM16 DETECTOR
%% =========================================================================

function xhat = local_qam16_detect(x)

levels = [-3 -1 1 3]/sqrt(10);

x = x(:);

xr = real(x);
xi = imag(x);

[~,ir] = min(abs(xr-levels),[],2);
[~,ii] = min(abs(xi-levels),[],2);

xhat = levels(ir)+1j*levels(ii);

xhat = xhat(:);

end


%% =========================================================================
%% LOCAL: VEHICULAR-A CHANNEL
%% =========================================================================

function [h,delayS] = local_vehicular_A(p,fd)

delays = [0 0.31 0.71 1.09 1.73 2.01];

pdB = [0 -1 -9 -10 -15 -20];

delayS = round(delays*1e-6/p.Ts);

amp = 10.^(pdB/20);
amp = amp / norm(amp);

L = max(delayS)+1;

h = zeros(L,1);

for k = 1:numel(delayS)

    g = (randn+1j*randn)/sqrt(2);

    h(delayS(k)+1) = amp(k)*g;

end

end


%% =========================================================================
%% LOCAL: TIME EVOLUTION USING J0 DOPPLER CORRELATION
%% =========================================================================

function h1 = local_time_evolve(h0,p,fd)

h0 = h0(:);

Tsym = (p.N+p.G)*p.Ts;

rho = besselj(0,2*pi*fd*Tsym);

rho = max(-1,min(1,rho));

h1 = zeros(size(h0));

for k = 1:length(h0)

    q = (randn+1j*randn)/sqrt(2);

    h1(k) = ...
        rho*h0(k) + ...
        sqrt(max(0,1-rho^2))*q;

end

end


%% =========================================================================
%% LOCAL: CHANNEL CFR
%% =========================================================================

function H = local_channel_cfr(h,N)

h = h(:);

if numel(h) > N
    error('Channel length exceeds FFT size.');
end

ht = [h;zeros(N-numel(h),1)];

H = fftshift(fft(ht,N));

H = H(:);

end