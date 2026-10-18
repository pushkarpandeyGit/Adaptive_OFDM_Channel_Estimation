%% ========================================================================
% DCS_PAPER_FIG6_FIG7_FIG8_MASTER.m
%
% ONE MASTER SCRIPT FOR THE PAPER-REPRODUCTION STAGE
%
% Paper:
%   H.-R. Park, "A Low-Complexity Channel Estimation for OFDM Systems
%   Based on CIR Length Adaptation," IEEE Access, 2022.
%
% PURPOSE
%   1) Reproduce Fig. 6 using the project's validated Fig. 6 result/code.
%   2) Reproduce Fig. 7 using the same estimator definitions.
%   3) Reproduce Fig. 8 using the paper's time-varying settings.
%
% IMPORTANT
%   - NO ML novelty is included.
%   - NO adaptive Hann/Hamming/Blackman novelty is included.
%   - The baseline is the paper only.
%
% FIG. 6:
%   SUI-2 and SUI-5, MSE, quasi-static.
%
% FIG. 7:
%   SUI-2 and SUI-5, uncoded 16-QAM SER, quasi-static.
%
% FIG. 8:
%   ITU-R Vehicular-A, fd=30 Hz, 2 decision-directed iterations.
%
% The project already contains a validated paper_fig6_reproduction.m and
% paper_fig6_results.mat. By default this master script loads those Fig. 6
% results instead of wasting time re-running the 1000-trial Fig. 6 simulation.
%
% Set RECOMPUTE_FIG6 = true if you want to regenerate Fig. 6.
%
% ========================================================================

clear; clc; close all;

%% ============================ USER SWITCHES ==============================
RUN_FIG6 = true;
RUN_FIG7 = true;
RUN_FIG8 = false;

RECOMPUTE_FIG6 = false;

% Monte-Carlo counts for the new Fig. 7 / Fig. 8 reproduction.
% Increase only after the curves are validated.
NUM_TRIALS_FIG7 = 400;
NUM_PACKETS_FIG8 = 200;

rng(42,'twister');

%% ============================= PROJECT ROOT ==============================
thisFile = mfilename('fullpath');

if isempty(thisFile)
    projectRoot = pwd;
else
    % This master file is intended to live in the project root or plots/.
    here = fileparts(thisFile);

    if exist(fullfile(here,'config','parameters.m'),'file')
        projectRoot = here;
    elseif exist(fullfile(fileparts(here),'config','parameters.m'),'file')
        projectRoot = fileparts(here);
    else
        error('Could not locate project root containing config/parameters.m.');
    end
end

addpath(genpath(projectRoot));

p = parameters();

fprintf('\n============================================================\n');
fprintf(' DCS PAPER REPRODUCTION MASTER\n');
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
    fprintf('Recomputing Fig. 6 with corrected paper SHT/CIR path.\n');
    fprintf('The old cached MAT file is NOT used because its adaptive-SVD\n');
    fprintf('curve was generated without the paper Hann-windowed SHT step.\n');

    SNR6 = 0:5:40;
    N6=256; Nu6=200; Np6=100; M6=128; G6=64;
    K6=56; T6=15; fixedRank6=G6+3; offset6=3;
    numMC6=1000;
    usedK6=(-100:-1).'; usedK6=[usedK6;(1:100).'];
    ltsK6=[(-100:2:-2) (2:2:100)].';

    % True second-order channel models used for Eq. (31).
    Rtrue2=local_physical_covariance(N6,usedK6,usedK6,[0 1.6 4.4],[0 -12 -15]);
    Rhp2=local_physical_covariance(N6,usedK6,ltsK6,[0 1.6 4.4],[0 -12 -15]);
    Rpp2=local_physical_covariance(N6,ltsK6,ltsK6,[0 1.6 4.4],[0 -12 -15]);

    Rhh5=local_physical_covariance(N6,usedK6,usedK6,[0 16 40],[0 -5 -10]);
    Rhp5=local_physical_covariance(N6,usedK6,ltsK6,[0 16 40],[0 -5 -10]);
    Rpp5=local_physical_covariance(N6,ltsK6,ltsK6,[0 16 40],[0 -5 -10]);

    mseL2=zeros(size(SNR6)); mseA2=zeros(size(SNR6));
    mseF2=zeros(size(SNR6)); mseW2=zeros(size(SNR6));
    mseL5=zeros(size(SNR6)); mseA5=zeros(size(SNR6));
    mseF5=zeros(size(SNR6)); mseW5=zeros(size(SNR6));

    for is=1:numel(SNR6)
        snrDb=SNR6(is);
        gammaP=4*10^(snrDb/10);  % 6-dB LTS training boost used in project baseline

        KL2=Rhp2/(Rpp2+(1/gammaP)*eye(Np6));
        KL5=Rhp5/(Rpp5+(1/gammaP)*eye(Np6));

        KF=local_svd_filter_from_L(G6,offset6,N6,usedK6,ltsK6,gammaP,G6+3);
        KW=local_build_wdft_matrix(N6,M6,ltsK6,usedK6);

        mseL2(is)=local_mse_eq31(Rtrue2,Rhp2,Rpp2,KL2,gammaP,Nu6);
        mseL5(is)=local_mse_eq31(Rhh5,Rhp5,Rpp5,KL5,gammaP,Nu6);
        mseF2(is)=local_mse_eq31(Rtrue2,Rhp2,Rpp2,KF,gammaP,Nu6);
        mseF5(is)=local_mse_eq31(Rhh5,Rhp5,Rpp5,KF,gammaP,Nu6);
        mseW2(is)=local_mse_eq31(Rtrue2,Rhp2,Rpp2,KW,gammaP,Nu6);
        mseW5(is)=local_mse_eq31(Rhh5,Rhp5,Rpp5,KW,gammaP,Nu6);

        [mseA2(is),meanL2,stdL2]=local_adaptive_mc_corrected( ...
            N6,M6,G6,K6,T6,offset6,gammaP,numMC6, ...
            usedK6,ltsK6,[0 1.6 4.4],[0 -12 -15]);

        [mseA5(is),meanL5,stdL5]=local_adaptive_mc_corrected( ...
            N6,M6,G6,K6,T6,offset6,gammaP,numMC6, ...
            usedK6,ltsK6,[0 16 40],[0 -5 -10]);

        fprintf('SNR=%2d dB | SUI2 LMMSE %.3e AD %.3e FIX %.3e WDFT %.3e | Lhat %.2f±%.2f | SUI5 AD %.3e Lhat %.2f±%.2f\n', ...
            snrDb,mseL2(is),mseA2(is),mseF2(is),mseW2(is),meanL2,stdL2,mseA5(is),meanL5,stdL5);
    end

    S=struct();
    S.SNR_values=SNR6;
    S.MSE_LMMSE_2=mseL2; S.MSE_ADAPT_2=mseA2;
    S.MSE_FIXED_2=mseF2; S.MSE_WDFT_2=mseW2;
    S.MSE_LMMSE_5=mseL5; S.MSE_ADAPT_5=mseA5;
    S.MSE_FIXED_5=mseF5; S.MSE_WDFT_5=mseW5;

    resultDir=fullfile(projectRoot,'evaluation');
    if ~exist(resultDir,'dir'), mkdir(resultDir); end
    save(fullfile(resultDir,'paper_fig6_results_corrected.mat'),'-struct','S');

    figure('Color','w','Position',[100 60 760 850]);
    subplot(2,1,1);
    semilogy(SNR6,mseL2,'k-*','LineWidth',1.2,'MarkerSize',7); hold on;
    semilogy(SNR6,mseA2,'b-s','LineWidth',1.2,'MarkerSize',5);
    semilogy(SNR6,mseF2,'r-o','LineWidth',1.2,'MarkerSize',5);
    semilogy(SNR6,mseW2,'k-v','LineWidth',1.2,'MarkerSize',6);
    grid on; box on; xlim([0 40]); ylim([1e-6 3e-1]); xticks(0:5:40);
    ylabel('Mean Square Error'); title('(a)');
    legend('LMMSE with L=5','Adaptive SVD (sim)','Fixed SVD with L=64','Windowed DFT with L=64','Location','northeast');

    subplot(2,1,2);
    semilogy(SNR6,mseL5,'k-*','LineWidth',1.2,'MarkerSize',7); hold on;
    semilogy(SNR6,mseA5,'b-s','LineWidth',1.2,'MarkerSize',5);
    semilogy(SNR6,mseF5,'r-o','LineWidth',1.2,'MarkerSize',5);
    semilogy(SNR6,mseW5,'k-v','LineWidth',1.2,'MarkerSize',6);
    grid on; box on; xlim([0 40]); ylim([1e-5 1]); xticks(0:5:40);
    xlabel('SNR (dB)'); ylabel('Mean Square Error'); title('(b)');
    legend('LMMSE with L=40','Adaptive SVD (sim)','Fixed SVD with L=64','Windowed DFT with L=64','Location','northeast');
    sgtitle('Comparison of MSEs for Quasi-Static Channels');
    exportgraphics(gcf,fullfile(resultDir,'FIG6_MASTER_CORRECTED.png'),'Resolution',300);
    fprintf('\nFig. 6 corrected run complete.\n');
end

%% =============================== FIG. 7 ================================
if RUN_FIG7

    fprintf('\n============================================================\n');
    fprintf(' FIGURE 7 - QUASI-STATIC 16-QAM SER\n');
    fprintf('============================================================\n');

    SNR7 = 10:5:40;

    % Paper: 16-QAM, 1 GHz, +/-20 ppm CFO, CPE tracking.
    % The estimator comparison itself is performed from the second LTS.
    %
    % For the reproduction baseline we use the equivalent frequency-domain
    % second-LTS observation and a pilot-based CPE tracker. The CFO estimate
    % is assumed corrected as stated in the paper; the remaining common phase
    % is tracked from the eight IEEE 802.16 pilot tones.

    mapping = subcarrier_mapping(p);
    usedK = mapping.usedIndex(:);
    dataK = mapping.dataLocations(:);
    pilotK = mapping.pilotLocations(:);

    % Exact paper LTS locations.
    ltsK = [-100:2:-2 2:2:100].';

    % Build fixed estimator bank at each SNR.
    ser = zeros(5,length(SNR7),2);

    fprintf('Trials/SNR/channel = %d\n',NUM_TRIALS_FIG7);

    for ch = 1:2

        if ch == 1
            chanName = 'SUI-2';
        else
            chanName = 'SUI-5';
        end

        fprintf('\n%s\n',chanName);

        for is = 1:length(SNR7)

            snrDb = SNR7(is);

            % IEEE 802.16 LTS training boost used by the project's
            % established reproduction.
            gammaP = 10^((snrDb+6)/10);

            % True/model covariance.
            if ch == 1
                Ltrue = 5;
            else
                Ltrue = 40;
            end

            Rtrue = local_uniform_covariance(p.N,Ltrue);

            Rhh = local_extract_covariance(Rtrue,usedK,usedK,p.N);
            Rhp = local_extract_covariance(Rtrue,usedK,ltsK,p.N);
            Rpp = local_extract_covariance(Rtrue,ltsK,ltsK,p.N);

            % LMMSE.
            K_lmmse = Rhp/(Rpp+(1/gammaP)*eye(100));

            % Fixed SVD, L=64, p=67.
            R64 = local_uniform_covariance(p.N,64);
            Rhp64 = local_extract_covariance(R64,usedK,ltsK,p.N);
            Rpp64 = local_extract_covariance(R64,ltsK,ltsK,p.N);

            [U64,D64] = eig((Rpp64+Rpp64')/2);
            [lam64,ord] = sort(real(diag(D64)),'descend');
            U64 = U64(:,ord);

            Dp64 = local_build_Dp(lam64,67,gammaP);
            K_fixed = Rhp64*U64*Dp64*U64';

            % WDFT estimator.
            K_wdft = local_build_wdft_matrix( ...
                p.N,128,ltsK,usedK);

            % Adaptive SVD filter bank.
            Kbank = cell(p.G,1);

            for L = 1:p.G
                RL = local_uniform_covariance(p.N,L);
                RhpL = local_extract_covariance(RL,usedK,ltsK,p.N);
                RppL = local_extract_covariance(RL,ltsK,ltsK,p.N);

                [UL,DL] = eig((RppL+RppL')/2);
                [lamL,ordL] = sort(real(diag(DL)),'descend');
                UL = UL(:,ordL);

                rankL = min(L+p.adaptiveSVDOffset,100);

                DpL = local_build_Dp(lamL,rankL,gammaP);

                Kbank{L} = RhpL*UL*DpL*UL';
            end

            err = zeros(5,NUM_TRIALS_FIG7);

            for tr = 1:NUM_TRIALS_FIG7

                % True stochastic channel.
                if ch == 1
                    h = local_sui2(p);
                else
                    h = local_sui5(p);
                end

                Hfull = local_cfr(h,p.N);

                Hused = Hfull(usedK+p.N/2+1);
                Hlts = Hfull(ltsK+p.N/2+1);

                % Second-LTS LS observation.
                nls = sqrt(1/(2*gammaP)) * ...
                    (randn(100,1)+1j*randn(100,1));

                HLS = Hlts+nls;

                % SHT.
                cir = local_ls_to_cir(HLS,128);
                sigma2 = mean(abs(cir(65:120)).^2);

                Lhat = local_sht(cir,sigma2,15,64);
                Lhat = max(1,min(64,Lhat));

                K_adapt = Kbank{Lhat};

                Hhat1 = K_lmmse*HLS;
                Hhat2 = K_adapt*HLS;
                Hhat3 = K_fixed*HLS;
                Hhat4 = K_wdft*HLS;

                % Generate normalized 16-QAM data.
                Ns = length(dataK);
                txSym = local_qam16(Ns);
                txSym = txSym(:);

                % Eight known BPSK pilots.
                pilotSym = ones(length(pilotK),1);

                % Construct transmitted used-subcarrier vector.
                Xused = zeros(length(usedK),1);

                [~,idD] = ismember(dataK,usedK);
                [~,idP] = ismember(pilotK,usedK);

                Xused(idD) = txSym;
                Xused(idP) = pilotSym;

                % Received data in frequency domain.
                nd = sqrt(1/(2*10^(snrDb/10))) * ...
                    (randn(length(usedK),1)+1j*randn(length(usedK),1));

                Y = Hused.*Xused + nd;

                % CPE tracking from known pilots.
                ests = {Hhat1,Hhat2,Hhat3,Hhat4};

                for ee = 1:4
                    % Each estimator already outputs the 200 used-subcarrier
                    % CFR. Do NOT interpolate it again.
                    HestUsed = ests{ee}(:);

                    [~,idP2] = ismember(pilotK,usedK);

                    z = Y(idP2)./(HestUsed(idP2).*pilotSym+eps);
                    phi = angle(sum(z));

                    Ye = Y.*exp(-1j*phi);

                    [~,idD2] = ismember(dataK,usedK);

                    Hdata = HestUsed(idD2);

                    shat = Ye(idD2)./(Hdata+eps);

                    rxSym = local_qam16_detect(shat);

                    err(ee,tr) = mean(rxSym ~= txSym);
                end

                % Perfect channel knowledge.
                z = Y(idP)./(Hused(idP).*pilotSym+eps);
                phi = angle(sum(z));

                Ye = Y.*exp(-1j*phi);
                shatPerfect = Ye(idD)./(Hused(idD)+eps);

                rxPerfect = local_qam16_detect(shatPerfect);

                err(5,tr) = mean(rxPerfect ~= txSym);
            end

            ser(:,is,ch) = mean(err,2);

            fprintf('%s  SNR=%2d dB: ',chanName,snrDb);
            fprintf('LMMSE %.3e | AD %.3e | FIX %.3e | WDFT %.3e | PERFECT %.3e\n', ...
                ser(1,is,ch),ser(2,is,ch),ser(3,is,ch), ...
                ser(4,is,ch),ser(5,is,ch));
        end
    end

    figure('Color','w','Position',[100 50 760 850]);

    subplot(2,1,1);
    semilogy(SNR7,ser(1,:,1),'k-*','LineWidth',1.0); hold on;
    semilogy(SNR7,ser(2,:,1),'b-s','LineWidth',1.0);
    semilogy(SNR7,ser(3,:,1),'r-o','LineWidth',1.0);
    semilogy(SNR7,ser(4,:,1),'k-v','LineWidth',1.0);
    semilogy(SNR7,ser(5,:,1),'m-+','LineWidth',1.0);
    grid on; box on;
    xlabel('SNR (dB)');
    ylabel('Symbol Error Rate');
    title('(a) SUI-2');
    legend('LMMSE with L=5','Adaptive SVD with Lhat', ...
        'Fixed SVD with L=64','Windowed DFT with L=64', ...
        'Perfect channel knowledge','Location','southwest');
    xlim([10 40]);

    subplot(2,1,2);
    semilogy(SNR7,ser(1,:,2),'k-*','LineWidth',1.0); hold on;
    semilogy(SNR7,ser(2,:,2),'b-s','LineWidth',1.0);
    semilogy(SNR7,ser(3,:,2),'r-o','LineWidth',1.0);
    semilogy(SNR7,ser(4,:,2),'k-v','LineWidth',1.0);
    semilogy(SNR7,ser(5,:,2),'m-+','LineWidth',1.0);
    grid on; box on;
    xlabel('SNR (dB)');
    ylabel('Symbol Error Rate');
    title('(b) SUI-5');
    legend('LMMSE with L=40','Adaptive SVD with Lhat', ...
        'Fixed SVD with L=64','Windowed DFT with L=64', ...
        'Perfect channel knowledge','Location','southwest');
    xlim([10 40]);

    sgtitle('Comparison of SERs for Quasi-Static Channels');

    save(fullfile(projectRoot,'evaluation','paper_fig7_master_results.mat'), ...
        'SNR7','ser');

    exportgraphics(gcf, ...
        fullfile(projectRoot,'evaluation','FIG7_MASTER_REPRODUCTION.png'), ...
        'Resolution',300);

    fprintf('\nFig. 7 complete.\n');
end

%% =============================== FIG. 8 ================================
if RUN_FIG8

    fprintf('\n============================================================\n');
    fprintf(' FIGURE 8 - TIME-VARYING VEHICULAR-A CHANNEL\n');
    fprintf('============================================================\n');

    fprintf(['This section uses the paper settings fd=30 Hz and two ', ...
             'decision-directed iterations.\n']);

    % The project did not contain a complete Fig. 8 simulator. This master
    % therefore reconstructs the paper's iterative decision-directed flow
    % using a correlated Vehicular-A tap process.
    %
    % It is deliberately separated from the already validated Fig. 6 code.
    %
    % IMPORTANT:
    % This is a reproduction implementation, not a claim of bit-identical
    % source-code recovery from the paper.

    SNR8 = 0:5:40;

    mse8 = zeros(4,length(SNR8));
    ser8 = zeros(4,length(SNR8));

    mapping = subcarrier_mapping(p);
    usedK = mapping.usedIndex(:);
    dataK = mapping.dataLocations(:);
    pilotK = mapping.pilotLocations(:);

    ltsK = [-100:2:-2 2:2:100].';

    fd = p.fd;
    numIter = p.numIterations;

    fprintf('fd         = %g Hz\n',fd);
    fprintf('iterations = %d\n',numIter);
    fprintf('packets/SNR = %d\n',NUM_PACKETS_FIG8);

    for is = 1:length(SNR8)

        snrDb = SNR8(is);
        gammaP = 10^((snrDb+6)/10);

        errMSE = zeros(4,NUM_PACKETS_FIG8);
        errSER = zeros(4,NUM_PACKETS_FIG8);

        for pkt = 1:NUM_PACKETS_FIG8

            % Initial channel and slowly time-varying taps.
            [h0,h1] = local_vehA_pair(p,fd);

            H0 = local_cfr(h0,p.N);
            H1 = local_cfr(h1,p.N);

            H0u = H0(usedK+p.N/2+1);
            H1u = H1(usedK+p.N/2+1);

            H0p = H0(ltsK+p.N/2+1);

            nls = sqrt(1/(2*gammaP))* ...
                (randn(100,1)+1j*randn(100,1));

            HLS = H0p+nls;

            % Start from LMMSE-like / adaptive estimators.
            R5 = local_uniform_covariance(p.N,5);
            Rhp5 = local_extract_covariance(R5,usedK,ltsK,p.N);
            Rpp5 = local_extract_covariance(R5,ltsK,ltsK,p.N);
            K1 = Rhp5/(Rpp5+(1/gammaP)*eye(100));

            R40 = local_uniform_covariance(p.N,40);
            Rhp40 = local_extract_covariance(R40,usedK,ltsK,p.N);
            Rpp40 = local_extract_covariance(R40,ltsK,ltsK,p.N);
            K2 = Rhp40/(Rpp40+(1/gammaP)*eye(100));

            R64 = local_uniform_covariance(p.N,64);
            Rhp64 = local_extract_covariance(R64,usedK,ltsK,p.N);
            Rpp64 = local_extract_covariance(R64,ltsK,ltsK,p.N);

            [U64,D64] = eig((Rpp64+Rpp64')/2);
            [lam64,ord] = sort(real(diag(D64)),'descend');
            U64 = U64(:,ord);
            K3 = Rhp64*U64*local_build_Dp(lam64,67,gammaP)*U64';

            cir = local_ls_to_cir(HLS,128);
            sigma2 = mean(abs(cir(65:120)).^2);
            Lhat = max(1,min(64,local_sht(cir,sigma2,15,64)));

            RL = local_uniform_covariance(p.N,Lhat);
            RhpL = local_extract_covariance(RL,usedK,ltsK,p.N);
            RppL = local_extract_covariance(RL,ltsK,ltsK,p.N);
            [UL,DL] = eig((RppL+RppL')/2);
            [lamL,ordL] = sort(real(diag(DL)),'descend');
            UL = UL(:,ordL);
            K4 = RhpL*UL*local_build_Dp( ...
                lamL,min(Lhat+3,100),gammaP)*UL';

            Kset = {K1,K2,K3,K4};

            % First data symbol.
            tx = local_qam16(length(dataK));
            pilotSym = ones(length(pilotK),1);

            X = zeros(length(usedK),1);
            [~,idD] = ismember(dataK,usedK);
            [~,idP] = ismember(pilotK,usedK);
            X(idD)=tx;
            X(idP)=pilotSym;

            % Use the current channel for the data symbol.
            Htrue = H1u;

            noise = sqrt(1/(2*10^(snrDb/10))) * ...
                (randn(length(usedK),1)+1j*randn(length(usedK),1));

            Y = Htrue.*X+noise;

            for ee=1:4

                Hest = Kset{ee}*HLS;
                HestU = Hest(:);

                for it=1:numIter

                    z = Y(idP)./(HestU(idP).*pilotSym+eps);
                    phi = angle(sum(z));
                    Ye = Y.*exp(-1j*phi);

                    shat = Ye(idD)./(HestU(idD)+eps);
                    det = local_qam16_detect(shat);

                    % Decision-directed update:
                    % replace data observations with detected symbols and
                    % combine them with the known pilot observations.
                    Hobs = zeros(length(usedK),1);
                    known = false(length(usedK),1);

                    Hobs(idP)=Ye(idP)./pilotSym;
                    known(idP)=true;

                    Hobs(idD)=Ye(idD)./det;
                    known(idD)=true;

                    % Least-squares channel update over all used tones,
                    % then smooth it through a low-rank physical-domain
                    % projection. This is the project-consistent
                    % decision-directed realization.
                    hDD = local_used_to_cir(Hobs,usedK,p.N);
                    hDD(65:end)=0;

                    Hdd = local_cir_to_used(hDD,usedK,p.N);

                    % Blend with previous estimate for stability.
                    HestU = 0.5*HestU + 0.5*Hdd;
                end

                mseNow = mean(abs(HestU-Htrue).^2);

                z = Y(idP)./(HestU(idP).*pilotSym+eps);
                phi = angle(sum(z));
                Ye = Y.*exp(-1j*phi);
                shat = Ye(idD)./(HestU(idD)+eps);

                rx = local_qam16_detect(shat);

                errMSE(ee,pkt)=mseNow;
                errSER(ee,pkt)=mean(rx~=tx);
            end
        end

        mse8(:,is)=mean(errMSE,2);
        ser8(:,is)=mean(errSER,2);

        fprintf('SNR=%2d dB: ',snrDb);
        fprintf('MSE AD=%.3e | FIX=%.3e | L5=%.3e | L40=%.3e\n', ...
            mse8(4,is),mse8(3,is),mse8(1,is),mse8(2,is));
    end

    figure('Color','w','Position',[100 50 760 850]);

    subplot(2,1,1);
    semilogy(SNR8,mse8(1,:),'k-*','LineWidth',1.0); hold on;
    semilogy(SNR8,mse8(2,:),'r-o','LineWidth',1.0);
    semilogy(SNR8,mse8(3,:),'b-s','LineWidth',1.0);
    semilogy(SNR8,mse8(4,:),'m-d','LineWidth',1.0);
    grid on; box on;
    xlabel('SNR (dB)');
    ylabel('Mean Square Error');
    title('(a) MSE');
    legend('LMMSE L=5','LMMSE L=40','Fixed SVD L=64', ...
        'Adaptive SVD','Location','southwest');

    subplot(2,1,2);
    semilogy(SNR8,ser8(1,:),'k-*','LineWidth',1.0); hold on;
    semilogy(SNR8,ser8(2,:),'r-o','LineWidth',1.0);
    semilogy(SNR8,ser8(3,:),'b-s','LineWidth',1.0);
    semilogy(SNR8,ser8(4,:),'m-d','LineWidth',1.0);
    grid on; box on;
    xlabel('SNR (dB)');
    ylabel('Symbol Error Rate');
    title('(b) SER');
    legend('LMMSE L=5','LMMSE L=40','Fixed SVD L=64', ...
        'Adaptive SVD','Location','southwest');

    sgtitle('Time-Varying Channel: f_d = 30 Hz, Two Iterations');

    save(fullfile(projectRoot,'evaluation','paper_fig8_master_results.mat'), ...
        'SNR8','mse8','ser8');

    exportgraphics(gcf, ...
        fullfile(projectRoot,'evaluation','FIG8_MASTER_REPRODUCTION.png'), ...
        'Resolution',300);

    fprintf('\nFig. 8 complete.\n');
end

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function R = local_uniform_covariance(N,L)

    k=(0:N-1).';
    d=k-k.';

    R=zeros(N,N);

    same=(d==0);
    R(same)=1;

    idx=~same;

    z=1j*2*pi*L.*d(idx)/N;
    R(idx)=(1-exp(-z))./z;

    R=(R+R')/2;
end

function Rsub = local_extract_covariance(R,rowK,colK,N)

    rowK=rowK(:).';
    colK=colK(:).';

    ri=mod(rowK+N/2,N)+1;
    ci=mod(colK+N/2,N)+1;

    Rsub=R(ri,ci);
end

function Dp=local_build_Dp(lambda,pRank,gamma)

    lambda=real(lambda(:));
    Np=length(lambda);

    pRank=max(1,min(round(pRank),Np));

    Dp=zeros(Np,Np);

    for k=1:pRank
        Dp(k,k)=1/(lambda(k)+1/gamma);
    end
end

function K=local_build_wdft_matrix(N,M,ltsK,usedK)

    Np=100;
    L=64;

    % Paper Eq. (23) transform.
    n=(-50:50).';

    w=0.5+0.5*cos(2*pi*n/M);

    m=(0:M-1).';
    q=(-50:50);

    G=(1/sqrt(M))*exp( ...
        1j*2*pi/M*((m-M/2)*q));

    % Retain the causal 64-sample part in the same convention used by
    % the project's validated WDFT implementation.
    Gc=G(M/2+1:M,:);

    D=zeros(101,100);

    for ii=1:101

        kk=-100+2*(ii-1);

        if kk==0
            D(ii,50)=0.5;
            D(ii,51)=0.5;
        else
            idx=find(ltsK==kk,1);
            D(ii,idx)=1;
        end
    end

    A=zeros(length(usedK),L);
    taps=(0:L-1);

    for ii=1:length(usedK)
        A(ii,:)=exp(-1j*2*pi*usedK(ii)*taps/N);
    end

    K=A*Gc*diag(w)*D;
end

function cir=local_ls_to_cir(HLS,M)

    HLS=HLS(:);

    H101=zeros(101,1);

    H101(1:50)=HLS(1:50);
    H101(51)=0.5*(HLS(50)+HLS(51));
    H101(52:101)=HLS(51:100);

    n=(-50:50).';

    w=0.5+0.5*cos(2*pi*n/M);

    Hwin=H101.*w;

    m=(0:M-1).';
    q=-50:50;

    G=(1/sqrt(M))*exp( ...
        1j*2*pi/M*((m-M/2)*q));

    cir=fftshift(G*Hwin);
    cir=cir(:);
end

function Lhat=local_sht(cir,sigma2,T,G)

    threshold=T*sigma2;

    e=abs(cir(1:G)).^2;

    Lhat=0;

    for n=1:G
        if e(n)>threshold
            Lhat=n;
        end
    end
end

function H=local_cfr(h,N)

    h=h(:);

    hp=[h;zeros(N-length(h),1)];

    H=fftshift(fft(hp,N));
    H=H(:);
end

function h=local_sui2(p)

    delays=round(p.SUI2.pathDelay_us*1e-6/p.Ts);
    powers=10.^(p.SUI2.pathPower_dB/20);

    K=10^(p.SUI2.RiceanK_dB/10);

    h=zeros(1,max(delays)+1);

    los=sqrt(K/(K+1));
    nlos=sqrt(1/(K+1))*(randn+1j*randn)/sqrt(2);

    h(delays(1)+1)=powers(1)*(los+nlos);

    for ii=2:length(delays)
        h(delays(ii)+1)=powers(ii)* ...
            (randn+1j*randn)/sqrt(2);
    end
end

function h=local_sui5(p)

    delays=round(p.SUI5.pathDelay_us*1e-6/p.Ts);
    powers=10.^(p.SUI5.pathPower_dB/20);

    h=zeros(1,max(delays)+1);

    for ii=1:length(delays)
        h(delays(ii)+1)=powers(ii)* ...
            (randn+1j*randn)/sqrt(2);
    end
end

function h=local_vehA(p,fd)

    % ITU-R Vehicular-A path set used by the paper/project.
    delays=round(p.VehA.pathDelay_us*1e-6/p.Ts);

    powers=10.^(p.VehA.pathPower_dB/20);

    h=zeros(1,max(delays)+1);

    % Independent Rayleigh paths.
    % The temporal correlation is handled by local_vehA_pair().
    for ii=1:length(delays)
        h(delays(ii)+1)=powers(ii)* ...
            (randn+1j*randn)/sqrt(2);
    end
end

function [h0,h1]=local_vehA_pair(p,fd)

    h0=local_vehA(p,fd);

    % OFDM symbol duration = N*Ts.
    Tsym=p.N*p.Ts;

    rho=besselj(0,2*pi*fd*Tsym);

    h1=zeros(size(h0));

    for ii=1:length(h0)

        innovation=(randn+1j*randn)/sqrt(2);

        % Normalize the innovation relative to the tap power.
        tapPower=abs(h0(ii))^2;

        if tapPower<eps
            tapPower=1;
        end

        h1(ii)=rho*h0(ii)+sqrt(max(0,1-rho^2))* ...
            sqrt(tapPower)*innovation;
    end
end

function y=local_qam16(N)

    levels=[-3 -1 1 3];

    a=levels(randi(4,N,1));
    b=levels(randi(4,N,1));

    y=(a+1j*b)/sqrt(10);
end

function rx=local_qam16_detect(y)

    y=y(:);

    levels=[-3 -1 1 3];

    re=real(y)*sqrt(10);
    im=imag(y)*sqrt(10);

    reOut=zeros(size(re));
    imOut=zeros(size(im));

    for ii=1:length(levels)
        mask=abs(re-levels(ii))==min(abs(re-levels),[],2);
        reOut(mask)=levels(ii);

        mask2=abs(im-levels(ii))==min(abs(im-levels),[],2);
        imOut(mask2)=levels(ii);
    end

    % Return normalized 16-QAM symbols, matching local_qam16().
    % The previous version returned raw levels (-3,-1,1,3), while the
    % transmitted symbols are normalized by sqrt(10). That caused every
    % symbol comparison to fail, including the PERFECT-channel reference.
    rx=(reOut+1j*imOut)/sqrt(10);
end

function h=local_used_to_cir(Hused,usedK,N)

    L=64;
    A=exp(-1j*2*pi/N*(usedK(:)*(0:L-1)));

    h=A\Hused;
    h=h(:);
end

function H=local_cir_to_used(h,usedK,N)

    h=h(:);

    taps=0:length(h)-1;

    A=exp(-1j*2*pi/N*(usedK(:)*taps));

    H=A*h;
end

function Kf=local_svd_filter_from_L(Lmodel,offset,N,usedK,ltsK,gammaP,rankValue)
    R=local_uniform_covariance(N,Lmodel);
    Rhp=local_extract_covariance(R,usedK,ltsK,N);
    Rpp=local_extract_covariance(R,ltsK,ltsK,N);
    [U,D]=eig((Rpp+Rpp')/2);
    [lam,ord]=sort(real(diag(D)),'descend');
    U=U(:,ord);
    Dp=local_build_Dp(lam,rankValue,gammaP);
    Kf=Rhp*U*Dp*U';
end

function [mseMean,Lmean,Lstd]=local_adaptive_mc_corrected( ...
    N,M,G,Ksht,T,offset,gammaP,numMC,usedK,ltsK,delays,pdB)

    vals=zeros(numMC,1);
    Lvals=zeros(numMC,1);

    % Use the actual paper-style SHT path: 100 LTS tones -> insert DC=0,
    % Hann window -> centered 128-point transform -> fftshift.
    for t=1:numMC
        [Hu,Hp]=local_fractional_channel(usedK,ltsK,N,delays,pdB,delays(1)==0 && numel(delays)==3 && pdB(2)==-12);

        noise=sqrt(1/(2*gammaP))*(randn(100,1)+1j*randn(100,1));
        HLS=Hp+noise;

        cir=local_ls_to_cir(HLS,M);
        sigma2=mean(abs(cir(G+1:G+Ksht)).^2);
        Lhat=local_sht(cir,sigma2,T,G);
        Lhat=max(1,min(G,Lhat));
        Lvals(t)=Lhat;

        Kf=local_svd_filter_from_L(Lhat,offset,N,usedK,ltsK,gammaP,min(Lhat+offset,100));
        Hhat=Kf*HLS;
        vals(t)=mean(abs(Hhat-Hu).^2);
    end

    mseMean=mean(vals);
    Lmean=mean(Lvals);
    Lstd=std(Lvals);
end

function [Hu,Hp]=local_fractional_channel(usedK,ltsK,N,delays,pdB,isSUI2)
    % Direct frequency-domain path model, retaining the fractional SUI-2
    % delays 1.6 and 4.4 samples rather than rounding them to integer taps.
    Pamp=10.^(pdB/20);
    c=zeros(numel(delays),1);

    if isSUI2
        Kr=10^(3/10);
        c(1)=Pamp(1)*(sqrt(Kr/(Kr+1)) + ...
            sqrt(1/(Kr+1))*(randn+1j*randn)/sqrt(2));
        for d=2:numel(delays)
            c(d)=Pamp(d)*(randn+1j*randn)/sqrt(2);
        end
    else
        for d=1:numel(delays)
            c(d)=Pamp(d)*(randn+1j*randn)/sqrt(2);
        end
    end

    Hu=zeros(numel(usedK),1);
    Hp=zeros(numel(ltsK),1);
    for d=1:numel(delays)
        Hu=Hu+c(d)*exp(-1j*2*pi*delays(d)*usedK/N);
        Hp=Hp+c(d)*exp(-1j*2*pi*delays(d)*ltsK/N);
    end
end

function R=local_physical_covariance(N,rows,cols,delays,pdB)
    P=10.^(pdB/10);
    P=P/sum(P);
    d=rows(:)-cols(:).';
    R=zeros(numel(rows),numel(cols));
    for q=1:numel(delays)
        R=R+P(q)*exp(-1j*2*pi*delays(q)*d/N);
    end
end

function mse=local_mse_eq31(Rhh,Rhp,Rpp,K,gammaP,Nu)
    C=Rpp+(1/gammaP)*eye(size(Rpp,1));
    E=Rhh-Rhp*K'-K*Rhp'+K*C*K';
    mse=max(real(trace(E))/Nu,eps);
end
