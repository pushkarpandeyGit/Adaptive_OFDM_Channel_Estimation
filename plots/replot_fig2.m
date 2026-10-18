function fig2 = replot_fig2()
% REPLOT_FIG2 - Redraws Figure 2 instantly from saved data without re-simulating

matFile = which('paper_fig2_results.mat');
if isempty(matFile)
    % Look in common folders
    if exist(fullfile(pwd,'evaluation','paper_fig2_results.mat'),'file')
        matFile = fullfile(pwd,'evaluation','paper_fig2_results.mat');
    elseif exist(fullfile(pwd,'plots','paper_fig2_results.mat'),'file')
        matFile = fullfile(pwd,'plots','paper_fig2_results.mat');
    elseif exist('paper_fig2_results.mat','file')
        matFile = 'paper_fig2_results.mat';
    else
        error('paper_fig2_results.mat not found. Please run paper_fig2_reproduction once.');
    end
end

data = load(matFile);
tap = data.tap;
E2_rect_noisy = data.E2_rect_noisy;
E2_hann_noisy = data.E2_hann_noisy;
E2_rect_clean = data.E2_rect_clean;
E2_hann_clean = data.E2_hann_clean;
E5_rect_noisy = data.E5_rect_noisy;
E5_hann_noisy = data.E5_hann_noisy;
E5_rect_clean = data.E5_rect_clean;
E5_hann_clean = data.E5_hann_clean;

fig2 = figure('Name','Figure 2 - Effect of Windowing','Color','w', ...
              'Units','pixels','Position',[100 50 650 800]);

% (a) SUI-2 model
subplot(2,1,1);
plot(tap, E2_rect_noisy, '-b', 'LineWidth', 1.0); hold on;
plot(tap, E2_hann_noisy, '-k', 'LineWidth', 1.0);
plot(tap, E2_rect_clean, '-.b', 'LineWidth', 0.9);
plot(tap, E2_hann_clean, '-.k', 'LineWidth', 0.9);
grid on; box on;
xlim([0 128]); ylim([-40 50]);
xlabel('Time Index (T_s)', 'FontSize', 10);
ylabel('Tap Energy (dB)', 'FontSize', 10);
legend('Rectangular window', 'Hann window', ...
       'Rectangular window (noise-free)', 'Hann window (noise-free)', ...
       'Location', 'northeast', 'FontSize', 8);
text(0.5, -0.22, '(a)', 'Units', 'normalized', ...
     'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'normal');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);

% (b) SUI-5 model
subplot(2,1,2);
plot(tap, E5_rect_noisy, '-b', 'LineWidth', 1.0); hold on;
plot(tap, E5_hann_noisy, '-k', 'LineWidth', 1.0);
plot(tap, E5_rect_clean, '-.b', 'LineWidth', 0.9);
plot(tap, E5_hann_clean, '-.k', 'LineWidth', 0.9);
grid on; box on;
xlim([0 128]); ylim([-40 50]);
xlabel('Time Index (T_s)', 'FontSize', 10);
ylabel('Tap Energy (dB)', 'FontSize', 10);
legend('Rectangular window', 'Hann window', ...
       'Rectangular window (noise-free)', 'Hann window (noise-free)', ...
       'Location', 'northeast', 'FontSize', 8);
text(0.5, -0.22, '(b)', 'Units', 'normalized', ...
     'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'normal');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);

fprintf('Figure 2 redrawn successfully.\n');
end
