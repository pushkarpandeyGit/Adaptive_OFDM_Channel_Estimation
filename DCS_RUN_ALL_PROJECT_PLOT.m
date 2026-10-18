%% ========================================================================
% DCS_RUN_ALL_PROJECT_PLOT.m
%
% ONE-CLICK DCS PROJECT PIPELINE - FULL IEEE ACCESS REPRODUCTION
%
% Project:
%   Adaptive OFDM Channel Estimation
%
% Paper:
%   H.-R. Park, "A Low-Complexity Channel Estimation for OFDM Systems
%   Based on CIR Length Adaptation," IEEE Access, 2022.
%
% It generates:
%   Fig. 2  - Effect of windowing on CIR length estimation (SUI-2 & SUI-5)
%   Fig. 3  - Distributions of CIR length estimates (SHT, SOT, GAIC)
%   Fig. 4  - SHT threshold comparison (T = 5, 15, 25)
%   Fig. 5  - Average CIR-length comparison between SHT and SOT
%   Fig. 6  - Quasi-static MSE comparison (LMMSE, Adaptive SVD, Fixed SVD, WDFT)
%   Fig. 7  - Quasi-static 16-QAM SER comparison with insets
%   Fig. 8  - Time-varying Vehicular-A MSE and SER (fd = 30 Hz, 2 DD iterations)
%   Table 2 - Complexity comparison
% ========================================================================

clear; clc; close all;

%% ----------------------------- SETTINGS ---------------------------------
RUN_FIG2   = true;
RUN_FIG3   = true;
RUN_FIG4   = true;
RUN_FIG5   = true;
RUN_FIG6   = true;
RUN_FIG7   = true;
RUN_FIG8   = true;
RUN_TABLE2 = true;
%% -------------------------------------------------------------------------

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
    error(['Put DCS_RUN_ALL_PROJECT_PLOT.m in the project root ' ...
           'or make sure it is run from the project root.']);
end

addpath(genpath(projectRoot));

evalDir  = fullfile(projectRoot,'evaluation');
plotsDir = fullfile(projectRoot,'plots');

if ~exist(evalDir,'dir'),  mkdir(evalDir);  end
if ~exist(plotsDir,'dir'), mkdir(plotsDir); end

fprintf('\n============================================================\n');
fprintf(' DCS PROJECT - ONE CLICK PLOT PIPELINE\n');
fprintf('============================================================\n');
fprintf('Project root: %s\n',projectRoot);
fprintf('Fig 2  : %d\n',RUN_FIG2);
fprintf('Fig 3  : %d\n',RUN_FIG3);
fprintf('Fig 4  : %d\n',RUN_FIG4);
fprintf('Fig 5  : %d\n',RUN_FIG5);
fprintf('Fig 6  : %d\n',RUN_FIG6);
fprintf('Fig 7  : %d\n',RUN_FIG7);
fprintf('Fig 8  : %d\n',RUN_FIG8);
fprintf('Table 2: %d\n',RUN_TABLE2);
fprintf('============================================================\n');

%% ------------------------------- FIG. 2 ---------------------------------
if RUN_FIG2
    fprintf('\n>>> FIGURE 2\n');
    run_and_save_fig2(projectRoot);
end

%% ------------------------------- FIG. 3 ---------------------------------
if RUN_FIG3
    fprintf('\n>>> FIGURE 3\n');
    run_and_save_fig3(projectRoot);
end

%% ------------------------------- FIG. 4 ---------------------------------
if RUN_FIG4
    fprintf('\n>>> FIGURE 4\n');
    run_and_save_fig4(projectRoot);
end

%% ------------------------------- FIG. 5 ---------------------------------
if RUN_FIG5
    fprintf('\n>>> FIGURE 5\n');
    run_and_save_fig5(projectRoot);
end

%% ------------------------------- FIG. 6 ---------------------------------
if RUN_FIG6
    fprintf('\n>>> FIGURE 6\n');
    run_and_save_fig6(projectRoot);
end

%% ------------------------------- FIG. 7/8 -------------------------------
if RUN_FIG7 || RUN_FIG8
    fprintf('\n>>> FIGURES 7/8\n');
    run_existing_fig7_fig8_master(projectRoot, RUN_FIG7, RUN_FIG8);
end

%% ------------------------------ TABLE 2 ---------------------------------
if RUN_TABLE2
    fprintf('\n>>> TABLE 2\n');
    run_table2(projectRoot);
end

fprintf('\n============================================================\n');
fprintf(' DCS PROJECT PIPELINE FINISHED\n');
fprintf('============================================================\n');
fprintf('Check:\n');
fprintf('  %s\n',evalDir);
fprintf('for generated PNG/MAT/CSV outputs.\n');
fprintf('============================================================\n');


%% ========================================================================
% LOCAL FUNCTIONS
%% ========================================================================

function run_and_save_fig2(projectRoot)

f = fullfile(projectRoot,'plots','paper_fig2_reproduction.m');
if ~isfile(f)
    f = fullfile(projectRoot,'evaluation','paper_fig2_reproduction.m');
end
assert(isfile(f),'Missing Fig. 2 script: %s',f);

run(f);

pngFile = fullfile(projectRoot,'evaluation','paper_fig2_reproduction.png');
assert(isfile(pngFile), 'Fig. 2 PNG was not generated: %s', pngFile);

fprintf('Fig. 2: official PNG result confirmed.\n');

end


function run_and_save_fig3(projectRoot)

f = fullfile(projectRoot,'evaluation','paper_fig3_reproduction.m');
if ~isfile(f)
    f = fullfile(projectRoot,'plots','paper_fig3_reproduction.m');
end
assert(isfile(f),'Missing Fig. 3 script: %s',f);

run(f);

pngFile = fullfile(projectRoot,'evaluation','paper_fig3_reproduction.png');
assert(isfile(pngFile), 'Fig. 3 PNG was not generated: %s', pngFile);

fprintf('Fig. 3: official PNG + MAT results confirmed.\n');

end


function run_and_save_fig4(projectRoot)

f = fullfile(projectRoot,'evaluation','paper_fig4_reproduction.m');
if ~isfile(f)
    f = fullfile(projectRoot,'plots','plot_fig4_reproduction.m');
end
assert(isfile(f),'Missing Fig. 4 script: %s',f);

run(f);

pngFile = fullfile(projectRoot,'evaluation','paper_fig4_reproduction.png');
assert(isfile(pngFile), 'Fig. 4 PNG was not generated: %s', pngFile);

fprintf('Fig. 4: official PNG + MAT results confirmed.\n');

end


function run_and_save_fig5(projectRoot)

pRoot = projectRoot;
f = fullfile(pRoot,'plots','plot_fig5_reproduction.m');
assert(isfile(f),'Missing Fig. 5 script: %s',f);

run(f);
if ~exist('projectRoot','var') || isempty(projectRoot), projectRoot = pRoot; end

required = { ...
    fullfile(projectRoot,'evaluation','paper_fig5_reproduction.png'), ...
    fullfile(projectRoot,'evaluation','paper_fig5_results.mat')};

for k = 1:numel(required)
    if ~isfile(required{k})
        error('Fig. 5 expected output was not created: %s',required{k});
    end
end

fprintf('Fig. 5: official PNG + MAT results confirmed.\n');

end


function run_and_save_fig6(projectRoot)

pRoot = projectRoot;
f = fullfile(pRoot,'plots','paper_fig6_reproduction.m');
assert(isfile(f),'Missing Fig. 6 script: %s',f);

run(f);
if ~exist('projectRoot','var') || isempty(projectRoot), projectRoot = pRoot; end

matFile = fullfile(projectRoot,'evaluation','paper_fig6_results.mat');
assert(isfile(matFile), 'Fig. 6 script did not create %s',matFile);

pngFile = fullfile(projectRoot,'evaluation','paper_fig6_reproduction.png');
if ~isfile(pngFile)
    figs = sort_figures(findall(0,'Type','figure'));
    if ~isempty(figs)
        exportgraphics(figs(end),pngFile,'Resolution',300);
    end
end

fprintf('Fig. 6: official MAT result confirmed.\n');

end


function run_existing_fig7_fig8_master(projectRoot,run7,run8)

src = fullfile(projectRoot,'plots', ...
    'DCS_PAPER_REPRODUCTION_MASTER_FINAL_CORRECTED.m');

assert(isfile(src),'Missing existing Fig. 7/8 master: %s',src);

txt = fileread(src);

% The outer runner handles Fig. 6.
txt = strrep(txt,'RUN_FIG6 = true;','RUN_FIG6 = false;');
txt = strrep(txt,'RUN_FIG6=true;','RUN_FIG6=false;');

if run7
    txt = strrep(txt,'RUN_FIG7 = false;','RUN_FIG7 = true;');
    txt = strrep(txt,'RUN_FIG7=false;','RUN_FIG7=true;');
else
    txt = strrep(txt,'RUN_FIG7 = true;','RUN_FIG7 = false;');
    txt = strrep(txt,'RUN_FIG7=true;','RUN_FIG7=false;');
end

if run8
    txt = strrep(txt,'RUN_FIG8 = false;','RUN_FIG8 = true;');
    txt = strrep(txt,'RUN_FIG8=false;','RUN_FIG8=true;');
else
    txt = strrep(txt,'RUN_FIG8 = true;','RUN_FIG8 = false;');
    txt = strrep(txt,'RUN_FIG8=true;','RUN_FIG8=false;');
end

% Convert common Unicode punctuation/symbols to ASCII.
pairs = { ...
    char(8211), '-', ...
    char(8212), '-', ...
    char(8722), '-', ...
    char(8594), '->', ...
    char(8592), '<-', ...
    char(8804), '<=', ...
    char(8805), '>=', ...
    char(177),  '+/-', ...
    char(215),  'x', ...
    char(176),  ' deg', ...
    char(956),  'u', ...
    char(8230), '...', ...
    char(8216), '''', ...
    char(8217), '''', ...
    char(8220), '"', ...
    char(8221), '"', ...
    char(160),  ' '};

for k = 1:2:numel(pairs)
    txt = strrep(txt,pairs{k},pairs{k+1});
end

% Remove any remaining non-ASCII characters.
d = double(txt);
txt = txt(d <= 127);

tmp = fullfile(projectRoot,'plots','DCS_RUNTIME_FIG78_MASTER.m');

fid = fopen(tmp,'w');
assert(fid>=0,'Could not create temporary runtime Fig. 7/8 master.');
cleanupFile = onCleanup(@() safe_delete(tmp)); %#ok<NASGU>

fwrite(fid,txt,'char');
fclose(fid);

fprintf('Running ASCII-safe Fig. 7/8 runtime copy...\n');
pRootBackup = projectRoot;
run7Backup  = run7;
run8Backup  = run8;

oldDir = pwd;
cd(fullfile(pRootBackup,'plots'));
cleanupCd = onCleanup(@() cd(oldDir));
run(tmp);

if ~exist('projectRoot','var') || isempty(projectRoot), projectRoot = pRootBackup; end
if ~exist('run7','var'), run7 = run7Backup; end
if ~exist('run8','var'), run8 = run8Backup; end

if run7
    assert(isfile(fullfile(projectRoot,'evaluation', ...
        'DCS_MASTER_FIG7.png')), ...
        'Fig. 7 PNG was not generated.');
    assert(isfile(fullfile(projectRoot,'evaluation', ...
        'DCS_MASTER_FIG7_RESULTS.mat')), ...
        'Fig. 7 MAT results were not generated.');
    fprintf('Fig. 7: PNG + MAT results confirmed.\n');
end

if run8
    assert(isfile(fullfile(projectRoot,'evaluation', ...
        'DCS_MASTER_FIG8.png')), ...
        'Fig. 8 PNG was not generated.');
    assert(isfile(fullfile(projectRoot,'evaluation', ...
        'DCS_MASTER_FIG8_RESULTS.mat')), ...
        'Fig. 8 MAT results were not generated.');
    fprintf('Fig. 8: PNG + MAT results confirmed.\n');
end

end

function run_table2(projectRoot)

f = fullfile(projectRoot,'plots','table2_complexity.m');
assert(isfile(f),'Missing Table 2 script: %s',f);

oldDir = pwd;
cleanup = onCleanup(@() cd(oldDir));

cd(fullfile(projectRoot,'plots'));
run(f);

srcMat = fullfile(projectRoot,'plots','paper_table2_complexity.mat');
srcCsv = fullfile(projectRoot,'plots','paper_table2_complexity.csv');

if isfile(srcMat)
    copyfile(srcMat,fullfile(projectRoot,'evaluation', ...
        'paper_table2_complexity.mat'));
end

if isfile(srcCsv)
    copyfile(srcCsv,fullfile(projectRoot,'evaluation', ...
        'paper_table2_complexity.csv'));
end

fprintf('Table 2: MAT + CSV results confirmed.\n');

end


function figs = sort_figures(figs)

if isempty(figs)
    return;
end

nums = zeros(numel(figs),1);
for k=1:numel(figs)
    nums(k) = double(figs(k).Number);
end

[~,idx] = sort(nums);
figs = figs(idx);

end


function safe_delete(f)

if isfile(f)
    try
        delete(f);
    catch
    end
end

end
