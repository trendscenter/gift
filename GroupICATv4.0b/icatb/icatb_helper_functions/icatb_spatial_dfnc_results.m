%% Spatial dFNC Results
function icatb_spatial_dfnc_results(param_file, htmlDir)

if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('title', 'Select spatial dFNC Parameter File', 'typeEntity', 'file', 'typeSelection', 'single', 'filter', '*_sdfnc.mat');
    drawnow;
    if (isempty(param_file))
        error('ICA/sdFNC parameter file is not selected');
    end
end

drawnow;

[outputDir, fN, extn] = fileparts(param_file);
if (isempty(outputDir))
    outputDir = pwd;
end
param_file = fullfile(outputDir, [fN, extn]);

load (param_file);

if (~exist('sdfncInfo', 'var'))
    error('Selected file is not a valid spatial dfnc file');
end

if (~isfield(sdfncInfo, 'resultFiles'))
    error('Please run spatial dfnc in order to do post-processing');
end

sdfncInfo.outputDir = outputDir;
if (~exist('htmlDir', 'var'))
    htmlDir = fullfile(outputDir, [sdfncInfo.userInput.prefix, '_sdfnc_html']);
end

if (exist(htmlDir, 'dir') ~= 7)
    mkdir(htmlDir);
end

sdfncInfo.htmlDir = htmlDir;

% Plot t-maps
sdFNCResults = plotComps(sdfncInfo);
writeHTML(htmlDir, sdFNCResults, 'tmaps.html', 'T-maps');

% Mean FC
sdFNCResults = plotMeanFC(sdfncInfo);
writeHTML(htmlDir, sdFNCResults, 'mean_fc.html', 'Mean FC');

try
    % Two sample t-test
    sdFNCResults = plotTtestWin(sdfncInfo);
    writeHTML(htmlDir, sdFNCResults, 'ttest_fc_win.html', 'Two sample t-test Of Windowed FC');
catch
end

% Plot KLD
sdFNCResults = plotKLD(sdfncInfo);
writeHTML(htmlDir, sdFNCResults, 'kld.html', 'Kullback Lieber Divergence');

try
    % Median test
    sdFNCResults = plotMedianFC(sdfncInfo);
    sdFNCResults(end + 1) = plotMannWhitney(sdfncInfo);
    writeHTML(htmlDir, sdFNCResults, 'median_test.html', 'Median test');
catch
end

% Cluster states
sdFNCResults = plotStates(sdfncInfo);
writeHTML(htmlDir, sdFNCResults(1), 'avg_states.html', 'Average FC Of States');
writeHTML(htmlDir, sdFNCResults(2), 'hist_fc_states.html', 'Histogram Of State Size');
writeHTML(htmlDir, sdFNCResults(3), 'tm_states.html', 'Transition Matrix');


% Copy help files
copyfile(fullfile(fileparts(which('gift')), 'icatb_helpManual', '*.gif'), htmlDir);
copyfile(fullfile(fileparts(which('icatb_spatial_dfnc_results.m')), 'icatb_spatial_dfnc_results.html'), htmlDir);
copyfile(fullfile(fileparts(which('icatb_spatial_dfnc_results.m')), 'icatb_sdfnc_toc.html'), htmlDir);


% Open help file
icatb_openHTMLHelpFile(fullfile(htmlDir, 'icatb_spatial_dfnc_results.html'));

function sdFNCResults = plotComps(sdfncInfo)
%% One sample t-test maps
% One sample t-test is computed on IVA components across all data-sets and windows

outputDir = sdfncInfo.outputDir;
post_process_params = sdfncInfo.postprocess.params;
compFiles = icatb_rename_4d_file(fullfile(outputDir, sdfncInfo.compFiles));
compFiles = compFiles(post_process_params.comps, :);
df = (sdfncInfo.numOfSub*sdfncInfo.numOfSess*sdfncInfo.windows - 1);
threshold = post_process_params.threshold;
threshold_type = post_process_params.threshold_type;
structFile = fullfile (fileparts(which('gift.m')), 'icatb_templates', 'ch2bet.nii');
htmlDir = sdfncInfo.htmlDir;

outFiles = cell(1, size(compFiles, 1));
for nF = 1:size(compFiles, 1)
    
    gH = figure('color', 'w', 'name', ['T-map ',  num2str(post_process_params.comps(nF))], 'visible', 'on', 'tag', ['sdfnc_figures', num2str(nF)]);
    
    if (nF == 1)
        cmap = getCmap('positive');
    end
    movegui(gH, 'center');
    cn = deblank(compFiles(nF, :));
    dat = icatb_loadData(cn);
    good_inds = find(abs(dat) > eps);
    pvals = icatb_get_pvalue(dat(good_inds), df, 0);
    if (strcmpi(threshold_type, 'fdr'))
        [p_masked, p_fdr] = icatb_fdr(pvals, threshold);
        p_masked = pvals < p_masked;
        str = ['p < ', num2str(threshold), ' FDR corrected '];
    else
        p_masked = pvals < threshold;
        str = ['p < ', num2str(threshold), ' uncorrected'];
    end
    good_inds = good_inds(p_masked);
    
    if (isempty(good_inds))
        warning(['No voxels found for component with the specified threshold p <', num2str(threshold), ' ', str]);
    else
        sh = axes('parent', gH, 'position', [0.1, 0.1, 0.8, 0.8], 'units', 'normalized');
        minICAVal = min(abs(dat(good_inds)));
        maxICAVal = max(abs(dat(good_inds)));
        plotStackedOrtho(cn, 'structfile', structFile, 'image_values', 'positive', 'convert_to_zscores', 'no', 'threshold',  min(abs(dat(good_inds))), 'set_to_max_voxel', 1, ...
            'get_interp_data', 1, 'cmap', cmap, 'axesh', sh, 'colorbar_label', cellstr(char(num2str(minICAVal, '%0.1f'), num2str(maxICAVal, '%0.1f'))), ...
            'colorbar', true, 'labels', ['Component ', icatb_returnFileIndex(nF)]);
        outFile = [sdfncInfo.userInput.prefix, '_sdfnc_tmap_', num2str(post_process_params.comps(nF)), '.png'];
        myFrame = getframe(gH);
        imwrite(myFrame.cdata, fullfile(htmlDir, outFile));
        outFiles{nF} = outFile;
        close(gH);
    end
    
end

outFiles = outFiles(icatb_good_cells(outFiles));
sdFNCResults(1).title = 'T-maps';
sdFNCResults(1).text = ['One sample t-test is computed on IVA components across all data-sets and windows (', str, ').'];
sdFNCResults(1).files = outFiles;

function sdFNCResults = plotMeanFC(sdfncInfo)
%% Mean Of Standard deviation Of FC
% Standard deviation is computed across windows and median of standard
% deviation is computed for each group and averaged across subjects.


outputDir = sdfncInfo.outputDir;
groupNames = cellstr(char(sdfncInfo.userInput.group.name));
post_process_params = sdfncInfo.postprocess.params;

prefix = sdfncInfo.userInput.prefix;
post_process_file = [prefix, '_sdfnc_post_process.mat'];
post_process_file = fullfile(outputDir, post_process_file);

htmlDir = sdfncInfo.htmlDir;

load(post_process_file, 'MIStdVals');

sz = get(0, 'ScreenSize');
sz = 0.8*sz;
fig_position = [50, 50, sz(3), sz(4)];

num_rows = ceil(sqrt(length(groupNames)));
num_cols = ceil(length(groupNames)/num_rows);

countA = 0;
CLIM = [];
gH = figure('color', 'w', 'name', 'Mean Of STD Of FC', 'position', fig_position, 'visible', 'on', 'tag', 'sdfnc_figures_fcstd');
for nrow = 1:num_rows
    for ncol = 1:num_cols
        countA = countA + 1;
        if (countA <= length(groupNames))
            M = squeeze(mean(MIStdVals{countA}, 3));
            sh(countA) = subplot(num_rows, num_cols, countA);
            [FH, AH, CH, IH] = icatb_plot_FNC(M, [], cellstr(num2str(post_process_params.comps(:))), (1:length(post_process_params.comps)), gH, ' ', sh(countA));
            colormap(jet(64));
            title(['Avg std of ', groupNames{countA}], 'parent', sh(countA), 'fontweight', 'bold', 'fontsize', 13);
            maxCLIM = [CLIM, max((M(:)))];
            minCLIM = [CLIM, min((M(:)))];
            delete(CH);
        end
    end
end

set(sh, 'CLIM', [minCLIM, maxCLIM]);
ch = colorbar;
pos = get(sh(end), 'position');
pos(1) = pos(1) + pos(3) + 0.05;
pos(3) = 0.013;
sPos = get(sh(1), 'position');
ePos = get(sh(end), 'position');
pos(2) = 0.5*(sPos(2) - ePos(2));
set(ch, 'position', pos);


outFile = [sdfncInfo.userInput.prefix, '_mean_std_fc.png'];
myFrame = getframe(gH);
imwrite(myFrame.cdata, fullfile(htmlDir, outFile));
close(gH);

sdFNCResults(1).title = 'Mean Of Standard deviation Of FC';
sdFNCResults(1).text = 'Standard deviation is computed across windows and median of standard deviation is computed for each group and averaged across subjects.';
sdFNCResults(1).files = {outFile};

function sdFNCResults = plotTtestWin(sdfncInfo)
%% Two sample t-test results
% Two sample t-test is computed between groups at each window.

num_windows = sdfncInfo.windows;
outputDir = sdfncInfo.outputDir;
prefix = sdfncInfo.userInput.prefix;
post_process_file = [prefix, '_sdfnc_post_process.mat'];
post_process_file = fullfile(outputDir, post_process_file);
sz = get(0, 'ScreenSize');
sz = 0.8*sz;
fig_position = [50, 50, sz(3), sz(4)];
post_process_params = sdfncInfo.postprocess.params;
htmlDir = sdfncInfo.htmlDir;

%try
%resultsInfo = load(post_process_file, 'groupCombNames', 'ttest_results');
resultsInfo = load(post_process_file);
if (~isfield(resultsInfo, 'ttest_results'))
    return;
end

ttest_results = resultsInfo.ttest_results;
groupCombNames = resultsInfo.groupCombNames;

num_rows = ceil(sqrt(num_windows));
num_cols = ceil(num_windows/num_rows);
outFiles = cell(1, length(groupCombNames));
for nC = 1:length(groupCombNames)
    tmp = squeeze(ttest_results(nC, :, :, :, 1));
    tmp(isnan(tmp))=0;
    CLIM = max(abs(tmp(:)));
    gH = figure('color', 'w', 'name', ['Two sample t-test of ', groupCombNames{nC}], 'position', fig_position, 'visible', 'on', 'tag', ['sdfnc_figures_ttest', num2str(nC)]);
    clear sh;
    
    sh = zeros(1, num_windows);
    
    for nWin = 1:num_windows
        sh(nWin) = subplot(num_rows, num_cols, nWin);
        [FH, AH, CH, IH] = icatb_plot_FNC(squeeze(tmp(nWin, :, :)), [-CLIM, CLIM], cellstr(num2str(post_process_params.comps(:))), (1:length(post_process_params.comps)), gH, ...
            ' ',  sh(nWin));
        colormap(jet(64));
        if ((nWin == round(num_cols/2)))
            title(char(groupCombNames{nC}, ' '), 'parent', sh(nWin), 'fontweight', 'bold', 'fontsize', 13);
        end
        ylabel(['Window ', num2str(nWin)], 'parent', sh(nWin));
        delete(CH);
    end
    
    ch = colorbar;
    pos = get(sh(end), 'position');
    pos(1) = pos(1) + pos(3) + 0.12;
    pos(3) = 0.013;
    set(ch, 'position', pos);
    ylabel('T-values', 'parent', ch);
    
    outFile = [sdfncInfo.userInput.prefix, '_ttest_fc_win.png'];
    myFrame = getframe(gH);
    imwrite(myFrame.cdata, fullfile(htmlDir, outFile));
    outFiles{nC} = outFile;
    close(gH);
    
end

if (strcmpi(post_process_params.threshold_type, 'fdr'))
    str = ['p < ', num2str(post_process_params.threshold), ' FDR corrected '];
else
    str = ['p < ', num2str(post_process_params.threshold), ' uncorrected'];
end

sdFNCResults(1).title = 'Two sample t-test Of Windowed FC';
sdFNCResults(1).text = ['Two sample t-test is computed between groups at each window (', str, ')'];
sdFNCResults(1).files = outFiles;


function sdFNCResults = plotKLD(sdfncInfo)
%% KL divergence between windows
% Average Kullback-Leibler divergence of functional connectivity is
% computed between sequential windows for each group

num_windows = sdfncInfo.windows;
outputDir = sdfncInfo.outputDir;
prefix = sdfncInfo.userInput.prefix;
post_process_file = [prefix, '_sdfnc_post_process.mat'];
post_process_file = fullfile(outputDir, post_process_file);
load(post_process_file, 'KL');
groupNames = cellstr(char(sdfncInfo.userInput.group.name));
htmlDir = sdfncInfo.htmlDir;

gH = figure('color', 'w', 'name', 'Average KLD for FNC', 'visible', 'on', 'tag', 'sdfnc_figures_kld');
values = cell(1, length(groupNames));
errors = cell(1, length(groupNames));

for nG = 1:length(groupNames)
    tmp = KL{nG};
    values{nG} = mean(tmp);
    errors{nG} = std(tmp)/sqrt(size(tmp, 1));
end

bw_legend = cell(1, num_windows - 1);
for nWin = 1:length(bw_legend)
    bw_legend{nWin} = ['w', num2str(nWin), ' & w', num2str(nWin + 1)];
end
values = cat(1, values{:});
errors = cat(1, errors{:});
plotBars(values, errors, 0.8, groupNames, winter, bw_legend);


outFile = [sdfncInfo.userInput.prefix, '_kld.png'];
myFrame = getframe(gH);
imwrite(myFrame.cdata, fullfile(htmlDir, outFile));
close(gH);

sdFNCResults(1).title = 'Kullback Liebler Divergence';
sdFNCResults(1).text = 'Average Kullback-Leibler divergence of functional connectivity is computed between sequential windows for each group';
sdFNCResults(1).files = {outFile};


function sdFNCResults = plotMedianFC(sdfncInfo)
%% Median Of Standard Deviation Of FC
% Median is computed on standard deviation of functional connectivity. Standard deviation is
% computed across windows.

groupNames = cellstr(char(sdfncInfo.userInput.group.name));
outputDir = sdfncInfo.outputDir;
prefix = sdfncInfo.userInput.prefix;
post_process_file = [prefix, '_sdfnc_post_process.mat'];
post_process_file = fullfile(outputDir, post_process_file);
sz = get(0, 'ScreenSize');
sz = 0.8*sz;
fig_position = [50, 50, sz(3), sz(4)];
post_process_params = sdfncInfo.postprocess.params;
htmlDir = sdfncInfo.htmlDir;

%resultsInfo = load(post_process_file, 'MIStdVals', 'groupCombNames');
resultsInfo = load(post_process_file);
if (~isfield(resultsInfo, 'MIStdVals'))
    return;
end

MIStdVals = resultsInfo.MIStdVals;
groupCombNames = resultsInfo.groupCombNames;


gH = figure('color', 'w', 'name', 'Median of FC STD', 'position', fig_position, 'visible', 'on', 'tag', 'sdfnc_figures_mdian_fc');
num_rows = ceil(sqrt(length(groupNames)));
num_cols = ceil(length(groupNames)/num_rows);
clear sh;
sh = zeros(1, sdfncInfo.numGroups);
CLIM = [];
for nG = 1:sdfncInfo.numGroups
    tmp = squeeze(median(MIStdVals{nG}, 3));
    sh(nG) = subplot(num_rows, num_cols, nG);
    [FH, AH, CH, IH] = icatb_plot_FNC(tmp, [], cellstr(num2str(post_process_params.comps(:))), (1:length(post_process_params.comps)), gH, ' ', sh(nG));
    title(['Median of FC STD Of ', groupNames{nG}], 'parent', sh(nG), 'fontweight', 'bold', 'fontsize', 13);
    colormap(jet(64));
    CLIM = max([CLIM, max(abs(tmp(:)))]);
    delete(CH);
end

set(sh, 'CLIM', [0, CLIM]);
ch = colorbar;
pos = get(sh(end), 'position');
pos(1) = pos(1) + pos(3) + 0.05;
pos(3) = 0.013;
sPos = get(sh(1), 'position');
ePos = get(sh(end), 'position');
pos(2) = 0.5*(sPos(2) - ePos(2));
set(ch, 'position', pos);


outFile = [sdfncInfo.userInput.prefix, '_median_test.png'];
myFrame = getframe(gH);
imwrite(myFrame.cdata, fullfile(htmlDir, outFile));
close(gH);

sdFNCResults(1).title = 'Median Of Standard Deviation Of FC';
sdFNCResults(1).text = 'Median is computed on standard deviation of functional connectivity. Standard deviation is computed across windows.';
sdFNCResults(1).files = {outFile};

function sdFNCResults = plotMannWhitney(sdfncInfo)
%% Mann-whitney U-test
% Mann-whitney U-test is computed for each standard deviation of functional connectivity between
% groups. The thicker the line the more significant is the median value
% over time.

outputDir = sdfncInfo.outputDir;
groupNames = cellstr(char(sdfncInfo.userInput.group.name));
group_combinations = nchoosek(1:sdfncInfo.numGroups, 2);
prefix = sdfncInfo.userInput.prefix;
post_process_file = [prefix, '_sdfnc_post_process.mat'];
post_process_file = fullfile(outputDir, post_process_file);
post_process_params = sdfncInfo.postprocess.params;
htmlDir = sdfncInfo.htmlDir;

%load(post_process_file, 'median_test_results', 'groupCombNames');
resultsInfo = load(post_process_file);
if (~isfield(resultsInfo, 'median_test_results'))
    return;
end

median_test_results = resultsInfo.median_test_results;
groupCombNames = resultsInfo.groupCombNames;

numComp = length(post_process_params.comps(:));
num_rows = ceil(sqrt(length(groupCombNames)));
num_cols = ceil(length(groupCombNames)/num_rows);

%try
icnames = cellstr([repmat('IC ', numComp, 1), num2str(post_process_params.comps(:))]);
sh = zeros(1, length(groupCombNames));
gH = figure('color', 'w', 'name', 'FC STD Between Groups', 'visible', 'on', 'tag', 'sdfnc_figures_mann_whiteny');
for nComb = 1:length(groupCombNames)
    g1Ind = group_combinations(nComb, 1);
    g2Ind = group_combinations(nComb, 2);
    sh(nComb) = subplot(num_rows, num_cols, nComb);
    a1 = squeeze(median_test_results(nComb, :, :, 1));
    a1(isnan(a1)) = 0;
    a1(abs(a1) ~= 0) = 1;
    a2 = squeeze(median_test_results(nComb, :, :, 2));
    a2(isnan(a2)) = 0;
    plot_AdjMatrix(a1, icnames, a2, groupNames([g1Ind, g2Ind]));
    title(groupCombNames{nComb}, 'parent', sh(nComb), 'fontweight', 'bold', 'fontsize', 13);
end

outFile = [sdfncInfo.userInput.prefix, '_mann_whitney.png'];
myFrame = getframe(gH);
imwrite(myFrame.cdata, fullfile(htmlDir, outFile));
close(gH);

sdFNCResults(1).title = 'Mann-whitney U-test';
sdFNCResults(1).text = 'Mann-whitney U-test is computed for each standard deviation of functional connectivity between groups. The thicker the line the more significant is the median value over time.';
sdFNCResults(1).files = {outFile};

% catch
% end


function sdFNCResults = plotStates(sdfncInfo)


outputDir = sdfncInfo.outputDir;
groupNames = cellstr(char(sdfncInfo.userInput.group.name));
prefix = sdfncInfo.userInput.prefix;
post_process_file = [prefix, '_sdfnc_post_process.mat'];
post_process_file = fullfile(outputDir, post_process_file);
htmlDir = sdfncInfo.htmlDir;

load(post_process_file, 'clusterInfo');

colors = {'-ob', '-or', '-og', '-ok', '-om', '-oy', '-ok'};

sz = get(0, 'ScreenSize');
sz = 0.7*sz;
fig_position = [50, 50, sz(3), sz(4)];
post_process_params = sdfncInfo.postprocess.params;

outFiles = cell(length(clusterInfo), 3);

for nComp = 1:length(clusterInfo)
    
    gH = figure('color', 'w', 'name', ['Component ', num2str(sdfncInfo.postprocess.params.comps(nComp))], 'visible', 'on', 'tag', ['sdfnc_figures_cluster', num2str(nComp)], ...
        'position', fig_position);
    
    num_rows = ceil(sqrt(sdfncInfo.postprocess.params.num_clusters));
    num_cols = ceil(sdfncInfo.postprocess.params.num_clusters/num_rows);
    
    
    %% Markov Modeling results
    % K-means is computed for each component functional connectivity values with all other components. Average functional connectivity is computed for within state (with standard error
    % mean)
    for nCluster = 1:sdfncInfo.postprocess.params.num_clusters
        sh = subplot(num_rows, num_cols, nCluster);
        errorH = zeros(1, length(groupNames));
        for nG = 1:length(groupNames)
            if (~isempty(clusterInfo{nComp}.state_mean{nG, nCluster}))
                errorH(nG) = errorbar(clusterInfo{nComp}.state_mean{nG, nCluster}, clusterInfo{nComp}.state_sem{nG, nCluster}, colors{nG});
            else
                continue;
            end
            hold on;
        end
        hold on;
        title(['State ', num2str(nCluster)], 'parent', sh);
        if (nCluster == 1)
            ylabel(['Component ', num2str(sdfncInfo.postprocess.params.comps(nComp))], 'parent', sh);
        end
        axis tight;
        if (nCluster == sdfncInfo.postprocess.params.num_clusters)
            legendStr = groupNames;
            legendStr = legendStr(errorH ~= 0);
            errorH = errorH(errorH ~= 0);
            if (~isempty(legendStr))
                legend(errorH, legendStr{:});
            end
            clear legendStr;
        end
    end
    
    outFile = [sdfncInfo.userInput.prefix, '_avg_states', num2str(sdfncInfo.postprocess.params.comps(nComp)), '.png'];
    myFrame = getframe(gH);
    imwrite(myFrame.cdata, fullfile(htmlDir, outFile));
    close(gH);
    
    outFiles{nComp, 1} = outFile;
    
    clear outFile;
    
    %% Histogram for state size
    % Counts for each state and group are plotted. Maximum count is windows
    % x subjects.
    legendStr = [repmat('State ', sdfncInfo.postprocess.params.num_clusters, 1), num2str((1:sdfncInfo.postprocess.params.num_clusters)')];
    gH = figure('color', 'w', 'name', ['Frequency for component ', num2str(sdfncInfo.postprocess.params.comps(nComp))], 'visible', 'on', 'tag', ['sdfnc_figures_hist', num2str(nComp)]);
    bh = bar(clusterInfo{nComp}.frequency_states);
    set(gca, 'XTickLabel', groupNames);
    colormap(winter);
    ylabel('Count (Max = windows x subjects)');
    legend(legendStr, 'location', 'BestOutSide');
    
    outFile = [sdfncInfo.userInput.prefix, '_hist_states', num2str(sdfncInfo.postprocess.params.comps(nComp)), '.png'];
    myFrame = getframe(gH);
    imwrite(myFrame.cdata, fullfile(htmlDir, outFile));
    close(gH);
    outFiles{nComp, 2} = outFile;
    
    clear outFile;
    
    
    %% Transition matrix
    % Transition matrix is computed for each group
    gH = figure('color', 'w', 'name', ['Transition Matrix for component ', num2str(sdfncInfo.postprocess.params.comps(nComp))], 'visible', 'on', 'tag', ['sdfnc_figures_trans', num2str(nComp)]);
    num_rows = ceil(sqrt(length(groupNames)));
    num_cols = ceil(length(groupNames)/num_rows);
    strs = cellstr(num2str((1:sdfncInfo.postprocess.params.num_clusters)'));
    sh = zeros(1, length(groupNames));
    for nG = 1:length(groupNames)
        sh(nG) = subplot(num_rows, num_cols, nG);
        icatb_plot_matrix(clusterInfo{nComp}.transition_matrix{nG}, strs, strs, 'axesh', sh(nG), 'title', groupNames{nG}, 'cmap', hot, 'clim', [0, 1], 'colorbar', 0);
        set(findobj(gH, 'type', 'axes'), 'Ycolor', 'k');
        set(findobj(gH, 'type', 'axes'), 'Xcolor', 'k');
        axis image;
        ylabel(['Component ', num2str(sdfncInfo.postprocess.params.comps(nComp))], 'parent', sh(nG));
    end
    
    ch = colorbar;
    if (length(groupNames) > 1)
        pos = get(sh(end), 'position');
        pos(1) = pos(1) + pos(3) + 0.12;
        pos(3) = 0.02;
        sPos = get(sh(1), 'position');
        ePos = get(sh(end), 'position');
        pos(2) = 0.5*(sPos(2) - ePos(2));
        set(ch, 'position', pos);
    end
    
    outFile = [sdfncInfo.userInput.prefix, '_tm_states', num2str(sdfncInfo.postprocess.params.comps(nComp)), '.png'];
    myFrame = getframe(gH);
    imwrite(myFrame.cdata, fullfile(htmlDir, outFile));
    close(gH);
    outFiles{nComp, 3} = outFile;
    
    clear outFile;
    
    
end



sdFNCResults(1).title = 'Markov Modeling results';
sdFNCResults(1).text = 'K-means is computed for each component functional connectivity values with all other components. Average functional connectivity is computed for within state (with standard error mean)';
outFile = outFiles(:, 1);
outFile = outFile(icatb_good_cells(outFile));
sdFNCResults(1).files = outFile;

clear outFile;

sdFNCResults(2).title = 'Histogram for state size';
sdFNCResults(2).text = 'Counts for each state and group are plotted. Maximum count is windows X subjects.';
outFile = outFiles(:, 2);
outFile = outFile(icatb_good_cells(outFile));
sdFNCResults(2).files = outFile;

clear outFile;


sdFNCResults(3).title = 'Transition matrix';
sdFNCResults(3).text = 'Transition matrix is computed for each group';
outFile = outFiles(:, 3);
outFile = outFile(icatb_good_cells(outFile));
sdFNCResults(3).files = outFile;

clear outFile;


function plot_AdjMatrix(A,text1, a2, gNames)

[m,n]=size(A);
if m~=n
    error('Adjacency matrix must be square.')
end

sign1 = sign(a2);

minLineWidth = min(abs(a2(a2 ~= 0)));


vertices=exp(j*2*pi/n.*(0:n-1));
vertices=vertices(:);
x=real(vertices);
y=imag(vertices);

scatter(x,y,'filled'); hold on;
if ~isempty(text1)
    for i = 1:length(text1)
        text(x(i)+0.02,y(i),text1{i},'FontSize',14)
    end
end

[row,col] = find(A~=0);

posHandles = [];
negHandles = [];
for i = 1:length(row)
    x_line = [x(row(i)),x(col(i))];
    y_line = [y(row(i)),y(col(i))];
    if ~isempty(sign1)
        if sign1(row(i),col(i)) == 1
            lineWidth = abs(a2(row(i), col(i)))/minLineWidth;
            posHandles(end + 1) = plot(x_line, y_line, '-r', 'linewidth', lineWidth); hold on;
        elseif sign1(row(i),col(i)) == -1
            lineWidth = abs(a2(row(i), col(i)))/minLineWidth;
            negHandles(end + 1) = plot(x_line, y_line, '-b', 'linewidth', lineWidth); hold on;
        elseif sign1(row(i),col(i)) == 0
            plot(x_line, y_line, '-k'); hold on;
        end
    else
        plot(x_line, y_line, '-'); hold on;
    end
    clear x_line y_line;
end

handles = [];
legendStr = {};
try
    handles = [handles, posHandles(end)];
    legendStr{end + 1} = [gNames{1}, ' - ', gNames{2}];
catch
end

try
    handles = [handles, negHandles(end)];
    legendStr{end + 1} = [gNames{2}, ' - ', gNames{1}];
catch
end

try
    lh = legend(handles, legendStr);
    set(lh, 'location', 'BestOutSide');
catch
end

axis([-1.3 1.3 -1.3 1.3]);
axis square;
grid off;
axis off;


function plotBars(barvalues, errors, barWidth, groupNames, cmap, bw_legend)

if ((icatb_get_matlab_version <= 2013) || strcmpi(version('-release'), '2014a'))
    icatb_barweb(barvalues, errors, barWidth, groupNames, [], [], [], cmap, 'y', bw_legend, 'plot');
else
    inds = round(linspace(1, 64, size(barvalues, 2)));
    
    ticks = [];
    
    count = 0;
    for j = 1:size(barvalues, 1)
        for i = 1:size(barvalues, 2)
            count = count + 1;
            bh(count)=bar(i + (j-1)*(size(barvalues, 2)) + 2*(j-1), barvalues(j,i), barWidth, 'facecolor', cmap(inds(i), :));
            if (i==ceil(size(barvalues, 2)/2))
                ticks = [ticks, i + (j-1)*(size(barvalues, 2)) + 2*(j-1)];
            end
            set(gca, 'xtick', '');
            set(gca, 'xticklabel', '');
            hold on;
            errorbar(i + (j-1)*(size(barvalues, 2)) + 2*(j-1), barvalues(j,i), errors(j,i), 'kx', 'linewidth', 1);
        end
    end
    
    set(gca, 'xtick', ticks);
    set(gca, 'xticklabel', groupNames);
    legend(bh(1: size(barvalues, 2)), bw_legend{:}, 'location', 'best');
end



function [sliceXY, sliceXZ, sliceYZ] = returnSlices(data, voxelcoords)

sliceXY = rot90(reshape(data(:, :, voxelcoords(end)), size(data, 1), size(data, 2)));
sliceXZ = rot90(reshape(data(:, voxelcoords(2), :), size(data, 1),size(data, 3)));
sliceYZ = rot90(reshape(data(voxelcoords(1), :, :), size(data, 2), size(data, 3)));

function data = stackData(slices, minVal)

[m1, n1] = size(slices{1});
[m2, n2] = size(slices{2});
[m3, n3] = size(slices{3});

maxSizeX = max([m1, m2, m3]);
maxSizeY = max([n1, n2, n3]);

data = minVal*ones(maxSizeX, [n1 + n2 + n3]);

e = 0;
for nS = 1:length(slices)
    tmp = slices{nS};
    s = e + 1;
    e = e + size(tmp, 2);
    inda = ceil(maxSizeX/2) - ceil(size(tmp, 1)/2) + 1;
    indb = inda + size(tmp, 1) - 1;
    data(inda:indb, s:e) = tmp;
end


function c = charC(a, b)

len = max([length(a), length(b)]);

c = repmat(' ', 2, len);

s = ceil(len/2 - length(a)/2) + 1;
c(1, s:s+length(a) - 1) = a;
s = ceil(len/2 - length(b)/2) + 1;
c(2, s:s + length(b) - 1) = b;


function cmap = getCmap(image_values)

load icatb_colors coldhot_sensitive;

returnValue = strmatch(lower(image_values), {'positive and negative', 'positive', 'absolute value', 'negative'}, 'exact');

if (returnValue == 1)
    cmap = coldhot_sensitive(1:4:end, :);
elseif (returnValue == 4)
    cmap = coldhot_sensitive(1:128, :);
    cmap = cmap(1:2:end, :);
else
    cmap = coldhot_sensitive(129:end, :);
    cmap = cmap(1:2:end, :);
end

function deleteFigures

try
    chkH = findobj('-regexp', 'tag', 'sdfnc_figures');
    delete(chkH);
catch
end


function writeHTML(htmlDir, sdFNCResults, outFile, titleStr)


html_file = fullfile(htmlDir, outFile);
start_string = ['<html><head><title>', titleStr, '</title></head>'];

results_string1 = [];
for nR = 1:length(sdFNCResults)
    results_string1 = [results_string1, '<h2 align = "center">', sdFNCResults(nR).title, '</h2>'];
    results_string1 = [results_string1, '<hr>'];
    results_string1 = [results_string1, '<p align = "center">', sdFNCResults(nR).text, '</p>'];
    files =  sdFNCResults(nR).files;
    num_cols = 2;
    num_rows = ceil(length(files)/num_cols);
    countF = 0;
    results_string1 = [results_string1, '<table>'];
    for nrows = 1:num_rows
        results_string1 = [results_string1, '</tr>'];
        for ncols = 1:num_cols
            countF = countF + 1;
            if (countF <= length(files))
                results_string1 = [results_string1,  '<td align="center"><img align="center" src = "', files{countF}, '" > </img></td>'];
            end
        end
        results_string1 = [results_string1, '</tr>'];
    end
    results_string1 = [results_string1, '</table>'];
    results_string1 = [results_string1, '<p>&nbsp;&nbsp;</p>'];
end

end_string =  '</html>';


results_string = [start_string,  results_string1, end_string];

dlmwrite(html_file, results_string, '');


function plotStackedOrtho(file_name, varargin)

useColorbar = 1;
threshold = 1;
convert_to_zscores = 'no';
image_values = 'positive';
labels = '';
cmap = hot(64);
structFile = fullfile (fileparts(which('gift.m')), 'icatb_templates', 'ch2bet.nii');

for n = 1:2:length(varargin)
    if (strcmpi(varargin{n}, 'structfile'))
        structFile = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'image_values'))
        image_values = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'cmap'))
        cmap = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'threshold'))
        threshold = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'convert_to_zscores'))
        convert_to_zscores = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'axesh'))
        sh = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'colorbar'))
        useColorbar = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'labels'))
        labels = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'colorbar_label'))
        xTicklabel = varargin{n + 1};
    end
end

hD = icatb_orth_views(file_name, 'structfile', structFile, 'image_values', image_values, 'convert_to_zscores', convert_to_zscores, 'threshold', threshold, 'set_to_max_voxel', 1, ...
    'get_interp_data', 1, 'cmap', cmap);
if (~exist('sh', 'var'))
    sh = gca;
end
[sliceXY, sliceXZ, sliceYZ] = returnSlices(hD.data, hD.maxVoxelPos);
tmp = stackData({sliceYZ, sliceXZ, sliceXY}, 0);
imagesc(tmp, [1, 200]);
colormap(hD.cmap);
axis(sh, 'image');
set(sh, 'Ytick', []);
set(sh, 'Xtick', [])
realCoords = (hD.maxVoxelPos - hD.voxelOrigin).*hD.VOX;

str = charC (labels, ['(', num2str(realCoords(1)), ',', num2str(realCoords(2)), ',', num2str(realCoords(3)), ')']);
title(str, 'parent', sh, 'horizontalalignment', 'center', 'fontweight', 'bold');
if (useColorbar)
    ch = colorbar('location', 'southoutside');
    set(ch, 'units', 'normalized');
    chPos = get(sh, 'position');
    chPos(1) = chPos(1) + 0.5*chPos(3) - 0.2;
    chPos(2) = chPos(2) - 0.05;
    chPos(3) = 0.4;
    chPos(4) = 0.05;
    set(ch, 'position', chPos);
    set(ch, 'xLim', [1, 100]);
    xTicks = get(ch, 'xTick');
    set(ch, 'yTick', []);
    if (exist('xTicklabel', 'var'))
        set(ch, 'xTick', [xTicks(1), xTicks(end)]);
        set(ch, 'xTicklabel', xTicklabel);
    end
end