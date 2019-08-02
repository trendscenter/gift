function icatb_gica_html_report(param_file, opts)

global GICA_PARAM_FILE;

icatb_defaults;
global EXPERIMENTAL_TR;
global PARAMETER_INFO_MAT_FILE;
global CONVERT_Z;
global IMAGE_VALUES;
global THRESHOLD_VALUE;
global ANATOMICAL_PLANE;
global TIMECOURSE_POSTPROCESS;
global FONT_COLOR;

if (~isempty(GICA_PARAM_FILE))
    param_file = GICA_PARAM_FILE;
end

if (~exist('param_file', 'var'))
    filterP = ['*', PARAMETER_INFO_MAT_FILE, '*.mat'];
    param_file = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Parameter File', 'filter', filterP);
end

outputDir = fileparts(param_file);
if (isempty(outputDir))
    outputDir = pwd;
end

load(param_file);


sesInfo.outputDir = outputDir;

structFile = fullfile(fileparts(which('groupica.m')), 'icatb_templates', 'ch2bet.nii');
try
    structFile = opts.anatomical_file;
catch
end

slices_in_mm = (-40:4:72);
try
    slices_in_mm = opts.slices_in_mm;
catch
end

convert_to_zscores = CONVERT_Z;
try
    convert_to_zscores = opts.convert_to_zscores;
catch
end


threshold = str2num(THRESHOLD_VALUE);
try
    threshold = opts.threshold;
catch
end

image_values = IMAGE_VALUES;
try
    image_values = opts.image_values;
catch
end

slice_plane = ANATOMICAL_PLANE;
try
    slice_plane = opts.anatomical_plane;
catch
end


groupsInfo = [];

try
    groupsInfo = opts.groupsInfo;
catch
end

%% Group ICA Parameters
sesInfo.structFile = structFile;
sesInfo.slice_plane = slice_plane;
sesInfo.image_values = image_values;
sesInfo.threshold = threshold;
sesInfo.groupsInfo = groupsInfo;
sesInfo.convert_to_zscores = convert_to_zscores;
newText = dispParams(sesInfo);
disp(newText);

drawnow;

%%

if (sesInfo.which_analysis == 2)
    
    %% ICASSO Plots
    icassoResultsFile = fullfile(outputDir, [sesInfo.userInput.prefix, '_icasso_results.mat']);
    try
        if (sesInfo.write_analysis_steps_in_dirs)
            icassoResultsFile = fullfile(outputDir, [sesInfo.userInput.prefix, '_ica_files', filesep, 'icasso_results.mat']);
        end
    catch
    end
    if (exist(icassoResultsFile, 'file'))
        load(icassoResultsFile, 'sR');
        try
            icassoShow(sR, 'L', sesInfo.numComp, 'colorlimit', [.8 .9]);
        catch
        end
    end
    
end

%% Mean Components
% Mean across all subjects and sessions is computed for each component
%
% * *a) Timecourse* - Mean timecourse is converted to z-scores.
% * *b) Spectra* - Timecourses spectra is computed for each data-set and averaged across sessions. Mean and standard error of mean is shown in the figure.
% * *c) Montage* - Axial slices are shown.
% * *d) Ortho slices* - Ortho plot is shown for the peak voxel and coordinates are reported.

compFileNaming = sesInfo.icaOutputFiles(1).ses(1).name;
currentFile = deblank(compFileNaming(1, :));

if ~exist(currentFile, 'file')
    [zipFileName, files_in_zip] = icatb_getViewingSet_zip(currentFile, [], 'real', sesInfo.zipContents);
    if (~isempty(zipFileName))
        icatb_unzip(regexprep(zipFileName, ['.*\', filesep], ''), fullfile(outputDir, fileparts(currentFile)));
    end
end

compFiles = icatb_rename_4d_file(icatb_fullFile('directory', outputDir, 'files', compFileNaming));
icaTimecourse = icatb_loadICATimeCourse(compFiles);

postProcessFile = fullfile(outputDir, [sesInfo.userInput.prefix, '_postprocess_results.mat']);
if (~exist(postProcessFile, 'file'))
    icatb_postprocess_timecourses(param_file);
end

load(postProcessFile);

if (sesInfo.numOfSess > 1)
    spectra_tc_all = reshape(mean(spectra_tc_all, 2), sesInfo.numOfSub, length(freq), sesInfo.numComp);
else
    spectra_tc_all = reshape(squeeze(spectra_tc_all), sesInfo.numOfSub, length(freq), sesInfo.numComp);
end


clear tc;

freq_limits = [0.1, 0.15];
try
    freq_limits = TIMECOURSE_POSTPROCESS.spectra.freq_limits;
catch
end



if (~exist('spatial_maps_MI', 'var'))
    countS = 0;
    fprintf('\n');
    for nSub = 1:sesInfo.numOfSub
        for nSess = 1:sesInfo.numOfSess
            countS = countS + 1;
            ic = icatb_loadComp(sesInfo, (1:sesInfo.numComp), 'subjects', nSub, 'sessions', nSess, 'vars_to_load', 'ic');
            tmp = icatb_compute_mi(ic');
            if (countS == 1)
                spatial_maps_MI = zeros(sesInfo.numOfSub, sesInfo.numOfSess, sesInfo.numComp, sesInfo.numComp);
            end
            spatial_maps_MI(nSub, nSess, :, :) = tmp;
        end
    end
end

cmap = getCmap(image_values);

fALFF = zeros(size(compFiles, 1), 1);
dynamic_range = zeros(size(compFiles, 1), 1);
components = (1:size(compFiles, 1))';

for nF = 1:size(compFiles, 1)
    
    cn = deblank(compFiles(nF, :));
    
    [dd, fN, extn] = fileparts(cn);
    
    gH = icatb_getGraphics([fN, extn], 'graphics', 'imviewer', 'on');
    set(gH, 'resize', 'on');
    
    xOffSet = 0.05;
    yOffSet = 0.05;
    
    
    width2 = 0.5;
    height2 = 0.5;
    sh = axes('parent', gH, 'units', 'normalized', 'position', [0.01, yOffSet, width2, height2]);
    icatb_image_viewer(cn, 'structfile', structFile, 'image_values', image_values, 'convert_to_zscores', convert_to_zscores, 'threshold', threshold, ...
        'slices_in_mm', slices_in_mm, 'anatomical_view', slice_plane, 'axesh', sh, 'colorbar', 0, 'labels', ' ');
    
    
    width = 0.4;
    height = 0.4;
    
    sh = axes('parent', gH, 'units', 'normalized', 'position', [width2 + 0.03, height2/2-(height/2), width, height]);
    plotStackedOrtho(cn, 'structfile', structFile, 'image_values', image_values, 'convert_to_zscores', convert_to_zscores, 'threshold', threshold, 'set_to_max_voxel', 1, ...
        'get_interp_data', 1, 'cmap', cmap, 'axesh', sh, 'colorbar', false, 'labels', 'Peak Coordinates (mm)', 'colorbar', true, 'colorbar_label', true);
    
    yPos = 0.65;
    width = 0.4;
    height = 0.25;
    sh = axes('parent', gH, 'units', 'normalized', 'position', [xOffSet, yPos, width, height]);
    tc =  icatb_detrend(icaTimecourse(:, nF));
    tc_label = 'Z-scores';
    if (strcmpi(convert_to_zscores, 'yes'))
        tc = tc/std(tc);
        tc_label = 'Data units';
    end
    icatb_plotTimecourse('data', tc, 'parent', sh, 'color', 'm', 'title', ['Component ', icatb_returnFileIndex(nF)], 'xlabelstr', 'Scans', 'ylabelstr', tc_label);
    
    sh = axes('parent', gH, 'units', 'normalized', 'position', [xOffSet+width+0.1, yPos, width, height]);
    clear tc;
    tc.data = squeeze(spectra_tc_all(:, :, nF));
    tmp_dynamicrange = zeros(1, size(tc.data, 1));
    tmp_fALFF = tmp_dynamicrange;
    for nS = 1:length(tmp_dynamicrange)
        [tmp_dynamicrange(nS), tmp_fALFF(nS)] = icatb_get_spec_stats(tc.data(nS, :), freq, freq_limits);
    end
    
    fALFF(nF) = mean(tmp_fALFF);
    dynamic_range(nF) = mean(tmp_dynamicrange);
    
    if (nF == 1)
        dynamic_range_all = zeros(length(tmp_dynamicrange), size(compFiles, 1));
        fALFF_all = dynamic_range_all;
    end
    
    dynamic_range_all(:, nF) = tmp_dynamicrange;
    fALFF_all(:, nF) = tmp_fALFF;
    
    tc.xAxis = freq;
    tc.isSpectra = 1;
    tc.xlabelStr = 'Frequency (Hz)';
    tc.ylabelStr = 'Power';
    tc.titleStr = sprintf('Dynamic range: %0.3f, Power_L_F/Power_H_F: %0.3f',  dynamic_range(nF), fALFF(nF));
    icatb_plot_spectra(sh, tc);
    
    clear tc;
    
end

% Append fALFF and dynamic range info to post process file
save(postProcessFile, 'dynamic_range', 'fALFF', 'dynamic_range_all', 'fALFF_all', '-append');

%% Spectral Summary
% * *a) dynamic_range* - Difference between the peak power and minimum power at frequencies to the right of the peak.
% * *b) fALFF* - Low frequency to high frequency power ratio.
%
varNames = {'ComponentNumber', 'DynamicRange', 'fALFF'};
try
    T = table(components, dynamic_range, fALFF, 'VariableNames', varNames);
    disp(T);
catch
    fprintf('%20s\t%20s\t%20s\n', varNames{:});
    fprintf('\n');
    for n = 1:length(components)
        fprintf('%20s\t%20s\t%20s\n', num2str(components(n), '%d'), num2str(dynamic_range(n), '%0.3f'), num2str(fALFF(n), '%0.3f'));
    end
    fprintf('\n');
end


%% Temporal Stats On Beta Weights
% Multiple regression is done using the timecourses from SPM design matrix as model and ICA timecourses as observations. R^2 values for each
% component are shown in bar plot. For each component, one sample t-test results of each session and condition are shown in the bar plots.
%
temporal_stats_betas = [];
try
    temporal_stats_betas = opts.temporal_stats_betas;
catch
end

if (~isempty(temporal_stats_betas))
    
    varNames = {'ComponentNumber', 'Rsquare'};
    
    try
        T = table(temporal_stats_betas.component_numbers(:), temporal_stats_betas.R2(:), 'VariableNames', varNames);
        disp(T);
    catch
        fprintf('%20s\t%20s\n', varNames{:});
        fprintf('\n');
        for n = 1:length(temporal_stats_betas.component_numbers)
            fprintf('%20s\t%20s\n', num2str(temporal_stats_betas.component_numbers(n), '%d'), num2str(temporal_stats_betas.R2(n), '%0.3f'));
        end
        fprintf('\n');
    end
    
    icatb_disp_temp_regress_results(temporal_stats_betas, 0);
    
    %     fH = icatb_getGraphics('Regress values', 'timecourse', 'regress_betas', 'on');
    %     set(fH, 'resize', 'on');
    %     %fH = figure('color', 'w');
    %     sh = axes('units', 'normalized', 'position', [0.1, 0.3, 0.85, 0.6]);
    %     bh = bar(temporal_stats_betas.R2(:));
    %     colormap(hsv(1));
    %     set(sh,'XTick',(1:sesInfo.numComp));
    %     set(sh, 'XTickLabel', num2str(temporal_stats_betas.component_numbers(:)));
    %     ylabel('R^2', 'parent', sh, 'interpreter', 'tex');
    %     xlabel('Components', 'parent', sh);
    %     title('R^2', 'parent', sh, 'interpreter', 'tex');
    %     axis(sh, 'tight');
    %     set(sh, 'Xcolor', FONT_COLOR);
    %     set(sh, 'Ycolor', FONT_COLOR);
    %
    %     selectedRegressors = temporal_stats_betas.selectedRegressors;
    %     numCond = length(selectedRegressors);
    %     numOfSess = sesInfo.numOfSess;
    %     ttest_betas = temporal_stats_betas.ttest.tstat;
    %     pval_betas = temporal_stats_betas.ttest.pval;
    %     comps_betas = temporal_stats_betas.component_numbers;
    %     pval_sig1 = 0.01;
    %     pval_sig2 = 0.001;
    %     maxDF = max(squeeze(temporal_stats_betas.ttest.df(:)));
    %
    %     tval_sig1 = abs(icatb_spm_invTcdf(pval_sig1/2, maxDF));
    %     tval_sig2 = abs(icatb_spm_invTcdf(pval_sig2/2, maxDF));
    %
    %     YLIMA = [min([ttest_betas(:);-tval_sig2-0.1]), max([ttest_betas(:);tval_sig2 + 0.1])];
    %
    %
    %     if ((numOfSess == 1) || (numCond == 1))
    %         ttest_betas = squeeze(ttest_betas);
    %         pval_betas = squeeze(pval_betas);
    %         ncolsT = ceil(sqrt(size(ttest_betas, 1)));
    %         nrowsT = ceil(size(ttest_betas, 1) / ncolsT);
    %         for nCond = 1:size(ttest_betas, 2)
    %             tmpTT = ttest_betas(:, nCond);
    %             tmpTP = pval_betas(:, nCond);
    %             %tmpTP(isfinite(tmpTP)==0) = Inf;
    %             %tmpTT(tmpTP > pval_sig) = 0;
    %             fH = icatb_getGraphics('Tvalues', 'timecourse', 'regress_betas_T', 'on');
    %             set(fH, 'resize', 'on');
    %             sh = axes('units', 'normalized', 'position', [0.1, 0.16, 0.85, 0.65]);
    %             bh = bar(tmpTT);
    %             set(sh, 'Ylim', YLIMA);
    %             %set(sh,'color', 'w');
    %             colormap(hsv(1));
    %             ylabel('T-values', 'parent', sh);
    %             xlabel('Components', 'parent', sh);
    %             title(selectedRegressors{nCond}, 'parent', sh);
    %             set(sh,'XTick',(1:sesInfo.numComp));
    %             set(sh, 'XTickLabel', num2str(temporal_stats_betas.component_numbers(:)));
    %
    %             hold on;
    %             xLIMA = [1, sesInfo.numComp];
    %             set(sh, 'XLim', xLIMA);
    %             plot(xLIMA, tval_sig1*ones(1, 2), '--y', 'linewidth', 1.5);
    %             hold on;
    %             plot(xLIMA, tval_sig2*ones(1, 2), '--g', 'linewidth', 1.5);
    %             hold on;
    %             plot(xLIMA, -tval_sig1*ones(1, 2), '--y', 'linewidth', 1.5);
    %             hold on;
    %             plot(xLIMA, -tval_sig2*ones(1, 2), '--g', 'linewidth', 1.5);
    %             hold off;
    %             axis(sh, 'tight');
    %
    %             set(sh, 'Xcolor', FONT_COLOR);
    %             set(sh, 'Ycolor', FONT_COLOR);
    %             legendStr = {'Comp'; ['p < ', num2str(pval_sig1, '%0.2f')]; ['p < ', num2str(pval_sig2, '%0.3f')]};
    %             legend(legendStr{:}, 'location', 'bestoutside');
    %
    %         end
    %
    %     else
    %
    %         legendStr = cellstr([repmat('Run ', sesInfo.numOfSess, 1),  num2str((1:sesInfo.numOfSess)')]);
    %         legendStr = [legendStr; ['p < ', num2str(pval_sig1, '%0.2f')]; ['p < ', num2str(pval_sig2, '%0.3f')]];
    %         for nComp = 1:sesInfo.numComp
    %             tmpTP = squeeze(pval_betas(:, nComp, :));
    %             tmpTT = squeeze(ttest_betas(:, nComp, :));
    %             %tmpTP(isfinite(tmpTP)==0) = Inf;
    %             %tmpTT(tmpTP > pval_sig) = 0;
    %             %fH = figure('color', 'w');
    %             fH = icatb_getGraphics('T values', 'timecourse', 'regress_betas_T', 'on');
    %             set(fH, 'resize', 'on');
    %             sh = axes('units', 'normalized', 'position', [0.1, 0.16, 0.85, 0.65]);
    %             bh = bar(tmpTT');
    %             set(sh, 'Ylim', YLIMA);
    %             %set(sh,'color', 'w');
    %             colormap(hsv(sesInfo.numOfSess));
    %             ylabel('T-values', 'parent', sh);
    %             %xlabel('Regressors', 'parent', sh);
    %             title(['Component ', num2str(temporal_stats_betas.component_numbers(nComp))], 'parent', sh);
    %             set(sh,'XTick',(1:length(selectedRegressors)));
    %             set(sh, 'XTickLabel', selectedRegressors, 'TickDir', 'out');
    %
    %             hold on;
    %             %xLIMA = get(sh, 'XLim');
    %             xLIMA = [1, length(selectedRegressors)];
    %             set(sh, 'XLim', xLIMA);
    %             plot(xLIMA, tval_sig1*ones(1, 2), '--y', 'linewidth', 1.5);
    %             hold on;
    %             plot(xLIMA, tval_sig2*ones(1, 2), '--g', 'linewidth', 1.5);
    %             hold on;
    %             plot(xLIMA, -tval_sig1*ones(1, 2), '--y', 'linewidth', 1.5);
    %             hold on;
    %             plot(xLIMA, -tval_sig2*ones(1, 2), '--g', 'linewidth', 1.5);
    %             hold off;
    %             try
    %                 set(sh, 'XTickLabelRotation', 90);
    %             catch
    %                 xticklabel_rotate([],90,[],'interpreter','none', 'Color', FONT_COLOR);
    %             end
    %
    %             legend(legendStr{:}, 'location', 'bestoutside');
    %             axis(sh, 'tight');
    %
    %             set(sh, 'Xcolor', FONT_COLOR);
    %             set(sh, 'Ycolor', FONT_COLOR);
    %
    %         end
    %
    %     end
    
end



%% Kurtosis of timecourses and spatial maps
% Mean across subjects is reported in table. Figure shows mean+/- SEM across subjects
if (exist('kurt_comp', 'var'))
    
    tmp_tc = squeeze(mean(kurt_comp.tc, 2));
    tmp_ic = squeeze(mean(kurt_comp.ic, 2));
    
    if (sesInfo.numOfSub > 1)
        values = [mean(tmp_tc); mean(tmp_ic)];
        errors = [std(tmp_tc)/sqrt(size(tmp_tc, 1)); std(tmp_ic)/sqrt(size(tmp_ic, 1))];
    else
        values = [tmp_tc(:)'; tmp_ic(:)'];
        errors = zeros(size(values));
    end
    
    varNames = {'ComponentNumber', 'Timecourses', 'SpatialMaps'};
    
    try
        T = table(components, values(1, :)',  values(2, :)', 'VariableNames', varNames);
        disp(T);
    catch
        fprintf('%20s\t%20s\t%20s\n', varNames{:});
        fprintf('\n');
        for n = 1:length(components)
            fprintf('%20s\t%20s\t%20s\n', num2str(components(n), '%d'), num2str(values(1, n), '%0.3f'), num2str(values(2, n), '%0.3f'));
        end
        fprintf('\n');
    end
    
    drawnow;
    bw_legend = cellstr([repmat('Comp ', length(components), 1), num2str(components)]);
    kurtH = figure('color', 'w', 'name', 'Kurtosis (Timecourses)', 'position', get(gH, 'position'));
    %plotBars(values, errors, 0.8, {'Timecourses', 'Spatial Maps'}, winter(64), bw_legend);
    plotBars(values(1, :), errors(1, :), 0.8, {'Timecourses'}, winter(64), bw_legend);
    title('Kurtosis (Mean +/- SEM)');
    set(gca, 'Ylim', [min(values(:)-abs(errors(:)))-0.1, max(values(:)+abs(errors(:))+0.1)]);
    
    pause(1);
    kurtH = figure('color', 'w', 'name', 'Kurtosis (SM)', 'position', get(gH, 'position'));
    plotBars(values(2, :), errors(2, :), 0.8, {'Spatial Maps'}, winter(64), bw_legend);
    title('Kurtosis (Mean +/- SEM)');
    set(gca, 'Ylim', [min(values(:)-abs(errors(:)))-0.1, max(values(:)+abs(errors(:))+0.1)]);
    
    clear values errors;
    drawnow;
    
end

%% FNC correlations
% Functional network connectivity correlations are computed for each
% data-set and averaged across sessions.
if (sesInfo.numOfSess > 1)
    fnc_corrs_all = reshape(mean(fnc_corrs_all, 2), sesInfo.numOfSub, sesInfo.numComp, sesInfo.numComp);
else
    fnc_corrs_all = reshape(squeeze(fnc_corrs_all), sesInfo.numOfSub, sesInfo.numComp, sesInfo.numComp);
end

fnc_grps = cell(1, length(groupsInfo));
grp_names = cell(1, length(groupsInfo));
CLIMG = [];

if (sesInfo.numOfSub > 1)
    for n = 1:length(groupsInfo)
        grp_names{n} = groupsInfo(n).name;
        tmp = icatb_z_to_r(squeeze(mean(fnc_corrs_all(groupsInfo(n).val, :, :))));
        fnc_grps{n} = tmp;
        CLIMG = max([CLIMG, max(abs(tmp(:)))]);
    end
    fnc_corrs_all = squeeze(mean(fnc_corrs_all));
else
    fnc_corrs_all = squeeze(fnc_corrs_all);
end
fnc_corrs_all = icatb_z_to_r(fnc_corrs_all);
comps = (1:sesInfo.numComp)';
CLIM = max(abs(fnc_corrs_all(:)));
gH = figure('color', 'w');
set(gH, 'resize', 'on');
axesH = axes('parent', gH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
icatb_plot_FNC(fnc_corrs_all, [-CLIM, CLIM], cellstr(num2str(comps)), (1:length(comps)), gH, [], axesH);
colormap(jet(64));
title('Average FNC Correlations', 'parent', axesH);

for n = 1:length(grp_names)
    gH = figure('color', 'w');
    set(gH, 'resize', 'on');
    axesH = axes('parent', gH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
    icatb_plot_FNC(fnc_grps{n}, [-CLIMG, CLIMG], cellstr(num2str(comps)), (1:length(comps)), gH, [], axesH);
    colormap(jet(64));
    title(['Average FNC Correlations of ', grp_names{n}], 'parent', axesH);
end


%% FNC metrics of component spatial maps
% Mutual information is computed between components spatially and averaged
% across data-sets.

if (sesInfo.numOfSess > 1)
    spatial_maps_MI = reshape(mean(spatial_maps_MI, 2), sesInfo.numOfSub, sesInfo.numComp, sesInfo.numComp);
else
    spatial_maps_MI = reshape(squeeze(spatial_maps_MI), sesInfo.numOfSub, sesInfo.numComp, sesInfo.numComp);
end

fnc_grps = cell(1, length(groupsInfo));
grp_names = cell(1, length(groupsInfo));
CLIMG = [];

if (sesInfo.numOfSub > 1)
    for n = 1:length(groupsInfo)
        grp_names{n} = groupsInfo(n).name;
        tmp = (squeeze(mean(spatial_maps_MI(groupsInfo(n).val, :, :))));
        fnc_grps{n} = tmp;
        CLIMG = max([CLIMG, max(abs(tmp(:)))]);
    end
end

if (sesInfo.numOfSub > 1)
    spatial_maps_MI = squeeze(mean(spatial_maps_MI));
else
    spatial_maps_MI = squeeze(spatial_maps_MI);
end


comps = (1:sesInfo.numComp)';
CLIM = [min(abs(spatial_maps_MI(:))), max(abs(spatial_maps_MI(:)))];
gH = figure('color', 'w');
set(gH, 'resize', 'on');
axesH = axes('parent', gH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
icatb_plot_FNC(spatial_maps_MI, CLIM, cellstr(num2str(comps)), (1:length(comps)), gH, ' ', axesH);
colormap(jet(64));
title('Average FNC metrics (Spatial maps)', 'parent', axesH);

for n = 1:length(grp_names)
    gH = figure('color', 'w');
    set(gH, 'resize', 'on');
    axesH = axes('parent', gH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
    icatb_plot_FNC(fnc_grps{n}, [0, CLIMG], cellstr(num2str(comps)), (1:length(comps)), gH, ' ', axesH);
    colormap(jet(64));
    title(['Average spatial FNC metrics of ', grp_names{n}], 'parent', axesH);
end

if (exist('files_in_zip', 'var') && ~isempty(files_in_zip))
    icatb_cleanupFiles(files_in_zip, outputDir);
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

data = minVal*ones(maxSizeX, n1 + n2 + n3);

e = 0;
for nS = 1:length(slices)
    tmp = slices{nS};
    s = e + 1;
    e = e + size(tmp, 2);
    inda = ceil(maxSizeX/2) - ceil(size(tmp, 1)/2) + 1;
    indb = inda + size(tmp, 1) - 1;
    data(inda:indb, s:e) = tmp;
end



function plotStackedOrtho(file_name, varargin)

icatb_defaults;
global FONT_COLOR;

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
    end
end

hD = icatb_orth_views(file_name, 'structfile', structFile, 'image_values', image_values, 'convert_to_zscores', convert_to_zscores, 'threshold', threshold, 'set_to_max_voxel', 1, ...
    'get_interp_data', 1, 'cmap', cmap);
if (~exist('sh', 'var'))
    sh = gca;
end
[sliceXY, sliceXZ, sliceYZ] = returnSlices(hD.data, hD.maxVoxelPos);
tmp = stackData({sliceYZ, sliceXZ, sliceXY}, 101);
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
    chPos(2) = chPos(2) - 0.02;
    chPos(3) = 0.4;
    chPos(4) = 0.025;
    set(ch, 'position', chPos);
    set(ch, 'xLim', [1, 100]);
    xTicks = get(ch, 'xTick');
    set(ch, 'yTick', []);
    good_inds = abs(hD.data) > eps;
    minICAVal = hD.minICAIM;
    maxICAVal = hD.maxICAIM;
    set(ch, 'xTick', [xTicks(1), xTicks(end)]);
    set(ch, 'xTicklabel', cellstr(char(num2str(minICAVal, '%0.1f'), num2str(maxICAVal, '%0.1f'))));
    set(ch, 'XColor', FONT_COLOR, 'YColor', FONT_COLOR);
end

function cmap = getCmap(image_values)

load icatb_colors coldhot coldhot_sensitive;

returnValue = strmatch(lower(image_values), {'positive and negative', 'positive', 'absolute value', 'negative'}, 'exact');

if (returnValue == 1)
    coldhot_sensitive = coldhot;
    cmap = coldhot_sensitive(1:4:end, :);
elseif (returnValue == 4)
    cmap = coldhot_sensitive(1:128, :);
    cmap = cmap(1:2:end, :);
else
    cmap = coldhot_sensitive(129:end, :);
    cmap = cmap(1:2:end, :);
end

function c = charC(a, b)

len = max([length(a), length(b)]);

c = repmat(' ', 2, len);

s = ceil(len/2 - length(a)/2) + 1;
c(1, s:s+length(a) - 1) = a;
s = ceil(len/2 - length(b)/2) + 1;
c(2, s:s + length(b) - 1) = b;


function newText = dispParams(sesInfo)

icaAlgo = icatb_icaAlgorithm;
selected_ica_algorithm = deblank(icaAlgo(sesInfo.algorithm, :));

D(1).string = '....................................................';

D(size(D,2)+1).string = '';

modalityType = 'fmri';
try
    modalityType = sesInfo.modality;
catch
end

if (strcmpi(modalityType, 'fmri'))
    D(size(D,2)+1).string = ['Number of Subjects : ', num2str(sesInfo.numOfSub)];
    D(size(D,2)+1).string = ['Number of Sessions : ', num2str(sesInfo.numOfSess)];
end

D(size(D,2)+1).string = ['Number of Independent Components : ', num2str(sesInfo.numComp)];
D(size(D,2)+1).string = ['ICA Algorithm : ', selected_ica_algorithm];

if (~strcmpi(modalityType, 'smri'))
    D(size(D,2)+1).string = ['Number Of Scans/Timepoints : ', num2str(sesInfo.numOfScans)];
else
    D(size(D,2)+1).string = ['Number Of Subjects : ', num2str(sesInfo.numOfScans)];
end

[modality, dataTitle] = icatb_get_modality;

if(isempty(sesInfo.userInput.maskFile))
    D(size(D,2)+1).string = ['Mask File : Default Mask Created From ', dataTitle, ' Data'];
else
    [pathstr,name] = fileparts(sesInfo.userInput.maskFile);
    D(size(D,2)+1).string = ['Mask File : ',name];
end


preproc_options = icatb_preproc_data;
preproc_options = cellstr(preproc_options);
preproc_type = 'remove mean per timepoint';
if (isfield(sesInfo, 'preproc_type'))
    preproc_type = sesInfo.preproc_type;
end

preprocInd = strmatch(lower(preproc_type), lower(preproc_options), 'exact');

D(size(D,2)+1).string = ['Data Pre-processing Type : ', preproc_options{preprocInd}];

pcaTypes = icatb_pca_options;
pcaTypes = cellstr(pcaTypes);
pcaType = 'standard';
if (isfield(sesInfo, 'pcaType'))
    pcaType = sesInfo.pcaType;
end

if (isnumeric(pcaType))
    pcaType = pcaTypes{pcaType};
end

pcaInd = strmatch(lower(pcaType), lower(pcaTypes), 'exact');

D(size(D,2)+1).string = ['PCA Type : ', pcaTypes{pcaInd}];

if (~strcmpi(modalityType, 'smri'))
    
    multiSubGroupPCA = ((sesInfo.numReductionSteps == 1) && (sesInfo.numOfSub*sesInfo.numOfSess > 1));
    
    if (~multiSubGroupPCA || strcmpi(selected_ica_algorithm, 'iva-gl'))
        groupPCAOpts = char('Subject Specific', 'Grand Mean');
        group_pca_type = 'Subject Specific';
        if (isfield(sesInfo, 'group_pca_type'))
            group_pca_type = sesInfo.group_pca_type;
        end
        groupPCAInd = strmatch(lower(group_pca_type), lower(groupPCAOpts), 'exact');
        group_pca_type = deblank(groupPCAOpts(groupPCAInd, :));
        
        D(size(D,2)+1).string = ['Group PCA Type : ', group_pca_type];
    end
end

useTemporalICA = 0;
if (strcmpi(modalityType, 'fmri'))
    group_ica_type = 'spatial';
    try
        group_ica_type = sesInfo.group_ica_type;
    catch
    end
    useTemporalICA = strcmpi(group_ica_type, 'temporal');
    D(size(D,2)+1).string = ['Group ICA Type : ', upper(group_ica_type(1)), group_ica_type(2:end)];
end

if strcmpi(selected_ica_algorithm, 'moo-icar')
    selected_ica_algorithm = 'gig-ica';
end


if (~(strcmpi(selected_ica_algorithm, 'iva-gl') || strcmpi(selected_ica_algorithm, 'gig-ica') || strcmpi(selected_ica_algorithm, 'constrained ica (spatial)') || useTemporalICA))
    backReconType = 'Regular';
    if (isfield(sesInfo, 'backReconType'))
        backReconType = sesInfo.backReconType;
    end
    
    backReconOptions = cellstr(icatb_backReconOptions);
    backReconInd = strmatch(lower(backReconType), lower(backReconOptions), 'exact');
    
    backReconType = backReconOptions{backReconInd};
    if ((sesInfo.numReductionSteps == 1) && (sesInfo.numOfSub*sesInfo.numOfSess > 1))
        backReconType = 'Spatial-temporal Regression';
    end
    
    D(size(D, 2)+1).string = ['Back Reconstruction Type : ', backReconType];
end

calibrateOptions = icatb_scaleICA;
D(size(D,2)+1).string = ['Scaling Components : ', deblank(calibrateOptions(sesInfo.scaleType + 1, :))];

stability_analysis = 'none';
whichAnalysis = 1;
try
    whichAnalysis = sesInfo.which_analysis;
catch
end

if (whichAnalysis == 2)
    stability_analysis = 'ICASSO';
end

if (whichAnalysis == 3)
    stability_analysis = 'MST';
end

D(size(D,2)+1).string = ['Stability analysis type : ', stability_analysis];

checkAnalysisMode = 'Serial';
try
    checkAnalysisMode = sesInfo.parallel_info.mode;
catch
end

D(size(D,2)+1).string = ['Group analysis mode: ', upper(checkAnalysisMode(1)), checkAnalysisMode(2:end)];
D(size(D,2)+1).string = ['Anatomical file: ', sesInfo.structFile];
D(size(D,2)+1).string = ['Slice Plane: ', upper(sesInfo.slice_plane(1)), sesInfo.slice_plane(2:end)];
D(size(D,2)+1).string = ['Image values: ', upper(sesInfo.image_values(1)), sesInfo.image_values(2:end)];
D(size(D,2)+1).string = ['Convert to Z-scores: ', sesInfo.convert_to_zscores];
D(size(D,2)+1).string = ['Threshold: ', num2str(sesInfo.threshold)];
if (~isempty(sesInfo.groupsInfo))
    D(size(D,2)+1).string = 'Group Info: ';
    for nG = 1:length(sesInfo.groupsInfo)
        D(size(D,2)+1).string = ['Name: ', sesInfo.groupsInfo(nG).name, ' Subjects: ', num2str(sesInfo.groupsInfo(nG).val)];
    end
end
D(size(D,2)+1).string = '';
D(size(D,2)+1).string = '....................................................';
D(size(D,2)+1).string = '';
newText = char(D.string);


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


function hText = xticklabel_rotate(XTick,rot,varargin)
%hText = xticklabel_rotate(XTick,rot,XTickLabel,varargin)     Rotate XTickLabel
%
% Syntax: xticklabel_rotate
%
% Input:
% {opt}     XTick       - vector array of XTick positions & values (numeric)
%                           uses current XTick values or XTickLabel cell array by
%                           default (if empty)
% {opt}     rot         - angle of rotation in degrees, 90? by default
% {opt}     XTickLabel  - cell array of label strings
% {opt}     [var]       - "Property-value" pairs passed to text generator
%                           ex: 'interpreter','none'
%                               'Color','m','Fontweight','bold'
%
% Output:   hText       - handle vector to text labels
%
% Example 1:  Rotate existing XTickLabels at their current position by 90?
%    xticklabel_rotate
%
% Example 2:  Rotate existing XTickLabels at their current position by 45? and change
% font size
%    xticklabel_rotate([],45,[],'Fontsize',14)
%
% Example 3:  Set the positions of the XTicks and rotate them 90?
%    figure;  plot([1960:2004],randn(45,1)); xlim([1960 2004]);
%    xticklabel_rotate([1960:2:2004]);
%
% Example 4:  Use text labels at XTick positions rotated 45? without tex interpreter
%    xticklabel_rotate(XTick,45,NameFields,'interpreter','none');
%
% Example 5:  Use text labels rotated 90? at current positions
%    xticklabel_rotate([],90,NameFields);
%
% Example 6:  Multiline labels
%    figure;plot([1:4],[1:4])
%    axis([0.5 4.5 1 4])
%    xticklabel_rotate([1:4],45,{{'aaa' 'AA'};{'bbb' 'AA'};{'ccc' 'BB'};{'ddd' 'BB'}})
%
% Note : you can not RE-RUN xticklabel_rotate on the same graph.
%



% This is a modified version of xticklabel_rotate90 by Denis Gilbert
% Modifications include Text labels (in the form of cell array)
%                       Arbitrary angle rotation
%                       Output of text handles
%                       Resizing of axes and title/xlabel/ylabel positions to maintain same overall size
%                           and keep text on plot
%                           (handles small window resizing after, but not well due to proportional placement with
%                           fixed font size. To fix this would require a serious resize function)
%                       Uses current XTick by default
%                       Uses current XTickLabel is different from XTick values (meaning has been already defined)

% Brian FG Katz
% bfgkatz@hotmail.com
% 23-05-03
% Modified 03-11-06 after user comment
%	Allow for exisiting XTickLabel cell array
% Modified 03-03-2006
%   Allow for labels top located (after user comment)
%   Allow case for single XTickLabelName (after user comment)
%   Reduced the degree of resizing
% Modified 11-jun-2010
%   Response to numerous suggestions on MatlabCentral to improve certain
%   errors.
% Modified 23-sep-2014
%   Allow for mutliline labels


% Other m-files required: cell2mat
% Subfunctions: none
% MAT-files required: none
%
% See also: xticklabel_rotate90, TEXT,  SET

% Based on xticklabel_rotate90
%   Author: Denis Gilbert, Ph.D., physical oceanography
%   Maurice Lamontagne Institute, Dept. of Fisheries and Oceans Canada
%   email: gilbertd@dfo-mpo.gc.ca  Web: http://www.qc.dfo-mpo.gc.ca/iml/
%   February 1998; Last revision: 24-Mar-2003

% check to see if xticklabel_rotate has already been here (no other reason for this to happen)
if isempty(get(gca,'XTickLabel')),
    error('xticklabel_rotate : can not process, either xticklabel_rotate has already been run or XTickLabel field has been erased')  ;
end

% if no XTickLabel AND no XTick are defined use the current XTickLabel
%if nargin < 3 & (~exist('XTick') | isempty(XTick)),
% Modified with forum comment by "Nathan Pust" allow the current text labels to be used and property value pairs to be changed for those labels
if (nargin < 3 || isempty(varargin{1})) & (~exist('XTick') | isempty(XTick)),
    xTickLabels = get(gca,'XTickLabel')  ; % use current XTickLabel
    if ~iscell(xTickLabels)
        % remove trailing spaces if exist (typical with auto generated XTickLabel)
        temp1 = num2cell(xTickLabels,2)         ;
        for loop = 1:length(temp1),
            temp1{loop} = deblank(temp1{loop})  ;
        end
        xTickLabels = temp1                     ;
    end
    varargin = varargin(2:length(varargin));
end

% if no XTick is defined use the current XTick
if (~exist('XTick') | isempty(XTick)),
    XTick = get(gca,'XTick')        ; % use current XTick
end

%Make XTick a column vector
XTick = XTick(:);

if ~exist('xTickLabels'),
    % Define the xtickLabels
    % If XtickLabel is passed as a cell array then use the text
    if (length(varargin)>0) & (iscell(varargin{1})),
        xTickLabels = varargin{1};
        varargin = varargin(2:length(varargin));
    else
        xTickLabels = num2str(XTick);
    end
end

if length(XTick) ~= length(xTickLabels),
    error('xticklabel_rotate : must have same number of elements in "XTick" and "XTickLabel"')  ;
end

%Set the Xtick locations and set XTicklabel to an empty string
set(gca,'XTick',XTick,'XTickLabel','')

if nargin < 2,
    rot = 90 ;
end

% Determine the location of the labels based on the position
% of the xlabel
hxLabel = get(gca,'XLabel');  % Handle to xlabel
xLabelString = get(hxLabel,'String');

% if ~isempty(xLabelString)
%    warning('You may need to manually reset the XLABEL vertical position')
% end

set(hxLabel,'Units','data');
xLabelPosition = get(hxLabel,'Position');
y = xLabelPosition(2);

%CODE below was modified following suggestions from Urs Schwarz
y=repmat(y,size(XTick,1),1);
% retrieve current axis' fontsize
fs = get(gca,'fontsize');

if ~iscell(xTickLabels)
    % Place the new xTickLabels by creating TEXT objects
    hText = text(XTick, y, xTickLabels,'fontsize',fs);
else
    % Place multi-line text approximately where tick labels belong
    for cnt=1:length(XTick),
        hText(cnt) = text(XTick(cnt),y(cnt),xTickLabels{cnt},...
            'VerticalAlignment','top', 'UserData','xtick');
    end
end

% Rotate the text objects by ROT degrees
%set(hText,'Rotation',rot,'HorizontalAlignment','right',varargin{:})
% Modified with modified forum comment by "Korey Y" to deal with labels at top
% Further edits added for axis position
xAxisLocation = get(gca, 'XAxisLocation');
if strcmp(xAxisLocation,'bottom')
    set(hText,'Rotation',rot,'HorizontalAlignment','right',varargin{:})
else
    set(hText,'Rotation',rot,'HorizontalAlignment','left',varargin{:})
end

% Adjust the size of the axis to accomodate for longest label (like if they are text ones)
% This approach keeps the top of the graph at the same place and tries to keep xlabel at the same place
% This approach keeps the right side of the graph at the same place

set(get(gca,'xlabel'),'units','data')           ;
labxorigpos_data = get(get(gca,'xlabel'),'position')  ;
set(get(gca,'ylabel'),'units','data')           ;
labyorigpos_data = get(get(gca,'ylabel'),'position')  ;
set(get(gca,'title'),'units','data')           ;
labtorigpos_data = get(get(gca,'title'),'position')  ;

set(gca,'units','pixel')                        ;
set(hText,'units','pixel')                      ;
set(get(gca,'xlabel'),'units','pixel')          ;
set(get(gca,'ylabel'),'units','pixel')          ;
% set(gca,'units','normalized')                        ;
% set(hText,'units','normalized')                      ;
% set(get(gca,'xlabel'),'units','normalized')          ;
% set(get(gca,'ylabel'),'units','normalized')          ;

origpos = get(gca,'position')                   ;

% textsizes = cell2mat(get(hText,'extent'))       ;
% Modified with forum comment from "Peter Pan" to deal with case when only one XTickLabelName is given.
x = get( hText, 'extent' );
if iscell( x ) == true
    textsizes = cell2mat( x ) ;
else
    textsizes = x;
end

largest =  max(textsizes(:,3))                  ;
longest =  max(textsizes(:,4))                  ;

laborigext = get(get(gca,'xlabel'),'extent')    ;
laborigpos = get(get(gca,'xlabel'),'position')  ;

labyorigext = get(get(gca,'ylabel'),'extent')   ;
labyorigpos = get(get(gca,'ylabel'),'position') ;
leftlabdist = labyorigpos(1) + labyorigext(1)   ;

% assume first entry is the farthest left
leftpos = get(hText(1),'position')              ;
leftext = get(hText(1),'extent')                ;
leftdist = leftpos(1) + leftext(1)              ;
if leftdist > 0,    leftdist = 0 ; end          % only correct for off screen problems

% botdist = origpos(2) + laborigpos(2)            ;
% newpos = [origpos(1)-leftdist longest+botdist origpos(3)+leftdist origpos(4)-longest+origpos(2)-botdist]
%
% Modified to allow for top axis labels and to minimize axis resizing
if strcmp(xAxisLocation,'bottom')
    newpos = [origpos(1)-(min(leftdist,labyorigpos(1)))+labyorigpos(1) ...
        origpos(2)+((longest+laborigpos(2))-get(gca,'FontSize')) ...
        origpos(3)-(min(leftdist,labyorigpos(1)))+labyorigpos(1)-largest ...
        origpos(4)-((longest+laborigpos(2))-get(gca,'FontSize'))]  ;
else
    newpos = [origpos(1)-(min(leftdist,labyorigpos(1)))+labyorigpos(1) ...
        origpos(2) ...
        origpos(3)-(min(leftdist,labyorigpos(1)))+labyorigpos(1)-largest ...
        origpos(4)-(longest)+get(gca,'FontSize')]  ;
end
set(gca,'position',newpos)                      ;

% readjust position of text labels after resize of plot
set(hText,'units','data')                       ;
for loop= 1:length(hText),
    set(hText(loop),'position',[XTick(loop), y(loop)])  ;
end

% adjust position of xlabel and ylabel
laborigpos = get(get(gca,'xlabel'),'position')  ;
set(get(gca,'xlabel'),'position',[laborigpos(1) laborigpos(2)-longest 0])   ;

% switch to data coord and fix it all
set(get(gca,'ylabel'),'units','data')                   ;
set(get(gca,'ylabel'),'position',labyorigpos_data)      ;
set(get(gca,'title'),'position',labtorigpos_data)       ;

set(get(gca,'xlabel'),'units','data')                   ;
labxorigpos_data_new = get(get(gca,'xlabel'),'position')  ;
set(get(gca,'xlabel'),'position',[labxorigpos_data(1) labxorigpos_data_new(2)])   ;


% Reset all units to normalized to allow future resizing
set(get(gca,'xlabel'),'units','normalized')          ;
set(get(gca,'ylabel'),'units','normalized')          ;
set(get(gca,'title'),'units','normalized')          ;
set(hText,'units','normalized')                      ;
set(gca,'units','normalized')                        ;

if nargout < 1,
    clear hText
end

