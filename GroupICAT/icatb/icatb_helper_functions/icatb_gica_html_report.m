function icatb_gica_html_report(param_file, opts)

global GICA_PARAM_FILE;

icatb_defaults;
global PARAMETER_INFO_MAT_FILE;
global CONVERT_Z;
global IMAGE_VALUES;
global THRESHOLD_VALUE;
global ANATOMICAL_PLANE;
global TIMECOURSE_POSTPROCESS;
global FONT_COLOR;
global GICA_RESULTS_SUMMARY;
global UI_FS;

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

structFile = fullfile(fileparts(which('groupica.m')), 'icatb_templates', 'ch2bet_3x3x3.nii');
try
    structFile = opts.anatomical_file;
catch
end

[~, structHInfo] = icatb_returnHInfo(structFile);

if ~(all(structHInfo.VOX == structHInfo.VOX(1)))
    disp('Spatial template is not isotropic. Reslicing spatial template ...');
    structFile = noisecloud_spm_coregister(fullfile(fileparts(which('groupica.m')), 'icatb_templates', 'ch2bet_3x3x3.nii'), ...
        deblank(structFile), '', outputDir);
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
spatial_maps_MI = [];
compute_mi = 1;
compute_kurtosis = 1;
compute_fnc = 1;
try
    compute_mi = GICA_RESULTS_SUMMARY.compute.mi;
    compute_kurtosis = GICA_RESULTS_SUMMARY.compute.kurtosis;
    compute_fnc = GICA_RESULTS_SUMMARY.compute.fnc;
catch
end

try
    groupsInfo = opts.groupsInfo;
catch
end

saveFigInfo = 0;
try
    saveFigInfo = opts.save_info;
catch
end

figVisible = 'on';
if (saveFigInfo)
    
    figVisible = 'off';
end

resultsFormat = 'html';
try
    resultsFormat = opts.format;
catch
end

resultsDir =  fullfile(outputDir, [sesInfo.userInput.prefix, '_gica_results']);
try
    resultsDir = opts.outputDir;
catch
end

printRes = '';

try
    printRes = GICA_RESULTS_SUMMARY.print_resolution;
catch
end

if (isempty(printRes))
    printRes = '-r72';
end

if (saveFigInfo)
    if (exist(resultsDir, 'dir') ~= 7)
        mkdir(resultsDir);
    end
    pdfPrefix = mfilename;
end

resultsInfo = [];

try
    print_opts = GICA_RESULTS_SUMMARY.print_opts;
catch
    print_opts = {'-noui'};
end
printOptions = [{'-dpdf', printRes}, print_opts];
if (icatb_get_matlab_version < 2016)
    printOptions(strcmpi(printOptions, '-bestfit')) = [];
    printOptions(strcmpi(printOptions, '-fillpage')) = [];
end


components = (1:sesInfo.numComp)';
subjects_to_use = (1:sesInfo.numOfSub);

use_parallel = 0;

try
    computeInfo = opts.compute;
    subjects_to_use = computeInfo.subjects;
    components = computeInfo.components;
catch
end

try
    use_parallel = computeInfo.use_parallel;
catch
end

try
    compute_mi = opts.compute.mi;
catch
end

try
    compute_kurtosis = opts.compute.kurtosis;
catch
end

try
    compute_fnc = opts.compute.fnc;
catch
end

components = components(:);

%% Group ICA Parameters
sesInfo.structFile = structFile;
sesInfo.slice_plane = slice_plane;
sesInfo.image_values = image_values;
sesInfo.threshold = threshold;
sesInfo.groupsInfo = groupsInfo;
sesInfo.convert_to_zscores = convert_to_zscores;
newText = dispParams(sesInfo, subjects_to_use, components);
disp(newText);

pdfFiles = {};

htmlSummaryStr = [];

if (saveFigInfo)
    if (~strcmpi(resultsFormat, 'pdf'))
        icatb_print_table(cellstr(newText), fullfile(resultsDir, 'params.txt'));
        resultsInfo(end + 1).title = 'Group ICA Parameters';
        resultsInfo(end).text = ' ';
        resultsInfo(end).files = {'params.txt'};
        htmlSummaryStr(end + 1).title = 'Group ICA Parameters';
        htmlSummaryStr(end).tag =  'group_ica_params';
        htmlSummaryStr(end).str = get_result_strings(resultsDir, resultsInfo(end), htmlSummaryStr(end).tag);
        %writeHTML(resultsDir, resultsInfo(end), 'introduction.html', 'Group ICA Parameters');
    else
        print_strs = cellstr(char('Group ICA parameters:', '', newText));
        yPos = 1.1;
        gH = icatb_getGraphics(print_strs{1}, 'graphics', 'contentspdf', figVisible);
        aH = axes('units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8], 'parent', gH, 'visible', 'off', 'color', get(gH, 'color'));
        for n = 1:length(print_strs)
            if (n == 1)
                tmpFont = UI_FS;
            else
                tmpFont = UI_FS - 2;
            end
            text(0.45, yPos, print_strs{n}, 'units', 'normalized', 'horizontalalignment', 'center', 'verticalAlignment', 'middle', 'color', FONT_COLOR, 'fontsize', tmpFont);
            yPos = yPos - 0.075;
        end
        
        tmpImFile = [pdfPrefix, '_', icatb_returnFileIndex(1), '.pdf'];
        set(gH, 'PaperPositionMode', 'auto');
        print(gH, fullfile(resultsDir, tmpImFile), printOptions{:});
        countPdfs = 1;
        pdfFiles{countPdfs} = fullfile(resultsDir, tmpImFile);
        delete(gH);
        
    end
end

resultsInfo =[];

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
    
    icassoFigs = [];
    outFiles = {''};
    if (exist(icassoResultsFile, 'file'))
        load(icassoResultsFile, 'sR');
        try
            icassoShow(sR, 'L', sesInfo.numComp, 'colorlimit', [.8 .9]);
            try
                icassoFigs = findobj(0, 'type', 'figure');
                icassoFigNames = get(icassoFigs, 'name');
                icassoFigInds = icatb_good_cells(regexpi(icassoFigNames, 'icasso:'));
                icassoFigs = icassoFigs(icassoFigInds);
                if (saveFigInfo)
                    if (~strcmpi(resultsFormat, 'pdf'))
                        outFiles = printFigs(icassoFigs, resultsDir, 'icasso');
                        resultsInfo(end+1).title = 'ICASSO Plots';
                        resultsInfo(end).text = ' ';
                        resultsInfo(end).files = outFiles;
                        %writeHTML(resultsDir, resultsInfo(end), 'icasso.html', 'ICASSO');
                        htmlSummaryStr(end + 1). title = 'ICASSO Plots';
                        htmlSummaryStr(end).tag =  'icasso';
                        htmlSummaryStr(end).str = get_result_strings(resultsDir, resultsInfo(end), htmlSummaryStr(end).tag);
                    else
                        
                        for nPdfs = 1:length(icassoFigs)
                            countPdfs = countPdfs + 1;
                            tmpImFile = [pdfPrefix, '_', num2str(countPdfs), '.pdf'];
                            set(icassoFigs(nPdfs), 'PaperPositionMode', 'auto');
                            print(icassoFigs(nPdfs), fullfile(resultsDir, tmpImFile), printOptions{:});
                            pdfFiles{countPdfs} = fullfile(resultsDir, tmpImFile);
                        end
                        
                        
                    end
                    delete(icassoFigs);
                end
            catch
            end
        catch
        end
    end
    
end

drawnow;

%% Mean Components
% Mean across all subjects and sessions is computed for each component
%
% * *a) Timecourse* - Mean timecourse is converted to z-scores.
% * *b) Spectra* - Timecourses spectra is computed for each data-set and averaged across sessions. Mean and standard error of mean is shown in the figure.
% * *c) Montage* - Axial slices are shown.
% * *d) Ortho slices* - Ortho plot is shown for the peak voxel and coordinates are reported.

resultsInfo = [];
outFiles = {''};
% if (saveFigInfo)
%
%     resultsInfo(end + 1).title = 'Mean Components';
%     resultsInfo(end).text = ['<ul> <li> <b> a) Timecourse </b> - Mean timecourse is converted to z-scores. </li>', ...
%         '<li> <b> b) Spectra </b> - Timecourses spectra is computed for each data-set and averaged across sessions. Mean and standard error of mean is shown in the figure. </li> ', ...
%         '<li> <b> c) Montage </b> - Axial slices are shown. </li>', ...
%         '<li> <b> a) Ortho slices </b> - Ortho plot is shown for the peak voxel and coordinates are reported. </li> </ul>'];
%
%
% end


compFileNaming = sesInfo.icaOutputFiles(1).ses(1).name;
currentFile = deblank(compFileNaming(1, :));


if (~exist('computeInfo', 'var'))
    
    if ~exist(fullfile(outputDir, currentFile), 'file')
        [zipFileName, files_in_zip] = icatb_getViewingSet_zip(currentFile, [], 'real', sesInfo.zipContents);
        if (~isempty(zipFileName))
            icatb_unzip(regexprep(zipFileName, ['.*\', filesep], ''), fullfile(outputDir, fileparts(currentFile)));
        end
    end
    
    compFiles = icatb_rename_4d_file(icatb_fullFile('directory', outputDir, 'files', compFileNaming));
    
    postProcessFile = fullfile(outputDir, [sesInfo.userInput.prefix, '_postprocess_results.mat']);
    
else
    
    if (use_parallel)
        icatb_par_postprocess_timecourses(param_file, 'subjects', subjects_to_use, 'components', components, 'outputDir', resultsDir, ...
            'compute_mi', compute_mi, 'compute_kurtosis', compute_kurtosis, 'compute_fnc', compute_fnc);
    else
        icatb_postprocess_timecourses(param_file, 'subjects', subjects_to_use, 'components', components, 'outputDir', resultsDir, ...
            'compute_mi', compute_mi, 'compute_kurtosis', compute_kurtosis, 'compute_fnc', compute_fnc);
    end
    
    postProcessFile = fullfile(resultsDir, [sesInfo.userInput.prefix, '_postprocess_results.mat']);
    
    writeComponentMaps (sesInfo, compFileNaming, subjects_to_use, components, resultsDir);
    
    compFiles = icatb_rename_4d_file(icatb_fullFile('directory', resultsDir, 'files', compFileNaming));
    
end

numComp = length(components);
numOfSub = length(subjects_to_use);

icaTimecourse = icatb_loadICATimeCourse(compFiles);

%postProcessFile = fullfile(outputDir, [sesInfo.userInput.prefix, '_postprocess_results.mat']);
if (~exist(postProcessFile, 'file'))
    icatb_par_postprocess_timecourses(param_file);
end

varsInMat = whos('-file', postProcessFile);
chkAgg = strmatch('aggregate', cellstr(char(varsInMat.name)), 'exact');
if (isempty(chkAgg))
    aggregate = getAggInfo(sesInfo, postProcessFile);
else
    load(postProcessFile, 'aggregate');
end

% matFileInfo = matfile(postProcessFile, 'Writable', true);
% freq = matFileInfo.freq;
%
% matFileInfo.fALFF = [];
% matFileInfo.dynamic_range = [];
% matFileInfo.dynamic_range_all = [];
% matFileInfo.fALFF_all = [];
%
% disp('Computing mean FNC of timecourses ...');
% [fnc_corrs_all, CLIM_FNC, fnc_grps, CLIMG, grp_names] = computeFNC(postProcessFile, 'fnc_corrs_all', 1, numComp, numOfSub, groupsInfo);
% disp('Done computing mean of FNC');
%
% if (~isempty(who(matFileInfo, 'spatial_maps_MI')))
%     disp('Computing mean FNC of spatial maps ...');
%     [spatial_maps_MI, CLIM_mi, fnc_grps_mi, CLIMG_mi] = computeFNC(postProcessFile, 'spatial_maps_MI', 0, numComp, numOfSub, groupsInfo);
%     disp('Done computing mean FNC of spatial maps');
% end
%
% closeGCP;

%load(postProcessFile);

numOfSess = sesInfo.numOfSess;

clear tc;

freq_limits = [0.1, 0.15];
try
    freq_limits = TIMECOURSE_POSTPROCESS.spectra.freq_limits;
catch
end

cmap = getCmap(image_values);

fALFF = aggregate.spectra.fALFF;
dynamic_range = aggregate.spectra.dynamic_range;

resultsInfo = [];
tmp_fnames = cell(1, size(compFiles, 1));


drawnow;

for nF = 1:size(compFiles, 1)
    
    cn = deblank(compFiles(nF, :));
    
    [dd, fN, extn] = fileparts(cn);
    
    gH = icatb_getGraphics([fN, extn, '(Comp ', icatb_returnFileIndex(components(nF)), ')'], 'graphics', 'imviewer', figVisible);
    set(gH, 'resize', 'on');
    %set(gH, 'colormap', cmap);
    colormap(cmap);
    
    xOffSet = 0.05;
    yOffSet = 0.05;
    
    
    width2 = 0.5;
    height2 = 0.5;
    sh = axes('parent', gH, 'units', 'normalized', 'position', [0.01, yOffSet, width2, height2]);
    icatb_image_viewer(cn, 'structfile', structFile, 'image_values', image_values, 'convert_to_zscores', convert_to_zscores, 'threshold', threshold, ...
        'slices_in_mm', slices_in_mm, 'anatomical_view', slice_plane, 'axesh', sh, 'colorbar', 0, 'labels', ...
        {['(Comp ', icatb_returnFileIndex(components(nF)), ')']});
    
    width = 0.4;
    height = 0.4;
    
    sh = axes('parent', gH, 'units', 'normalized', 'position', [width2 + 0.03, height2/2-(height/2), width, height]);
    plotStackedOrtho(cn, 'structfile', structFile, 'image_values', image_values, 'convert_to_zscores', convert_to_zscores, 'threshold', threshold, 'set_to_max_voxel', 1, ...
        'get_interp_data', 1, 'cmap', cmap, 'axesh', sh, 'colorbar', false, 'labels', 'Peak Coordinates (mm)', 'colorbar', true, 'colorbar_label', true);
    
    drawnow;
    
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
    icatb_plotTimecourse('data', tc, 'parent', sh, 'color', 'm', 'title', ['Component ', icatb_returnFileIndex(components(nF))], ...
        'xlabelstr', 'Scans', 'ylabelstr', tc_label);
    
    sh = axes('parent', gH, 'units', 'normalized', 'position', [xOffSet+width+0.1, yPos, width, height]);
    clear tc;
    
    % chkSize = size(matFileInfo, 'spectra_tc_all');
    
    %pause(0.1);
    
    %     if (numel(chkSize) == 4)
    %         tc.data = reshape(squeeze(mean(matFileInfo.spectra_tc_all(:, :, :, nF), 2)), numOfSub, length(freq));
    %     else
    %         if (numel(chkSize) == 3)
    %             tc.data = matFileInfo.spectra_tc_all(:, :, nF);
    %         else
    %             tc.data = matFileInfo.spectra_tc_all(:, nF);
    %             tc.data = tc.data';
    %         end
    %     end
    %
    %     tc.data = reshape(tc.data, numOfSub, length(freq));
    %
    %     tmp_dynamicrange = zeros(1, size(tc.data, 1));
    %     tmp_fALFF = tmp_dynamicrange;
    %     for nS = 1:length(tmp_dynamicrange)
    %         [tmp_dynamicrange(nS), tmp_fALFF(nS)] = icatb_get_spec_stats(tc.data(nS, :), freq, freq_limits);
    %     end
    %
    %     if (~isempty(matFileInfo.fALFF))
    %         matFileInfo.fALFF(1, nF) = mean(tmp_fALFF);
    %         matFileInfo.dynamic_range(1, nF) = mean(tmp_dynamicrange);
    %     else
    %         matFileInfo.fALFF = reshape(mean(tmp_fALFF), 1, 1);
    %         matFileInfo.dynamic_range = reshape(mean(tmp_dynamicrange), 1, 1);
    %     end
    %
    %     tmp_dynamicrange = tmp_dynamicrange(:);
    %     tmp_fALFF = tmp_fALFF(:);
    %
    %     if (~isempty(matFileInfo.fALFF_all))
    %         matFileInfo.dynamic_range_all(:, nF) = tmp_dynamicrange;
    %         matFileInfo.fALFF_all(:, nF) = tmp_fALFF;
    %     else
    %         matFileInfo.dynamic_range_all = tmp_dynamicrange;
    %         matFileInfo.fALFF_all = tmp_fALFF;
    %     end
    
    tc.xAxis = aggregate.spectra.freq;
    tc.isSpectra = 1;
    tc.xlabelStr = 'Frequency (Hz)';
    tc.ylabelStr = 'Power';
    tc.titleStr = sprintf('Dynamic range: %0.3f, Power_L_F/Power_H_F: %0.3f', dynamic_range(nF), fALFF(nF));
    %icatb_plot_spectra(sh, tc);
    icatb_plot_with_ste_area(sh, aggregate.spectra.freq, aggregate.spectra.mean(:, nF)', aggregate.spectra.sem(:, nF)', 'm', [.5 .5 1]);
    xlabel(tc.xlabelStr, 'parent', sh, 'Interpreter', 'tex');
    ylabel(tc.ylabelStr, 'parent', sh, 'Interpreter', 'tex');
    title(tc.titleStr, 'parent', sh, 'Interpreter', 'tex');
    
    set(sh, 'XColor', FONT_COLOR, 'YColor', FONT_COLOR);
    
    axis(sh, 'tight');
    
    clear tc;
    
    drawnow;
    
    %meanH(end + 1) = gH;
    
    if (nF == 1)
        figurePos = get(gH, 'position');
    end
    
    
    if (saveFigInfo)
        if (~strcmpi(resultsFormat, 'pdf'))
            if (nF == size(compFiles, 1))
                resultsInfo(end + 1).title = 'Mean Components';
                resultsInfo(end).text = ['<ul> <li> <b> a) Timecourse </b> - Mean timecourse is converted to z-scores. </li>', ...
                    '<li> <b> b) Spectra </b> - Timecourses spectra is computed for each data-set and averaged across sessions. Mean and standard error of mean is shown in the figure. </li> ', ...
                    '<li> <b> c) Montage </b> - Axial slices are shown. </li>', ...
                    '<li> <b> a) Ortho slices </b> - Ortho plot is shown for the peak voxel and coordinates are reported. </li> </ul>'];
            end
            
            tmp_fname = printFigs(gH, resultsDir, ['mean_comp', '_', icatb_returnFileIndex(nF)], 'no');
            tmp_fnames{nF} = tmp_fname{1};
            
            %writeHTML(resultsDir, resultsInfo(end), 'mean_comp.html', 'Mean Components');
            if (nF == size(compFiles, 1))
                resultsInfo(end).files = tmp_fnames;
                htmlSummaryStr(end + 1). title = 'Mean Components';
                htmlSummaryStr(end).tag = 'mean_comp';
                htmlSummaryStr(end).str = get_result_strings(resultsDir, resultsInfo(end), htmlSummaryStr(end).tag);
            end
            
        else
            
            %for nPdfs = 1:length(meanH)
            countPdfs = countPdfs + 1;
            tmpImFile = [pdfPrefix, '_', num2str(nF), '.pdf'];
            set(gH, 'PaperPositionMode', 'auto');
            print(gH, fullfile(resultsDir, tmpImFile), printOptions{:});
            pdfFiles{countPdfs} = fullfile(resultsDir, tmpImFile);
            %end
        end
        
        delete(gH);
    end
    
end


% Append fALFF and dynamic range info to post process file
%save(postProcessFile, 'dynamic_range', 'fALFF', 'dynamic_range_all', 'fALFF_all', '-append');

%% Spectral Summary
% * *a) dynamic_range* - Difference between the peak power and minimum power at frequencies to the right of the peak.
% * *b) fALFF* - Low frequency to high frequency power ratio.
%
%fALFF = matFileInfo.fALFF;
%dynamic_range = matFileInfo.dynamic_range;
fALFF = aggregate.spectra.fALFF;
dynamic_range = aggregate.spectra.dynamic_range;

outFiles = {''};
resultsInfo = [];
if (saveFigInfo)
    if (~strcmpi(resultsFormat, 'pdf'))
        
        resultsInfo(end + 1).title = 'Spectral Summary';
        resultsInfo(end).text = ['<ul> <li> <b> a) dynamic range </b> - Difference between the peak power and minimum power at frequencies to the right of the peak. </li>', ...
            '<li> <b> b) fALFF </b> - Low frequency to high frequency power ratio. </li> ', ...
            '</ul>'];
        
        numPara = 1;
        varStruct(numPara).tag = 'Component Number';
        varStruct(numPara).value = components(:);
        
        numPara = numPara + 1;
        varStruct(numPara).tag = 'DynamicRange';
        varStruct(numPara).value = num2str(dynamic_range(:), '%0.3f');
        
        numPara = numPara + 1;
        varStruct(numPara).tag = 'fALFF';
        varStruct(numPara).value = num2str(fALFF(:), '%0.3f');
        resultsInfo(end).files = {'spectra_summary.txt'};
        icatb_printToFile(fullfile(resultsDir,  resultsInfo(end).files{1}), varStruct, '', 'column_wise');
        
        % writeHTML(resultsDir, resultsInfo(end), 'spectral_summary.html', 'Spectral Summary');
        htmlSummaryStr(end + 1). title = 'Spectral Summary';
        htmlSummaryStr(end).tag = 'spectral_summary';
        htmlSummaryStr(end).str = get_result_strings(resultsDir, resultsInfo(end), htmlSummaryStr(end).tag);
        
    end
    
else
    
    
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
    
end


%% Temporal Stats On Beta Weights
% Multiple regression is done using the timecourses from SPM design matrix as model and ICA timecourses as observations. R^2 values for each
% component are shown in bar plot. For each component, one sample t-test results of each session and condition are shown in the bar plots.
%

drawnow;

temporal_stats_betas = [];
try
    temporal_stats_betas = opts.temporal_stats_betas;
catch
end

resultsInfo = [];

if (~isempty(temporal_stats_betas))
    
    varNames = {'ComponentNumber', 'Rsquare'};
    
    
    outFiles = {''};
    if (saveFigInfo)
        
        if (~strcmpi(resultsFormat, 'pdf'))
            
            resultsInfo(end + 1).title = 'Temporal Stats On Beta Weights';
            resultsInfo(end).text = ['Multiple regression is done using the timecourses from SPM design matrix as model and ICA timecourses as observations. R^2 values for each', ...
                'component are shown in bar plot. For each component, one sample t-test results of each session and condition are shown in the bar plots.'];
            
            clear varStruct;
            
            numPara = 1;
            varStruct(numPara).tag = varNames{1};
            varStruct(numPara).value = temporal_stats_betas.component_numbers(:);
            
            numPara = numPara + 1;
            varStruct(numPara).tag = varNames{2};
            varStruct(numPara).value = num2str( temporal_stats_betas.R2(:), '%0.3f');
            
            resultsInfo(end).files = {'temporal_rsq_info.txt'};
            icatb_printToFile(fullfile(resultsDir,  resultsInfo(end).files{1}), varStruct, '', 'column_wise');
        end
        
        %writeHTML(resultsDir, resultsInfo(end), 'spectral_summary.html', 'Spectral Summary');
        
    else
        
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
        
    end
    
    icatb_disp_temp_regress_results(temporal_stats_betas, 0);
    
    if (saveFigInfo)
        
        regress_figsH = findobj(0, 'tag', 'regress_betas');
        regress_betas_T = findobj(0, 'tag', 'regress_betas_T');
        regress_figsH(end + 1:end + length(regress_betas_T)) = regress_betas_T;
        if (~strcmpi(resultsFormat, 'pdf'))
            resultsInfo(end + 1).title = ' ';
            resultsInfo(end).text = ' ';
            regressPNGFiles = printFigs(regress_figsH, resultsDir, 'temporal_sort');
            resultsInfo(end).files = regressPNGFiles;
            %writeHTML(resultsDir, resultsInfo(end-1:end), 'temporal_sort.html', 'Temporal sort on beta weights');
            htmlSummaryStr(end + 1). title = 'Temporal sort on beta weights';
            htmlSummaryStr(end).tag = 'temporal_sort';
            htmlSummaryStr(end).str = get_result_strings(resultsDir, resultsInfo(end-1:end), htmlSummaryStr(end).tag);
            
            
        else
            
            for nPdfs = 1:length(regress_figsH)
                countPdfs = countPdfs + 1;
                tmpImFile = [pdfPrefix, '_', num2str(countPdfs), '.pdf'];
                set(regress_figsH(nPdfs), 'PaperPositionMode', 'auto');
                print(regress_figsH(nPdfs), fullfile(resultsDir, tmpImFile), printOptions{:});
                pdfFiles{countPdfs} = fullfile(resultsDir, tmpImFile);
            end
            
        end
        delete(regress_figsH);
        
    end
    
end


resultsInfo = [];

drawnow;

%% Kurtosis of timecourses and spatial maps
% Mean across subjects is reported in table. Figure shows mean+/- SEM across subjects
%if (isfield(matFileInfo, 'kurt_comp') && ~isempty(matFileInfo.kurt_comp))
if (isfield(aggregate, 'kurt'))
    
    %     size_kurt_tc = size(getfield(matFileInfo.kurt_comp, 'tc'));
    %
    %     if (numel(size_kurt_tc) == 4)
    %         % average across sessions
    %         tmp_tc = squeeze(mean(getfield(matFileInfo.kurt_comp, 'tc'), 2));
    %         tmp_ic = squeeze(mean(getfield(matFileInfo.kurt_comp, 'ic'), 2));
    %     else
    %         tmp_tc = getfield(matFileInfo.kurt_comp, 'tc');
    %         tmp_ic = getfield(matFileInfo.kurt_comp, 'ic');
    %     end
    %
    %     if (numOfSub > 1)
    %         values = [mean(tmp_tc); mean(tmp_ic)];
    %         errors = [std(tmp_tc)/sqrt(size(tmp_tc, 1)); std(tmp_ic)/sqrt(size(tmp_ic, 1))];
    %     else
    %         values = [tmp_tc(:)'; tmp_ic(:)'];
    %         errors = zeros(size(values));
    %     end
    
    values = [aggregate.kurt.tc.mean; aggregate.kurt.ic.mean];
    errors = [aggregate.kurt.tc.sem; aggregate.kurt.ic.sem];
    
    varNames = {'ComponentNumber', 'Timecourses', 'SpatialMaps'};
    
    if (saveFigInfo)
        
        if (~strcmpi(resultsFormat, 'pdf'))
            
            resultsInfo(end + 1).title = 'Kurtosis of timecourses and spatial maps';
            resultsInfo(end).text = 'Mean across subjects is reported in table. Figure shows mean+/- SEM across subjects';
            
            clear varStruct;
            
            numPara = 1;
            varStruct(numPara).tag = varNames{1};
            varStruct(numPara).value = components(:);
            
            numPara = numPara + 1;
            varStruct(numPara).tag = varNames{2};
            varStruct(numPara).value = num2str(values(1, :)', '%0.3f');
            
            numPara = numPara + 1;
            varStruct(numPara).tag = varNames{3};
            varStruct(numPara).value = num2str(values(2, :)', '%0.3f');
            
            resultsInfo(end).files = {'kurtosis_info.txt'};
            icatb_printToFile(fullfile(resultsDir,  resultsInfo(end).files{1}), varStruct, '', 'column_wise');
        end
        
    else
        
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
        
    end
    
    drawnow;
    tmp_values = values(1, :);
    tmp_errors = errors(1, :);
    bw_legend = cellstr([repmat('Comp ', length(components), 1), num2str(components(:))]);
    kurtH(1) = figure('color', 'w', 'name', 'Kurtosis (Timecourses)', 'position', figurePos, 'visible', figVisible);
    %plotBars(values, errors, 0.8, {'Timecourses', 'Spatial Maps'}, winter(64), bw_legend);
    plotBars(tmp_values, tmp_errors, 0.8, {'Timecourses'}, winter(64), bw_legend);
    title('Kurtosis (Mean +/- SEM)');
    set(gca, 'Ylim', [min(tmp_values(:) -abs(tmp_errors(:))) - 0.1, max(tmp_values(:) + abs(tmp_errors(:)) + 0.1)]);
    
    %pause(1);
    tmp_values = values(2, :);
    tmp_errors = errors(2, :);
    kurtH(2) = figure('color', 'w', 'name', 'Kurtosis (SM)', 'position', figurePos, 'visible', figVisible);
    plotBars(tmp_values, tmp_errors, 0.8, {'Spatial Maps'}, winter(64), bw_legend);
    title('Kurtosis (Mean +/- SEM)');
    set(gca, 'Ylim', [min(tmp_values(:) -abs(tmp_errors(:))) - 0.1, max(tmp_values(:) + abs(tmp_errors(:)) + 0.1)]);
    
    clear values errors;
    drawnow;
    
    if (saveFigInfo)
        
        if (~strcmpi(resultsFormat, 'pdf'))
            kurtosisPNGFiles = printFigs(kurtH, resultsDir, 'kurtosis');
            resultsInfo(end + 1).title = ' ';
            resultsInfo(end).text = ' ';
            resultsInfo(end).files = kurtosisPNGFiles;
            %resultsInfo(end).files(2:2+length(kurtosisPNGFiles)-1) = kurtosisPNGFiles;
            %writeHTML(resultsDir, resultsInfo(end-1:end), 'kurtosis_info.html', 'Kurtosis of timecourses and spatial maps');
            htmlSummaryStr(end + 1). title = 'Kurtosis of timecourses and spatial maps';
            htmlSummaryStr(end).tag = 'kurtosis';
            htmlSummaryStr(end).str = get_result_strings(resultsDir, resultsInfo(end-1:end), htmlSummaryStr(end).tag);
        else
            for nPdfs = 1:length(kurtH)
                countPdfs = countPdfs + 1;
                tmpImFile = [pdfPrefix, '_', num2str(countPdfs), '.pdf'];
                set(kurtH(nPdfs), 'PaperPositionMode', 'auto');
                print(kurtH(nPdfs), fullfile(resultsDir, tmpImFile), printOptions{:});
                pdfFiles{countPdfs} = fullfile(resultsDir, tmpImFile);
            end
        end
        delete(kurtH);
    end
    
end

drawnow;

%% FNC correlations
% Functional network connectivity correlations are computed for each
% data-set and averaged across sessions.
if (isfield(aggregate, 'fnc'))
    
    fnc_corrs_all = aggregate.fnc.mean;
    fnc_corrs_all = icatb_z_to_r(fnc_corrs_all);
    CLIM = max(abs(fnc_corrs_all(:)));
    CLIM = [-CLIM, CLIM];
    gH = figure('color', 'w', 'visible', figVisible);
    set(gH, 'resize', 'on');
    axesH = axes('parent', gH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
    icatb_plot_FNC(fnc_corrs_all, CLIM, cellstr(num2str(components(:))), (1:length(components(:))), gH, [], axesH);
    colormap(jet(64));
    title('Average FNC Correlations', 'parent', axesH);
    
    fncHandles(1) = gH;
    
    % for n = 1:length(grp_names)
    %     gH = figure('color', 'w', 'visible', figVisible);
    %     set(gH, 'resize', 'on');
    %     axesH = axes('parent', gH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
    %     icatb_plot_FNC(fnc_grps{n}, CLIMG, cellstr(num2str(components(:))), (1:length(components)), gH, [], axesH);
    %     colormap(jet(64));
    %     title(['Average FNC Correlations of ', grp_names{n}], 'parent', axesH);
    %     fncHandles(end + 1) = gH;
    % end
    
    
    resultsInfo = [];
    if (saveFigInfo)
        
        if (~strcmpi(resultsFormat, 'pdf'))
            resultsInfo(end + 1).title = 'FNC correlations on timecourses';
            resultsInfo(end).text = 'Functional network connectivity correlations are computed for each data-set using timecourses and averaged across sessions.';
            FNCTPNGFiles = printFigs(fncHandles, resultsDir, 'FNC_time');
            resultsInfo(end).files = FNCTPNGFiles;
            %writeHTML(resultsDir, resultsInfo(end), 'fnc_timecourses.html', 'FNC correlations on timecourses');
            htmlSummaryStr(end + 1).title = 'FNC correlations on timecourses';
            htmlSummaryStr(end).tag  = 'fnc_corrs_tc';
            htmlSummaryStr(end).str = get_result_strings(resultsDir, resultsInfo(end), htmlSummaryStr(end).tag);
            
        else
            
            for nPdfs = 1:length(fncHandles)
                countPdfs = countPdfs + 1;
                tmpImFile = [pdfPrefix, '_', num2str(countPdfs), '.pdf'];
                set(fncHandles(nPdfs), 'PaperPositionMode', 'auto');
                print(fncHandles(nPdfs), fullfile(resultsDir, tmpImFile), printOptions{:});
                pdfFiles{countPdfs} = fullfile(resultsDir, tmpImFile);
            end
            
        end
        
        delete(fncHandles);
    end
    
end

drawnow;

%% FNC metrics of component spatial maps
% Mutual information is computed between components spatially and averaged
% across data-sets.

if (isfield(aggregate, 'mi'))
    spatial_maps_MI =  aggregate.mi.mean;
    CLIM = [min(abs(spatial_maps_MI(:))), max(abs(spatial_maps_MI(:)))];
    %     fnc_grps = fnc_grps_mi;
    %     CLIMG = CLIMG_mi;
    %     CLIM = CLIM_mi;
    
    gH = figure('color', 'w', 'visible', figVisible);
    set(gH, 'resize', 'on');
    axesH = axes('parent', gH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
    icatb_plot_FNC(spatial_maps_MI, CLIM, cellstr(num2str(components(:))), (1:length(components)), gH, ' ', axesH);
    colormap(jet(64));
    title('Average FNC metrics (Spatial maps)', 'parent', axesH);
    
    clear fncHandles;
    fncHandles(1) = gH;
    
    %     for n = 1:length(grp_names)
    %         gH = figure('color', 'w', 'visible', figVisible);
    %         set(gH, 'resize', 'on');
    %         axesH = axes('parent', gH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
    %         icatb_plot_FNC(fnc_grps{n}, CLIMG, cellstr(num2str(components(:))), (1:length(components)), gH, ' ', axesH);
    %         colormap(jet(64));
    %         title(['Average spatial FNC metrics of ', grp_names{n}], 'parent', axesH);
    %         fncHandles(end + 1) = gH;
    %     end
    
    resultsInfo = [];
    if (saveFigInfo)
        
        if (~strcmpi(resultsFormat, 'pdf'))
            resultsInfo(end + 1).title = 'FNC metrics of component spatial maps';
            resultsInfo(end).text = 'Mutual information is computed between components spatially and averaged across data-sets.';
            FNCSPNGFiles = printFigs(fncHandles, resultsDir, 'FNC_maps');
            resultsInfo(end).files = FNCSPNGFiles;
            %writeHTML(resultsDir, resultsInfo(end), 'fnc_spatial_maps.html', 'FNC metrics of component spatial maps');
            htmlSummaryStr(end + 1). title = 'FNC metrics of component spatial maps';
            htmlSummaryStr(end).tag = 'fnc_corrs_sm';
            htmlSummaryStr(end).str = get_result_strings(resultsDir, resultsInfo(end), htmlSummaryStr(end).tag);
        else
            for nPdfs = 1:length(fncHandles)
                countPdfs = countPdfs + 1;
                tmpImFile = [pdfPrefix, '_', num2str(countPdfs), '.pdf'];
                set(fncHandles(nPdfs), 'PaperPositionMode', 'auto');
                print(fncHandles(nPdfs), fullfile(resultsDir, tmpImFile), printOptions{:});
                pdfFiles{countPdfs} = fullfile(resultsDir, tmpImFile);
            end
        end
        delete(fncHandles);
    end
end

if (saveFigInfo)
    if (strcmpi(resultsFormat, 'pdf'))
        append_pdfs2(fullfile(resultsDir, [pdfPrefix, '.pdf']), pdfFiles{:}); % new append_pdfs2 function
        for nF = 1:length(pdfFiles)
            try
                delete(pdfFiles{nF});
            catch
            end
        end
    else
        writeHTML2(fullfile(resultsDir, [mfilename, '.html']), htmlSummaryStr);
    end
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
structFile = fullfile (fileparts(which('gift.m')), 'icatb_templates', 'ch2bet_3x3x3.nii');

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


function newText = dispParams(sesInfo, subjects_to_use, all_components)


dasets_in_use = zeros(1, length(subjects_to_use)*sesInfo.numOfSess);
endInd = 0;
for tmpS = 1:length(subjects_to_use)
    ia = (subjects_to_use(tmpS) - 1)*sesInfo.numOfSess + 1;
    ib = subjects_to_use(tmpS)*sesInfo.numOfSess;
    startInd = endInd + 1;
    endInd = endInd + sesInfo.numOfSess;
    dasets_in_use(startInd:endInd) = (ia:ib);
end

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
    D(size(D,2)+1).string = ['Number of Subjects : ', num2str(length(subjects_to_use))];
    D(size(D,2)+1).string = ['Number of Sessions : ', num2str(sesInfo.numOfSess)];
    if(isfield(sesInfo.userInput, 'TR'))
        Trs = sesInfo.userInput.TR;
        if (length(Trs) == 1)
            Trs = repmat(Trs, 1, length(subjects_to_use));
        else
            Trs = Trs(subjects_to_use);
        end
        
        D(size(D,2)+1).string = ['TR in seconds (Mean, Standard deviation): ', '(', num2str(mean(Trs)), ', ', num2str(std(Trs)), ')'];
    end
end



D(size(D,2)+1).string = ['Number of Independent Components : ', num2str(length(all_components))];
D(size(D,2)+1).string = ['ICA Algorithm : ', selected_ica_algorithm];

if (~strcmpi(modalityType, 'smri'))
    D(size(D,2)+1).string = ['Number Of Scans/Timepoints (Mean, Standard deviation): ', '(', num2str(mean(sesInfo.diffTimePoints(dasets_in_use))), ', ', num2str(std(sesInfo.diffTimePoints(dasets_in_use))), ')'];
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

stability_opts = {'None', 'ICASSO', 'MST', 'Cross ISI'};

D(size(D,2)+1).string = ['Stability analysis type : ', stability_opts{whichAnalysis}];

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
    legend(bh(1: size(barvalues, 2)), bw_legend{:}, 'location', 'bestoutside');
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


function writeHTML(htmlDir, sdFNCResults, outFile, titleStr)

html_file = fullfile(htmlDir, outFile);
start_string = ['<html><head><title>', titleStr, '</title>', ...
    '<meta http-equiv="Content-Type" content="HTML, DHTML, CSS, XML, XHTML, JavaScript, VBScript; charset=iso-8859-1">', ...
    '<link rel="stylesheet" type="text/css" href="./style.css"></head>'];

results_string1 = [];
for nR = 1:length(sdFNCResults)
    results_string1 = [results_string1, '<h2 align = "center">', sdFNCResults(nR).title, '</h2>'];
    results_string1 = [results_string1, '<hr>'];
    results_string1 = [results_string1, '<p align = "center">', sdFNCResults(nR).text, '</p>'];
    files =  sdFNCResults(nR).files;
    
    [pp,ff,extn]= fileparts(files{1});
    
    if (~strcmpi(extn, '.txt'))
        
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
        
    else
        
        for nrows = 1:length(files)
            allLines = icatb_textscan(fullfile(htmlDir, files{nrows}));
            allLines = allLines(icatb_good_cells(allLines));
            for nLine = 1:length(allLines)
                
                allCols = textscan(allLines{nLine}, '%s', 'whitespace', '\t');
                allCols = allCols{1};
                
                if (nLine == 1)
                    tableStr = '<table align = "center" border = "1" width = "35%">';
                    tDStart = '<th align = "center">';
                    tDEnd = '</th>';
                else
                    tDStart = '<td align = "center">';
                    tDEnd = '</td>';
                end
                
                tableStr = [tableStr, '<tr>'];
                
                for nCols = 1:length(allCols)
                    tableStr = [tableStr, tDStart, allCols{nCols}, tDEnd];
                end
                
                tableStr = [tableStr, '</tr>'];
                
            end
            results_string1 = [tableStr, '</table>'];
        end
    end
    
    results_string1 = [results_string1, '<p>&nbsp;&nbsp;</p>'];
end

end_string =  '</html>';


results_string = [start_string,  results_string1, end_string];

dlmwrite(html_file, results_string, '');



function results_string = get_result_strings(htmlDir, sdFNCResults, tagStr)

%html_file = fullfile(htmlDir, outFile);
% start_string = ['<html><head><title>', titleStr, '</title>', ...
%     '<meta http-equiv="Content-Type" content="HTML, DHTML, CSS, XML, XHTML, JavaScript, VBScript; charset=iso-8859-1">', ...
%     '<link rel="stylesheet" type="text/css" href="./style.css"></head>'];

results_string1 = [];
for nR = 1:length(sdFNCResults)
    results_string1 = [results_string1, '<h2 align = "center"><a name = "', tagStr, '">', sdFNCResults(nR).title, '</a></h2>'];
    %results_string1 = [results_string1, '<hr>'];
    results_string1 = [results_string1, '<p align = "center">', sdFNCResults(nR).text, '</p>'];
    files =  sdFNCResults(nR).files;
    
    [pp,ff,extn]= fileparts(files{1});
    
    if (~strcmpi(extn, '.txt'))
        
        num_cols = 1;
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
        
    else
        
        for nrows = 1:length(files)
            allLines = icatb_textscan(fullfile(htmlDir, files{nrows}));
            allLines = allLines(icatb_good_cells(allLines));
            for nLine = 1:length(allLines)
                
                allCols = textscan(allLines{nLine}, '%s', 'whitespace', '\t');
                allCols = allCols{1};
                
                if (nLine == 1)
                    tableStr = '<table align = "center" border = "1" width = "35%">';
                    tDStart = '<th align = "center">';
                    tDEnd = '</th>';
                else
                    tDStart = '<td align = "center">';
                    tDEnd = '</td>';
                end
                
                tableStr = [tableStr, '<tr>'];
                
                for nCols = 1:length(allCols)
                    tableStr = [tableStr, tDStart, allCols{nCols}, tDEnd];
                end
                
                tableStr = [tableStr, '</tr>'];
                
            end
            results_string1 = [results_string1, tableStr, '</table>'];
        end
    end
    
    results_string1 = [results_string1, '<p>&nbsp;&nbsp;</p>'];
end

end_string =  '</html>';


results_string = ['<p> </p> <hr>', results_string1, ' <p> </p>'];

%dlmwrite(html_file, results_string, '');



function writeHTML2(fileName, resultsStr)


start_string = '<html><head><title> Group ICA Results </title></head>';
start_string = [start_string, '<p> </p><h2> Contents </h2><p><ul>'];


titles = cellstr(char(resultsStr.title));

for nR = 1:length(titles)
    start_string = [start_string, '<li><h3><a href="#', resultsStr(nR).tag, '">', titles{nR}, '</a></h3></li>'];
end

start_string = [start_string, '</ul></p><p> </p>'];
end_string =  '</html>';
results{1} = start_string;
results{2} = char(resultsStr.str);
results{3} = end_string;


dlmwrite(fileName, char(results), '');




function imNames = printFigs(Figs, outdir, filename, append_prefix)

icatb_defaults;
global GICA_RESULTS_SUMMARY;

if (~exist('append_prefix', 'var'))
    append_prefix = 'yes';
end

PRINT_RESOLUTION = '';
try
    PRINT_RESOLUTION = GICA_RESULTS_SUMMARY.print_resolution;
catch
end

if (isempty(PRINT_RESOLUTION))
    PRINT_RESOLUTION = '-r72';
end

try
    print_opts = GICA_RESULTS_SUMMARY.print_opts;
catch
    print_opts = {'-noui'};
end


printOptions = [{'-dpng', PRINT_RESOLUTION}, print_opts];

printOptions(strcmpi(printOptions, '-bestfit')) = [];
printOptions(strcmpi(printOptions, '-fillpage')) = [];

imNames = cell(1, length(Figs));
for nF = 1:length(Figs)
    if (strcmpi(append_prefix, 'yes'))
        imNames{nF} = [filename, '_', icatb_returnFileIndex(nF), '.png'];
    else
        imNames{nF} = [filename, '.png'];
    end
    print(Figs(nF), fullfile(outdir, imNames{nF}), printOptions{:});
end



function writeComponentMaps (sesInfo, compFileNaming, subjects_to_use, components, resultsDir)


tp = [];
dasets_in_use = zeros(1, length(subjects_to_use)*sesInfo.numOfSess);
endInd = 0;
for tmpS = 1:length(subjects_to_use)
    ia = (subjects_to_use(tmpS) - 1)*sesInfo.numOfSess + 1;
    ib = subjects_to_use(tmpS)*sesInfo.numOfSess;
    tp = min([tp, min(sesInfo.diffTimePoints([ia, ib]))]);
    startInd = endInd + 1;
    endInd = endInd + sesInfo.numOfSess;
    dasets_in_use(startInd:endInd) = (ia:ib);
end

meanIm = zeros(length(components), length(sesInfo.mask_ind));
meanTC = zeros(tp, length(components));


sesInfo.userInput.files = [];
sesInfo.inputFiles = [];
try
    sesInfo.userInput.ICA_Options = {};
    sesInfo.ICA_Options = {};
catch
end

numOfSess = sesInfo.numOfSess;
parfor nDataset = 1:length(dasets_in_use)
    nD = dasets_in_use(nDataset);
    nSub = ceil(nD/numOfSess);
    nSess = mod(nD - 1, numOfSess) + 1;
    [tc, ic]  = icatb_loadComp(sesInfo, components, 'subjects', nSub, 'sessions', nSess);
    meanIm = meanIm + ic';
    meanTC =  meanTC + tc(1:tp, :);
end


meanIm = meanIm / length(dasets_in_use);
meanTC = meanTC / length(dasets_in_use);

try
    if exist(fullfile(resultsDir, compFileNaming), 'file')
        delete (fullfile(resultsDir, compFileNaming));
    end
catch
end

icatb_saveICAData(compFileNaming, meanIm, meanTC, sesInfo.mask_ind, length(components), sesInfo.HInfo, 'real', [], resultsDir);



function closeGCP

if (~isempty(which('parpool')))
    try
        poolobj = gcp('nocreate');
        delete(poolobj);
    catch
    end
else
    try
        matlabpool close;
    catch
    end
end


function aggregate = getAggInfo(sesInfo, postProcessFile)
% Get aggregate info

icatb_defaults
global TIMECOURSE_POSTPROCESS;

freq_limits = [0.1, 0.15];
try
    freq_limits = TIMECOURSE_POSTPROCESS.spectra.freq_limits;
catch
end

matFileInfo = matfile(postProcessFile);
freq = matFileInfo.freq;
spectra_size = size(matFileInfo, 'spectra_tc_all');

try
    subjects = matFileInfo.subjects;
catch
end

try
    components = matFileInfo.components;
catch
end

if (~exist('subjects', 'var'))
    subjects = (1:sesInfo.numOfSub);
end

subjects (subjects > max(sesInfo.numOfSub)) = [];

if (isempty(subjects))
    error('Please check the subjects variable passed');
end


if (~exist('components', 'var'))
    components = (1:sesInfo.numComp);
end

components (components > max(sesInfo.numComp)) = [];

if (isempty(components))
    error('Please check the components variable passed');
end

numComp = length(components);
numOfSub = length(subjects);


if (~isempty(who(matFileInfo, 'fnc_corrs_all')))
    disp('Computing mean FNC of timecourses ...');
    aggregate.fnc.mean = computeFNC(postProcessFile, 'fnc_corrs_all', numComp, numOfSub);
    disp('Done computing mean of FNC');
end

if (~isempty(who(matFileInfo, 'spatial_maps_MI')))
    disp('Computing mean FNC of spatial maps ...');
    aggregate.mi.mean = computeFNC(postProcessFile, 'spatial_maps_MI', numComp, numOfSub);
    disp('Done computing mean FNC of spatial maps');
end

closeGCP;

dynamic_range = zeros(1, numComp);
fALFF = zeros(1, numComp);
spectra_mean = zeros(length(freq), numComp);
spectra_sem = zeros(length(freq), numComp);

for nF = 1:numComp
    
    if (numel(spectra_size) == 4)
        tc = reshape(squeeze(mean(matFileInfo.spectra_tc_all(:, :, :, nF), 2)), numOfSub, length(freq));
    else
        if (numel(spectra_size) == 3)
            tc = matFileInfo.spectra_tc_all(:, :, nF);
        else
            tc = matFileInfo.spectra_tc_all(:, nF);
            tc = tc';
        end
    end
    
    numOfSub = size(tc, 1);
    
    tc = reshape(tc, numOfSub, length(freq));
    
    tmp_dynamicrange = zeros(1, size(tc, 1));
    tmp_fALFF = tmp_dynamicrange;
    for nS = 1:length(tmp_dynamicrange)
        [tmp_dynamicrange(nS), tmp_fALFF(nS)] = icatb_get_spec_stats(tc(nS, :), freq, freq_limits);
    end
    
    dynamic_range(nF) = mean(tmp_dynamicrange);
    fALFF(nF) = mean(tmp_fALFF);
    spectra_mean(:, nF) = mean(tc, 1);
    spectra_sem(:, nF) = (std(tc, 0, 1) + eps) / sqrt(size(tc, 1));
    
    clear tc
    
    
end

aggregate.spectra.dynamic_range = dynamic_range;
aggregate.spectra.fALFF = fALFF;
aggregate.spectra.mean = spectra_mean;
aggregate.spectra.sem = spectra_sem;
aggregate.spectra.freq = freq;

clear dynamic_range fALFF spectra_mean spectra_sem;

if (~isempty(who(matFileInfo, 'kurt_comp')) && ~isempty(matFileInfo.kurt_comp))
    size_kurt_tc = size(getfield(matFileInfo.kurt_comp, 'tc'));
    if (numel(size_kurt_tc) == 4)
        % average across sessions
        tmp_tc = squeeze(mean(getfield(matFileInfo.kurt_comp, 'tc'), 2));
        tmp_ic = squeeze(mean(getfield(matFileInfo.kurt_comp, 'ic'), 2));
    else
        tmp_tc = getfield(matFileInfo.kurt_comp, 'tc');
        tmp_ic = getfield(matFileInfo.kurt_comp, 'ic');
    end
    
    if (numOfSub > 1)
        values = [mean(tmp_tc); mean(tmp_ic)];
        errors = [std(tmp_tc)/sqrt(size(tmp_tc, 1)); std(tmp_ic)/sqrt(size(tmp_ic, 1))];
    else
        values = [tmp_tc(:)'; tmp_ic(:)'];
        errors = zeros(size(values));
    end
    aggregate.kurt.tc.mean = values(1, :);
    aggregate.kurt.tc.sem = errors(1, :);
    aggregate.kurt.ic.mean = values(2, :);
    aggregate.kurt.ic.sem = errors(2, :);
end

function fnc_corrs_all = computeFNC(postProcessFile, varName, numComp, numOfSub)
% Compute FNC

matFileInfo = matfile (postProcessFile);
if (numOfSub > 1)
    fnc_corrs_all = zeros(numComp, numComp);
    parfor nSubIn = 1:numOfSub
        matFileInfo = matfile (postProcessFile);
        tmp = matFileInfo.(varName)(nSubIn, :, :, :);
        if (size(tmp, 2) > 1)
            tmp = reshape(mean(tmp, 2), numComp, numComp);
        else
            tmp = reshape(squeeze(tmp), numComp, numComp);
        end
        fnc_corrs_all = fnc_corrs_all + tmp;
    end
    fnc_corrs_all = fnc_corrs_all ./ numOfSub;
    
else
    fnc_corrs_all = squeeze(matFileInfo.fnc_corrs_all);
end
