function icatb_sbm_html_report(param_file, opts)

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
global UI_FS;
global GICA_RESULTS_SUMMARY;

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

if (~strcmpi(sesInfo.modality, 'smri'))
    error('You need to use sbm toolbox to display results');
end


sesInfo.outputDir = outputDir;


if (~strcmpi(sesInfo.modality, 'smri'))
    error('You need to use gift toolbox to display results');
end


isGIFTI = 0;
try
    isGIFTI = isfield(sesInfo.userInput.HInfo.V, 'gifti');
catch
end

if (isGIFTI)
    disp('GIFTI format is not supported');
    return;
end


if (~isGIFTI)
    structFile = fullfile(fileparts(which('groupica.m')), 'icatb_templates', 'ch2bet_3x3x3.nii');
    try
        structFile = opts.anatomical_file;
    catch
    end
else
    structFile = '';
    %structFile = sesInfo.HInfo.V.gifti;
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


resultsFormat = 'html';
try
    resultsFormat = opts.format;
catch
end

resultsDir =  fullfile(outputDir, [sesInfo.userInput.prefix, '_sbm_results']);
try
    resultsDir = opts.outputDir;
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

%% SBM ICA Parameters
sesInfo.structFile = structFile;
sesInfo.slice_plane = slice_plane;
sesInfo.image_values = image_values;
sesInfo.threshold = threshold;
sesInfo.groupsInfo = groupsInfo;
sesInfo.convert_to_zscores = convert_to_zscores;
newText = dispParams(sesInfo);
disp(newText);



pdfFiles = {};

htmlSummaryStr = [];

if (saveFigInfo)
    if (~strcmpi(resultsFormat, 'pdf'))
        icatb_print_table(cellstr(newText), fullfile(resultsDir, 'params.txt'));
        resultsInfo(end + 1).title = 'SBM ICA Parameters';
        resultsInfo(end).text = ' ';
        resultsInfo(end).files = {'params.txt'};
        htmlSummaryStr(end + 1).title = 'SBM ICA Parameters';
        htmlSummaryStr(end).tag =  'sbmica_params';
        htmlSummaryStr(end).str = get_result_strings(resultsDir, resultsInfo(end), htmlSummaryStr(end).tag);
        %writeHTML(resultsDir, resultsInfo(end), 'introduction.html', 'Group ICA Parameters');
    else
        print_strs = cellstr(char('SBM ICA parameters:', '', newText));
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
        print(gH, '-dpdf', printRes, '-noui', fullfile(resultsDir, tmpImFile));
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
                            print(icassoFigs(nPdfs), '-dpdf', printRes, '-noui', fullfile(resultsDir, tmpImFile));
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

resultsInfo = [];
outFiles = {''};

%% SBM components
% Group spatial maps and subject loading coefficients
%
% * *a) Subject loading coefficients*
% * *b) Group spatial maps*

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
icaTimecourse = detrend(icaTimecourse, 0);
if (strcmpi(convert_to_zscores, 'yes'))
    icaTimecourse = icaTimecourse*diag(1./std(icaTimecourse));
end

kurt_tc = kurt(icaTimecourse, size(icaTimecourse, 2));


fprintf('\n');
ic = icatb_loadComp(sesInfo, (1:sesInfo.numComp), 'subjects', 1, 'sessions', 1, 'vars_to_load', 'ic');
kurt_ic = kurt(ic);
spatial_maps_MI = icatb_compute_mi(ic');

components = (1:size(compFiles, 1))';

returnValue = strmatch(lower(image_values), {'positive and negative', 'positive', 'absolute value', 'negative'}, 'exact');
cmap = icatb_getColormap(1, returnValue, 0);
if (strcmpi(convert_to_zscores, 'yes'))
    cbar_title='Z-scores';
else
    scale_opts = icatb_scaleICA;
    if ~ischar(sesInfo.scaleType)
        cbar_title = deblank(scale_opts(sesInfo.scaleType+1, :));
    else
        cbar_title = sesInfo.scaleType;
    end
end

meanH = [];
resultsInfo = [];
tmp_fnames = cell(1, size(compFiles, 1));
for nF = 1:size(compFiles, 1)
    
    cn = deblank(compFiles(nF, :));
    
    [dd, fN, extn] = fileparts(cn);
    
    if (~isGIFTI)
        gH = icatb_getGraphics([fN, extn], 'graphics', 'imviewer', 'on');
        set(gH, 'resize', 'on');
        sh = subplot(4, 1, 1);
    else
        gH = icatb_getGraphics([fN, extn], 'timecourse', 'imviewer', 'on');
        sh = gca;
    end
    
    tc =  icaTimecourse(:, nF);
    tc_label = 'Data units';
    if (strcmpi(convert_to_zscores, 'yes'))
        tc_label = 'Z-scores';
    end
    icatb_plotTimecourse('data', tc, 'parent', sh, 'color', 'm', 'title', ['Component ', icatb_returnFileIndex(nF)], 'xlabelstr', 'Subjects', 'ylabelstr', tc_label);
    
    if (~isGIFTI)
        sh = subplot(4, 1, 2:4);
        icatb_image_viewer(cn, 'structfile', structFile, 'image_values', image_values, 'convert_to_zscores', convert_to_zscores, 'threshold', threshold, ...
            'slices_in_mm', slices_in_mm, 'anatomical_view', slice_plane, 'axesh', sh, 'colorbar', 1, 'labels', ' ');
    else
        gC = gifti(cn);
        tmp = icatb_applyDispParameters_comp(gC.cdata(:)', strcmpi(convert_to_zscores, 'yes'), returnValue, threshold);
        colormap(cmap);
        gC.cdata = reshape(tmp, size(gC.cdata));
        cn2=[cn(1:end-4),'_tmp','.gii'];
        save(gC, cn2);
        icatb_display_gifti(cn2, 'colorbar_title', cbar_title, 'title', ['Component ', icatb_returnFileIndex(nF)]);
        
        delete(cn2);
        
    end
    
    if (saveFigInfo)
        if (~strcmpi(resultsFormat, 'pdf'))
            
            if (nF == size(compFiles, 1))
                resultsInfo(end + 1).title = 'SBM Components';
                resultsInfo(end).text = ['<ul> <li> <b> a) Subject loading coefficients </b> </li>', ...
                    '<li> <b> b) Group spatial maps </b> </li> </ul>'];
            end
            
            tmp_fname = printFigs(gH, resultsDir, ['group_comp', '_', icatb_returnFileIndex(nF)], 'no');
            tmp_fnames{nF} = tmp_fname{1};
            
            %resultsInfo(end).files = printFigs(meanH, resultsDir, 'group_comp');
            %writeHTML(resultsDir, resultsInfo(end), 'mean_comp.html', 'Mean Components');
            if (nF == size(compFiles, 1))
                resultsInfo(end).files = tmp_fnames;
                htmlSummaryStr(end + 1). title = 'Group Components';
                htmlSummaryStr(end).tag = 'group_comp';
                htmlSummaryStr(end).str = get_result_strings(resultsDir, resultsInfo(end), htmlSummaryStr(end).tag);
            end
            
        else
            
            %for nPdfs = 1:length(meanH)
            countPdfs = countPdfs + 1;
            tmpImFile = [pdfPrefix, '_', num2str(countPdfs), '.pdf'];
            set(gH, 'PaperPositionMode', 'auto');
            print(gH, '-dpdf', printRes, '-noui', fullfile(resultsDir, tmpImFile));
            pdfFiles{countPdfs} = fullfile(resultsDir, tmpImFile);
            %end
        end
        
        delete(gH);
    end
    
    
end

%resultsInfo = [];
%figurePos = get(gH, 'position');

% if (saveFigInfo)
%     if (~strcmpi(resultsFormat, 'pdf'))
%         resultsInfo(end + 1).title = 'SBM Components';
%         resultsInfo(end).text = ['<ul> <li> <b> a) Subject loading coefficients </b> </li>', ...
%             '<li> <b> b) Group spatial maps </b> </li> </ul>'];
%
%         resultsInfo(end).files = printFigs(meanH, resultsDir, 'group_comp');
%         %writeHTML(resultsDir, resultsInfo(end), 'mean_comp.html', 'Mean Components');
%         htmlSummaryStr(end + 1). title = 'Group Components';
%         htmlSummaryStr(end).tag = 'group_comp';
%         htmlSummaryStr(end).str = get_result_strings(resultsDir, resultsInfo(end), htmlSummaryStr(end).tag);
%
%     else
%
%         for nPdfs = 1:length(meanH)
%             countPdfs = countPdfs + 1;
%             tmpImFile = [pdfPrefix, '_', num2str(countPdfs), '.pdf'];
%             set(meanH(nPdfs), 'PaperPositionMode', 'auto');
%             print(meanH(nPdfs), '-dpdf', printRes, '-noui', fullfile(resultsDir, tmpImFile));
%             pdfFiles{countPdfs} = fullfile(resultsDir, tmpImFile);
%         end
%     end
%
%     delete(meanH);
% end

outFiles = {''};
resultsInfo = [];

%% PCA components variance summary
% * *Percent variance* - Percent variance is reported for each component.
%

if (isfield(sesInfo, 'pca_variances'))
    pca_variances = sesInfo.pca_variances;
    pca_variances = pca_variances(:);
    varNames = {'ComponentNumber', 'PercentVariance'};
    
    
    if (saveFigInfo)
        if (~strcmpi(resultsFormat, 'pdf'))
            
            resultsInfo(end + 1).title = 'PCA Percent Variance Summary';
            resultsInfo(end).text = '<ul> <li> Percent variance - Percent variance is reported for each PCA component </li> </ul>';
            
            numPara = 1;
            varStruct(numPara).tag = varNames{1};
            varStruct(numPara).value = components(:);
            
            numPara = numPara + 1;
            varStruct(numPara).tag = varNames{2};
            varStruct(numPara).value = num2str(pca_variances(:), '%0.3f');
            
            resultsInfo(end).files = {'pca_summary.txt'};
            icatb_printToFile(fullfile(resultsDir,  resultsInfo(end).files{1}), varStruct, '', 'column_wise');
            
            htmlSummaryStr(end + 1). title = 'PCA Percent Variance Summary';
            htmlSummaryStr(end).tag = 'pca_summary';
            htmlSummaryStr(end).str = get_result_strings(resultsDir, resultsInfo(end), htmlSummaryStr(end).tag);
            
        end
        
    else
        
        try
            T = table(components, pca_variances, 'VariableNames', varNames);
            disp(T);
        catch
            fprintf('%20s\t%20s\n', varNames{:});
            fprintf('\n');
            for n = 1:length(components)
                fprintf('%20s\t%20s\n', num2str(components(n), '%d'), num2str(pca_variances(n), '%0.3f'));
            end
            fprintf('\n');
        end
        
    end
end

outFiles = {''};
resultsInfo = [];

%% ICA components variance summary
% * *Percent variance* - Percent variance is reported for each component.
%
if (isfield(sesInfo, 'ica_variances'))
    ica_variances = sesInfo.ica_variances;
    ica_variances = ica_variances(:);
    varNames = {'ComponentNumber', 'PercentVariance'};
    
    
    if (saveFigInfo)
        if (~strcmpi(resultsFormat, 'pdf'))
            
            resultsInfo(end + 1).title = 'ICA Percent Variance Summary';
            resultsInfo(end).text = '<ul> <li> Percent variance - Percent variance is reported for each ICA component </li> </ul>';
            
            numPara = 1;
            varStruct(numPara).tag = varNames{1};
            varStruct(numPara).value = components(:);
            
            numPara = numPara + 1;
            varStruct(numPara).tag = varNames{2};
            varStruct(numPara).value = num2str(ica_variances(:), '%0.3f');
            
            resultsInfo(end).files = {'ica_summary.txt'};
            icatb_printToFile(fullfile(resultsDir,  resultsInfo(end).files{1}), varStruct, '', 'column_wise');
            
            htmlSummaryStr(end + 1). title = 'ICA Percent Variance Summary';
            htmlSummaryStr(end).tag = 'ica_summary';
            htmlSummaryStr(end).str = get_result_strings(resultsDir, resultsInfo(end), htmlSummaryStr(end).tag);
            
        end
        
    else
        
        try
            T = table(components, ica_variances, 'VariableNames', varNames);
            disp(T);
        catch
            fprintf('%20s\t%20s\n', varNames{:});
            fprintf('\n');
            for n = 1:length(components)
                fprintf('%20s\t%20s\n', num2str(components(n), '%d'), num2str(ica_variances(n), '%0.3f'));
            end
            fprintf('\n');
        end
        
    end
end


outFiles = {''};
resultsInfo = [];

%% Kurtosis of loading coefficients and spatial maps
%

varNames = {'ComponentNumber', 'Timecourses', 'SpatialMaps'};

if (saveFigInfo)
    
    if (~strcmpi(resultsFormat, 'pdf'))
        
        resultsInfo(end + 1).title = 'Kurtosis of timecourses and spatial maps';
        resultsInfo(end).text = 'Kurtosis is reported for each component.';
        
        clear varStruct;
        
        numPara = 1;
        varStruct(numPara).tag = varNames{1};
        varStruct(numPara).value = components(:);
        
        numPara = numPara + 1;
        varStruct(numPara).tag = varNames{2};
        varStruct(numPara).value = num2str(kurt_tc(:), '%0.3f');
        
        numPara = numPara + 1;
        varStruct(numPara).tag = varNames{3};
        varStruct(numPara).value = num2str(kurt_ic(:), '%0.3f');
        
        resultsInfo(end).files = {'kurtosis_info.txt'};
        icatb_printToFile(fullfile(resultsDir,  resultsInfo(end).files{1}), varStruct, '', 'column_wise');
    end
    
else
    
    try
        T = table(components, kurt_tc(:),  kurt_ic(:), 'VariableNames', varNames);
        disp(T);
    catch
        fprintf('%20s\t%20s\t%20s\n', varNames{:});
        fprintf('\n');
        for n = 1:length(components)
            fprintf('%20s\t%20s\t%20s\n', num2str(components(n), '%d'), num2str(kurt_tc(n), '%0.3f'), num2str(kurt_ic(n), '%0.3f'));
        end
        fprintf('\n');
    end
    
end

gH = figure('color', 'w');
set(gH, 'resize', 'on');
sh = axes('parent', gH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
bar(kurt_tc, 'parent', sh);
xlabel('Components', 'parent', sh);
ylabel('Kurtosis', 'parent', sh);
title('Loading coefficents', 'parent', sh);
axis(sh, 'tight');

kurtH(1) = gH;

gH = figure('color', 'w');
set(gH, 'resize', 'on');
sh = axes('parent', gH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
bar(kurt_ic, 'parent', sh);
xlabel('Components', 'parent', sh);
ylabel('Kurtosis', 'parent', sh);
title('Spatial maps', 'parent', sh);
axis(sh, 'tight');
kurtH(2) = gH;

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
            print(kurtH(nPdfs), '-dpdf', printRes, '-noui', fullfile(resultsDir, tmpImFile));
            pdfFiles{countPdfs} = fullfile(resultsDir, tmpImFile);
        end
    end
    delete(kurtH);
end


outFiles = {''};
resultsInfo = [];

%% FNC correlations of subject loadings
% Functional network connectivity correlations
fnc_corrs_all = icatb_corr(icaTimecourse);
fnc_corrs_all(1:size(fnc_corrs_all, 1) + 1:end) = 0;
fnc_grps = cell(1, length(groupsInfo));
grp_names = cell(1, length(groupsInfo));
CLIMG = [];


for n = 1:length(groupsInfo)
    grp_names{n} = groupsInfo(n).name;
    tmp = icatb_corr(icaTimecourse(groupsInfo(n).val, :));
    tmp(1:size(tmp, 1) + 1:end) = 0;
    fnc_grps{n} = tmp;
    CLIMG = max([CLIMG, max(abs(tmp(:)))]);
end

fncHandles = [];
comps = (1:sesInfo.numComp)';
CLIM = max(abs(fnc_corrs_all(:)));
gH = figure('color', 'w');
set(gH, 'resize', 'on');
axesH = axes('parent', gH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
icatb_plot_FNC(fnc_corrs_all, [-CLIM, CLIM], cellstr(num2str(comps)), (1:length(comps)), gH, [], axesH);
colormap(jet(64));
title('FNC Correlations', 'parent', axesH);
fncHandles(end + 1) = gH;

for n = 1:length(grp_names)
    gH = figure('color', 'w');
    set(gH, 'resize', 'on');
    axesH = axes('parent', gH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
    icatb_plot_FNC(fnc_grps{n}, [-CLIMG, CLIMG], cellstr(num2str(comps)), (1:length(comps)), gH, [], axesH);
    colormap(jet(64));
    title(['FNC Correlations of ', grp_names{n}], 'parent', axesH);
    fncHandles(end + 1) = gH;
end


%fncHandles(1) = gH;
resultsInfo = [];
if (saveFigInfo)
    
    if (~strcmpi(resultsFormat, 'pdf'))
        resultsInfo(end + 1).title = 'FNC correlations on loadings';
        resultsInfo(end).text = 'Functional network connectivity correlations on subject loading coefficients.';
        FNCTPNGFiles = printFigs(fncHandles, resultsDir, 'FNC_loadings');
        resultsInfo(end).files = FNCTPNGFiles;
        %writeHTML(resultsDir, resultsInfo(end), 'fnc_timecourses.html', 'FNC correlations on timecourses');
        htmlSummaryStr(end + 1).title = 'FNC correlations on subject loadings';
        htmlSummaryStr(end).tag  = 'fnc_corrs_loadings';
        htmlSummaryStr(end).str = get_result_strings(resultsDir, resultsInfo(end), htmlSummaryStr(end).tag);
        
    else
        
        for nPdfs = 1:length(fncHandles)
            countPdfs = countPdfs + 1;
            tmpImFile = [pdfPrefix, '_', num2str(countPdfs), '.pdf'];
            set(fncHandles(nPdfs), 'PaperPositionMode', 'auto');
            print(fncHandles(nPdfs), '-dpdf', printRes, '-noui', fullfile(resultsDir, tmpImFile));
            pdfFiles{countPdfs} = fullfile(resultsDir, tmpImFile);
        end
        
    end
    
    delete(fncHandles);
end



%% FNC metrics of component spatial maps
% Mutual information is computed between components spatially
comps = (1:sesInfo.numComp)';
CLIM = [min(abs(spatial_maps_MI(:))), max(abs(spatial_maps_MI(:)))];
gH = figure('color', 'w');
set(gH, 'resize', 'on');
axesH = axes('parent', gH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
icatb_plot_FNC(spatial_maps_MI, CLIM, cellstr(num2str(comps)), (1:length(comps)), gH, ' ', axesH);
colormap(jet(64));
title('FNC metrics (Spatial maps)', 'parent', axesH);

resultsInfo = [];
fncHandles = gH;
if (saveFigInfo)
    
    if (~strcmpi(resultsFormat, 'pdf'))
        resultsInfo(end + 1).title = 'FNC metrics of component spatial maps';
        resultsInfo(end).text = 'Mutual information is computed between components spatially.';
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
            print(fncHandles(nPdfs), '-dpdf', printRes, '-noui', fullfile(resultsDir, tmpImFile));
            pdfFiles{countPdfs} = fullfile(resultsDir, tmpImFile);
        end
    end
    delete(fncHandles);
end


if (saveFigInfo)
    if (strcmpi(resultsFormat, 'pdf'))
        append_pdfs(fullfile(resultsDir, [pdfPrefix, '.pdf']), pdfFiles{:});
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


if (strcmpi(modalityType, 'fmri'))
    preproc_options = icatb_preproc_data;
    preproc_options = cellstr(preproc_options);
    preproc_type = 'remove mean per timepoint';
    if (isfield(sesInfo, 'preproc_type'))
        preproc_type = sesInfo.preproc_type;
    end
    
    preprocInd = strmatch(lower(preproc_type), lower(preproc_options), 'exact');
    
    D(size(D,2)+1).string = ['Data Pre-processing Type : ', preproc_options{preprocInd}];
end

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

if (strcmpi(modalityType, 'fmri'))
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

if (~strcmpi(modalityType, 'smri'))
    checkAnalysisMode = 'Serial';
    try
        checkAnalysisMode = sesInfo.parallel_info.mode;
    catch
    end
    D(size(D,2)+1).string = ['Group analysis mode: ', upper(checkAnalysisMode(1)), checkAnalysisMode(2:end)];
end

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


function k = kurt(x, numComp)
% Kurtosis

if (size(x, 1) == 1)
    if (numel(x) == length(x))
        x = repmat(x, numComp, 1);
    end
end

x = icatb_remove_mean(x);
s2 = mean(x.^2);
m4 = mean(x.^4);
k = (m4 ./ (s2.^2));


function writeHTML2(fileName, resultsStr)


start_string = '<html><head><title> SBM ICA Results </title></head>';
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
            results_string1 = [results_string1, tableStr, '</table>'];
        end
    end
    
    results_string1 = [results_string1, '<p>&nbsp;&nbsp;</p>'];
end

end_string =  '</html>';


results_string = ['<p> </p> <hr>', results_string1, ' <p> </p>'];

%dlmwrite(html_file, results_string, '');


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

imNames = cell(1, length(Figs));
for nF = 1:length(Figs)
    if (strcmpi(append_prefix, 'yes'))
        imNames{nF} = [filename, '_', icatb_returnFileIndex(nF), '.png'];
    else
        imNames{nF} = [filename, '.png'];
    end
    print(Figs(nF), '-dpng', PRINT_RESOLUTION, '-noui', fullfile(outdir, imNames{nF}));
end