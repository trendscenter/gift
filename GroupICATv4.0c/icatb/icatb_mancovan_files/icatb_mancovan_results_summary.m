function icatb_mancovan_results_summary(mancovanInfo)

try
    resultsDir = mancovanInfo.display.dirName;
catch
    resultsDir = [mancovanInfo.prefix, '_results_summary'];
end


try
    mancovanInfo.userInput.compFiles = mancovanInfo.display.compFiles;
catch
end

display_steps = {'features', 'mult', 'uni'};

try
    display_steps = mancovanInfo.display.steps;
catch
end

t_threshold = 1.5;
try
    t_threshold = mancovanInfo.display.t_threshold;
catch
end

mancovanInfo.display.t_threshold = t_threshold;

resultsDir = fullfile(mancovanInfo.outputDir, resultsDir);

%resultsDir = fullfile(mancovanInfo.outputDir, [mancovanInfo.prefix, '_results_summary']);
htmlFile = fullfile(resultsDir, 'icatb_mancovan_results_summary.html');

if (exist(resultsDir, 'dir') ~= 7)
    mkdir (resultsDir);
end

saveFigInfo = 0;
try
    saveFigInfo = mancovanInfo.display.save_output;
catch
end

if (isdeployed)
    saveFigInfo = 1;
end

if (saveFigInfo)
    disp('Generating reults summary. Please wait ....');
end

try
    structFile = mancovanInfo.display.structFile;
catch
    structFile = fullfile(fileparts(which('gift.m')), 'icatb_templates', 'ch2bet.nii');
end

try
    imStr = mancovanInfo.display.image_values;
catch
    imStr = 'Positive';
end

try
    threshdesc = mancovanInfo.display.threshdesc;
catch
    threshdesc = 'fdr';
end

mancovanInfo.display.threshdesc = threshdesc;

mancovanInfo.display.image_values = imStr;
mancovanInfo.display.structFile = structFile;
htmlSummaryStr = [];
resultsInfo = [];

%% Features display (Spatial maps, Timecourses spectra, FNC correlations)
% * *a) Spatial maps and spectra* - Orthogonal slices of T-maps for each component and timecourses spectra is shown.
% * *b) FNC correlations* -  FNC correlations are averaged across subjects.
%
if (~isempty(strmatch(lower('features'), lower(display_steps), 'exact')))
    gH = icatb_plot_mancova_features(mancovanInfo);
    
    
    if (saveFigInfo)
        
        resultsInfo(end + 1).title = 'Features';
        resultsInfo(end).text = ['<ul> <li> <b> a) Spatial maps and spectra </b> - Orthogonal slices of T-maps for each component and timecourses spectra is shown. </li>', ...
            '<li> <b> b) FNC correlations </b> - FNC correlations are averaged across subjects. </li> ', ...
            ' </li> </ul>'];
        
        resultsInfo(end).files = printFigs([gH.H], resultsDir, 'features_comp');
        htmlSummaryStr(end + 1). title = 'Features';
        htmlSummaryStr(end).tag = 'features_comp';
        htmlSummaryStr(end).str = get_result_strings(resultsDir, resultsInfo(end), htmlSummaryStr(end).tag);
    end
end

drawnow;

%% Multivariate results
%  Multivariate tests are done on the features to determine the significant covariates which are later used in the univariate tests on each feature.
designCriteria = 'mancova';
try
    designCriteria = mancovanInfo.designCriteria;
catch
end

if (isfield(mancovanInfo, 'univInfo') && ~isempty(mancovanInfo.univInfo))
    designCriteria = 'none';
end

clear gH;

if (~isempty(strmatch(lower('mult'), lower(display_steps), 'exact')))
    if (strcmpi(designCriteria, 'mancova'))
        gH(1).H = icatb_plot_mult_mancovan(mancovanInfo);
        if (isfield(mancovanInfo, 'time'))
            % Time 1
            gH(end + 1).H = icatb_plot_mult_mancovan(mancovanInfo, 1);
            % Time 2
            gH(end + 1).H  = icatb_plot_mult_mancovan(mancovanInfo, 2);
        end
        
        if (saveFigInfo)
            resultsInfo(end + 1).title = 'Multivariate results';
            resultsInfo(end).text = ['<ul> <li>', ...
                'Multivariate tests are done on the features to determine the significant covariates which are later used in the univariate tests on each feature.</li> </ul>'];
            
            resultsInfo(end).files = printFigs([gH.H], resultsDir, 'multivariate_results');
            htmlSummaryStr(end + 1). title = 'Multivariate results';
            htmlSummaryStr(end).tag = 'multivariate';
            htmlSummaryStr(end).str = get_result_strings(resultsDir, resultsInfo(end), htmlSummaryStr(end).tag);
        end
        
    end
end

clear gH;
drawnow;

%% Univariate results
%
% * *a) Spatial maps* - T-maps of the significant covariate are shown as composite t-maps. Beta-values are averaged over significant clusters.
% * *b) Spectra* - Univariate t-tests are done using the significant covariates on the spectra. Beta-values ate averaged over frequency bands.
% * *c) FNC* - Univariate t-tests are done using the significant covariates on the FNC correlations.
% Connectogram of FNC correlations is also shown. Thumbnails of mean component maps are also plotted.
outDir = mancovanInfo.outputDir;
outputFiles = mancovanInfo.outputFiles;

if (~isempty(strmatch(lower('uni'), lower(display_steps), 'exact')))
    
    start_terms = {};
    for nO = 1:length(outputFiles)
        for nR = 1:length(outputFiles(nO).filesInfo.result_files)
            load(fullfile(outDir, outputFiles(nO).filesInfo.result_files{nR}), 'UNI');
            for nT = 1:length(UNI.tests)
                for nC = 1:length(UNI.stats{nT}.Contrast)
                    if (~isempty(UNI.stats{nT}.Contrast{nC}))
                        start_terms{end + 1} = UNI.stats{nT}.Contrast{nC};
                    else
                        start_terms{end + 1} = UNI.tests{nT};
                    end
                end
            end
        end
    end
    
    if (~isempty(start_terms))
        [ddd, indbb] = unique(start_terms);
        start_terms = start_terms(sort(indbb));
    end
    
    load(fullfile(mancovanInfo.outputDir, mancovanInfo.outputFiles(1).filesInfo.result_files{1}), 'UNI');
    if (~exist('UNI', 'var') || isempty(UNI))
        error('Please run mancova in order to view results');
    end
    
    gH = [];
    for nF = 1:length(start_terms)
        
        mancovanInfo.covariatesToPlot = start_terms(nF);
        figs  = icatb_plot_univariate_results(mancovanInfo);
        if (~isempty(figs))
            gH = [gH, [figs.H]];
        end
        drawnow;
        
        if (isfield(mancovanInfo, 'time'))
            % Time 1
            figs = icatb_plot_univariate_results(mancovanInfo, 1);
            if (~isempty(figs))
                gH = [gH, [figs.H]];
            end
            drawnow;
            % Time 2
            figs  =  icatb_plot_univariate_results(mancovanInfo, 2);
            if (~isempty(figs))
                gH = [gH, [figs.H]];
            end
            drawnow;
        end
    end
    
    if( ~isempty(gH))
        if (saveFigInfo)
            resultsInfo(end + 1).title = 'Univariate results';
            resultsInfo(end).text = ['<ul> <li> <b> a) Spatial Maps </b>  - T-maps of the significant covariate are shown as composite t-maps. Beta-values are averaged over significant clusters. </li>', ...
                '<li> <b> b) FNC Spectra </b> - Univariate t-tests are done using the significant covariates on the spectra. Beta-values ate averaged over frequency bands. </li> ', ...
                '<li> <b> c) FNC </b> - Univariate t-tests are done using the significant covariates on the FNC correlations. Connectogram of FNC correlations is also shown. Thumbnails of mean component maps are also plotted.</li></ul>'];
            resultsInfo(end).files = printFigs(gH, resultsDir, 'univariate_results');
            htmlSummaryStr(end + 1). title = 'Univariate results';
            htmlSummaryStr(end).tag = 'univariate';
            htmlSummaryStr(end).str = get_result_strings(resultsDir, resultsInfo(end), htmlSummaryStr(end).tag);
        end
    end
    
end

if (saveFigInfo)
    writeHTML2(htmlFile, htmlSummaryStr);
    disp('Done');
    close all;
end


function applyDefs(figH)

for nH = 1:length(figH)
    
    H = figH(nH).H;
    
    pos = get(0, 'defaultFigurePosition');
    set(H, 'color', 'w');
    %set(H, 'position', pos);
    set(H, 'resize', 'on');
    set(H, 'visible', 'on');
    titleH = get(findobj(H, 'type', 'axes'), 'title');
    if (iscell(titleH))
        % titleH = cell2mat(titleH);
        for nT = 1:length(titleH)
            set(titleH{nT}, 'color', 'k');
        end
    else
        set(titleH, 'color', 'k');
    end
    axesH = findobj(H, 'type', 'axes');
    set(axesH, 'YColor', 'k', 'XColor', 'k');
    set(axesH, 'fontname', 'times');
    set(axesH, 'fontsize', 12);
    textH = findobj(H, 'type', 'text');
    set(textH, 'color', 'k');
    set(textH, 'fontname', 'times');
    set(textH, 'fontsize', 12);
    try
        C = findobj(H, 'type', 'colorbar');
        set(C, 'color', 'k');
        set(get(C, 'label'), 'color', 'k');
    catch
    end
    
end


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


function imNames = printFigs(Figs, outdir, filename)

icatb_defaults;
global GICA_RESULTS_SUMMARY;

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
    imNames{nF} = [filename, '_', icatb_returnFileIndex(nF), '.png'];
    print(Figs(nF), '-dpng', PRINT_RESOLUTION, '-noui', fullfile(outdir, imNames{nF}));
end

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
