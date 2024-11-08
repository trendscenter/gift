function icatb_display_fnc_ica(param_file, opts)
% Display FNC ICA components
%

icatb_defaults;
global GICA_RESULTS_SUMMARY;
global UI_FS;
global FONT_COLOR;

if ~exist('param_file', 'var')
    param_file = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Parameter File', 'filter', '*ica*param*.mat');
    if isempty(param_file)
        error('Parameter file is not selected for display');
    end
end

if (ischar(param_file))
    load(param_file);
else
    sesInfo = param_file;
end

outputDir = sesInfo.outputDir;

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

resultsFormat = 'html';
try
    resultsFormat = opts.format;
catch
end

resultsDir =  fullfile(outputDir, [sesInfo.userInput.prefix, '_fnc_results']);
try
    resultsDir = opts.outputDir;
catch
end


if (saveFigInfo)
    if (exist(resultsDir, 'dir') ~= 7)
        mkdir(resultsDir);
    end
    pdfPrefix = mfilename;
end

colorbar_label = '';
if strcmpi(sesInfo.scaleType, 'z-scores')
    colorbar_label = 'z-scores';
end


reference_file = '';
try
    reference_file = sesInfo.userInput.reference_file;
catch
end


ica_mat_file = fullfile(sesInfo.outputDir, [sesInfo.prefix, '_ica_c1-1.mat']);
load(ica_mat_file);

fnc_components = 1:size(ic, 2);

if (isempty(reference_file))
    comp_network_names = {'All', fnc_components};
else
    comp_network_names = sesInfo.userInput.network_summary_opts;
end

load icatb_colors coldhot;
cmap_corr = coldhot(1:4:end, :);

fncComps = comp_network_names(:, 2);
fncComps = [fncComps{:}];
fncComps = fncComps(:);

newText = dispParams(sesInfo);

disp(newText);

%% FNC ICA Parameters
pdfFiles = {};
htmlSummaryStr = [];
resultsInfo =[];

if (saveFigInfo)
    if (~strcmpi(resultsFormat, 'pdf'))
        icatb_print_table(cellstr(newText), fullfile(resultsDir, 'params.txt'));
        resultsInfo(end + 1).title = 'FNC ICA Parameters';
        resultsInfo(end).text = ' ';
        resultsInfo(end).files = {'params.txt'};
        htmlSummaryStr(end + 1).title = 'FNC ICA Parameters';
        htmlSummaryStr(end).tag =  'fncica_params';
        htmlSummaryStr(end).str = get_result_strings(resultsDir, resultsInfo(end), htmlSummaryStr(end).tag);
        %writeHTML(resultsDir, resultsInfo(end), 'introduction.html', 'Group ICA Parameters');
    else
        print_strs = cellstr(char('FNC ICA parameters:', '', newText));
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


%% FNC Components
fncHandles = [];
for nComp = 1:sesInfo.numComp
    
    FNCM = squeeze(ic(nComp, :, :));
    FNCM(isfinite(FNCM) == 0) = 0;
    CLIM = max(abs(FNCM(:)));
    network_values = cellfun(@length, comp_network_names(:,2), 'UniformOutput', false);
    network_values = [network_values{:}];
    %fH = figure('color', 'w');
    gH = icatb_getGraphics('FNC Correlations', 'graphics', 'fnc_corrs', 'on');
    set(gH, 'resize', 'on');
    pos = get(gH, 'position');
    %set(gH, 'position', pos);
    CLIM_corr = [-CLIM, CLIM];
    icatb_plot_FNC(FNCM, CLIM_corr, cellstr(num2str(fncComps)), (1:length(fncComps)), gH, colorbar_label, gca, ...
        network_values, comp_network_names(:,1));
    title(['FNC Comp ', icatb_returnFileIndex(nComp)]);
    %     set(findobj(gH, 'type', 'line'), 'linewidth', 1.5);
    %     set(findobj(gH, 'type', 'line'), 'color', [0, 0, 0]);
    colormap(cmap_corr);
    
    fncHandles(end + 1) = gH;
    
    
    drawnow;
    
    if (~isempty(reference_file))
        
        % Connectogram view - FNC correlations are shown using bezier curves and thumbnails of spatial maps are shown in a circle. Components within the same network are shown in the same color.
        fH = icatb_plot_connectogram([], comp_network_names, 'C', FNCM, 'image_file_names', reference_file, 'cmap', cmap_corr, 'CLIM', CLIM_corr, ...
            'title', ['FNC Comp ', icatb_returnFileIndex(nComp)], 'colorbar_label', colorbar_label);
        fncHandles(end + 1) = fH;
    end
    
    drawnow;
    
    
end



resultsInfo = [];
if (saveFigInfo)
    
    if (~strcmpi(resultsFormat, 'pdf'))
        resultsInfo(end + 1).title = 'FNC correlations on loadings';
        resultsInfo(end).text = 'Functional network connectivity correlations on subject loading coefficients.';
        FNCTPNGFiles = printFigs(fncHandles, resultsDir, 'FNC_loadings');
        resultsInfo(end).files = FNCTPNGFiles;
        %writeHTML(resultsDir, resultsInfo(end), 'fnc_timecourses.html', 'FNC correlations on timecourses');
        htmlSummaryStr(end + 1).title = 'FNC Components';
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

function writeHTML2(fileName, resultsStr)


start_string = '<html><head><title> FNC ICA Results </title></head>';
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


function newText = dispParams(sesInfo)

icaAlgo = icatb_icaAlgorithm;
selected_ica_algorithm = deblank(icaAlgo(sesInfo.algorithm, :));

D(1).string = '....................................................';

D(size(D,2)+1).string = '';

modalityType = 'FNC';

D(size(D,2)+1).string = ['Number of Subjects : ', num2str(sesInfo.numOfSub)];
D(size(D,2)+1).string = ['Modality : ', modalityType];
D(size(D,2)+1).string = ['Number of Independent Components : ', num2str(sesInfo.numComp)];
D(size(D,2)+1).string = ['ICA Algorithm : ', selected_ica_algorithm];

D(size(D,2)+1).string = ['Scaling Components : ', sesInfo.scaleType];

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

D(size(D,2)+1).string = '';
D(size(D,2)+1).string = '....................................................';
D(size(D,2)+1).string = '';
newText = char(D.string);