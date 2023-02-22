function icatb_network_summary(network_opts)


icatb_defaults;
global UI_FS;
global IMAGE_VALUES;
global CONVERT_Z;
global THRESHOLD_VALUE;
global FONT_COLOR;
global GICA_RESULTS_SUMMARY;


image_values = 'positive';
thresholds = 1;
structFile = fullfile(fileparts(which('groupica.m')), 'icatb_templates', 'ch2bet_3x3x3.nii');
convert_to_z = 'yes';
prefix = 'IC';
fnc_colorbar_label = 'Corr';
fnc_matrix_file = [];

outputDir = pwd;
try
    outputDir = network_opts.outputDir;
catch
end

if (~exist(outputDir, 'dir'))
    mkdir(outputDir);
end

try
    file_names = network_opts.file_names;
catch
end

try
    fnc_matrix_file = network_opts.fnc_matrix_file;
catch
end

try
    image_values = network_opts.image_values;
catch
end

try
    thresholds = network_opts.threshold;
catch
end

try
    structFile = network_opts.structFile;
catch
end

try
    convert_to_z = network_opts.convert_to_z;
catch
end

try
    prefix = network_opts.prefix;
catch
end

try
    fnc_colorbar_label = network_opts.fnc_colorbar_label;
catch
end

conn_threshold = [];
try
    conn_threshold = network_opts.conn_threshold;
catch
end

if (ischar(conn_threshold))
    conn_threshold = str2num(conn_threshold);
end

display_type = 'slices';
try
    display_type = network_opts.display_type;
catch
end

slice_plane = 'sagittal';
try
    slice_plane = network_opts.slice_plane;
catch
end

imWidth = [];
try
    imWidth = network_opts.imWidth;
catch
end

if (ischar(imWidth))
    imWidth = str2num(imWidth);
end

load icatb_colors coldhot;
cmap_corr = coldhot(1:4:end, :);
try
    cmap_corr = network_opts.cmap;
catch
end

CLIM_corr = [];
try
    CLIM_corr = network_opts.CLIM;
catch
end

if (ischar(CLIM_corr))
    CLIM_corr = str2num(CLIM_corr);
end

comp_network_names = network_opts.comp_network_names;

save_info = 0;

try
    save_info = network_opts.save_info;
catch
end

if (isdeployed)
    save_info = 1;
end

try
    resultsFormat = GICA_RESULTS_SUMMARY.format;
catch
    resultsFormat = 'html';
end

try
    resultsFormat = network_opts.format;
catch
end


pdfPrefix = [prefix, '_network_summary'];


printRes = '';
try
    print_res = GICA_RESULTS_SUMMARY.print_resolution;
catch
end

if (isempty(printRes))
    printRes = '-r72';
end

% for i = 1:2:length(varargin)
%     if (strcmpi(varargin{i}, 'image_values'))
%         image_values = varargin{i + 1};
%     elseif (strcmpi(varargin{i}, 'convert_to_z'))
%         convert_to_z = varargin{i + 1};
%     elseif (strcmpi(varargin{i}, 'threshold'))
%         thresholds = varargin{i + 1};
%     elseif (strcmpi(varargin{i}, 'anatomical_file'))
%         structFile = varargin{i + 1};
%     elseif (strcmpi(varargin{i}, 'prefix'))
%         prefix = varargin{i + 1};
%         elseif (strcmpi(varargin{i}, 'fnc_colorbar_label'))
%             fnc_colorbar_label = varargin{i + 1};
%     end
%     end
% end
%
% if (~exist('file_names', 'var') || isempty(file_names))
%     file_names = icatb_selectEntry('typeSelection', 'multiple', 'typeEntity', 'file', 'filter', '*.img;*.nii', 'title', 'Select component nifti images ...', 'fileType', 'image');
% end
%
%
% useGUI = 0;
% if (~exist('comp_network_names', 'var'))
%     useGUI = 1;
% end
%
% if (isempty(file_names))
%     error('Files are not selected');
% end
%
% if (~exist('fnc_matrix_file', 'var'))
%     fnc_matrix_file = icatb_selectEntry('typeSelection', 'single', 'typeEntity', 'file', 'filter', '*.mat;*.txt;*.dat', 'title', 'Select FNC correlations file');
% end

drawnow;

file_names = icatb_rename_4d_file(file_names);

% Load fnc matrix
FNCM = loadFNC(fnc_matrix_file);


drawnow;


fncComps = comp_network_names(:, 2);
fncComps = [fncComps{:}];
fncComps = fncComps(:);

if (~isempty(FNCM))
    if ((size(FNCM, 1) ~= length(fncComps)) && (size(FNCM, 1) ~= size(file_names, 1)))
        error('FNC matrix dimensions does not match the number of components entered');
    end
    
    if (size(FNCM, 1) == size(file_names, 1))
        % Truncate
        FNCM = FNCM(fncComps, fncComps);
    end
end

if (length(thresholds) == 1)
    thresholds = thresholds.*ones(1, size(comp_network_names, 1));
end

lengths = cellfun(@length, comp_network_names(:,2), 'UniformOutput', false);
lengths = cellfun(@length, comp_network_names(:,2), 'UniformOutput', false);
[dd, maxRInds] = max([lengths{:}]);

printInfo = [];

%% Rendering: Multiple components are rendered on surface of a 'standard' brain
nrows = ceil(size(comp_network_names, 1)/sqrt(size(comp_network_names, 1)));
ncols = ceil(size(comp_network_names, 1)/nrows);
render_pngs = {};
countPdfs = 0;
for nF = 1:size(comp_network_names,1)
    fH = icatb_getGraphics('Render', 'normal', 'render_view', 'on');
    set(fH, 'resize', 'on');
    set(fH,'menubar','none');
    
    rData = icatb_display_composite(file_names(comp_network_names{nF, 2}, :), 'convert_to_zscores', convert_to_z, 'threshold', thresholds(1), 'image_values', image_values, ...
        'display_type', 'render');
    
    nLabels = comp_network_names{nF, 2};
    nLabels = cellstr(num2str(nLabels(:)));
    sh = axes('parent', fH, 'units', 'normalized', 'position', [0.16, 0.16, 0.72, 0.72]);
    image(rData.slices, 'parent', sh);
    colormap(rData.cmap);
    title(comp_network_names{nF, 1}, 'parent', sh, 'color', FONT_COLOR);
    try
        aspRatio = size(rData.slices, 1)/size(rData.slices, 2);
        set(sh, 'PlotBoxAspectRatio', [1, aspRatio, 1]);
    catch
        axis(sh, 'image;');
    end
    axis(sh, 'off');
    ch = colorbar('peer', sh);
    set(sh, 'clim', [1, size(rData.cmap, 1)]);
    pos = get(ch, 'position');
    origPos = pos(4);
    pos(4) = 0.7*pos(4);
    pos(2) = pos(2) + 0.5*(origPos - pos(4));
    pos(1) = pos(1) + pos(3) + 0.03;
    colorLen = rData.colorLen;
    colorMid = ceil(colorLen/2) + 1;
    xticks = colorMid:colorLen:(lengths{nF}*colorLen);
    set(ch, 'ylim', [1, colorLen*lengths{nF}]);
    set(ch, 'position', pos);
    set(ch, 'ytick', xticks);
    set(ch, 'YtickLabel', nLabels);
    set(ch, 'color', FONT_COLOR);
    tmp_f_name = [prefix, '_render_', num2str(nF)];
    outFNames = fullfile(outputDir, [tmp_f_name, '.fig']);
    saveas(fH, outFNames);
    pngFNames = printFigToPNG(fH, outputDir, tmp_f_name);
    
    render_pngs(end + 1) = pngFNames;
    
    if (strcmpi(resultsFormat, 'pdf'))
        countPdfs = countPdfs + 1;
        tmpImFile = [pdfPrefix, '_', icatb_returnFileIndex(countPdfs), '.pdf'];
        set(fH, 'PaperPositionMode', 'auto');
        print(fH, '-dpdf', print_res, '-noui', fullfile(outputDir, tmpImFile));
        pdfFiles{countPdfs} = fullfile(outputDir, tmpImFile);
    end
    
    if (save_info)
        delete(fH);
    end
    
end


printInfo(end + 1).title = 'Rendering';
printInfo(end).text = 'Rendering: Multiple components are rendered on surface of a ''standard'' brain';
printInfo(end).files = render_pngs;
printInfo(end).tag = 'render';
printInfo(end).str = get_result_strings(outputDir, printInfo(end), printInfo(end).tag);


%% FNC correlations: Correlations are visualized in a matrix plot
if (~isempty(FNCM))
    %     load icatb_colors coldhot;
    %     cmap = coldhot(1:4:end, :);
    FNCM(isfinite(FNCM) == 0) = 0;
    CLIM = max(abs(FNCM(:)));
    network_values = cellfun(@length, comp_network_names(:,2), 'UniformOutput', false);
    network_values = [network_values{:}];
    %fH = figure('color', 'w');
    gH = icatb_getGraphics('FNC Correlations', 'graphics', 'fnc_corrs', 'on');
    pos = get(gH, 'position');
    %set(gH, 'position', pos);
    icatb_plot_FNC(FNCM, [-CLIM, CLIM], cellstr(num2str(fncComps)), (1:length(fncComps)), gH, fnc_colorbar_label, gca, ...
        network_values, comp_network_names(:,1));
    %     set(findobj(gH, 'type', 'line'), 'linewidth', 1.5);
    %     set(findobj(gH, 'type', 'line'), 'color', [0, 0, 0]);
    colormap(jet);
    SetProps(gH);
    tmp_f_name = [prefix, '_FNC'];
    outFNames = fullfile(outputDir, [tmp_f_name, '.fig']);
    saveas(gH, outFNames);
    pngFNames = printFigToPNG(gH, outputDir, tmp_f_name);
    
    if (strcmpi(resultsFormat, 'pdf'))
        countPdfs = countPdfs + 1;
        tmpImFile = [pdfPrefix, '_', icatb_returnFileIndex(countPdfs), '.pdf'];
        set(gH, 'PaperPositionMode', 'auto');
        print(gH, '-dpdf', print_res, '-noui', fullfile(outputDir, tmpImFile));
        pdfFiles{countPdfs} = fullfile(outputDir, tmpImFile);
    end
    
    
    if (save_info)
        delete(gH);
    end
    
    printInfo(end + 1).title = 'FNC correlations';
    printInfo(end).text = 'FNC correlations: Correlations are visualized in a matrix plot';
    printInfo(end).files = pngFNames;
    printInfo(end).tag = 'fnc_corrs';
    printInfo(end).str = get_result_strings(outputDir, printInfo(end), printInfo(end).tag);
    
    drawnow;
    
    %% Connectogram view - FNC correlations are shown using bezier curves and thumbnails of spatial maps are shown in a circle. Components within the same network are shown in the same color.
    fH = icatb_plot_connectogram([], comp_network_names, 'C', FNCM, 'threshold', thresholds(1), 'image_file_names', file_names, 'colorbar_label', fnc_colorbar_label, 'cmap', cmap_corr, ...
        'slice_plane', slice_plane, 'conn_threshold', conn_threshold, 'imwidth', imWidth, 'display_type', display_type, 'CLIM', CLIM_corr, ...
        'template_file', structFile);
    tmp_f_name = [prefix, '_connectogram'];
    outFNames = fullfile(outputDir, [tmp_f_name, '.fig']);
    saveas(fH, outFNames);
    pngFNames = printFigToPNG(fH, outputDir, tmp_f_name);
    
    
    if (strcmpi(resultsFormat, 'pdf'))
        countPdfs = countPdfs + 1;
        tmpImFile = [pdfPrefix, '_', icatb_returnFileIndex(countPdfs), '.pdf'];
        set(fH, 'PaperPositionMode', 'auto');
        print(fH, '-dpdf', print_res, '-noui', '-bestfit', fullfile(outputDir, tmpImFile));
        pdfFiles{countPdfs} = fullfile(outputDir, tmpImFile);
    end
    
    if (save_info)
        delete(fH);
    end
    
    printInfo(end + 1).title = 'Connectogram';
    printInfo(end).text = 'Connectogram view - FNC correlations are shown using bezier curves and thumbnails of spatial maps are shown in a circle. Components within the same network are shown in the same color.';
    printInfo(end).files = pngFNames;
    printInfo(end).tag = 'connectogram';
    printInfo(end).str = get_result_strings(outputDir, printInfo(end), printInfo(end).tag);
    
end

pause(1);

drawnow;

%% Multiple components are displayed in a composite plot. Orthogonal slices are shown.
gH = icatb_plot_composite_orth_views(file_names, comp_network_names, 'convert_to_zscores', convert_to_z, 'threshold', thresholds(1), 'image_values', image_values, ...
    'anatomical_file', structFile);
tmp_f_name = [prefix, '_composite_orth_views'];
outFNames = fullfile(outputDir, [tmp_f_name, '.fig']);
saveas(gH, outFNames);
pngFNames = printFigToPNG(gH, outputDir, tmp_f_name);

printInfo(end + 1).title = 'Composite Orthogonal Slices';
printInfo(end).text = 'Multiple components are displayed in a composite plot. Orthogonal slices are shown.';
printInfo(end).files = pngFNames;
printInfo(end).tag = 'composite_ortho';
printInfo(end).str = get_result_strings(outputDir, printInfo(end), printInfo(end).tag);

if (strcmpi(resultsFormat, 'pdf'))
    countPdfs = countPdfs + 1;
    tmpImFile = [pdfPrefix, '_', icatb_returnFileIndex(countPdfs), '.pdf'];
    set(gH, 'PaperPositionMode', 'auto');
    print(gH, '-dpdf', print_res, '-noui', '-bestfit', fullfile(outputDir, tmpImFile));
    pdfFiles{countPdfs} = fullfile(outputDir, tmpImFile);
end

if (save_info)
    delete(gH);
end

pause(1);

drawnow;

if (0)
    % skip this part for now as it involves cropping images and text
    % automatically
    
    %% Stacked ortho slices are shown for each component in the network. Title shows component numbers plotted from top to bottom
    icatb_groupNetworks(file_names, comp_network_names(:, 1), comp_network_names(:, 2), 'convert_to_z', convert_to_z, 'threshold', thresholds(1), 'image_values', image_values, ...
        'structfile', structFile, 'interMediatePlot', 1);
    tmp_f_name = [prefix, '_orth_views'];
    %outFNames = fullfile(outputDir, [tmp_f_name, '.fig']);
    gH = findobj(0, 'tag', 'group_networks');
    
    for nG = 1:length(gH)
        
        outFNames = fullfile(outputDir, [tmp_f_name, '_', num2str(nG), '.fig']);
        
        savefig(gH(nG), outFNames);
        
        if (strcmpi(resultsFormat, 'pdf'))
            countPdfs = countPdfs + 1;
            tmpImFile = [pdfPrefix, '_', icatb_returnFileIndex(countPdfs), '.pdf'];
            set(gH(nG), 'PaperPositionMode', 'auto');
            print(gH(nG), '-dpdf', print_res, '-noui', '-bestfit', fullfile(outputDir, tmpImFile));
            pdfFiles{countPdfs} = fullfile(outputDir, tmpImFile);
        end
        
    end
    
    pngFNames = printFigToPNG(gH, outputDir, tmp_f_name);
    printInfo(end + 1).title = 'Stacked Ortho Slices';
    printInfo(end).text = 'Stacked ortho slices are shown for each component in the network. Title shows component numbers plotted from top to bottom';
    printInfo(end).files = pngFNames;
    printInfo(end).tag = 'stacked_ortho';
    printInfo(end).str = get_result_strings(outputDir, printInfo(end), printInfo(end).tag);
    
    if (save_info)
        delete(gH);
    end
    
end


if (save_info)
    if (strcmpi(resultsFormat, 'pdf'))
        append_pdfs(fullfile(outputDir, [pdfPrefix, '.pdf']), pdfFiles{:});
        for nF = 1:length(pdfFiles)
            try
                delete(pdfFiles{nF});
            catch
            end
        end
    else
        writeHTML2(fullfile(outputDir, [prefix, '_network_summary.html']), printInfo);
    end
end


function FNCM = loadFNC(file_name)

if (isempty(file_name))
    FNCM = [];
end

if (ischar(file_name))
    [pp, fn, extn] = fileparts(file_name);
    
    if (strcmpi(extn, '.nii') || strcmpi(extn, '.img'))
        FNCM = icatb_loadData(file_name);
    else
        FNCM = icatb_load_ascii_or_mat(file_name);
    end
    
else
    FNCM = file_name;
end

if (numel(FNCM) == length(FNCM))
    FNCM = icatb_vec2mat(FNCM);
end


function imNames = printFigToPNG(Figs, outdir, filename)

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
            results_string1 = [results_string1, tableStr, '</table>'];
        end
    end
    
    results_string1 = [results_string1, '<p>&nbsp;&nbsp;</p>'];
end

end_string =  '</html>';


results_string = ['<p> </p> ', results_string1, ' <p> </p>'];



function writeHTML2(fileName, resultsStr)


start_string = '<html><head><title> Network Summary </title></head>';
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


function SetProps(H)

%pos = get(0, 'defaultFigurePosition');
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
%set(axesH, 'fontsize', 8);
textH = findobj(H, 'type', 'text');
set(textH, 'color', 'k');
set(textH, 'fontname', 'times');
%set(textH, 'fontsize', 8);
try
    C = findobj(H, 'type', 'colorbar');
    set(C, 'color', 'k');
    set(get(C, 'label'), 'color', 'k');
catch
end
