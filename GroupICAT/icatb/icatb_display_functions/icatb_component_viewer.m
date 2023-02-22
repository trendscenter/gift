function icatb_component_viewer(paramFile, varargin)
% Spectral viewer tool
%


icatb_defaults;
global DETRENDNUMBER;
global PARAMETER_INFO_MAT_FILE;
global EXPERIMENTAL_TR;
global FONT_COLOR;
global UI_FS;

filterP = ['*', PARAMETER_INFO_MAT_FILE, '*.mat'];

if (~exist('paramFile', 'var'))
    paramFile = icatb_selectEntry('typeEntity', 'file', 'title', 'Select ICA parameter file', 'filter', filterP, 'typeselection', 'single');
end

drawnow;

load(paramFile);

if (~exist('sesInfo', 'var'))
    error('Selected file is not a valid parameter file');
end


for i = 1:2:length(varargin)
    if (strcmpi(varargin{i}, 'tr'))
        TR = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'comps'))
        compNumbers = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'detrend_no'))
        detrendNumber = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'image_values'))
        image_values = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'convert_to_zscores'))
        convert_to_zscores = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'threshold'))
        threshold = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'anatomical_file'))
        anatomical_file = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'compFiles'))
        compFiles = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'freq_limits'))
        freq_limits = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'average_timecourses'))
        average_timecourses = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'output_format'))
        output_format = varargin{i + 1};
    end
end

outputDir = fileparts(paramFile);
if (isempty(outputDir))
    outputDir = pwd;
end

cd (outputDir);

sesInfo.outputDir = outputDir;

if (~exist('detrendNumber', 'var'))
    detrendNumber = DETRENDNUMBER;
end

startPath = fileparts(which('gift.m'));
startPath = fullfile(startPath, 'icatb_templates');

% Select anatomical file
if (~exist('anatomical_file', 'var'))
    anatomical_file = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Structural Image', 'filter', '*.img;*.nii', ...
        'typeselection', 'single', 'fileType', 'image', 'filenumbers', 1, 'startpath', startPath);
end

drawnow;

[startPath, f_name, extn] = fileparts(icatb_parseExtn(anatomical_file));

% check the image extension
if ~strcmpi(extn, '.nii') && ~strcmpi(extn, '.img')
    error('Structural image should be in NIFTI or Analyze format');
end

cd (outputDir);

%% Uncompress component files
if (sesInfo.numOfSub*sesInfo.numOfSess > 1)
    listStr = getListStr(sesInfo.icaOutputFiles, sesInfo.outputDir);
    index = icatb_listdlg('PromptString', 'Select Components Set', 'SelectionMode','single',...
        'ListString', listStr, 'movegui', 'center', 'windowStyle', 'modal');
    if (isempty(index))
        error('Components set is not selected');
    end
    viewing_set = deblank(listStr(index, :));
    % load components
    dashIndex = icatb_findstr(viewing_set, '-');
    spaceIndex = icatb_findstr(viewing_set, ' ');
    a = viewing_set(1 : dashIndex(1) - 1);
    b = viewing_set(dashIndex(1) + 1 : spaceIndex(1));
    compFiles = sesInfo.icaOutputFiles(str2num(a)).ses(str2num(b)).name;
    clear a b dashIndex spaceIndex;
else
    compFiles = sesInfo.icaOutputFiles(1).ses(1).name;
    viewing_set = getListStr(sesInfo.icaOutputFiles, sesInfo.outputDir);
    viewing_set = deblank(viewing_set(1, :));
end

% Select components
listStr = num2str((1:sesInfo.numComp)');

if ~exist('compNumbers', 'var')
    title_fig = 'Select component/components';
    compNumbers = icatb_listdlg('PromptString', title_fig, 'SelectionMode', 'multiple', 'ListString', listStr, ...
        'movegui', 'center', 'windowStyle', 'modal', 'title_fig', title_fig);
end

compNumbers = compNumbers(:)';

if (isempty(compNumbers))
    error('Component/Components are selected');
end

if (max(compNumbers) > sesInfo.numComp)
    error('Max value of component numbers selected exceed the no. of components');
end


% Select component parameters
numParameters = 1;

if (~exist('image_values', 'var'))
    opts = char('Positive', 'Positive and Negative', 'Absolute Value', 'Negative');
    inputText(numParameters).promptString = 'Select image values';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = opts;
    inputText(numParameters).dataType = 'string';
    inputText(numParameters).tag = 'image_values';
    inputText(numParameters).enable = 'on';
end

if (~exist('convert_to_zscores', 'var'))
    numParameters = numParameters + 1;
    inputText(numParameters).promptString = 'Do you want to convert to z-scores?';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = char('Yes', 'No');
    inputText(numParameters).dataType = 'string';
    inputText(numParameters).tag = 'convert_to_zscores';
    inputText(numParameters).enable = 'on';
end

if (~exist('threshold', 'var'))
    numParameters = numParameters + 1;
    inputText(numParameters).promptString = 'Enter threshold';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = '1';
    inputText(numParameters).dataType = 'numeric';
    inputText(numParameters).tag = 'threshold';
    inputText(numParameters).enable = 'on';
end

if (~exist('TR', 'var'))
    if (isfield(sesInfo, 'TR'))
        TR = sesInfo.TR;
    else
        TR = EXPERIMENTAL_TR;
    end
    numParameters = numParameters + 1;
    inputText(numParameters).promptString = 'Enter TR in seconds';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(TR);
    inputText(numParameters).dataType = 'numeric';
    inputText(numParameters).tag = 'TR';
    inputText(numParameters).enable = 'on';
end

if (~exist('freq_limits', 'var'))
    numParameters = numParameters + 1;
    inputText(numParameters).promptString = 'Enter low frequency and high frequency to compute fALFF (Hz)';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = '[0.1, 0.15]';
    inputText(numParameters).dataType = 'numeric';
    inputText(numParameters).tag = 'freq_limits';
    inputText(numParameters).enable = 'on';
    
end


if (~exist('average_timecourses', 'var'))
    numParameters = numParameters + 1;
    inputText(numParameters).promptString = 'Do you want to average timecourses across sessions?';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = char('Yes', 'No');
    inputText(numParameters).dataType = 'string';
    inputText(numParameters).tag = 'average_timecourses';
    inputText(numParameters).enable = 'on';
    
end

if (~exist('output_format', 'var'))
    numParameters = numParameters + 1;
    inputText(numParameters).promptString = 'Select results format';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = char('HTML', 'PDF', 'None');
    inputText(numParameters).dataType = 'string';
    inputText(numParameters).tag = 'output_format';
    inputText(numParameters).enable = 'on';
    
end

if (exist('inputText', 'var'))
    answers = icatb_inputDialog('inputtext', inputText, 'Title', 'Component Parameters', 'handle_visibility',  'on', 'windowStyle', 'modal');
    if (isempty(answers))
        error('Component parameters are not selected');
    end
    for nA = 1:length(answers)
        val = answers{nA};
        eval([inputText(nA).tag, '=val;']);
    end
end

if (~exist('output_format', 'var'))
    output_format = 'none';
end

if (length(TR) == 1)
    TR = repmat(TR, 1, sesInfo.numOfSub);
else
    if (length(TR) ~= sesInfo.numOfSub)
        error('Length of TR must match the number of subjects');
    end
end

% Compute spectra
disp('Loading subject timecourses ...');
minTpLength = min(sesInfo.diffTimePoints);
doInterp = 0;

if strcmpi(average_timecourses, 'yes')
    
    if (~all(TR == min(TR)))
        ratiosTR = (TR(:)')./min(TR);
        [numN, denN] = rat(ratiosTR);
        chkTp = ceil((sesInfo.diffTimePoints(1:sesInfo.numOfSess:end).*numN)./denN);
        minTpLength = min(chkTp);
        doInterp = 1;
    end
    
    for nSub = 1:sesInfo.numOfSub
        timecourses = icatb_loadComp(sesInfo, compNumbers, 'subjects', nSub, 'vars_to_load', 'tc', 'average_runs', 1, 'detrend_no', detrendNumber);
        timecourses = interpTC(timecourses, TR, nSub, minTpLength, doInterp);
        [temp_spectra, freq] = icatb_get_spectra(timecourses', min(TR));
        temp_spectra = temp_spectra./repmat(sum(temp_spectra, 2), [1, size(temp_spectra, 2)]);
        temp_spectra = temp_spectra';
        if (nSub == 1)
            power_spectra = zeros(sesInfo.numOfSub, size(temp_spectra, 1), length(compNumbers));
        end
        power_spectra(nSub, :, :) = temp_spectra;
    end
    
else
    
    if (~all(TR == min(TR)))
        ratiosTR = (TR(:)')./min(TR);
        [numN, denN] = rat(ratiosTR);
        chkTp = ceil((sesInfo.diffTimePoints(:)'.*numN)./denN);
        minTpLength = min(chkTp);
        doInterp = 1;
    end
    
    countS = 0;
    for nSub = 1:sesInfo.numOfSub
        for nSess = 1:sesInfo.numOfSess
            countS = countS + 1;
            timecourses = icatb_loadComp(sesInfo, compNumbers, 'subjects', nSub, 'sessions', nSess, 'vars_to_load', 'tc', 'truncate_tp', 1, 'detrend_no', ...
                detrendNumber);
            timecourses = interpTC(timecourses, TR, nSub, minTpLength, doInterp);
            %timecourses = timecourses(1:min(sesInfo.diffTimePoints), :);
            [temp_spectra, freq] = icatb_get_spectra(timecourses', min(TR));
            temp_spectra = temp_spectra./repmat(sum(temp_spectra, 2), [1, size(temp_spectra, 2)]);
            temp_spectra = temp_spectra';
            if (countS == 1)
                power_spectra = zeros(sesInfo.numOfSub, sesInfo.numOfSess, size(temp_spectra, 1), length(compNumbers));
            end
            power_spectra(nSub, nSess, :, :) = temp_spectra;
        end
    end
    if (sesInfo.numOfSess > 1)
        power_spectra = mean(power_spectra, 2);
    end
end


power_spectra = reshape(power_spectra, sesInfo.numOfSub, length(freq), length(compNumbers));

zipFileName = {};
currentFile = deblank(compFiles(1, :));
if ~exist(currentFile, 'file')
    [zipFileName, files_in_zip] = icatb_getViewingSet_zip(currentFile, [], 'real', sesInfo.zipContents);
    if (~isempty(zipFileName))
        icatb_unzip(regexprep(zipFileName, ['.*\', filesep], ''), fullfile(outputDir, fileparts(currentFile)));
    end
end

if (~exist('anatomical_file', 'var'))
    anatomical_file = fullfile(outputDir, currentFile);
end

compFiles = icatb_fullFile('files', compFiles, 'directory', outputDir);
compFiles = icatb_rename_4d_file(compFiles);

% Plot results
graphicsH = [];
count = 0;
dynamicrange = zeros(1, length(compNumbers));
fALFF = zeros(1, length(compNumbers));
imfFiles = cell(1, length(compNumbers) + 1);
htmlDir = fullfile(outputDir, [sesInfo.userInput.prefix, '_component_viewer_html']);

if (~strcmpi(output_format, 'none'))
    
    gH = icatb_getGraphics('Component viewer parameters', 'graphics', 'orthoviews', 'on');
    aH = axes('units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8], 'parent', gH, 'visible', 'off', 'color', get(gH, 'color'));
    
    print_strs{1} = 'Component viewer parameters:';
    print_strs{end + 1} = '';
    if (length(unique(TR)) == 1)
        print_strs{end + 1} = ['TR: ', num2str(TR(1)), ' sec'];
    else
        print_strs{end + 1} = ['TR: ', num2str(TR), ' sec'];
    end
    
    print_strs{end + 1} = ['Viewing set: ', viewing_set];
    print_strs{end + 1} = ['Components: ', num2str(compNumbers)];
    print_strs{end + 1} = ['Anatomical file: ', anatomical_file];
    print_strs{end + 1} = ['Image values: ', image_values];
    print_strs{end + 1} = ['Convert to Z-scores: ', convert_to_zscores];
    print_strs{end + 1} = ['Convert to threshold: ', num2str(threshold)];
    print_strs{end + 1} = ['Frequency limits for fALFF: ', num2str(freq_limits(1)), ', ', num2str(freq_limits(2)), ' Hz'];
    print_strs{end + 1} = ['Average timecourses: ', average_timecourses];
    
    yPos = 1.1;
    for n = 1:length(print_strs)
        if (n == 1)
            tmpFont = UI_FS;
        else
            tmpFont = UI_FS - 2;
        end
        text(0.45, yPos, print_strs{n}, 'units', 'normalized', 'horizontalalignment', 'center', 'verticalAlignment', 'middle', 'color', FONT_COLOR, 'fontsize', tmpFont);
        yPos = yPos - 0.075;
    end
    
    if (strcmpi(output_format, 'pdf'))
        tmpImFile = 'contents.pdf';
        print(gH, '-dpdf', '-noui', fullfile(outputDir, tmpImFile));
        imfFiles{1} = fullfile(outputDir, tmpImFile);
    elseif (strcmpi(output_format, 'html'))
        if (~exist(htmlDir, 'dir'))
            mkdir(htmlDir);
        end
        tmpImFile = 'contents.png';
        print(gH, '-dpng', '-noui', fullfile(htmlDir, tmpImFile));
        imfFiles{1} = tmpImFile;
    end
    
    try
        delete(gH);
    catch
    end
end

for ncomp = compNumbers
    currentCompFile = deblank(compFiles(ncomp, :));
    disp(['Plotting Component ', num2str(ncomp), ' ...']);
    count = count + 1;
    [tc, tmpd, tmpf, dynamicrange(count), fALFF(count)] = getSpecStats(squeeze(power_spectra(:, :, count)), freq_limits, freq);
    if (count == 1)
        dynamicrange_all = zeros(length(tmpd), length(compNumbers));
        fALFF_all = dynamicrange_all;
    end
    dynamicrange_all(:, count) = tmpd;
    fALFF_all(:, count) = tmpf;
    freq = tc.xAxis;
    H = icatb_orth_views(currentCompFile, 'structfile', anatomical_file, 'image_values', image_values, 'convert_to_zscores', convert_to_zscores, 'threshold', ...
        threshold, 'set_to_max_voxel', 1, 'tc', tc, 'labels', '', 'fig_title', [viewing_set, ' ', icatb_returnFileIndex(ncomp)]);
    graphicsH(length(graphicsH) + 1).H = H;
    fprintf('\n');
    [dddd, ffN, extn] = fileparts(currentCompFile);
    [extn, compFileN] = icatb_parseExtn(extn);
    if (strcmpi(output_format, 'pdf'))
        tmpImFile = [ffN, '_', num2str(compFileN), '.pdf'];
        print(H, '-dpdf', '-r300', '-noui', fullfile(outputDir, tmpImFile));
        imfFiles{count+1} = fullfile(outputDir, tmpImFile);
    elseif (strcmpi(output_format, 'html'))
        if (~exist(htmlDir, 'dir'))
            mkdir(htmlDir);
        end
        tmpImFile = [ffN, '_', num2str(compFileN), '.png'];
        print(H, '-dpng', '-r300', '-noui', fullfile(htmlDir, tmpImFile));
        imfFiles{count+1} = tmpImFile;
    end
end

icatb_plotNextPreviousExitButtons(graphicsH);

spectra_file = fullfile(outputDir, [sesInfo.userInput.prefix, '_spectra_results.mat']);
icatb_save(spectra_file, 'power_spectra', 'freq', 'dynamicrange', 'fALFF', 'compNumbers', 'freq_limits', 'TR', 'dynamicrange_all', 'fALFF_all');

disp(['Power spectra results are saved in ', spectra_file]);
disp('The variables are as follows:');
disp('a. power_spectra - Power spectra of dimensions subjects x spectral length x components');
disp('b. freq - Frequency in HZ');
disp('c. dynamicrange - Mean dynamic range over subjects ');
disp('d. fALFF - Mean ratio of low frequency power to the high frequency power over subjects');
disp('e. compNumbers - Component numbers.');
disp('f. freq_limits - Low frequeny (l < LF) and High frequency (f > HF) used to compute fALFF.');
fprintf('\n');

% Cleanup files
if (~isempty(zipFileName))
    icatb_delete_file_pattern(files_in_zip, outputDir);
end


if (strcmpi(output_format, 'pdf'))
    outImgFile = fullfile(outputDir, [sesInfo.userInput.prefix, '_component_viewer.pdf']);
    disp(['Writing PDF file ', outImgFile, ' ...']);
    append_pdfs(outImgFile, imfFiles{:});
    for nF = 1:length(imfFiles)
        try
            delete(imfFiles{nF});
        catch
        end
    end
    open(outImgFile);
elseif (strcmpi(output_format, 'html'))
    results_html(htmlDir, imfFiles, [sesInfo.userInput.prefix, '_component_viewer.html']);
end

fprintf('Done\n\n');


function [tc, dynamicrange, fALFF, mean_dyn, mfALFF] = getSpecStats(spectra_tc, freq_limits, freq)
%% Compute spectra

tc.data = spectra_tc;
tc.xAxis = freq;
tc.isSpectra = 1;
tc.xlabelStr = 'Frequency (Hz)';
tc.ylabelStr = 'Power';
dynamicrange = zeros(1, size(tc.data, 1));
fALFF = dynamicrange;
for nS = 1:size(tc.data, 1)
    [dynamicrange(nS), fALFF(nS)] = icatb_get_spec_stats(tc.data(nS, :), tc.xAxis, freq_limits);
end
mean_dyn = mean(dynamicrange);
mfALFF = mean(fALFF);
tc.titleStr = sprintf('Dynamic range: %0.3f, Power_L_F/Power_H_F: %0.3f', mean_dyn, mfALFF);


function listStr = getListStr(outputFiles, outputDir)

outputFiles = chkOutputFiles(outputFiles, outputDir);
counter = 1;
numOfSets = length(outputFiles);
% loop over number of sets
for jj = 1 : numOfSets
    % loop over sessions
    for kk = 1 : length(outputFiles(jj).ses)
        str = deblank(outputFiles(jj).ses(kk).name(1, :));
        [pathstr fileName] = fileparts(str);
        underScoreIndex = icatb_findstr(fileName, '_');
        str = fileName(1:underScoreIndex(end) - 1);
        str = [num2str(jj),'-', num2str(kk), ' ', str];
        options(counter).str = str;
        counter = counter + 1;
    end
end

listStr = char(options.str);

function outputFiles = chkOutputFiles(outputFiles, outDir)
%% Truncate output files
%

setsToInclude = [];
for jj = 1:length(outputFiles)
    file = deblank(outputFiles(jj).ses(1).name(1, :));
    file = fullfile(outDir, file);
    if (~chkFile(file))
        continue;
    end
    setsToInclude = [setsToInclude, jj];
end

outputFiles = outputFiles(setsToInclude);

function status = chkFile(file)

[pathstr, fN, extn] = fileparts(deblank(file));

lastPos = icatb_findstr(fN, '_');
fN2 = fN(1:lastPos(end));

status = 0;
if (exist(fullfile(pathstr, [fN, extn]), 'file') || exist(fullfile(pathstr, [fN2, '.zip']), 'file'))
    status = 1;
end


function timecourses = interpTC(timecourses, TR, nSub, minTpLength, doInterp)

if (doInterp)
    % Interpolate timecourses if needed for variable TRs across
    % subjects
    if (~all(TR == min(TR)))
        interpFactor = TR(nSub)/min(TR);
        [num, denom] = rat(interpFactor);
        timecourses = resample(timecourses, num, denom);
    end
end
timecourses = timecourses(1:minTpLength, :);


function results_html(outputDir, files, outFile)
%% HTML report
%

helpStr = 'Creating HTML file. Please wait ...';
helpH = helpdlg(helpStr);
disp(helpStr);

html_file = fullfile(outputDir, outFile);
if (~exist(outputDir, 'dir'))
    mkdir(outputDir);
end
start_string = '<html><head><title> Component viewer </title></head>';
start_string = [start_string, '<p> </p><h1> Component viewer </h1><p><ul>'];

results_string1 = '';
for nR = 1:length(files);
    results_string1 = [results_string1, '<p align = "center"> <img src = "', files{nR}, '" height=600> </img>'];
end

end_string =  '</html>';

results_string = [start_string,  results_string1, end_string];

dlmwrite(html_file, results_string, '');

icatb_openHTMLHelpFile(html_file);

try
    delete(helpH);
catch
end

fprintf('Done\n');
