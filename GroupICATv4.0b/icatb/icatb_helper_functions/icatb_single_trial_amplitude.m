function icatb_single_trial_amplitude(param_file)
%% Calculate and plot single trial amplitudes. To use this utility
% you need to temporally sort components using the regressor of interest in Sorting
% GUI.
%
% Inputs:
% 1. param_file - Parameter file
%

%% Load defaults
icatb_defaults;
global PARAMETER_INFO_MAT_FILE;
global SMOOTHPARA;
global SMOOTHINGVALUE;
global EVENTAVG_WINDOW_SIZE;
global FONT_COLOR;
% FONT DEFAULTS
global UI_FONTNAME;
global UI_FONTUNITS;
global UI_FS;

%% Load parameter and validate parameter file
filterP = ['*', PARAMETER_INFO_MAT_FILE, '*.mat'];

if ~exist('param_file', 'var')
    param_file = icatb_selectEntry('typeEntity', 'file', 'title', 'Select parameter file for computing single trial amplitudes', 'filter', ...
        filterP);
end

if isempty(param_file)
    error('Error:ParameterFile', 'Parameter file is not selected\n');
end

drawnow;

load(param_file);

outputDir = fileparts(param_file);

if isempty(outputDir)
    outputDir = pwd;
end

cd(outputDir);

if ~exist('sesInfo', 'var')
    error('Error:ParameterFile', 'Selected file %s is not a valid parameter file\n', param_file);
end

if (~sesInfo.isInitialized)
    error('Please run the analysis in order to compute single trial amplitudes');
end

modalityType = icatb_get_modality;

if ~strcmpi(modalityType, 'fmri')
    error('This utility works with fMRI modality only');
end

%% Select structural image
structFile = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Structural Image To Overlay Components', 'filter', '*.img;*.nii', ...
    'typeselection', 'single', 'fileType', 'image', 'filenumbers', 1);

if isempty(structFile)
    error('Structural file is not selected');
end

drawnow;

dispParmeters.structFile = structFile;
dispParmeters = icatb_displayGUI_defaults('init', [], dispParmeters , 'off');

returnValue = dispParmeters.returnValue;
convertToZ = dispParmeters.convertToZ;
threshold = dispParmeters.thresholdvalue;

%% Get fields from sesInfo
output_prefix = sesInfo.userInput.prefix;
numOfSub = sesInfo.numOfSub;
numOfSess = sesInfo.numOfSess;
numComp = sesInfo.numComp;
calibrateMATFileName = sesInfo.calibrate_components_mat_file;
compFiles = str2mat(sesInfo.icaOutputFiles(1).ses.name);

compStr = [repmat('Component ', numComp, 1), num2str((1:numComp)')];

% Select components
selComp = icatb_listdlg('PromptString', 'Select components ...', 'SelectionMode', 'multiple', 'ListString', ...
    compStr, 'movegui', 'center', 'title_fig', 'Select components');
if (isempty(selComp))
    error('Components are not selected to view the single trial amplitudes');
end
% End for selecting components

[zipFileName, files_in_zip] = icatb_getViewingSet_zip(compFiles, [], 'real', sesInfo.zipContents);

compFiles = icatb_fullFile('directory', outputDir, 'files', compFiles);

% Temporal regression files
regressMatFile = fullfile(outputDir, [output_prefix, '_temporal_regression.mat']);
regressTxtFile = fullfile(outputDir, [output_prefix, '_temporal_regression.txt']);

% Single trial results file
singleTrialResultsFile = fullfile(outputDir, [output_prefix, '_single_trial_results.mat']);

if (~exist('calculate_single_trial_amp', 'var'))
    calculate_single_trial_amp = 1;
end

if (~exist(regressMatFile, 'file'))
    error('Error:RegressionParam', 'Regression parameters MAT file %s doesn''t exist.\nPlease do temporal sorting on all data-sets', regressMatFile);
end

drawnow;

load(regressMatFile);

if (~exist('regressInfo', 'var')) && (~isfield(regressInfo, 'regressionParameters'))
    error('Error:RegressionParam', 'Selected file %s\nis not a valid regression parameters file', regressMatFile);
end


if (~isfield(regressInfo, 'regressor_all_datasets'))
    %% Check if regressor names exist in MAT file
    
    %% Get regression names from text file
    
    if (~exist(regressTxtFile, 'file'))
        error('Error:RegressFile', 'Temporal regression parameters file %s doesn''t exist', regressTxtFile);
    end
    
    regressorNames = getRegressorsTextFile(regressTxtFile, numOfSub, numOfSess, numComp);
    
    regressInfo.regressor_all_datasets = regressorNames;
    
    icatb_save(regressMatFile, 'regressInfo');
    
end

regressorNames = regressInfo.regressor_all_datasets;

%% Number of conditions
numCond = length(regressorNames{1, 1});

checkNumDataSets = ceil(size(regressInfo.regressionParameters, 1)/numCond);

if (numOfSub*numOfSess ~= checkNumDataSets)
    error('Error:RegressionParam', 'You need to run temporal sorting on all data-sets to use this utility');
end

%% Sort components based on R-square information
componentNumbers = regressInfo.componentNumbers;
[componentNumbers, sort_ind] = sort(componentNumbers);
regressInfo.componentNumbers = componentNumbers;
regressInfo.regressionValues = regressInfo.regressionValues(sort_ind);
regressInfo.regressionParameters = regressInfo.regressionParameters(:, sort_ind);

%% Truncate the component numbers
regressInfo.componentNumbers = regressInfo.componentNumbers(selComp);
regressInfo.regressionValues = regressInfo.regressionValues(selComp);
regressInfo.regressionParameters = regressInfo.regressionParameters(:, selComp);
numSelComp = length(selComp);

if (calculate_single_trial_amp)
    
    subjectICAFiles = icatb_parseOutputFiles('icaOutputFiles', sesInfo.icaOutputFiles, 'numOfSub', sesInfo.numOfSub, 'numOfSess', sesInfo.numOfSess);
    
    disp('Calculating single trial amplitudes ...');
    fprintf('\n');
    
    %% Get best regressor for each dataset
    best_regress_names = get_best_regressor(regressInfo, numOfSub, numOfSess);
    
    %% Get SPM design matrix details
    spmMatFlag = sesInfo.userInput.spmMatFlag;
    spmFiles = cellstr(str2mat(sesInfo.userInput.designMatrix.name));
    
    loadSPMEachDataSet = 1;
    if strcmpi(spmMatFlag, 'same_sub_same_sess') || strcmpi(spmMatFlag, 'same_sub_diff_sess')
        loadSPMEachDataSet = 0;
        spmInfo = load(spmFiles{1});
    end
    
    % Initialise vars
    countDataSet = 0;
    
    %% Load Timecourses and get single trial amplitude and event averages
    % Loop over subjects
    for nSub = 1:numOfSub
        
        if (loadSPMEachDataSet)
            clear spmInfo;
            % Load for each subjec
            spmInfo = load(spmFiles{nSub});
        end
        
        % Time response
        TR = spmInfo.SPM.xY.RT;
        
        % Loop over sessions
        for nSess = 1:numOfSess
            
            countDataSet = countDataSet + 1;
            
            % Load calibrated component files
            tc = icatb_loadComp(sesInfo, selComp, 'subjects', nSub, 'sessions', nSess, 'vars_to_load', 'tc', 'subject_ica_files', subjectICAFiles);
            
            % smooth ica time course
            if strcmpi(SMOOTHPARA, 'yes')
                tc = icatb_gauss_smooth1D(tc, SMOOTHINGVALUE);
            end
            
            % Detrend timecourse
            tc = icatb_detrend(tc, 1, size(tc, 1));
            
            % Loop over components
            for nComp = 1:numSelComp
                selectedRegressor = best_regress_names{nSub, nSess, nComp};
                % get onsets
                [data_sessionNumber, onset_number] = icatb_get_onsets(selectedRegressor, spmInfo.SPM);
                if (data_sessionNumber > length(spmInfo.SPM.Sess))
                    data_sessionNumber = length(spmInfo.SPM.Sess);
                end
                %%%%%% Modifying selected onset
                if strcmpi(spmInfo.SPM.xBF.UNITS, 'scans')
                    selectedOnset = round(spmInfo.SPM.Sess(data_sessionNumber).U(onset_number).ons);
                else
                    selectedOnset = round((spmInfo.SPM.Sess(data_sessionNumber).U(onset_number).ons) / TR);
                end
                
                if (selectedOnset(1) == 0)
                    selectedOnset(1) = 1;
                end
                
                %% Compute single trial amplitudes and event average
                [temp_beta_weights, temp_eventAvg] = icatb_compute_single_trial_amp(tc(:, nComp), TR, EVENTAVG_WINDOW_SIZE, ...
                    selectedOnset);
                
                if ((countDataSet == 1) && (nComp == 1))
                    % Initialise event average and single trial amplitude terms
                    eventAvg = zeros(length(temp_eventAvg), numOfSub*numOfSess, numSelComp);
                    tmp_single_trial_amps = cell(numOfSub*numOfSess, numSelComp);
                end
                
                %% Store Single trial amplitude terms and event average
                tmp_single_trial_amps{countDataSet, nComp} = temp_beta_weights(:);
                eventAvg(:, countDataSet, nComp) = temp_eventAvg(:);
                clear temp_beta_weights temp_eventAvg;
                
            end
            % End loop over components
        end
        % End loop over sessions
    end
    % End loop over subjects
    
    
    %% Handle unequal onset numbers
    size_ons = cellfun('size', tmp_single_trial_amps, 1);
    min_length = min(size_ons(:));
    single_trial_amps = zeros(min_length, numOfSub*numOfSess, numSelComp);
    countD = 0;
    % Loop over subjects
    for nSub = 1:numOfSub
        % Loop over sessions
        for nSess = 1:numOfSess
            countD = countD + 1;
            % Loop over components
            for nComp = 1:numSelComp
                temp = tmp_single_trial_amps{countD, nComp};
                single_trial_amps(:, countD, nComp) = temp(1:min_length);
                clear temp
            end
            % End loop over components
        end
        % End loop over sessions
    end
    % End loop over subjects
    
    
    % Mean of single trial amplitudes
    if (numOfSub*numOfSess > 1)
        singleTrialResults.sem_single_amps = squeeze(std(single_trial_amps, 0, 2)./sqrt(numOfSub*numOfSess - 1));
        singleTrialResults.sem_event_avg = squeeze(std(eventAvg, 0, 2)./sqrt(numOfSub*numOfSess - 1));
    end
    
    %% Mean of event average and single trial amplitudes
    singleTrialResults.best_regressors = best_regress_names;
    singleTrialResults.single_trial_amps = squeeze(mean(single_trial_amps, 2));
    singleTrialResults.eventAvg = squeeze(mean(eventAvg, 2));
    singleTrialResults.sel_comp = selComp;
    
    clear single_trial_amps eventAvg;
    
    %% Save single trial results
    icatb_save(singleTrialResultsFile, 'singleTrialResults');
    
else
    
    %% Load existing results file
    load(singleTrialResultsFile, 'singleTrialResults');
    
end


%% Display results
% 1. Get colormap
% 2. Loop over each component
% 3. Resize component image to the anatomical image
% 4. Plot the slices at the maximum voxel
% 5. Plot event average and single trial amplitudes
%

msgStr = 'Interpolating components. Please wait ...';
disp(msgStr);
fprintf('\n');
drawnow;

%% Interpolate components to structural dimension
[icasig, icaTimecourse, structuralImage, coords, HInfo] = icatb_loadICAData('structFile', structFile, 'compFiles', compFiles, ...
    'slicePlane', 'axial', 'comp_numbers', selComp, ...
    'convertToZ', convertToZ, 'threshValue', threshold, 'zipfile', zipFileName, ...
    'files_in_zip_file', files_in_zip, 'open_dialog', 'no');

structDIM = HInfo.DIM;

%% Get the maximum voxel position
maxLoc = zeros(numSelComp, 3);
for nComp = 1:size(icasig, 1)
    temp = icasig(nComp, :);
    [max_value, tmp_max_loc] = max(abs(temp(:)));
    [maxLoc(nComp, 1), maxLoc(nComp, 2), maxLoc(nComp, 3)] = ind2sub(structDIM, tmp_max_loc);
    clear temp;
end

%% Overlay components
[icasig, maxICAIM, minICAIM, maxInterval, minInterval] = icatb_overlayComp(icasig, structuralImage, returnValue, numSelComp);

maxICAIM = round(maxICAIM*10)/10;
minICAIM = round(minICAIM*10)/10;

CLIM = [minInterval, 2*maxInterval];

%% Get colormap associated with the image values
cm = icatb_getColormap(1, returnValue, 1);

numRows = 3;

numFigures = ceil(numSelComp/numRows);

xOffset = 0.08; yOffset = 0.08;
mainFigAxesWidth = (1 - 4*xOffset)/3;
mainFigAxesHeight = (1 - (numRows + 1)*yOffset)/numRows;
mainFigAxesHeight = min([mainFigAxesWidth, mainFigAxesHeight]);
mainFigAxesWidth = mainFigAxesHeight;

startComp = 1;
GraphicsH = repmat(struct('H', []), 1, numFigures);
%% Loop over figures
for nFig = 1:numFigures
    endComp = numRows*nFig;
    
    if (endComp > numSelComp)
        endComp = numSelComp;
    end
    
    comps = (startComp:endComp);
    GraphicsH(nFig).H = icatb_getGraphics(['Single Trial Amplitudes ', num2str(nFig)], 'Graphics', ['single_trial_fig', num2str(nFig)], 'on');
    countC = 0;
    
    mainFigAxesPos = [xOffset, 1 - yOffset - mainFigAxesHeight, mainFigAxesWidth, mainFigAxesHeight];
    
    %% Loop over components
    for nComp = comps
        countC = countC + 1;
        mainFigAxesPos(1) = 0.01;
        im = reshape(icasig(nComp, :), structDIM);
        
        %% Slices (XY, XZ, YZ)
        sliceXY = reshape(im(:, :, maxLoc(nComp, 3)), size(im, 1), size(im, 2));
        sliceXZ = reshape(im(:, maxLoc(nComp, 2), :), size(im, 1),size(im, 3));
        sliceYZ = reshape(im(maxLoc(nComp, 1), :, :), size(im, 2), size(im, 3));
        
        %% Plot Image first
        subH = axes('parent', GraphicsH(nFig).H, 'units', 'normalized', 'position', mainFigAxesPos, 'fontname', UI_FONTNAME, 'fontsize', ...
            UI_FS, 'fontunits', UI_FONTUNITS);
        compStr = icatb_returnFileIndex(selComp(nComp));
        title(['Comp ', compStr], 'parent', subH);
        axis(subH, 'off');
        
        imageAxesWidth = mainFigAxesPos(3)/2;
        imageAxisHeight = mainFigAxesPos(4)/2;
        
        %% Slice XZ
        imagePos = [mainFigAxesPos(1), mainFigAxesPos(2) + imageAxesWidth, imageAxesWidth, imageAxisHeight];
        plotImage(GraphicsH(nFig).H, sliceXZ, CLIM, imagePos);
        
        %% Slice YZ
        imagePos(1) = imagePos(1) + 0.01 + imagePos(3);
        plotImage(GraphicsH(nFig).H, sliceYZ, CLIM, imagePos);
        
        %% Slice XY
        imagePos(1) = mainFigAxesPos(1);
        imagePos(2) = imagePos(2) - 0.015 - imageAxesWidth;
        subH = plotImage(GraphicsH(nFig).H, sliceXY, CLIM, imagePos);
        colormap(cm);
        
        %% Plot colorbar
        colorbarPos = get(subH, 'position');
        colorbarPos(1) = colorbarPos(1) + colorbarPos(3) + 0.015;
        colorbarPos(3) = 0.02;
        colorbarH = colorbar('peer', subH);
        if ((icatb_get_matlab_version <= 2013) || strcmpi(version('-release'), '2014a'))
            set(colorbarH, 'fontname', UI_FONTNAME, 'fontsize', UI_FS - 2, 'fontunits', UI_FONTUNITS);
        else
            labelsH = get(colorbarH, 'Label');
            set(labelsH, 'fontname', UI_FONTNAME, 'fontsize', UI_FS - 2, 'fontunits', UI_FONTUNITS);
        end
        set(colorbarH, 'units', 'normalized');
        set(colorbarH, 'position', colorbarPos);
        set(colorbarH, 'YLim', [minInterval, maxInterval]);
        YTicks = get(colorbarH, 'YTick');
        set(colorbarH, 'XTick', []);
        set(colorbarH, 'YTick', []);
        set(colorbarH, 'XTickLabel', []);
        set(colorbarH, 'YTickLabel', []);
        if ((icatb_get_matlab_version <= 2013) || strcmpi(version('-release'), '2014a'))
            text(1.2, -0.01, num2str(minICAIM(nComp)), 'units', 'normalized', 'parent', colorbarH);
            text(1.2, 0.99, num2str(maxICAIM(nComp)), 'units', 'normalized', 'parent', colorbarH);
        else
            set(colorbarH, 'color', FONT_COLOR);
            set(colorbarH, 'YTick', [YTicks(1), YTicks(end)]);
            set(colorbarH, 'YTickLabel', char(num2str(minICAIM(nComp)), num2str(maxICAIM(nComp))));
        end
        
        %% Second column (Event Average)
        countC = countC + 1;
        mainFigAxesPos(1) = mainFigAxesPos(1) + xOffset + mainFigAxesPos(3);
        subH = axes('parent', GraphicsH(nFig).H, 'units', 'normalized', 'position', mainFigAxesPos, 'fontname', UI_FONTNAME, 'fontsize', ...
            UI_FS - 1);
        plot(singleTrialResults.eventAvg(:, nComp), 'm', 'parent', subH);
        if isfield(singleTrialResults, 'sem_event_avg')
            hold on;
            plot(singleTrialResults.eventAvg(:, nComp) + singleTrialResults.sem_event_avg(:, nComp), 'g:', 'parent', subH);
            hold on;
            plot(singleTrialResults.eventAvg(:, nComp) - singleTrialResults.sem_event_avg(:, nComp), 'g:', 'parent', subH);
            hold off;
            title('Event Average (Mean +/- SEM)', 'Parent', subH);
        else
            title('Event Average', 'Parent', subH);
        end
        xlabel('Window Size (Seconds)', 'Parent', subH);
        ylabel('Z-scores', 'Parent', subH);
        set(subH, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR);
        axis(subH, 'tight');
        
        %% Third column (Single Trial amps)
        countC = countC + 1;
        mainFigAxesPos(1) = mainFigAxesPos(1) + xOffset + mainFigAxesPos(3);
        subH = axes('parent', GraphicsH(nFig).H, 'units', 'normalized', 'position', mainFigAxesPos, 'fontname', UI_FONTNAME, 'fontsize', UI_FS - 1);
        plot(singleTrialResults.single_trial_amps(:, nComp), 'm', 'parent', subH);
        if isfield(singleTrialResults, 'sem_single_amps')
            hold on;
            plot(singleTrialResults.single_trial_amps(:, nComp) + singleTrialResults.sem_single_amps(:, nComp), 'g:', 'parent', subH);
            hold on;
            plot(singleTrialResults.single_trial_amps(:, nComp) - singleTrialResults.sem_single_amps(:, nComp), 'g:', 'parent', subH);
            hold off;
            title('Single Trial Amplitudes (Mean +/- SEM)', 'Parent', subH);
        else
            title('Single Trial Amplitudes', 'Parent', subH);
        end
        ylabel('Amplitude (\beta)', 'Parent', subH, 'interpreter', 'tex');
        xlabel('Generic Trial No.', 'parent', subH);
        set(subH, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR);
        axis(subH, 'tight');
        
        % Decrement y position
        mainFigAxesPos(2) = mainFigAxesPos(2) - yOffset - mainFigAxesHeight;
        
    end
    %% End loop over components
    startComp = endComp + 1;
end
%% End loop over figures


%% Plot prev, next and exit buttons
icatb_plotNextPreviousExitButtons(GraphicsH);

fprintf('Done plotting single trial amplitudes\n');


function best_regress_names = get_best_regressor(regressInfo, numOfSub, numOfSess)
%% Return the best regressor for each dataset and component
%

numCond = length(regressInfo.regressor_all_datasets{1, 1}); % Number of conditions
numComp = size(regressInfo.regressionParameters, 2); % Number of components

component_num = regressInfo.componentNumbers;
[component_num, sortInd] = sort(component_num);
regressionParameters = regressInfo.regressionParameters(:, sortInd);

% Initialise variable for storing best regressor for each data-set
best_regress_names = cell(numOfSub, numOfSess, numComp);

% Initialise vars
startT = 1;
endT = 0;

% Loop over subjects
for nSub = 1:numOfSub
    % Loop over sessions
    for nSess = 1:numOfSess
        endT = endT + numCond;
        indices = (startT:endT);
        % Loop over components
        for nComp = 1:numComp
            values = regressionParameters(indices, nComp);
            [val, ind] = max(values);
            temp = regressInfo.regressor_all_datasets{nSub, nSess};
            best_regress_names{nSub, nSess, nComp} = temp{ind};
        end
        % End loop over components
        startT = endT + 1;
    end
    % End loop over sessions
end
% End loop over subjects


function regressorNames = getRegressorsTextFile(regressTxtFile, numOfSub, numOfSess, numComp)
%% Get regressors from temporal sort text file
%

% Get only the first column
formatStr = ['%s', repmat('%*s', 1, numComp)];

all_lines = textread(regressTxtFile, formatStr, 'delimiter', '\t', 'whitespace', '');

%% Remove blank spaces
good_ind = icatb_good_cells(all_lines);
all_lines = all_lines(good_ind);

%% Check if the temporal sort is done on all data-sets
startInd = regexpi(all_lines, '^temporal sorting');

startInd = find(icatb_good_cells(startInd) ~= 0);

if isempty(startInd)
    error('You need to run temporal sorting on all data-sets');
end

checkAll = findstr(lower(all_lines{startInd(end)}), 'all data-sets');
if isempty(checkAll)
    error('You need to run temporal sorting on all data-sets');
end

%% Store only necessary lines
all_lines = all_lines(startInd(end) + 1:end);

all_lines = regexprep(all_lines, 'Subject \d+ Session \d+\s|Subject \d+\s', '');

% Match regressors
startInd = regexpi(all_lines, '.*Sn(');
good_ind = icatb_good_cells(startInd);

%% Truncate lines
all_lines = deblank(all_lines(good_ind));

numCond = ceil(length(all_lines)/(numOfSub*numOfSess));

% Initialise regressor names
regressorNames = cell(numOfSub, numOfSess);

startT = 1;
endT = 0;
count = 0;
% Loop over subjects
for nSub = 1:numOfSub
    % Loop over sessions
    for nSess = 1:numOfSess
        count = count + 1;
        endT = endT + numCond;
        regressorNames{nSub, nSess} = all_lines(startT:endT);
        startT = endT + 1;
    end
    % End loop over sessions
end
% End loop over subjects

function setImagePos(subHandle, imageAxisPos)
%% Set image axes position
%

yAxisRatio = imageAxisPos(4)/imageAxisPos(3);
xAxisRatio = imageAxisPos(3)/imageAxisPos(4);
if(yAxisRatio>1)
    yAxisRatio = 1;
else
    xAxisRatio = 1;
end

imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*yAxisRatio imageAxisPos(4)*xAxisRatio];

set(subHandle, 'position', imageAxisPos);


function subH = plotImage(graphicsH, data, CLIM, imagePos)
%% Function to plot the image at the specified position
%

icatb_defaults;
global UI_FONTNAME;
global UI_FONTUNITS;
global UI_FS;


subH = axes('parent', graphicsH, 'units', 'normalized', 'position', imagePos, 'fontname', UI_FONTNAME, 'fontsize', UI_FS, ...
    'fontunits', UI_FONTUNITS);
image(rot90(data), 'parent', subH, 'CDataMapping', 'scaled');
setImagePos(subH, imagePos);
set(subH, 'clim', CLIM); % set the axis positions to the specified
axis(subH, 'off');
