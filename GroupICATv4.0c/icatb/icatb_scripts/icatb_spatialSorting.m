function icatb_spatialSorting(param_file, selectedStr, templateFiles, dispParameters)
% Function to do spatial sorting of all data-sets. See
% icatb_example_spatial_sorting for running through batch mode

try
    icatb_defaults;
    global PARAMETER_INFO_MAT_FILE;
    %% Pass parameters for display
    global DETRENDNUMBER;
    global SMOOTHPARA;
    global SMOOTHINGVALUE;

    if ~exist('param_file', 'var')
        param_file = icatb_selectEntry('typeSelection', 'single', 'typeEntity', 'file', 'title', 'Select a valid ICA parameter file', ...
            'filter', ['*', PARAMETER_INFO_MAT_FILE, '*']);
    end

    load(param_file);

    if ~exist('sesInfo', 'var')
        error(['file: ', param_file, ' is not a valid ica parameter file']);
    end

    % Get the output directory
    [outputDir, fileName, extn] = fileparts(param_file);

    if isempty(outputDir)
        outputDir = pwd;
    end

    cd(outputDir);


    %%%%%%%%%% Get the required variables from sesInfo structure %%%%%%%%%%
    % Number of subjects
    numOfSub = sesInfo.numOfSub;
    numOfSess = sesInfo.numOfSess;

    % Number of components
    numComp = sesInfo.numComp;

    dataType = sesInfo.dataType;

    mask_ind = sesInfo.mask_ind;

    % First scan
    structFile = deblank(sesInfo.inputFiles(1).name(1, :));

    structVol = icatb_get_vol_nifti(structFile);

    DIM = structVol(1).dim(1:3);

    % slices in mm
    [parameters] = icatb_get_slice_def(structVol, 'axial');
    slices_in_mm = parameters.slices; clear parameters;

    zipContents.zipFiles = {};
    zipContents.files_in_zip(1).name = {};

    if isfield(sesInfo, 'zipContents')
        zipContents = sesInfo.zipContents;
    end

    %%%%%%%% End for getting the required vars from sesInfo %%%%%%%%%%


    %%%%%%%%%%%%% Get template data %%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('selectedStr', 'var')

        selectedStr = 'Same set of spatial templates for all data-sets';

        questionString = 'How do you want to select the spatial templates?';
        D(1).string = 'Same set of spatial templates for all data-sets';
        D(length(D) + 1).string = 'Different set of spatial templates for sessions';
        D(length(D) + 1).string = 'Different set of spatial templates for subjects and sessions';
        choiceString = str2mat(D.string);
        clear D;

        % Modify the choice string
        if numOfSub ==  1 & numOfSess > 1
            choiceString = choiceString(1:2, :);
        end

        if numOfSub*numOfSess > 1
            InputHandle = icatb_getGraphics('Selecting Spatial Template', 'normal', 'Spatial Template'); % figure handle
            answerSpatialTemplate = icatb_promptUI('popup', questionString, choiceString, 'numeric', InputHandle);
            close(InputHandle);
            selectedStr = deblank(choiceString(answerSpatialTemplate, :));
        end

        % Select spatial templates
        if strcmpi(selectedStr, 'different set of spatial templates for sessions')

            templateFiles = icatb_select_data('title', 'Select spatial templates', 'num_subjects', ...
                1, 'num_sessions', numOfSess, 'files_specification', 'equal', 'spm_check', 'no', ...
                'filter_string', '*.img', 'type_file_selection', 'multiple', 'figure_menu', 'data', ...
                'datatype', 'real');

            titlePrint = ['Spatial sorting of all data-sets using different templates for sessions'];

        elseif strcmpi(selectedStr, 'different set of spatial templates for subjects and sessions')

            templateFiles = icatb_select_data('title', 'Select spatial templates', 'num_subjects', ...
                numOfSub, 'num_sessions', numOfSess, 'files_specification', 'equal', 'spm_check', 'no', ...
                'filter_string', '*.img', 'type_file_selection', 'multiple', 'figure_menu', 'data', ...
                'datatype', 'real');

            titlePrint = ['Spatial sorting of all data-sets using different templates for subjects and sessions'];

        else

            % Select image files
            templateFiles.name = icatb_selectEntry('typeSelection', 'multiple', 'typeEntity', 'file', 'title', 'Select template/templates for spatial sorting', ...
                'filter', '*.img', 'filetype', 'image', 'filenumbers', 1);
            titlePrint = ['Spatial sorting of all data-sets using same templates for all subjects and sessions'];

        end

    else

        % do checking for template files
        if ~exist('templateFiles', 'var')
            error(['templateFiles variable is not passed']);
        end

        if ~isstruct(templateFiles)
            error('templateFiles variable must be a structure');
        end

        if strcmpi(selectedStr, 'different set of spatial templates for sessions')
            if length(templateFiles) < numOfSess
                error(['Check the templateFiles variable. It should have the length of the number of sessions']);
            else
                templateFiles = templateFiles(1:numOfSess);
            end
            titlePrint = ['Spatial sorting of all data-sets using different templates for sessions'];

        elseif  strcmpi(selectedStr, 'different set of spatial templates for subjects and sessions')
            if length(templateFiles) < numOfSub*numOfSess
                error(['Check the templateFiles variable. It should have the length of the number of data-sets']);
            else
                templateFiles = templateFiles(1:numOfSub*numOfSess);
            end
            titlePrint = ['Spatial sorting of all data-sets using different templates for subjects and sessions'];
        else

            templateFiles = templateFiles(1);
            titlePrint = ['Spatial sorting of all data-sets using same templates for all subjects and sessions'];
        end

    end

    drawnow;

    helpMsg = ['Loading template files ...'];
    disp(helpMsg);

    [helpHandle] = helpdlg(helpMsg, helpMsg);

    startTp = 1; endTp = 0;
    for nTemplates = 1:length(templateFiles)
        endTp = endTp + length(mask_ind);
        % current template
        currentTemplate = deblank(templateFiles(nTemplates).name);
        tempData = icatb_loadData(currentTemplate);
        tempData = permute(tempData, [4 1 2 3]);
        clear structuralImage;
        for nn = 1:size(currentTemplate, 1)
            currentData = squeeze(tempData(nn, :, :, :));
            currentData = currentData(:);
            spatialTemplates(startTp:endTp, nn) = currentData(mask_ind);
            clear currentData;
        end
        startTp = endTp + 1;
        clear tempData;
    end

    % Number of regressors
    numRegressors = size(spatialTemplates, 2);

    if strcmpi(selectedStr, 'same set of spatial templates for all data-sets')
        % Replicate the template over subjects and sessions
        spatialTemplates = repmat(spatialTemplates, numOfSub*numOfSess, 1);
    elseif strcmpi(selectedStr, 'different set of spatial templates for sessions')
        % Replicate only over session
        spatialTemplates = repmat(spatialTemplates, numOfSub, 1);
    end

    %%%%%%%%%%%%% End for getting template data %%%%%%%%%%%%%%%%%%%%%%%%%%


    drawnow;

    %%%%%%%%%%%%%%%%%%%%%%% Get component data %%%%%%%%%%%%%%%%%

    % Get the ICA Output files
    icaOutputFiles = sesInfo.icaOutputFiles;

    [subjectICAFiles, meanICAFiles, tmapICAFiles, meanALL_ICAFile] = icatb_parseOutputFiles('icaOutputFiles', icaOutputFiles, 'numOfSub', ...
        numOfSub, 'numOfSess', numOfSess, 'flagTimePoints', sesInfo.flagTimePoints);

    try
        delete(helpHandle);
    catch
    end

    helpMsg = 'Loading components...';
    [helpHandle] = helpdlg(helpMsg, helpMsg);

    % component data
    compData = zeros(numComp, length(mask_ind)*numOfSub*numOfSess);
    startPoint = 1; endPoint = 0;
    for nSub = 1:numOfSub

        for nSess = 1:numOfSess
            endPoint = endPoint + length(mask_ind);
            disp(['Loading Subject ', num2str(nSub), ' Session ', num2str(nSess), '...']);
            compFiles = subjectICAFiles(nSub).ses(nSess).name;
            % component files
            [zipFileName, files_in_zip] = icatb_getViewingSet_zip(compFiles, [], 'real', zipContents);

            compFiles = icatb_fullFile('directory', outputDir, 'files', compFiles);

            % Unzip files
            if ~isempty(zipFileName)
                icatb_unzip(zipFileName, outputDir);
            end

            drawnow;

            [icasig, HInfo, real_world_coords] = icatb_loadData(compFiles, 'real', [], [], [1:numComp]);

            if ~isempty(zipFileName)
                icatb_delete_file_pattern(files_in_zip, outputDir);
            end

            % Reshape icasig to components by voxels
            icasig = permute(icasig, [4 1 2 3]);

            % Structural volume
            HInfo.V = HInfo.V(1);

            % Reshape to 2d
            icasig = reshape(icasig, [numComp, prod(HInfo.DIM(1:3))]);

            compData(:, startPoint:endPoint) = icasig(:, mask_ind);
            clear icasig;
            startPoint = endPoint + 1;
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%% End for getting component data %%%%%%%%%%%%%%%%%%


    drawnow;

    %%%%%%%% Calculate Multiple Regression %%%%%%%%
    try
        delete(helpHandle);
    catch
    end

    helpMsg = 'Calculating Regression ...';
    [helpHandle] = helpdlg(helpMsg, helpMsg);

    componentValues = zeros(1, numComp);
    compDIMS = repmat(length(mask_ind), 1, numOfSub*numOfSess);
    % Multiple regression
    for nComp = 1:numComp
        disp(['Calculating Multiple regression for component ', num2str(nComp)]);
        % calculate regression
        [comparison, regressCoeff, ModelIndices] = icatb_multipleRegression(spatialTemplates, compData(nComp, :), ...
            numRegressors, numOfSub*numOfSess, compDIMS, 0);
        % store values from comparison
        componentValues(nComp) = comparison;
        regressionCoeff{nComp} = regressCoeff(ModelIndices);
        clear  regressCoeff;
    end

    clear compData; clear spatialTemplates;

    %%%%%%%% End for calculating Multiple Regression %%%%%%%%%%%


    %%%%%%%%%%%%% Sort Components %%%%%%%%%%%%%%%%%%%%

    % sort the components based on the criteria selected
    [sortedValues, sortedIndices] = sort(componentValues);

    clear componentValues;

    % transpose to column vector
    if size(sortedValues, 1) == 1
        sortedValues = sortedValues';
        sortedIndices = sortedIndices';
    end

    sortedValues = flipud(sortedValues); sortedIndices = flipud(sortedIndices); % flip the vector along rows

    temp = cell(1, length(regressionCoeff));
    for ii = 1:length(temp)
        temp{ii} = regressionCoeff{sortedIndices(ii)};
    end
    clear regressionCoeff;
    regressionCoeff = temp;
    clear temp;

    %%%%%%%% End for sorting %%%%%%%%%%%%


    %%%%%%%%%%%%%%%% Print regression parameters to a file %%%%%%%%%%%%%

    for nTemp = 1:length(templateFiles)

        for nn = 1:size(templateFiles(nTemp).name, 1)
            [pathstr, tempF(nn).name, extn] = fileparts(deblank(templateFiles(nTemp).name(nn, :)));
        end

        templateFiles(nTemp).name = str2mat(tempF.name);
        clear tempF;
    end


    numPara = 1;
    varStruct(numPara).tag = 'Component Number';
    varStruct(numPara).value = sortedIndices;

    numPara = numPara + 1;
    varStruct(numPara).tag = 'Multiple Regression';
    varStruct(numPara).value = sortedValues;

    countDataSets = 0;
    regressionParameters = cell2mat(regressionCoeff);

    % loop over all subjects
    for numgetSub = 1:numOfSub
        % Loop over the sessions involved
        for numgetSess = 1:numOfSess
            countDataSets = countDataSets + 1;
            tagPrefix = ['Subject ', num2str(numgetSub), ' Session ', num2str(numgetSess)];

            nVarStruct = length(varStruct);

            tempRegress = regressionParameters((countDataSets - 1)*numRegressors + 1:countDataSets*numRegressors, :)';

            % Get the appropriate regressors
            if strcmpi(selectedStr, 'same set of spatial templates for all data-sets')
                regressorNames = str2mat(templateFiles.name);
            elseif strcmpi(selectedStr, 'different set of spatial templates for sessions')
                regressorNames = str2mat(templateFiles(numgetSess).name);
            else
                regressorNames = str2mat(templateFiles(countDataSets).name);
            end


            for nPrint = 1:size(regressorNames, 1)
                nVarStruct = nVarStruct + 1;
                varStruct(nVarStruct).tag = [tagPrefix, ' ', deblank(regressorNames(nPrint, :))];
                varStruct(nVarStruct).value = tempRegress(:, nPrint);
            end

            clear tempRegress;
        end
    end

    outFileName = fullfile(outputDir, ['spatial_regression_sessions.txt']);

    msgStr = ['Regression parameters with the names of the spatial templates for the corresponding data set/sets are stored in ', ...
        outFileName];

    % print all the information to a text file
    icatb_printToFile(outFileName, varStruct, titlePrint, 'row_wise', 'append');
    disp(msgStr);

    %%%%%%%% End for printing regression parameters to a file %%%%%%%%%%%%%

    try
        delete(helpHandle);
    catch
    end

    %%%%%%%%%%%%%% Display Parameters %%%%%%%%%%%%%%%

    % 1 means positive and negative
    % 2 means positive
    % 3 means Absolute
    % 4 means Negative
    parameters.imagevalues = 1;
    % Anatomical plane like axial, sagital, coronal
    parameters.anatomicalplane = 'axial';
    % slices in mm (vector of real numbers)
    parameters.slicerange = slices_in_mm;
    % Number of images per figure
    parameters.imagesperfigure = 4;
    % convert to z scores
    parameters.convertToZ = 1;
    % Z Threshold
    parameters.thresholdvalue = 1;
    % Anatomical file required to display
    parameters.structFile = structFile;

    % Get input from the user: else display by default
    if exist('dispParameters', 'var')
        % display parameters
        varsNeeded = {'anatomicalplane', 'imagevalues', 'slicerange', 'imagesperfigure', 'convertToZ', 'thresholdvalue', ...
            'structFile'};
        for ii = 1:length(varsNeeded)
            if isfield(dispParameters, varsNeeded{ii})
                parameters = setfield(parameters, varsNeeded{ii}, getfield(dispParameters, varsNeeded{ii}));
            end
        end

    else
        return;

    end

    %%%%%%%%%% End for defining display parameters %%%%%%%



    %%%%%%%%%%%%%%%% Display Components %%%%%%%%%%%%%%%%%%%%%

    parameters.numComp = numComp; % Number of components
    parameters.compFiles = meanALL_ICAFile(1).name; % component files
    parameters.compNumbers = sortedIndices; % component numbers

    [zipFileName, files_in_zip] = icatb_getViewingSet_zip(parameters.compFiles, [], 'real', zipContents);
    parameters.compFiles = icatb_fullFile('directory', outputDir, 'files', parameters.compFiles);
    % load the component files here
    [icasig, icaTimecourse, structuralImage, coords, HInfo, parameters.text_left_right] = icatb_loadICAData('structFile', ...
        parameters.structFile, 'compFiles', parameters.compFiles, 'slicePlane', parameters.anatomicalplane, 'sliceRange', ...
        parameters.slicerange, 'comp_numbers', sortedIndices, 'convertToZ', parameters.convertToZ, 'returnValue', parameters.imagevalues, ...
        'threshValue', parameters.thresholdvalue, 'dataType', 'real', 'complexInfo', [], 'zipfile', zipFileName, 'files_in_zip_file', ...
        files_in_zip);
    %%%%%%%%% Apply structural image parameters %%%%%%%%


    % smooth ica time course
    if strcmpi(SMOOTHPARA, 'yes')
        icaTimecourse = icatb_gauss_smooth1D(icaTimecourse, SMOOTHINGVALUE);
    end

    structHInfo = HInfo;
    parameters.structHInfo = HInfo;
    structDIM = HInfo.DIM;
    parameters.structuralImage = structuralImage;
    parameters.filesOutputDir = outputDir;
    parameters.icasig = icasig;
    parameters.undetrendICA = icaTimecourse; % show the user detrend and undetrend
    parameters.icaTimecourse = icatb_detrend(icaTimecourse, 1, size(icaTimecourse, 1), DETRENDNUMBER); % detrended time course
    parameters.modelTimecourse = [];

    underScorePos = icatb_findstr(deblank(meanALL_ICAFile(1).name(1, :)), '_');
    parameters.figLabel = meanALL_ICAFile(1).name(1, 1:underScorePos(end));
    parameters.htmlFile = 'icatb_component_explorer.htm';
    for ii = 1:numComp
        compLabels(ii).string = ['Component ', num2str(sortedIndices(ii)) , ' ', 'Spatial regression = ', ...
            num2str(sortedValues(ii))];
    end
    parameters.compLabels = compLabels;
    clear prefix;
    % specify a flag here that the component explorer is accessed
    % independently
    flagComponentExplorer = 'no_displayGUI';

    icatb_componentExplore(parameters);

    %%%%%%%%%% End for displaying components %%%%%%%%%%%%%%%%%%%%

catch

    if exist('helpHandle', 'var')
        if ishandle(helpHandle)
            delete(helpHandle);
        end
    end

    icatb_displayErrorMsg;

end
