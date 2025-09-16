function [sortParameters] = icatb_sortComponents(varargin)
% sorts the components based on the time course or spatial
% template selected
%
% Input: inputParameters structure containing the input information. The
% fields are given below:
% 1. sortingCriteria - Currently there are four options like 'multiple
% regression', 'correlation', 'kurtosis', 'maximum voxel', add other
% criteria as needed
% 2. sortingType - options are 'temporal' and 'spatial'
% 3. icaTimecourse - concatenated ICA time course
% 4. modelTimecourse - model time course is the selected regressor during
% sorting. Model time course is not needed for kurtosis sorting criteria
% and therefore, modelTimecourse for kurtosis is [].
% 5. icasig - component images of dimension components by voxels to be
% spatially sorted. In case of temporal sorting icasig = [];
% 6. spmFileNames - SPM.mat file names selected during temporal sorting.
% 7. numSubjects - number of subjects involved in concatenation
% 8. numSessions - number of sessions involved in concatenation
% 9. numComp - number of components
% 10. diffTimePoints - store the lengths of the time courses
% 11. num_Regress - number of regressors used for sorting temporally
% 12. num_DataSets - number of data sets used to concatenate time courses
%
% Output:
% sortParameters structure containing sorted indices, sorted criteria values;
% regression coefficients are passed only for Multiple Regression
% detrended ICA time course for all components if specified, regressor
% names, SPM.mat file name, component images
% 1. sortedIndices - sorted indices
% 2. sortedValues - Values in descending order
% 3. sortingCriteria - sorting criteria as passed in input Parameter
% 4. sortingType - spatial or temporal
% 4. regressCoeff - regression coefficients only for multiple regression
% 5. icaTimecourse - sorted icaTimecourse
% 6. modelTimecourse - model time course selected. In case no model
% timecourse is selected modelTimecourse = [];
% 7. refInfo - reference information selected includes the number of
% subjects, sessions, regressor indices , regressor names selected. In case
% no regressor is selected refInfo = [];
% 8. icasig - component images (components by voxels). In case of temporal
% sorting icasig = [];
% 9. diffTimePoints - length of the time points for each data set used to
% concatenate. In case of spatial sorting diffTimePoints = [];

structHInfo = [];
spatialImage = '';

% loop over the n arguments
for ii = 1:2:nargin
    if strcmp(lower(varargin{ii}), 'sortingcriteria')
        sortingCriteria = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'sortingtype')
        sortingType = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'icatimecourse')
        icaTimecourse = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'modeltimecourse')
        modelTimecourse = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'icasig')
        icasig = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'spmfilenames')
        spmFileNames = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'numcomp')
        numComp = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'spatialtemplate')
        spatialTemplate = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'refinfo')
        refInfo = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'difftimepoints')
        diffTimePoints = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'num_regress')
        num_Regress = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'num_datasets')
        num_DataSets = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'structhinfo')
        structHInfo = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'viewingset')
        viewingSet = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'spatialimage')
        spatialImage = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'num_sort_subjects')
        getSubjects = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'num_sort_sessions')
        getSessions = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'output_dir')
        outputDir = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'input_prefix')
        inputPrefix = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'real_world_coords')
        real_world_coords = varargin{ii + 1};
    end
end

% load defaults
icatb_defaults;
global SMOOTHINGVALUE;
global SMOOTHPARA;
global PRINTTYPE_REGRESSORS;

maxVoxelPositions = [];

% smooth the time courses
icaTimecourse = icatb_smoothTimecourse(icaTimecourse, diffTimePoints, num_DataSets, numComp, ...
    SMOOTHPARA, SMOOTHINGVALUE);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Temporal Sorting %%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(lower(sortingType), 'temporal')
    switch lower(sortingCriteria)
        case 'multiple regression'
            t_begin = cputime;
            % loop over components
            for nComp = 1:numComp
                disp(['Calculating regression for component  ', num2str(nComp)]);
                [comparison, regressCoeff, ModelIndices, otherIndices, linearRegress, removeTrend, ...
                    icaTimecourse(:, nComp), sub_partial_corr, partialCorrSlopes] = ...
                    icatb_multipleRegression(modelTimecourse, icaTimecourse(:, nComp), num_Regress, num_DataSets, diffTimePoints);
                componentValues(nComp) = comparison;
                regressionCoeff{nComp} = regressCoeff;
                regressionLineFit{nComp} = linearRegress;
                subject_partial_corr{nComp} = sub_partial_corr;
                subject_partial_slopes{nComp} = partialCorrSlopes;
                clear sub_partial_corr;
                clear group_par_corr;
                clear regressCoeff;
                clear linearRegress;
                clear removeTrend;
            end
            t_end = cputime;
            %disp(['Time taken to compute regression is ', num2str(t_end - t_begin), ' seconds']);
            %sortParameters.regressionCoeff = regressionCoeff; % regression coefficients
        case 'correlation'
            % loop over components
            for nComp = 1:numComp
                [comparison, icaTimecourse(:, nComp), regressCoeff, ModelIndices] = icatb_correlateFunctions(modelTimecourse, icaTimecourse(:, nComp), ...
                    num_DataSets, diffTimePoints);
                componentValues(nComp) = comparison;
                regressionCoeff{nComp} = regressCoeff;
                clear  regressCoeff;
            end
        case 'kurtosis'
            % loop over components
            for nComp = 1:numComp
                [comparison, icaTimecourse(:, nComp)] = icatb_kurtosis(icaTimecourse(:, nComp), num_DataSets, diffTimePoints);
                componentValues(nComp) = comparison;
            end
    end

    %%%%%%%%%%%%%%%%%% spatial sorting %%%%%%%%%%%%%%%%%%
elseif strcmp(lower(sortingType), 'spatial')

    switch sortingCriteria
        case 'multiple regression'

            % Multiple regression
            for nComp = 1:numComp
                compMap = icasig(nComp, :);
                [comparison, regressCoeff, ModelIndices] = icatb_multipleRegression(spatialTemplate, compMap);
                % store values from comparison
                componentValues(nComp) = comparison;
                regressionCoeff{nComp} = regressCoeff;
                clear  regressCoeff;
            end

        case 'correlation'

            % Spatial Comparison
            for nComp = 1:numComp
                compMap = icasig(nComp, :);
                [comparison, tempmap, regressCoeff, ModelIndices] = icatb_correlateFunctions(spatialTemplate, compMap);
                clear tempmap;
                %store values from comparison
                componentValues(nComp) = comparison;
                regressionCoeff{nComp} = regressCoeff;
                clear  regressCoeff;

            end

        case 'kurtosis'

            for nComp = 1:numComp
                compMap = icasig(nComp, :);
                [comparison]= icatb_kurtosis(compMap);
                componentValues(nComp) = comparison;
            end

        case 'maximum voxel'
            % sort the components based on the maximum voxel criteria
            % Find the non-zero values and take the maximum of all the
            % voxels in the component

            maxVoxelPositions = zeros(numComp, 1);
            % Maximum voxel criteria
            for nComp = 1:numComp
                compMap = icasig(nComp, :);
                % Replace with the maximum voxel criteria
                [comparison, maxVoxelPos]= icatb_maxVoxelCriteria(spatialTemplate, compMap);
                componentValues(nComp) = comparison;
                maxVoxelPositions(nComp) = maxVoxelPos;
                %maxVoxelPositions(nComp) = real_world_coords(maxVoxelPos);
            end
    end
else
    error('select temporal or spatial sorting option');
end

% sort the components based on the criteria selected
[sortedValues sortedIndices] = sort(componentValues);

clear componentValues;

% transpose to column vector
if size(sortedValues, 1) == 1
    sortedValues = sortedValues';
    sortedIndices = sortedIndices';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% flip the corresponding sorting parameters %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables to flipped include sorting values, sorting indices, regression
% coeff, ICA timecourse
sortedValues = flipud(sortedValues); sortedIndices = flipud(sortedIndices); % flip the vector along rows

% regression coefficients
if exist('regressionCoeff', 'var')
    temp = cell(1, length(regressionCoeff));
    for ii = 1:length(temp)
        temp{ii} = regressionCoeff{sortedIndices(ii)};
    end
    clear regressionCoeff;
    regressionCoeff = temp;
    clear temp;
    sortParameters.regressionCoeff = regressionCoeff; % flipped regression coefficients
    clear regressionCoeff;
end

% regression line fit
if exist('regressionLineFit', 'var')
    temp = cell(1, length(regressionLineFit));
    for ii = 1:length(temp)
        temp{ii} = regressionLineFit{sortedIndices(ii)};
    end
    clear regressionLineFit;
    regressionLineFit = temp;
    clear temp;
    sortParameters.regressionLineFit = regressionLineFit; % flipped line fit
    clear regressionLineFit;
end

% Model indices
if exist('ModelIndices', 'var')
    sortParameters.comp_modelIndices = ModelIndices; % model indices information
    clear ModelIndices;
end

% Nuisance indices
if exist('otherIndices', 'var')
    sortParameters.comp_nuisanceIndices = otherIndices; % nusiance indices information
    clear otherIndices;
end


% Partial correlations
if exist('subject_partial_corr', 'var')
    clear temp;
    temp = subject_partial_corr; clear subject_partial_corr;
    temp1 = subject_partial_slopes; clear subject_partial_slopes;
    % Loop over components
    for nn = 1:numComp
        subject_partial_corr{nn} = temp{sortedIndices(nn)};
        subject_partial_slopes{nn} = temp1{sortedIndices(nn)};
    end
    % End loop over components
    clear temp; clear temp1;
end

% flip ICA time course
% loop over components
temp = zeros(size(icaTimecourse));

for ii = 1:numComp
    temp(:, ii) = icaTimecourse(:, sortedIndices(ii));
end

icaTimecourse = temp;
clear temp;

% sorting parameters
sortParameters.sortedIndices = sortedIndices;
clear sortedIndices;
sortParameters.sortedValues = sortedValues;
clear sortedValues;
sortParameters.sortingCriteria = sortingCriteria;
sortParameters.sortingType = sortingType;
sortParameters.icaTimecourse = icaTimecourse;
clear icaTimecourse;
% model time course may be detrended depending upon the sorting criteria
% used
sortParameters.modelTimecourse = modelTimecourse;
clear modelTimecourse;
% reference information contains the spm file names used, regressor
% information like names, indices, time courses, number of subjects and
% number of sessions used to concatenate
sortParameters.refInfo = refInfo;
clear refInfo;

% spatial maps order not reversed (taking into account the memory problems)
sortParameters.icasig = icasig;

% count for the time points
sortParameters.diffTimePoints = diffTimePoints;

if ~isempty(maxVoxelPositions)
    maxVoxelPositions = maxVoxelPositions(sortParameters.sortedIndices);
end

titlePrint = [sortingType, ' sorting of components of ', viewingSet, ' using ', sortingCriteria, ' criteria:'];
titlePrint = upper(titlePrint);
printType = PRINTTYPE_REGRESSORS;
writeOptional = 'append';

if strcmpi(sortingCriteria, 'correlation')
    valueTitle = 'Correlation Value';
elseif strcmpi(sortingCriteria, 'multiple regression')
    valueTitle = 'Multiple Regression Value';
elseif strcmpi(sortingCriteria, 'kurtosis')
    valueTitle = 'Kurtosis Value';
elseif strcmpi(sortingCriteria, 'maximum voxel')
    valueTitle = 'Maximum Voxel Value';
else
    valueTitle  = [sortingCriteria, ' Value'];
end

numPara = 1;
varStruct(numPara).tag = 'Component Number';
varStruct(numPara).value = sortParameters.sortedIndices;

numPara = numPara + 1;
varStruct(numPara).tag = valueTitle;
varStruct(numPara).value = sortParameters.sortedValues;

maxVoxelStr = [];
clear maxVoxelPos;

% regression parameters
if strcmpi(sortingCriteria, 'multiple regression') | strcmpi(sortingCriteria, 'correlation')

    regressMatFileName = fullfile(outputDir, [inputPrefix, '_', sortingType, '_regression.mat']);
    % regression values
    regressionValues = sortParameters.sortedValues;

    % component numbers
    componentNumbers = sortParameters.sortedIndices;


    % split the regressors corresponding to the data set
    regressionParameters = cell2mat(sortParameters.regressionCoeff);
    regressionParameters = regressionParameters(sortParameters.comp_modelIndices, :);

    if size(componentNumbers, 2) == 1
        componentNumbers = componentNumbers;
    end

    % get the time course printing parameters
    if strcmpi(sortingType, 'temporal')

        [varStruct, regressor_all_datasets] = form_name_regressors(getSubjects, getSessions, sortParameters, varStruct, regressionParameters);

    else
        nVarStruct = length(varStruct);
        tempRegress = regressionParameters';
        % loop over number of images
        for nPrint = 1:size(spatialImage, 1)
            nVarStruct = nVarStruct + 1;
            [pathstr, fName, extn] = fileparts(deblank(spatialImage(nPrint, :)));
            varStruct(nVarStruct).tag = [fName];
            varStruct(nVarStruct).value = tempRegress(:, nPrint);
        end


    end
    % end for checking spatial or temporal

    if strcmpi(sortingCriteria, 'multiple regression')
        regressInfo.regressionParameters = regressionParameters;
        regressInfo.regressionValues = sortParameters.sortedValues';
        clear regressionParameters;
        regressInfo.componentNumbers = componentNumbers';
        regressInfo.comments = ['Concatenated regression parameters are in sorted order.', ...
            ' Columns indicate component numbers as indicated in field componentNumbers.'];

        if exist('regressor_all_datasets', 'var')
            regressInfo.regressor_all_datasets = regressor_all_datasets;
        end

        icatb_save(regressMatFileName, 'regressInfo');
        disp(['The concatenated matrix of regression parameters in sorted order is stored in ', ...
            regressMatFileName]);
        clear regressInfo;
    end


elseif strcmpi(sortingCriteria, 'maximum voxel')
    nVarStruct = length(varStruct);
    for ii = 1:length(maxVoxelPositions)
        %[tempX, tempY, tempZ] = ind2sub(structHInfo.DIM, maxVoxelPositions(ii));

        coords = squeeze(real_world_coords(:, maxVoxelPositions(ii)));

        % Real world coordinates
        tempX = coords(1);
        tempY = coords(2);
        tempZ = coords(3);

        maxVoxelPos(ii).name = ['[', num2str(tempX), ', ', num2str(tempY), ', ', num2str(tempZ), ']'];
    end
    maxVoxelStr = str2mat(maxVoxelPos.name);
    varStruct(nVarStruct).tag = 'Real World Coordinates (mm)';
    varStruct(nVarStruct).value = maxVoxelStr;
    clear maxVoxelPos;
    %titlePrint = [titlePrint, '(structural file: ', structHInfo.V(1).fname, ')'];

end
% end for checking sorting criteria

% form text file name
if strcmpi(sortParameters.sortingCriteria, 'multiple regression')
    outFileName = fullfile(outputDir, [inputPrefix, '_', sortingType, '_regression.txt']);
    msgStr = ['Regression parameters with the names of the regressors for the corresponding data set/sets are stored in ', ...
        outFileName];
elseif strcmpi(sortParameters.sortingCriteria, 'maximum voxel')
    outFileName = fullfile(outputDir, [inputPrefix, '_max_voxel.txt']);
    msgStr = ['Maximum voxel results are saved in ', outFileName];
    printType = 'column_wise';
else
    outFileName = fullfile(outputDir, [inputPrefix, '_', sortingType, '_', sortParameters.sortingCriteria, '.txt']);
    msgStr = [sortingType, ' ', sortParameters.sortingCriteria, ' results are saved in ', outFileName];
end

% print all the information to a text file
icatb_printToFile(outFileName, varStruct, titlePrint, printType, writeOptional);
disp(msgStr);


clear varStruct;


numPara = 1;
varStruct(numPara).tag = 'Component Number';
varStruct(numPara).value = sortParameters.sortedIndices;
% get the time course printing parameters
if strcmpi(sortingType, 'temporal') & strcmpi(sortParameters.sortingCriteria, 'multiple regression')
    subject_partial_corr = cell2mat(subject_partial_corr);
    subject_partial_slopes = cell2mat(subject_partial_slopes);
    [varStruct] = form_name_regressors(getSubjects, getSessions, sortParameters, varStruct, subject_partial_corr);
    titlePrint = ['Partial correlation values of viewing set: ', viewingSet];
    outFileName = fullfile(outputDir, [inputPrefix, '_', sortingType, '_partial_corr.txt']);
    msgStr = ['Partial correlation values with the names of the regressors for the corresponding data set/sets are stored in ', ...
        outFileName];
    icatb_printToFile(outFileName, varStruct, titlePrint, printType, writeOptional);

    % Include slope information as well
    titlePrint = ['Beta Weights (Slopes) of the regressors of viewing set: ', viewingSet];
    varStruct = varStruct(1);
    [varStruct] = form_name_regressors(getSubjects, getSessions, sortParameters, varStruct, subject_partial_slopes);
    icatb_printToFile(outFileName, varStruct, titlePrint, printType, writeOptional);
    fprintf('\n');
    disp(msgStr);
end





% Sub function to get the regressor namings
function  [varStruct, regressor_all_datasets] = form_name_regressors(getSubjects, getSessions, sortParameters, varStruct, dataMat)


spmMatFlag = sortParameters.refInfo.spmMatFlag;

countDataSets = 0;

% Store the regressors for each dataset
regressor_all_datasets = cell(getSubjects, getSessions);

% loop over all subjects
for numgetSub = 1:getSubjects
    % Loop over the sessions involved
    for numgetSess = 1:getSessions
        countDataSets = countDataSets + 1;
        % check the size of the SPMFile
        if size(sortParameters.refInfo.SPMFile, 2) > 1
            % Pull the Regressor names here
            regressorNames = str2mat(sortParameters.refInfo.selectedRegressors(countDataSets).name);
        else
            % added the special case here
            if strcmp(spmMatFlag, 'same_sub_diff_sess')
                % Pull the session specific Regressor names here
                regressorNames = str2mat(sortParameters.refInfo.selectedRegressors(numgetSess).name);
            else
                % Pull the Regressor names here
                regressorNames = str2mat(sortParameters.refInfo.selectedRegressors.name);
            end
        end

        % Number of Regressors
        numRegressors = size(sortParameters.refInfo.modelIndex, 2);

        if getSubjects == 1 & getSessions == 1
            tagPrefix = [''];
        elseif getSubjects == 1 & getSessions > 1
            tagPrefix = ['Session ', num2str(numgetSess)];
        elseif getSubjects > 1 & getSessions == 1
            tagPrefix = ['Subject ', num2str(numgetSub)];
        else
            tagPrefix = ['Subject ', num2str(numgetSub), ' Session ', num2str(numgetSess)];
        end

        nVarStruct = length(varStruct);

        %tempRegress = regressionParameters((countDataSets - 1)*numRegressors + 1:countDataSets*numRegressors, :)';
        tempData = dataMat((countDataSets - 1)*numRegressors + 1:countDataSets*numRegressors, :)';

        for nPrint = 1:size(regressorNames, 1)
            nVarStruct = nVarStruct + 1;
            varStruct(nVarStruct).tag = [tagPrefix, ' ', deblank(regressorNames(nPrint, :))];
            varStruct(nVarStruct).value = tempData(:, nPrint);
        end

        clear tempRegress;

        regressorNames = cellstr(regressorNames);

        try
            regressorNames = strtrim(regressorNames);
        catch
            regressorNames = regexprep(regressorNames, '^\s', '');
        end

        % Store regressors
        regressor_all_datasets{numgetSub, numgetSess} = regressorNames;

    end

end