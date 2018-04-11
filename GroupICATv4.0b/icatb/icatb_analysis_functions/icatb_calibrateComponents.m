function sesInfo = icatb_calibrateComponents(sesInfo, statusHandle)
% Calibrated image maps are computed using the individual image maps and
% time courses.
% Input: sesInfo - structure containing all parameters necessary for group
% ica analysis

if(~exist('sesInfo','var'))
    [P] = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Parameter File', 'filter', '*param*.mat');
    if isempty(P)
        error('Parameter file is not selected for analysis');
    end
    [pathstr,fileName]=fileparts(P);
    outputDir = pathstr;
    % Make sure parameter file exists
    load(P);
    if ~exist('sesInfo', 'var')
        error(['The selected file ', P, ' does not contain the sesInfo variable']);
    end
else
    outputDir = sesInfo.outputDir;
end

sesInfo.outputDir = outputDir;

if sesInfo.isInitialized == 0
    error('Parameter file has not been initialized');
end

if ~isfield(sesInfo, 'dataType')
    sesInfo.dataType = 'real';
end

if ~exist('statusHandle', 'var')
    statusHandle = [];
end


[pp, paramFile] = fileparts(sesInfo.param_file);
clear pp;
paramFile = fullfile(outputDir, [paramFile, '.mat']);

appDataName = 'gica_waitbar_app_data';

if ~isempty(statusHandle)
    
    % get the status handles
    statusData = getappdata(statusHandle, appDataName);
    statusData.perCompleted = statusData.perCompleted + statusData.unitPerCompleted;
    setappdata(statusHandle, appDataName, statusData);
    set(statusHandle, 'name', [num2str(round(statusData.perCompleted*100)), '% analysis done']); waitbar(statusData.perCompleted, statusHandle);
    
end


mask_ind = sesInfo.mask_ind;

[modalityType, dataTitle, compSetFields] = icatb_get_modality;

disp(' ');
disp('---------------------------------------------------------------------');
disp('STARTING TO SCALE COMPONENT SETS');
disp('---------------------------------------------------------------------');

%load defaults
icatb_defaults;
global MAX_SUBJECTS_IN_MEMORY;
global SUBJECT_ICA_INDEX;
global WRITE_COMPLEX_IMAGES;
global FUNCTIONAL_DATA_FILTER;
global CENTER_IMAGES;
global DETRENDNUMBER;

if (isempty(CENTER_IMAGES))
    CENTER_IMAGES = 0;
end


[dd1, dd2, imExtn] = fileparts(FUNCTIONAL_DATA_FILTER);

if strcmpi(sesInfo.dataType, 'complex')
    % tag for identifying complex images
    WRITE_COMPLEX_IMAGES = sesInfo.write_complex_images;
end


icaAlgo = icatb_icaAlgorithm; % available ICA algorithms

algoVal = sesInfo.algorithm; % algorithm index

% selected ICA algorithm
algorithmName = deblank(icaAlgo(algoVal, :));

componentSigns = ones(1, sesInfo.numComp);

useTemporalICA = 0;
if (strcmpi(modalityType, 'fmri'))
    try
        useTemporalICA = strcmpi(sesInfo.group_ica_type, 'temporal');
    catch
    end
end

if strcmpi(algorithmName, 'moo-icar')
    algorithmName = 'gig-ica';
end

if (~useTemporalICA)
    if (~strcmpi(algorithmName, 'iva-gl') && ~strcmpi(algorithmName, 'iva-l') && ~strcmpi(algorithmName, 'constrained ica (spatial)') ...
            && ~strcmpi(algorithmName, 'gig-ica'))
        %--load reference image
        %---------------------------------------------------------
        icain=sesInfo.ica_mat_file;
        ref = load(fullfile(outputDir, [icain, '.mat']));
        refImage = ref.icasig;
        W = ref.W;
        
        pcain = [sesInfo.data_reduction_mat_file,num2str(sesInfo.numReductionSteps),'-',num2str(1), '.mat'];
        load(fullfile(outputDir, pcain), 'pcasig');
        
        if size(refImage, 2) == prod(sesInfo.HInfo.DIM(1:3))
            refImage = refImage(:, mask_ind);
        end
        
        componentSigns = zeros(1, size(refImage, 1));
        pcasig = W*pcasig';
        
        for nC = 1:length(componentSigns)
            componentSigns(nC) = sign(icatb_corr2(pcasig(nC, :), refImage(nC, :)));
        end
        
    end
end

clear ref;

conserve_disk_space = 0;
if (isfield(sesInfo, 'conserve_disk_space'))
    conserve_disk_space = sesInfo.conserve_disk_space;
end

if (strcmpi(algorithmName, 'gig-ica') || strcmpi(algorithmName, 'constrained ica (spatial)'))
    conserve_disk_space = 0;
end

backReconType = 'regular';
if (isfield(sesInfo, 'backReconType'))
    backReconType = sesInfo.backReconType;
end

if (~strcmpi(algorithmName, 'iva-gl') && ~strcmpi(algorithmName, 'iva-l') && ~strcmpi(algorithmName, 'gig-ica') && ~strcmpi(algorithmName, 'constrained ica (spatial)'))
    if (conserve_disk_space == 1)
        if (~strcmpi(backReconType, 'spatial-temporal regression'))
            [sesInfo.tcInfo, sesInfo.icInfo] = icatb_groupBackReconInfo(sesInfo, W);
        end
    end
end

dataSetsToRun = (1:sesInfo.numOfSub*sesInfo.numOfSess);

if (isfield(sesInfo, 'dataSetNo'))
    dataSetsToRun(dataSetsToRun < sesInfo.dataSetNo) = [];
    sesInfo = rmfield(sesInfo, 'dataSetNo');
end

sesInfo.center_images = 0;
sesInfo.detrendNumber = DETRENDNUMBER;

[subjectICAFiles, meanICAFiles, tmapICAFiles, meanALL_ICAFile] = icatb_parseOutputFiles('icaOutputFiles', sesInfo.icaOutputFiles, 'numOfSub', sesInfo.numOfSub, 'numOfSess', sesInfo.numOfSess, ...
    'flagTimePoints', sesInfo.flagTimePoints);

%% Center mean image distribution to zero
meanImOffsets = zeros(1, sesInfo.numComp);

if (strcmpi(modalityType, 'fmri') && (CENTER_IMAGES == 1))
    sesInfo.center_images = 1;
end

%if (~strcmpi(algorithmName, 'iva-gl') && ~strcmpi(algorithmName, 'moo-icar') && ~strcmpi(algorithmName, 'constrained ica (spatial)'))
if (conserve_disk_space ~= 1 && sesInfo.center_images)
    disp('Computing offset using the mean component maps which will be subtracted to the subject component maps ...');
    
    % Compute mean
    for nDataSet = 1:sesInfo.numOfSub*sesInfo.numOfSess
        load(fullfile(outputDir, [sesInfo.back_reconstruction_mat_file, num2str(nDataSet), '.mat']));
        tmp = getfield(compSet, compSetFields{1});
        clear compSet;
        if (nDataSet == 1)
            meanIm = zeros(size(tmp));
        end
        meanIm = meanIm + tmp;
        clear tmp;
    end
    
    meanIm = meanIm / (sesInfo.numOfSub*sesInfo.numOfSess);
    
    % Compute offsets based on mean image
    meanImOffsets = zeros(1, size(meanIm, 1));
    
    % Loop over components
    for nIm = 1:length(meanImOffsets)
        [dd, meanImOffsets(nIm)] = icatb_recenter_image(meanIm(nIm, :));
    end
    % End of loop over components
    
    clear meanIm;
    fprintf('Done\n');
end
%end

%%%%%%%% Test if the data sets have different number of time points %%%
% Count the number of files in each of the datasets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(sesInfo, 'flagTimePoints')
    flagTimePoints = sesInfo.flagTimePoints; % get the flag for time points
    diffTimePoints = sesInfo.diffTimePoints;
else
    % get the count of the files
    [diffTimePoints] = icatb_get_countTimePoints(sesInfo.userInput.files);
    
    % check time points
    checkTimePoints = find(diffTimePoints ~= diffTimePoints(1));
    
    if ~isempty(checkTimePoints)
        flagTimePoints = 'different_time_points';
    else
        flagTimePoints = 'same_time_points';
    end
    sesInfo.diffTimePoints = diffTimePoints;
    sesInfo.flagTimePoints = flagTimePoints;
end
% end for checking time points

if ~isfield(sesInfo, 'zipContents')
    sesInfo.zipContents.zipFiles = {};
    sesInfo.zipContents.files_in_zip.name = {};
end

%outfile = strrep(outfile,'.img','');
% naming of complex images
% [sesInfo, complexInfo] = icatb_name_complex_images(sesInfo, 'write');
%
% % global variable necessary for defining the real&imag or magnitude&phase
% WRITE_COMPLEX_IMAGES = sesInfo.userInput.write_complex_images;

parallelMode = 'serial';

try
    parallelMode = sesInfo.parallel_info.mode;
catch
end

toolboxNames = ver;
parallelCluster= ~isempty(find(strcmpi(cellstr(char(toolboxNames.Name)), 'parallel computing toolbox') == 1));

if (strcmpi(parallelMode, 'parallel') && parallelCluster)
    % Parallel
    statusHandle = [];
    sesInfo.scale_opts.componentSigns = componentSigns;
    sesInfo.scale_opts.meanImOffsets = meanImOffsets;
    sesInfo.scale_opts.compSetFields = compSetFields;
    sesInfo = icatb_parCalibrateComponents_cluster(sesInfo, dataSetsToRun, conserve_disk_space);
    
else
    % Serial
    for i = dataSetsToRun
        
        subjectNumber = ceil(i / sesInfo.numOfSess);
        sessNumber = mod(i - 1, sesInfo.numOfSess) + 1;
        disp(['--Subject ', num2str(subjectNumber),' Session ', num2str(sessNumber), '''s Component Set']);
        if (conserve_disk_space ~= 1)
            subFile = [sesInfo.back_reconstruction_mat_file, num2str(i), '.mat'];
            load(fullfile(outputDir, subFile));
        else
            compSet = getBackReconSet(sesInfo, i);
        end
        
        % IC
        tempIC = getfield(compSet, compSetFields{1});
        % TC
        tempTC = getfield(compSet, compSetFields{2});
        
        clear compSet;
        
        
        if (strcmpi(modalityType, 'smri'))
            sesInfo.ica_variances = compute_var(sesInfo, tempTC', tempIC');
        end
        
        % Since multi-subject BR comps are written as Time by components.
        % Reshape it as components by time
        if (sesInfo.numOfSub*sesInfo.numOfSess > 1)
            tempTC = tempTC';
        end
        
        
        
        % get sign of correlation between subject components and reference
        % component
        for compNum = 1:sesInfo.numComp
            tempIC(compNum, :) = tempIC(compNum, :) - meanImOffsets(compNum);
            %         subjectImage = detrend(tempIC(compNum, :), 0);
            %         if length(subjectImage) ~= size(refImage, 2)
            %             subjectImage = subjectImage(mask_ind);
            %         end
            %corSign = sign(icatb_corr(refImage(compNum, :), subjectImage));
            % check the sign function
            if sign(componentSigns(compNum)) == -1
                string = ['  Difference in Sign Between Reference Image and Subject ', num2str(subjectNumber), ' Session ', num2str(sessNumber), ...
                    ' Component ', num2str(compNum)];
                disp([string,' -> Changing Sign of component and timecourse']);
                tempIC(compNum, :) = tempIC(compNum, :) * -1;
                tempTC(compNum, :) = tempTC(compNum, :) * -1;
            end
        end
        drawnow;
        
        
        clear compSet;
        %[sesInfo, complexInfoRead] = icatb_name_complex_images(sesInfo, 'read');
        
        if (strcmpi(modalityType, 'fmri'))
            
            if (~sesInfo.center_images && (sesInfo.scaleType == 2))
                % Remove mean of component images
                disp('Removing mean of component images ...');
                tempIC = icatb_remove_mean(tempIC')';
                fprintf('Done\n');
            end
        end
        
        
        %scale timecourses and components
        %     [tempIC, tempTC] = icatb_scaleICA(tempIC, tempTC, sesInfo.inputFiles(i).name, ...
        %         sesInfo.scaleType, sesInfo.dataType, complexInfoRead, mask_ind);
        
        [tempIC, tempTC] = icatb_scaleICA(tempIC, tempTC, sesInfo.inputFiles(i).name, ...
            sesInfo.scaleType, sesInfo.dataType, [], mask_ind);
        
        if ~isempty(statusHandle)
            
            % get the status handles
            statusData = getappdata(statusHandle, appDataName);
            statusData.perCompleted = statusData.perCompleted + (statusData.unitPerCompleted / length(dataSetsToRun));
            setappdata(statusHandle, appDataName, statusData);
            set(statusHandle, 'name', [num2str(round(statusData.perCompleted*100)), '% analysis done']); waitbar(statusData.perCompleted, statusHandle);
            
        end
        
        if ~strcmpi(modalityType, 'eeg')
            
            if (conserve_disk_space ~= 1)
                %--- save data
                if strcmpi(imExtn, '.nii')
                    dispStr = ['...saving scaled ica data for subject ', num2str(subjectNumber), ' session ', num2str(sessNumber),...
                        ' in nifti format and as matlab file'];
                else
                    dispStr = ['...saving scaled ica data for subject ', num2str(subjectNumber), ' session ', num2str(sessNumber),...
                        ' in analyze format and as matlab file'];
                end
            else
                if strcmpi(imExtn, '.nii')
                    dispStr = ['...saving scaled ica data for subject ', num2str(subjectNumber), ' session ', num2str(sessNumber),...
                        ' in nifti format'];
                else
                    dispStr = ['...saving scaled ica data for subject ', num2str(subjectNumber), ' session ', num2str(sessNumber),...
                        ' in analyze format'];
                end
                
            end
            
        else
            
            dispStr = ['...saving scaled ica data for subject ', num2str(subjectNumber), ' session ', num2str(sessNumber),...
                ' as matlab file'];
            
        end
        
        disp(dispStr);
        
        if ((conserve_disk_space ~= 1) || strcmpi(modalityType, 'eeg'))
            
            % Get var names from compSetfields to store in calibration MAT files
            eval([compSetFields{1}, ' = tempIC;']);
            eval([compSetFields{2}, ' = tempTC;']);
            
            % save data in matlab file
            subFileOut = [sesInfo.calibrate_components_mat_file, num2str(subjectNumber), '-', num2str(sessNumber), '.mat'];
            subFileOut = fullfile(outputDir, subFileOut);
            drawnow;
            
            icatb_save(subFileOut, compSetFields{:});
            
            eval(['clear ', compSetFields{1}, ' ', compSetFields{2}]);
            
        end
        
        if ~strcmpi(modalityType, 'eeg')
            if ((conserve_disk_space == 1) && sesInfo.center_images)
                saveImgFiles(sesInfo, subjectICAFiles(subjectNumber).ses(sessNumber).name(1, :), tempIC, tempTC, 'no');
            else
                sesInfo = saveImgFiles(sesInfo, subjectICAFiles(subjectNumber).ses(sessNumber).name(1, :), tempIC, tempTC);
            end
        end
        
    end
    
end


if (isfield(sesInfo, 'tcInfo'))
    sesInfo = rmfield(sesInfo, 'tcInfo');
end

if (isfield(sesInfo, 'icInfo'))
    sesInfo = rmfield(sesInfo, 'icInfo');
end


tp = min(sesInfo.diffTimePoints);

if (conserve_disk_space == 1)
    
    %% Compute mean over datasets
    if (sesInfo.numOfSub*sesInfo.numOfSess > 1)
        fprintf('\nComputing mean components across all data-sets\n');
        if ~strcmpi(modalityType, 'eeg')
            filesToDelete = icatb_unZipSubjectMaps(sesInfo, subjectICAFiles);
            XYZ = icatb_get_voxel_coords(sesInfo.HInfo.DIM(1:3));
        end
        countDataSet = 0;
        % Loop over subjects
        for nSub = 1:sesInfo.numOfSub
            % Loop over sessions
            for nSess = 1:sesInfo.numOfSess
                countDataSet = countDataSet + 1;
                if ~strcmpi(modalityType, 'eeg')
                    tmpFile = icatb_fullFile('files', subjectICAFiles(nSub).ses(nSess).name, 'directory', outputDir);
                    ic = read_data(tmpFile, XYZ(:, mask_ind));
                    tc = icatb_loadICATimeCourse(tmpFile);
                else
                    tmpFile = [sesInfo.calibrate_components_mat_file, num2str(nSub), '-', num2str(nSess), '.mat'];
                    tmp = load(fullfile(outputDir, tmpFile));
                    ic = getfield(tmp, compSetFields{1});
                    tc = getfield(tmp, compSetFields{2});
                    clear tmp;
                end
                if (countDataSet == 1)
                    meanIC = zeros(size(ic));
                    meanTC = zeros(tp, size(tc, 2));
                end
                meanIC = meanIC + ic;
                meanTC = meanTC + tc(1:tp, :);
                clear ic tc;
            end
        end
        
        meanIC = meanIC / (sesInfo.numOfSub*sesInfo.numOfSess);
        meanTC = meanTC / (sesInfo.numOfSub*sesInfo.numOfSess);
    else
        meanIC = tempIC;
        meanTC = tempTC;
    end
    
    
    if (strcmpi(modalityType, 'fmri') && sesInfo.center_images)
        
        disp('Computing offset using the mean component maps which will be subtracted to the subject component maps ...');
        % Compute offsets based on mean image
        meanImOffsets = zeros(1, size(meanIC, 1));
        
        % Loop over components
        for nIm = 1:length(meanImOffsets)
            [meanIC(nIm, :), meanImOffsets(nIm)] = icatb_recenter_image(meanIC(nIm, :));
        end
        % End of loop over components
        fprintf('Done\n');
        
        fprintf('\n');
        disp('Applying the offset to subject component images ...');
        if (sesInfo.numOfSub*sesInfo.numOfSess > 1)
            filesToDelete = {};
            % Loop over subjects
            for subjectNumber = 1:sesInfo.numOfSub
                % Loop over sessions
                for sessNumber = 1:sesInfo.numOfSess
                    tmpFile = icatb_fullFile('files', subjectICAFiles(subjectNumber).ses(sessNumber).name, 'directory', outputDir);
                    tempIC = read_data(tmpFile, XYZ(:, mask_ind));
                    for compNum = 1:sesInfo.numComp
                        tempIC(compNum, :) = tempIC(compNum, :) - meanImOffsets(compNum);
                    end
                    tempTC = icatb_loadICATimeCourse(tmpFile);
                    if strcmpi(imExtn, '.nii')
                        dispStr = ['...saving scaled ica data for subject ', num2str(subjectNumber), ' session ', num2str(sessNumber),...
                            ' in nifti format'];
                    else
                        dispStr = ['...saving scaled ica data for subject ', num2str(subjectNumber), ' session ', num2str(sessNumber),...
                            ' in analyze format'];
                    end
                    disp(dispStr);
                    [sesInfo, filesInZip] = saveImgFiles(sesInfo, subjectICAFiles(subjectNumber).ses(sessNumber).name(1, :), tempIC, tempTC, [], 0);
                    clear tempIC tempTC;
                    if (~isempty(filesInZip))
                        filesToDelete{length(filesToDelete) + 1} = filesInZip;
                    end
                    clear filesInZip;
                end
                % End of loop over sessions
            end
            % End of loop over subjects
            icatb_cleanupFiles(filesToDelete, outputDir);
        else
            sesInfo = saveImgFiles(sesInfo, subjectICAFiles(subjectNumber).ses(sessNumber).name(1, :), meanIC, meanTC);
        end
        fprintf('Done\n');
        
    end
    
    
    if (sesInfo.numOfSub*sesInfo.numOfSess > 1)
        
        outfile = deblank(meanALL_ICAFile(1).name(1, :));
        [meanFiles, zipFileName, files_in_zip] = icatb_saveICAData(outfile, meanIC, meanTC, mask_ind, sesInfo.numComp,sesInfo.HInfo, sesInfo.dataType, [], outputDir);
        [ddd, outfile] = fileparts(outfile);
        fprintf('Done\n');
        currentZipFile = regexprep(outfile, '\_\d*$', '_.zip');
        
        if (ispc)
            chkCells = strcmpi(sesInfo.zipContents.zipFiles, currentZipFile);
        else
            chkCells = strcmp(sesInfo.zipContents.zipFiles, currentZipFile);
        end
        
        sesInfo.zipContents.zipFiles(chkCells) = [];
        sesInfo.zipContents.files_in_zip(chkCells) = [];
        
        % store the zip filenames to a structure
        countZip = length(sesInfo.zipContents.zipFiles) + 1;
        sesInfo.zipContents.zipFiles{countZip} = zipFileName;
        sesInfo.zipContents.files_in_zip(countZip).name = files_in_zip;
        
        
        if (~strcmpi(algorithmName, 'iva-gl') && ~strcmpi(algorithmName, 'iva-l')  && ~strcmpi(algorithmName, 'constrained ica (spatial)') ...
                && ~strcmpi(algorithmName, 'gig-ica'))
            
            disp('Comparing mean image with the aggregate ...');
            disp('Value shows how much the mean component is close w.r.t aggregate component');
            % Compare first component of aggregate data set with first component of
            % mean data set
            correlationValue = icatb_corr(refImage(1, :), meanIC(1, :));
            disp(['The comparison value is found to be ', num2str(correlationValue)]);
            fprintf('\n');
            
        end
        
        clear meanIC meanTC;
        
        
        if (exist('filesToDelete', 'var') && ~isempty(filesToDelete))
            icatb_cleanupFiles(filesToDelete, outputDir);
        end
        
    end
    
end


drawnow;
icatb_save(paramFile, 'sesInfo');

disp('---------------------------------------------------------------------');
disp('DONE SCALING COMPONENTS');
disp('---------------------------------------------------------------------');
disp(' ');

if ~isempty(statusHandle)
    
    % get the status handles
    statusData = getappdata(statusHandle, appDataName);
    statusData.perCompleted = statusData.perCompleted + statusData.unitPerCompleted;
    setappdata(statusHandle, appDataName, statusData);
    set(statusHandle, 'name', [num2str(round(statusData.perCompleted*100)), '% analysis done']); waitbar(statusData.perCompleted, statusHandle);
    
end


function [sesInfo, files_in_zip] = saveImgFiles(sesInfo, outfile, tempIC, tempTC, zipFiles, deleteFiles)
%% Save images and compress images if necessary
%

icatb_defaults;
global ZIP_IMAGE_FILES;

if (~exist('zipFiles', 'var') || isempty(zipFiles))
    zipFiles = ZIP_IMAGE_FILES;
end

if (~exist('deleteFiles', 'var'))
    deleteFiles = 1;
end

[fileNames, zipFileName, files_in_zip] = icatb_saveICAData(outfile, tempIC, tempTC, sesInfo.mask_ind, sesInfo.numComp, ...
    sesInfo.HInfo, sesInfo.dataType, [], sesInfo.outputDir, zipFiles, deleteFiles);


[ddd, outfile] = fileparts(deblank(outfile));
if (strcmpi(zipFiles, 'yes'))
    currentZipFile = regexprep(outfile, '\_\d*$', '_.zip');
    
    if (ispc)
        chkCells = strcmpi(sesInfo.zipContents.zipFiles, currentZipFile);
    else
        chkCells = strcmp(sesInfo.zipContents.zipFiles, currentZipFile);
    end
    
    sesInfo.zipContents.zipFiles(chkCells) = [];
    sesInfo.zipContents.files_in_zip(chkCells) = [];
    
    % store the zip filenames to a structure
    countZip = length(sesInfo.zipContents.zipFiles) + 1;
    sesInfo.zipContents.zipFiles{countZip} = zipFileName;
    sesInfo.zipContents.files_in_zip(countZip).name = files_in_zip;
    
    % Keep track of files added
    checkZipFiles = sesInfo.zipContents.zipFiles(cellfun('isempty', sesInfo.zipContents.zipFiles) == 0);
    if (~isempty(checkZipFiles))
        icatb_save(fullfile(sesInfo.outputDir, sesInfo.param_file), 'sesInfo');
    end
end

function Y = read_data(files, XYZ)
%% Read data
%

Y = icatb_spm_get_data(icatb_spm_vol(files), XYZ, 0);
Y(isfinite(Y) == 0) = 0;

function compSet = getBackReconSet(sesInfo, nDataSet)
%% Get back-reconstructed component of the subject
%

sesInfo.dataSetNo = nDataSet;
[sesInfo, compSet] = icatb_backReconstruct(sesInfo);


function pf = compute_var(sesInfo, A, S)
% Compute percent variance

data = icatb_remove_mean(icatb_read_data(sesInfo.inputFiles(1).name, [], sesInfo.mask_ind));
pf = icatb_compute_var_evd(data, A, S);