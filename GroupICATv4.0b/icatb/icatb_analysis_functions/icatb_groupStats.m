function [sesInfo] = icatb_groupStats(sesInfo, statusHandle)
% Group stats is performed to calculate the mean, t-map, std
% Input: sesInfo - structure containing all parameters necessary for group
% ica analysis.

if(~exist('sesInfo','var'))
    [P] = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Parameter File', 'filter', '*param*.mat');
    if isempty(P)
        error('Parameter file is not selected for analysis');
    end
    [pathstr, fileName] = fileparts(P);
    outputDir = pathstr;
    % Make sure parameter file exists
    load(P);
    if ~exist('sesInfo', 'var')
        error(['The selected file ', P, ' does not contain the sesInfo variable']);
    end
else
    outputDir = sesInfo.outputDir;
end

if ~exist('statusHandle', 'var')
    statusHandle = [];
end

if sesInfo.isInitialized == 0
    error('Parameter file has not been initialized');
end


conserve_disk_space = 0;
if (isfield(sesInfo, 'conserve_disk_space'))
    conserve_disk_space = sesInfo.conserve_disk_space;
end

icaAlgo = icatb_icaAlgorithm; % available ICA algorithms

algoVal = sesInfo.algorithm; % algorithm index

% selected ICA algorithm
algorithmName = deblank(icaAlgo(algoVal, :));


if (strcmpi(algorithmName, 'moo-icar'))
    algorithmName = 'gig-ica';
end

parallelMode = 'serial';
try
    parallelMode = sesInfo.parallel_info.mode;
catch
end

toolboxNames = ver;
parallelCluster = ~isempty(find(strcmpi(cellstr(char(toolboxNames.Name)), 'parallel computing toolbox') == 1));

runParallel = 0;
if (strcmpi(parallelMode, 'parallel') && parallelCluster)
    runParallel = 1;
end

if (conserve_disk_space ~= 1)
    
    [modalityType, dataTitle, compSetFields] = icatb_get_modality;
    
    % naming of complex images
    [sesInfo, complexInfo] = icatb_name_complex_images(sesInfo, 'write');
    
    statusHandle = findobj('tag', [sesInfo.userInput.prefix, 'waitbar']);
    
    appDataName = 'gica_waitbar_app_data';
    
    if ~isempty(statusHandle)
        
        % get the status handles
        statusData = getappdata(statusHandle, appDataName);
        statusData.perCompleted = statusData.perCompleted + statusData.unitPerCompleted;
        setappdata(statusHandle, appDataName, statusData);
        set(statusHandle, 'name', [num2str(round(statusData.perCompleted*100)), '% analysis done']); waitbar(statusData.perCompleted, statusHandle);
        
    end
    
    
    mask_ind = sesInfo.mask_ind;
    
    disp(' ');
    disp('---------------------------------------------------------------------');
    disp('STARTING GROUP STATS STEP');
    disp('---------------------------------------------------------------------');
    
    %load defaults
    icatb_defaults;
    global MAX_SUBJECTS_IN_MEMORY;
    maxSubjectsInMemory =MAX_SUBJECTS_IN_MEMORY;
    
    global GROUP_ICA_INDEX;
    global INDIVIDUAL_ICA_INDEX;
    global MEAN_INDEX;
    global TMAP_INDEX;
    global STD_INDEX;
    global MEAN_ALL_INDEX;
    global ZIP_IMAGE_FILES;
    
    %%%%%%%% Test if the data sets have different number of time points %%%
    % Count the number of files in each of the datasets
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(sesInfo, 'flagTimePoints')
        flagTimePoints = sesInfo.flagTimePoints; % get the flag for time points
        diffTimePoints = sesInfo.diffTimePoints; % different time points
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
    
    numTimePoints = min(diffTimePoints);
    
    subjectICAFiles = icatb_parseOutputFiles('icaOutputFiles', sesInfo.icaOutputFiles, 'numOfSub', sesInfo.numOfSub, 'numOfSess', sesInfo.numOfSess, 'flagTimePoints', sesInfo.flagTimePoints);
    
    checkMATFiles = 0;
    if ~strcmpi(modalityType, 'eeg')
        checkMATFiles = (length(dir(fullfile(outputDir, [sesInfo.calibrate_components_mat_file, '*.mat']))) ~= (sesInfo.numOfSub*sesInfo.numOfSess));
    end
    
    if (checkMATFiles)
        filesToDelete = icatb_unZipSubjectMaps(sesInfo, subjectICAFiles);
    end
    
    if(sesInfo.numOfSub ==1 && sesInfo.numOfSess ==1)
        
        disp('Only one subject, not calculating stats');
        
        if ~isempty(statusHandle)
            
            % get the status handles
            statusData = getappdata(statusHandle, appDataName);
            statusData.perCompleted = statusData.perCompleted + statusData.unitPerCompleted;
            setappdata(statusHandle, appDataName, statusData);
            set(statusHandle, 'name', [num2str(round(statusData.perCompleted*100)), '% analysis done']); waitbar(statusData.perCompleted, statusHandle);
            
        end
        
    elseif(sesInfo.numOfSub == 1 && sesInfo.numOfSess > 1)
        
        %     %--get mean of each subject's back reconstructed ica component
        %     %--------------------------------------------------------------
        disp('--calculating mean ica component and timecourse for all subjects and sessions');
        numSub = sesInfo.numOfSub;
        numSess = sesInfo.numOfSess;
        
        meanICASig = zeros(length(mask_ind), sesInfo.numComp, numSess);
        meanA = zeros(numTimePoints, sesInfo.numComp, numSess);
        
        %for i=1:numSub
        for ses=1:numSess
            
            [tc, ic] = icatb_loadComp(sesInfo, (1:sesInfo.numComp), 'subjects', 1, 'sessions', ses, 'vars_to_load', compSetFields, 'subject_ica_files', subjectICAFiles);
            ic = ic';
            if size(ic, 2) ~= length(mask_ind)
                ic = ic(:, mask_ind);
            end
            %get mean
            meanICASig(:, :, ses) = squeeze(meanICASig(:, :, ses)) + ic';
            meanA(:, :, ses) = squeeze(meanA(:, :, ses)) + tc(1:numTimePoints, :);
            
            clear ic tc;
            
        end
        %end
        %
        % Initialise mean for different sessions
        meanICASig_all = zeros(length(mask_ind), sesInfo.numComp, 1);
        meanA_all = zeros(numTimePoints, sesInfo.numComp, 1);
        %
        %     % sum the mean of each session
        for ses = 1:numSess
            %         % mean for different sessions
            meanICASig_all = meanICASig(:, :, ses) + meanICASig_all;
            meanA_all = meanA(:, :, ses) + meanA_all;
        end
        %
        %     % mean for different sessions
        meanICASig_all = meanICASig_all ./ numSess;
        meanA_all = meanA_all ./ numSess;
        disp('done calculating mean for different sessions');
        
        % Transpose all of these so that the dimensions are components
        % by voxels
        meanICASig = permute(meanICASig, [3, 2, 1]);
        meanA = permute(meanA, [3, 1, 2]);
        
        meanICASig_all = permute(meanICASig_all, [3, 2, 1]);
        meanA_all = permute(meanA_all, [3, 1, 2]);
        
        if ~isempty(statusHandle)
            
            % get the status handles
            statusData = getappdata(statusHandle, appDataName);
            statusData.perCompleted = statusData.perCompleted + statusData.unitPerCompleted;
            setappdata(statusHandle, appDataName, statusData);
            set(statusHandle, 'name', [num2str(round(statusData.perCompleted*100)), '% analysis done']); waitbar(statusData.perCompleted, statusHandle);
            
        end
        
        %
        %     %%% add mean for different sessions
        %     %save mean for different sessions images
        outfile= sesInfo.icaOutputFiles(1).ses(1).name(1,:);
        [meanFiles, zipFileName, files_in_zip] = icatb_saveICAData(outfile, squeeze(meanICASig_all(1,:,:)), ...
            squeeze(meanA_all(1,:,:)), mask_ind, sesInfo.numComp,sesInfo.HInfo, sesInfo.dataType, complexInfo, outputDir);
        
        %if ~isempty(sesInfo.zipContents.zipFiles)
        [ddd, outfile] = fileparts(outfile);
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
        %end
        
    else
        %--get mean of each subject's back reconstructed ica component
        %--------------------------------------------------------------
        disp('--calculating mean ica component and timecourse');
        numSub = sesInfo.numOfSub;
        numSess = sesInfo.numOfSess;
        meanICASig = zeros(length(mask_ind), sesInfo.numComp, numSess);
        meanA = zeros(numTimePoints, sesInfo.numComp, numSess);
        
        for i=1:numSub
            for ses=1:numSess
                
                %load next file
                [tc, ic] = icatb_loadComp(sesInfo, (1:sesInfo.numComp), 'subjects', i, 'sessions', ses, 'vars_to_load', compSetFields, 'subject_ica_files', subjectICAFiles);
                ic = ic';
                if size(ic, 2) ~= length(mask_ind)
                    ic = ic(:, mask_ind);
                end
                %get mean
                
                meanICASig(:, :, ses) = squeeze(meanICASig(:, :, ses)) + ic';
                meanA(:, :, ses) = squeeze(meanA(:, :, ses)) + tc(1:numTimePoints, :);
                clear ic tc;
            end
        end
        
        for ses=1:numSess
            meanICASig(:, :, ses) = meanICASig(:, :, ses)./numSub;
            meanA(:, :, ses) = meanA(:, :, ses)./numSub;
            disp(['done calculating mean for session ', num2str(ses)]);
        end
        
        meanICASig_all = zeros(length(mask_ind), sesInfo.numComp, 1);
        meanA_all = zeros(numTimePoints, sesInfo.numComp, 1);
        
        % sum the mean of each session
        for ses = 1:numSess
            % mean for different sessions
            meanICASig_all = meanICASig(:, :, ses) + meanICASig_all;
            meanA_all = meanA(:, :, ses) + meanA_all;
        end
        
        % mean for different sessions
        meanICASig_all = meanICASig_all ./ numSess;
        meanA_all = meanA_all ./ numSess;
        disp('done calculating mean for different sessions');
        
        %--get variance and standard deviation
        %-----------------------------------------------------
        disp('--calculating variance and standard deviation of components');
        varICASig = zeros(size(meanICASig));
        varA = zeros(size(meanA));
        
        numSub = sesInfo.numOfSub;
        
        
        for i=1:numSub
            
            for ses=1:numSess
                %load next file
                [tc, ic] = icatb_loadComp(sesInfo, (1:sesInfo.numComp), 'subjects', i, 'sessions', ses, 'vars_to_load', compSetFields, 'subject_ica_files', subjectICAFiles);
                ic = ic';
                if size(ic, 2) ~= length(mask_ind)
                    ic = ic(:, mask_ind);
                end
                varICASig(:, :, ses) = squeeze(varICASig(:, :, ses)) + (ic' - squeeze(meanICASig(:, :, ses))).^2;
                varA(:, :, ses) = squeeze(varA(:, :, ses)) + (tc(1:numTimePoints, :) - squeeze(meanA(:, :, ses))).^2;
                clear ic tc;
            end
            
            
        end
        
        % implement std for different timepoints
        
        for ses =1:numSess
            %variance and standard deviation for ica components
            varICASig(:, :, ses) = squeeze(varICASig(:, :, ses))./(numSub-1);
            stdICASig(:, :, ses) = sqrt(squeeze(varICASig(:, :, ses)));
            %variance and standard deviation for timecourse
            varA(:, :, ses) = squeeze(varA(:, :, ses))./(numSub-1);
            stdA(:, :, ses) = sqrt(squeeze(varA(:, :, ses)));
        end
        disp('done calculating variance and standard deviation');
        
        
        % implement tmaps for different timepoints
        
        %--calculate tmaps
        %-----------------------------------------------------
        disp('--calculating tmaps');
        tmapICASig = zeros(size(meanICASig));
        for ses=1:numSess
            for comp=1:sesInfo.numComp
                tmap_ind = find(stdICASig(:, comp, ses) > eps);
                divisor = squeeze(stdICASig(tmap_ind, comp, ses))./sqrt(numSub-1);
                tmapICASig(tmap_ind, comp, ses) = squeeze(meanICASig(tmap_ind, comp, ses)) ./ divisor;
            end
        end
        disp('done calculating tmaps');
        
        
        % Transpose all of these so that the dimensions are components
        % by voxels
        meanICASig = permute(meanICASig, [3, 2, 1]);
        meanA = permute(meanA, [3, 1, 2]);
        
        meanICASig_all = permute(meanICASig_all, [3, 2, 1]);
        meanA_all = permute(meanA_all, [3, 1, 2]);
        
        stdICASig = permute(stdICASig, [3, 2, 1]);
        stdA = permute(stdA, [3, 1, 2]);
        
        varICASig  = permute(varICASig , [3, 2, 1]);
        varA = permute(varA, [3, 1, 2]);
        
        tmapICASig = permute(tmapICASig, [3, 2, 1]);
        
        
        %-SAVE RESULTS
        %--------------------------------------------------
        disp('...saving group stats data...');
        
        
        for ses=1:numSess
            %save mean images
            outfile=sesInfo.icaOutputFiles(MEAN_INDEX).ses(ses).name(1,:);
            [fileNames, zipFileName, files_in_zip] = icatb_saveICAData(outfile, squeeze(meanICASig(ses,:,:)), ...
                squeeze(meanA(ses,:,:)), mask_ind, sesInfo.numComp,sesInfo.HInfo, sesInfo.dataType, complexInfo, outputDir);
            
            %if ~isempty(sesInfo.zipContents.zipFiles)
            [ddd, outfile] = fileparts(outfile);
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
            %end
            
            %save standard deviation images
            outfile=sesInfo.icaOutputFiles(STD_INDEX).ses(ses).name(1,:);
            [fileNames, zipFileName, files_in_zip] = icatb_saveICAData(outfile, ...
                squeeze(stdICASig(ses,:,:)), squeeze(stdA(ses,:,:)), mask_ind, sesInfo.numComp,sesInfo.HInfo, ...
                sesInfo.dataType, complexInfo, outputDir);
            
            %if ~isempty(sesInfo.zipContents.zipFiles)
            [ddd, outfile] = fileparts(outfile);
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
            %end
            
            %save tmap
            outfile=sesInfo.icaOutputFiles(TMAP_INDEX).ses(ses).name(1,:);
            [fileNames, zipFileName, files_in_zip] = icatb_saveICAData(outfile, ...
                squeeze(tmapICASig(ses,:,:)), squeeze(meanA(ses,:,:)), mask_ind, sesInfo.numComp,sesInfo.HInfo, ...
                sesInfo.dataType, complexInfo, outputDir);
            
            %if ~isempty(sesInfo.zipContents.zipFiles)
            [ddd, outfile] = fileparts(outfile);
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
            %end
            
        end
        
        if ~isempty(statusHandle)
            
            % get the status handles
            statusData = getappdata(statusHandle, appDataName);
            statusData.perCompleted = statusData.perCompleted + statusData.unitPerCompleted;
            setappdata(statusHandle, appDataName, statusData);
            set(statusHandle, 'name', [num2str(round(statusData.perCompleted*100)), '% analysis done']); waitbar(statusData.perCompleted, statusHandle);
            
        end
        
        
        %%% add mean for different sessions
        %save mean for different sessions images
        outfile = sesInfo.icaOutputFiles(MEAN_ALL_INDEX).ses(1).name(1,:);
        [meanFiles, zipFileName, files_in_zip] = icatb_saveICAData(outfile, squeeze(meanICASig_all(1,:,:)), ...
            squeeze(meanA_all(1,:,:)), mask_ind, sesInfo.numComp,sesInfo.HInfo, sesInfo.dataType, ...
            complexInfo, outputDir);
        
        %if ~isempty(sesInfo.zipContents.zipFiles)
        [ddd, outfile] = fileparts(outfile);
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
        %end
        
        
    end
    % end for switch
    
    
    if (sesInfo.numOfSub*sesInfo.numOfSess > 1 && (~strcmpi(algorithmName, 'iva-gl') && ~strcmpi(algorithmName, 'iva-l') ...
            && ~strcmpi(algorithmName, 'constrained ica (spatial)') && ~strcmpi(algorithmName, 'gig-ica')))
        disp('Comparing mean image with the aggregate ...');
        disp('Value shows how much the mean component is close w.r.t aggregate component');
        %load ica data
        icain =[sesInfo.ica_mat_file, '.mat'];
        load(fullfile(outputDir, icain), 'icasig');
        % Compare first component of aggregate data set with first component of
        % mean data set
        [correlationValue] = icatb_corr(squeeze(icasig(1, :)), squeeze(meanICASig_all(1, 1, :)));
        disp(['The comparison value is found to be ', num2str(correlationValue)]);
    end
    
    
    sesInfo = removeDupCells(sesInfo);
    
    [pp, fileName] = fileparts(sesInfo.userInput.param_file);
    drawnow;
    icatb_save(fullfile(outputDir, [fileName, '.mat']), 'sesInfo');
    
    %% Write FNC and spectra info by default in file *_postprocess_results.mat
    if (strcmpi(modalityType, 'fmri'))
        if (~runParallel)
            icatb_postprocess_timecourses(sesInfo);
        else
            icatb_par_postprocess_timecourses(sesInfo);
        end
    end
    
    if (exist('filesToDelete', 'var') && ~isempty(filesToDelete))
        icatb_cleanupFiles(filesToDelete, outputDir);
    end
    
    disp('---------------------------------------------------------------------');
    disp('ENDING GROUP STATS STEP');
    disp('---------------------------------------------------------------------');
    disp(' ');
    
    if ~isempty(statusHandle)
        % get the status handles
        statusData = getappdata(statusHandle, appDataName);
        statusData.perCompleted = statusData.perCompleted + statusData.unitPerCompleted;
        setappdata(statusHandle, appDataName, statusData);
        set(statusHandle, 'name', [num2str(round(statusData.perCompleted*100)), '% analysis done']); waitbar(statusData.perCompleted, statusHandle);
        
    end
    
    
end


function sesInfo = removeDupCells(sesInfo)
%% Remove duplicate cells

if ((isfield(sesInfo, 'zipContents'))  && (~isempty(sesInfo.zipContents.zipFiles)))
    chkCells = cellfun('isempty', sesInfo.zipContents.zipFiles);
    sesInfo.zipContents.zipFiles(chkCells) = [];
    sesInfo.zipContents.files_in_zip(chkCells) = [];
    if (ispc)
        [dd, inds] = unique(lower(sesInfo.zipContents.zipFiles));
    else
        [dd, inds] = unique(sesInfo.zipContents.zipFiles);
    end
    inds = sort(inds);
    sesInfo.zipContents.zipFiles = sesInfo.zipContents.zipFiles(inds);
    sesInfo.zipContents.files_in_zip = sesInfo.zipContents.files_in_zip(inds);
end