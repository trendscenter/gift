function [resume_info, sesInfo] = icatb_get_resume_info(sesInfo)
%% Check if the analysis can be resumed or not
%
% Inputs:
% sesInfo - Session information
%
% Outputs:
% resume_info - Resume information
%

icatb_defaults;
global NUM_RUNS_GICA;
global CENTER_IMAGES;
global DETRENDNUMBER;

if (isempty(CENTER_IMAGES))
    CENTER_IMAGES = 0;
end

if (isempty(NUM_RUNS_GICA))
    NUM_RUNS_GICA = 1;
end

resume_info = [];
groupNo = 1;
dataSetNo = 1;

parallelMode = 'serial';
try
    parallelMode = sesInfo.parallel_info.mode;
catch
end

if (strcmpi(parallelMode, 'parallel'))
    disp('!!!! Resume option is not available when group ICA is run in parallel mode. Try running the analysis steps individually in case the analysis is interrupted');
    disp('');
    return;
end


ica_types = cellstr(icatb_icaAlgorithm);

conserve_disk_space = 0;
if (isfield(sesInfo, 'conserve_disk_space'))
    conserve_disk_space = sesInfo.conserve_disk_space;
end

sesInfo.conserve_disk_space = conserve_disk_space;

clear conserve_disk_space;

if ~isfield(sesInfo.userInput, 'group_pca_type')
    sesInfo.userInput.group_pca_type = 'subject specific';
end

if ~isfield(sesInfo, 'group_pca_type')
    sesInfo.group_pca_type = 'subject specific';
end

oldSessInfo = sesInfo;

fprintf('\n');
fprintf('By default, running parameter initialization step to track changes in user input ...');
fprintf('\n');

%% Run parameter initialization by default
sesInfo.saveParamFile = 0;
sesInfo = icatb_parameterInitialization(sesInfo);


[modalityType, dataTitle, compSetFields] = icatb_get_modality;

cd(sesInfo.outputDir);

subjectFile = fullfile(sesInfo.outputDir, [sesInfo.userInput.prefix, 'Subject.mat']);

if (~exist(subjectFile, 'file'))
    files = sesInfo.userInput.files;
    numOfSub = sesInfo.userInput.numOfSub;
    numOfSess = sesInfo.userInput.numOfSess;
    SPMFiles = sesInfo.userInput.designMatrix;
    icatb_save(subjectFile, 'files', 'numOfSub', 'numOfSess', 'SPMFiles', 'modalityType');
    clear files SPMFiles;
end

zipContents.zipFiles = {};
if (isfield(sesInfo, 'zipContents'))
    zipContents = sesInfo.zipContents;
end

chkCell = cellfun('isempty', zipContents.zipFiles);
zipContents.zipFiles(chkCell) = [];
if (~isempty(zipContents.zipFiles))
    if (ispc)
        [dd, uniq_cells] = unique(lower(zipContents.zipFiles));
    else
        [dd, uniq_cells] = unique(zipContents.zipFiles);
    end
    zipContents.zipFiles = zipContents.zipFiles(sort(uniq_cells));
end

[subjectICAFiles, meanICAFiles, tmapICAFiles, meanALL_ICAFile, groupICAFiles, stdICAFiles] = icatb_parseOutputFiles('icaOutputFiles', sesInfo.icaOutputFiles, 'numOfSub', ...
    sesInfo.numOfSub, 'numOfSess', sesInfo.numOfSess);

try
    [old_subjectICAFiles, old_meanICAFiles, old_tmapICAFiles, old_meanALL_ICAFile, old_groupICAFiles, old_stdICAFiles] = icatb_parseOutputFiles('icaOutputFiles', ...
        oldSessInfo.icaOutputFiles, 'numOfSub', oldSessInfo.numOfSub, 'numOfSess', oldSessInfo.numOfSess);
catch
    
end

fprintf('\n');

disp('Checking missing analysis files ....');

if (sesInfo.conserve_disk_space ~= 1)
    stepsToCheck = (3:7);
else
    stepsToCheck = [3, 4, 6];
end

if strcmpi(ica_types{sesInfo.algorithm}, 'moo-icar')
    ica_types{sesInfo.algorithm} = 'gig-ica';
end

if (strcmpi(ica_types{sesInfo.algorithm}, 'gig-ica') || strcmpi(ica_types{sesInfo.algorithm}, 'constrained ica (spatial)'))
    % No data reduction
    stepsToCheck(stepsToCheck == 3) = [];
    % No back-reconstruction
    stepsToCheck(stepsToCheck == 5) = [];
end


filesInfo = repmat(struct('step', [], 'value', [], 'files', []), 1, length(stepsToCheck));

status = 1;

countStep = 0;
for nStep = stepsToCheck
    
    countStep = countStep + 1;
    
    %% Data reduction
    if (nStep == 3)
        
        fileNaming = sesInfo.data_reduction_mat_file;
        
        % Loop over groups
        for j = 1:sesInfo.numReductionSteps
            
            tmpFileN = [fileNaming, num2str(j), '-'];
            
            pos = findstr(filesep, tmpFileN);
            relativePath = '';
            if (~isempty(pos))
                relativePath = tmpFileN(1:pos(end));
            end
            
            chkOneStepPCA = strcmpi(ica_types{sesInfo.algorithm}, 'iva-gl') || strcmpi(ica_types{sesInfo.algorithm}, 'iva-l') || strcmpi(ica_types{sesInfo.algorithm}, 'gig-ica') || ...
                strcmpi(ica_types{sesInfo.algorithm}, 'constrained ica (spatial)');
            
            totalPCAs = sesInfo.reduction(j).numOfGroupsAfterCAT;
            
            if (~chkOneStepPCA && sesInfo.numReductionSteps == 1)
                totalPCAs = 1;
            end
            
            filesInfo(countStep).step = [[filesInfo(countStep).step], [repmat(nStep, 1, totalPCAs); repmat(j, 1, totalPCAs); (1:totalPCAs)]];
            allFiles =  repmat({''}, 1, totalPCAs);
            for nTPC = 1:totalPCAs
                allFiles{nTPC} = fullfile(sesInfo.outputDir, [tmpFileN, num2str(nTPC), '.mat']);
            end
            filesInfo(countStep).files = [filesInfo(countStep).files, allFiles];
            
            tmpFiles = dir(fullfile(sesInfo.outputDir, [tmpFileN, '*.mat']));
            
            if (~strcmpi(ica_types{sesInfo.algorithm}, 'iva-gl') && ~strcmpi(ica_types{sesInfo.algorithm}, 'iva-l'))
                if (j == sesInfo.numReductionSteps)
                    try
                        tmpFiles = tmpFiles(1);
                    catch
                    end
                end
            else
                tmpFiles((1:length(tmpFiles)) > sesInfo.numOfSub*sesInfo.numOfSess) = [];
            end
            
            
            dateInfo = char(tmpFiles.date);
            tmpFiles2 = char(tmpFiles.name);
            tmpFiles = tmpFiles2;
            tmpFiles = [repmat(relativePath, size(tmpFiles, 1), 1), tmpFiles];
            clear tmpFiles2;
            
            if (isempty(tmpFiles))
                status = 0;
                groupNo = j;
                dataSetNo = 1;
                break;
            else
                
                chkFirstPCs = 0;
                for nAllFiles = 1:length(allFiles)
                    if (~exist(allFiles{nAllFiles}, 'file'))
                        dataSetNo = nAllFiles;
                        chkFirstPCs = 1;
                        break;
                    end
                end
                
                if (chkFirstPCs)
                    status = 0;
                    groupNo = j;
                    break;
                end
                
                try
                    %% Detect change in pre-processing options
                    if (j == 1)
                        
                        preprocType = 'remove mean per timepoint';
                        old_preprocType = preprocType;
                        preproc_opts = icatb_preproc_data;
                        try
                            
                            old_preprocType = oldSessInfo.preproc_type;
                            preprocType = sesInfo.preproc_type;
                            
                            if (isnumeric(preprocType))
                                preprocType = lower(preproc_opts{preprocType});
                            end
                            
                            if (isnumeric(old_preprocType))
                                old_preprocType = lower(preproc_opts{old_preprocType});
                            end
                            
                        catch
                        end
                        
                        if (~strcmpi(preprocType, old_preprocType))
                            error('Preprocessing options changed');
                        end
                        
                        % Group PCA Type
                        if (sesInfo.numOfSub*sesInfo.numOfSess > 1)
                            if (~strcmpi(sesInfo.group_pca_type, oldSessInfo.group_pca_type))
                                error('Group PCA Type is changed.');
                            end
                        end
                        
                        
                        mask = zeros(prod(sesInfo.HInfo.DIM(1:3)), 1);
                        old_mask = mask;
                        mask(sesInfo.mask_ind) = 1;
                        old_mask(oldSessInfo.mask_ind) = 1;
                        
                        checkMaskInd = find((mask == old_mask) ~= 1);
                        
                        if (~isempty(checkMaskInd))
                            error('Mask changed');
                        end
                        
                    end
                    
                    %% No. of components
                    if (j == 1)
                        load(fullfile(sesInfo.outputDir, deblank(tmpFiles(1, :))), 'V');
                        if (~exist('V', 'var'))
                            error('V variable not found');
                        end
                        
                        compDims = size(V, 2);
                        clear V;
                    else
                        load(fullfile(sesInfo.outputDir, deblank(tmpFiles(1, :))), 'pcasig');
                        if (~exist('pcasig', 'var'))
                            error('pcasig variable not found');
                        end
                        compDims = size(pcasig, 2);
                        clear pcasig;
                    end
                    
                    if (compDims ~= sesInfo.reduction(j).numOfPCAfterReduction)
                        error('Dimensions of pcasig changed');
                    end
                    
                    %% Handle the case when no. of reduction steps is changed
                    if (j > 1)
                        if (oldSessInfo.numReductionSteps ~= sesInfo.numReductionSteps)
                            error('No of reduction steps changed');
                        end
                    end
                    
                catch
                    status = 0;
                    groupNo = j;
                    dataSetNo = 1;
                end
                
                
                if (~status)
                    break;
                end
                
                if (strcmpi(sesInfo.perfType, 'user specified settings'))
                    
                    [isPCAChanged, optionsChanged] = checkPCAOpts(oldSessInfo, sesInfo);
                    
                    if ((j == 1) && isPCAChanged)
                        status = 0;
                    end
                    
                    if ((j > 1) && optionsChanged)
                        status = 0;
                        groupNo = j;
                        dataSetNo = 1;
                    end
                    
                    if (~status)
                        break;
                    end
                end
                
                
                try
                    if (~strcmpi(oldSessInfo.group_ica_type, sesInfo.group_ica_type))
                        status = 0;
                        groupNo = j;
                        dataSetNo = 1;
                        break;
                    end
                catch
                end
                
                
                chkPcasig = 0;
                
                if (j == 1)
                    fS = whos('-file', fullfile(sesInfo.outputDir, [sesInfo.data_reduction_mat_file, '1-1.mat']));
                    fNames = cellstr(char(fS.name));
                    chkPcasig = isempty(strmatch('pcasig', fNames, 'exact'));
                    clear fNames fS;
                end
                
                chkSub = [];
                if (sesInfo.numReductionSteps ~= 1)
                    subjectsProcessed  = str2num(char(strrep(cellstr(tmpFiles(:, 1:end-4)), tmpFileN, '')));
                    [subjectsProcessed, subInds] = sort(subjectsProcessed(:)');
                    tmpFiles = tmpFiles(subInds, :);
                    dateInfo = dateInfo(subInds, :);
                    chkSub = (1:sesInfo.reduction(j).numOfGroupsAfterCAT);
                    chkSub(subjectsProcessed) = [];
                end
                
                if (chkPcasig)
                    chkSub = [];
                    load(fullfile(sesInfo.outputDir, [sesInfo.data_reduction_mat_file, '1-1.mat']), 'V', 'Lambda');
                    filesInfo(countStep).files = repmat(filesInfo(countStep).files(1), 1, length(filesInfo(countStep).files));
                    dateInfo = repmat(deblank(dateInfo(1, :)), length(filesInfo(countStep).files), 1);
                    if (size(V, 1) ~= sum(sesInfo.diffTimePoints))
                        status = 0;
                        groupNo = j;
                        dataSetNo = size(Lambda, 1) + 1;
                        break;
                    end
                end
                
                if (~isempty(chkSub))
                    status = 0;
                    groupNo = j;
                    dataSetNo = chkSub(1);
                    break;
                else
                    tmp = datenum(dateInfo);
                    tmp = tmp(:)';
                    filesInfo(countStep).value = [[filesInfo(countStep).value], tmp];
                    clear tmp;
                end
            end
            
            clear dateInfo tmpFiles tmpFileN subjectsProcessed;
            
        end
        % End of loop over groups
        
    end
    
    
    %% Calculate ICA step
    if (nStep == 4)
        
        fileNaming = sesInfo.ica_mat_file;
        
        filesInfo(countStep).step = [nStep; 1; 1];
        fileN = fullfile(sesInfo.outputDir, [fileNaming, '.mat']);
        filesInfo(countStep).files = {fileN};
        
        d = dir(fileN);
        
        if (isempty(d))
            status = 0;
            break;
        else
            tmp = datenum(datenum(d(1).date));
            tmp = tmp(:)';
            filesInfo(countStep).value = tmp;
            clear tmp;
        end
        
        
        %% Algorithm change and ICA options change
        try
            if (oldSessInfo.algorithm ~= sesInfo.algorithm)
                status = 0;
            else
                
                oldOptions = oldSessInfo.ICA_Options;
                newOptions = sesInfo.ICA_Options;
                
                if (length(newOptions) ~= length(oldOptions))
                    status = 0;
                end
                
                if (status)
                    % Loop over ICA options
                    for nOpts = 1:2:length(newOptions)
                        
                        tmpA = newOptions{nOpts + 1};
                        tmpOptsInds = strmatch(lower(newOptions{nOpts}), lower(oldOptions(1:2:end)), 'exact');
                        
                        if (isempty(tmpOptsInds))
                            status = 0;
                            break;
                        end
                        
                        tmpB = oldOptions{2*tmpOptsInds};
                        
                        if (ischar(tmpA))
                            tmpA = cellstr(tmpA);
                            tmpB = cellstr(tmpB);
                        end
                        
                        if (numel(tmpA) == length(tmpA))
                            tmpA = tmpA(:);
                            tmpB = tmpB(:);
                        end
                        
                        if (~isnumeric(tmpA))
                            tmpOptsInds = intersect(tmpA, tmpB);
                        else
                            if (length(tmpA) == numel(tmpA))
                                tmpOptsInds = find(prod(double(tmpA == tmpB), 2) == 1);
                            else
                                tmpOptsInds = find(abs(tmpA - tmpB) > eps);
                                if (~isempty(tmpOptsInds))
                                    break;
                                else
                                    continue;
                                end
                            end
                        end
                        
                        if (length(tmpOptsInds) ~= length(tmpA))
                            status = 0;
                            break;
                        end
                        
                    end
                    % End of loop over ICA options
                end
                
            end
            
        catch
        end
        
        if (~status)
            break;
        end
        
        if (~strcmpi(ica_types{sesInfo.algorithm}, 'iva-gl') && ~strcmpi(ica_types{sesInfo.algorithm}, 'iva-l') && ~strcmpi(ica_types{sesInfo.algorithm}, 'gig-ica') && ...
                ~strcmpi(ica_types{sesInfo.algorithm}, 'constrained ica (spatial)'))
            
            try
                
                if (oldSessInfo.which_analysis ~= sesInfo.which_analysis)
                    status = 0;
                end
                
                if (status)
                    
                    if (sesInfo.which_analysis == 1)
                        
                        if (NUM_RUNS_GICA ~= oldSessInfo.num_runs_gica)
                            status = 0;
                        end
                        
                    else
                        icasso_opts = struct('sel_mode', 'randinit', 'num_ica_runs', max([2, NUM_RUNS_GICA]));
                        old_icasso_opts = icasso_opts;
                        
                        if isfield(sesInfo, 'icasso_opts')
                            icasso_opts = sesInfo.icasso_opts;
                        end
                        
                        if isfield(oldSessInfo, 'icasso_opts')
                            old_icasso_opts = oldSessInfo.icasso_opts;
                        end
                        
                        if (~strcmpi(old_icasso_opts.sel_mode, icasso_opts.sel_mode) || (old_icasso_opts.num_ica_runs ~= icasso_opts.num_ica_runs))
                            status = 0;
                        end
                    end
                    
                end
                
            catch
            end
            
            
            if (~status)
                break;
            end
            
        end
        
    end
    
    %% Back-reconstruction step
    if (nStep == 5)
        
        
        fileNaming = sesInfo.back_reconstruction_mat_file;
        
        numBrs = sesInfo.numOfSub*sesInfo.numOfSess;
        filesInfo(countStep).step = [repmat(nStep, 1, numBrs); repmat(1, 1, numBrs); (1:numBrs)];
        allFiles =  repmat({''}, 1, numBrs);
        for nBrFile = 1:numBrs
            allFiles{nBrFile} = fullfile(sesInfo.outputDir, [fileNaming, num2str(nBrFile), '.mat']);
        end
        filesInfo(countStep).files = allFiles;
        
        
        %% Backrecon change
        backReconType = 'regular';
        old_backReconType = backReconType;
        try
            old_backReconType = oldSessInfo.backReconType;
            backReconType = sesInfo.backReconType;
        catch
        end
        
        if (~strcmpi(backReconType, old_backReconType))
            status = 0;
            break;
        end
        
        tmp = zeros(1, sesInfo.numOfSub*sesInfo.numOfSess);
        tmpFileN = fileNaming;
        tmpFiles = dir(fullfile(sesInfo.outputDir, [tmpFileN, '*.mat']));
        dateInfo = char(tmpFiles.date);
        tmpFiles = char(tmpFiles.name);
        
        pos = findstr(filesep, tmpFileN);
        relativePath = '';
        if (~isempty(pos))
            relativePath = tmpFileN(1:pos(end));
        end
        tmpFiles = [repmat(relativePath, size(tmpFiles, 1), 1), tmpFiles];
        
        if (isempty(tmpFiles))
            status = 0;
            dataSetNo = 1;
            break;
        else
            subjectsProcessed  = str2num(char(strrep(cellstr(tmpFiles(:, 1:end-4)), tmpFileN, '')));
            [subjectsProcessed, subInds] = sort(subjectsProcessed(:)');
            tmpFiles = tmpFiles(subInds, :);
            dateInfo = dateInfo(subInds, :);
            chkSub = (1:sesInfo.numOfSub*sesInfo.numOfSess);
            chkSub(subjectsProcessed) = [];
            if (~isempty(chkSub))
                status = 0;
                dataSetNo = chkSub(1);
                break;
            else
                tmp = datenum(dateInfo);
                tmp = tmp(:)';
                filesInfo(countStep).value = tmp;
            end
        end
        
        clear dateInfo tmpFiles tmpFileN subjectsProcessed tmp;
        
    end
    
    %% Calibrate components
    if (nStep == 6)
        
        fileNaming = sesInfo.calibrate_components_mat_file;
        tmp = zeros(1, sesInfo.numOfSub*sesInfo.numOfSess);
        tmpFiles = repmat({''}, 1, sesInfo.numOfSub*sesInfo.numOfSess);
        filesInfo(countStep).step = [repmat(nStep, 1, length(tmpFiles)); repmat(1, 1, length(tmpFiles)); (1:length(tmpFiles))];
        
        
        % Loop over datasets
        for k = 1:sesInfo.numOfSub*sesInfo.numOfSess
            subNum = ceil(k / sesInfo.numOfSess);
            sessNum = mod(k - 1, sesInfo.numOfSess) + 1;
            if ~(strcmpi(modalityType, 'eeg'))
                fileN = deblank(subjectICAFiles(subNum).ses(sessNum).name(end, :));
                zipFile = regexprep(fileN, '\_\d*\.\w{3}$', '_.zip');
                isZip = checkZip(zipContents, zipFile);
                if (isZip)
                    fileN = zipFile;
                end
                fileN = fullfile(sesInfo.outputDir, fileN);
            else
                fileN = fullfile(sesInfo.outputDir, [fileNaming, num2str(subNum), '-', num2str(sessNum), '.mat']);
            end
            
            d = dir(fileN);
            
            tmpFiles{k} = fileN;
            if (isempty(d))
                dataSetNo = k;
                status = 0;
                break;
            else
                tmp(k) = datenum(d(1).date);
                % Check detrend number and centering images
                if (strcmpi(modalityType, 'fmri'))
                    
                    if (k == 1)
                        
                        try
                            oldDetrendNo = oldSessInfo.detrendNumber;
                            if (oldDetrendNo ~= DETRENDNUMBER)
                                status = 0;
                                dataSetNo = 1;
                            end
                        catch
                        end
                        
                        if (status)
                            try
                                oldCenterIm = oldSessInfo.center_images;
                                
                                if (oldCenterIm ~= CENTER_IMAGES)
                                    status = 0;
                                    dataSetNo = 1;
                                end
                                
                            catch
                            end
                        end
                        
                    end
                    
                    if (~status)
                        break;
                    end
                    
                end
                % End for checking detrending and centering images
                
            end
            
        end
        % End of loop over data-sets
        
        filesInfo(countStep).files = tmpFiles;
        filesInfo(countStep).value = tmp;
        
        clear tmp tmpFiles;
        
        %% Scaling step change
        try
            if (oldSessInfo.scaleType ~= sesInfo.scaleType)
                status = 0;
            end
        catch
        end
        
        if (~status)
            break;
        end
        
        if (~strcmpi(modalityType, 'eeg'))
            
            %% Check if naming of files changed or not
            if (exist('old_subjectICAFiles', 'var'))
                if ispc
                    status = strcmpi(deblank(old_subjectICAFiles(1).ses(1).name(1, :)), deblank(subjectICAFiles(1).ses(1).name(1, :)));
                else
                    status = strcmp(deblank(old_subjectICAFiles(1).ses(1).name(1, :)), deblank(subjectICAFiles(1).ses(1).name(1, :)));
                end
                if (~status)
                    break;
                end
            end
            
        end
        
    end
    
    %% Group stats
    if ((nStep == 7) && (sesInfo.numOfSub*sesInfo.numOfSess > 1))
        
        if (sesInfo.numOfSub*sesInfo.numOfSess > 1)
            fileN = deblank(meanALL_ICAFile(1).name(end, :));
        else
            fileN = deblank(subjectICAFiles(1).ses(1).name(end, :));
        end
        
        if (~strcmpi(modalityType, 'eeg'))
            zipFile = regexprep(fileN, '\_\d*\.\w{3}$', '_.zip');
            isZip = checkZip(zipContents, zipFile);
            if (isZip)
                fileN = zipFile;
            end
            fileN = fullfile(sesInfo.outputDir, fileN);
        else
            fileN = fullfile(sesInfo.outputDir, [fileN(1:end-3), 'mat']);
        end
        
        filesInfo(countStep).step = [nStep; 1; 1];
        filesInfo(countStep).files = {fileN};
        
        d = dir(fileN);
        if (isempty(d))
            status = 0;
            break;
        else
            tmp = datenum(datenum(d(1).date));
            tmp = tmp(:)';
            filesInfo(countStep).value = tmp;
            clear tmp;
        end
        
        
        if (~strcmpi(modalityType, 'eeg'))
            
            %% Check if naming of files changed or not
            if (~isempty(meanALL_ICAFile(1).name))
                if (exist('old_meanALL_ICAFile', 'var'))
                    if (ispc)
                        status = strcmpi(deblank(old_meanALL_ICAFile(1).name(1, :)), deblank(meanALL_ICAFile(1).name(1, :)));
                    else
                        status = strcmp(deblank(old_meanALL_ICAFile(1).name(1, :)), deblank(meanALL_ICAFile(1).name(1, :)));
                    end
                    if (~status)
                        break;
                    end
                end
            end
            
            if (~isempty(meanICAFiles(1).name))
                if (exist('old_meanICAFiles', 'var'))
                    if (ispc)
                        status = strcmpi(deblank(old_meanICAFiles(1).name(1, :)), deblank(meanICAFiles(1).name(1, :)));
                    else
                        status = strcmp(deblank(old_meanICAFiles(1).name(1, :)), deblank(meanICAFiles(1).name(1, :)));
                    end
                    if (~status)
                        break;
                    end
                end
            end
            
            
            if (~isempty(tmapICAFiles(1).name))
                if (exist('old_tmapICAFiles', 'var'))
                    if (ispc)
                        status = strcmpi(deblank(old_tmapICAFiles(1).name(1, :)), deblank(tmapICAFiles(1).name(1, :)));
                    else
                        status = strcmp(deblank(old_tmapICAFiles(1).name(1, :)), deblank(tmapICAFiles(1).name(1, :)));
                    end
                    if (~status)
                        break;
                    end
                end
            end
            
            
            if (~isempty(stdICAFiles(1).name))
                if (exist('old_stdICAFiles', 'var'))
                    if (ispc)
                        status = strcmpi(deblank(old_stdICAFiles(1).name(1, :)), deblank(stdICAFiles(1).name(1, :)));
                    else
                        status = strcmp(deblank(old_stdICAFiles(1).name(1, :)), deblank(stdICAFiles(1).name(1, :)));
                    end
                    if (~status)
                        break;
                    end
                end
            end
            
        end
        
    end
    
    if (~status)
        break;
    end
    
end


steps = [filesInfo.step];
files = [filesInfo.files];

if (~status)
    
    inds = double(steps(1, :) == nStep) .* double(steps(2, :) == groupNo) .* double(steps(3, :) == dataSetNo);
    inds = find(inds == 1);
    
else
    
    disp('All files exist. Checking timestamp of files ...');
    
    inds = find(sign([0, diff([filesInfo.value])]) == -1);
    
end


if (isempty(inds))
    status = 1;
    inds = size(steps, 2);
else
    status = 0;
end

inds = inds(1);

inds2 = [inds - 1, inds];
inds2(inds2 == 0) = [];

for inds = inds2
    
    varsToCheck = {'pcasig', 'icasig', 'compSet', compSetFields{1}, compSetFields{1}};
    
    file = files{inds};
    
    if (strcmpi(file(end-3:end), '.zip') || strcmpi(file(end-3:end), '.img') || strcmpi(file(end-3:end), '.nii'))
        
        tmpNo = steps(3, inds);
        subNum = ceil(tmpNo / sesInfo.numOfSess);
        sessNum = mod(tmpNo - 1, sesInfo.numOfSess) + 1;
        
        if (sesInfo.numOfSub*sesInfo.numOfSess > 1)
            varsToCheck = {'pcasig', 'icasig', 'compSet', subjectICAFiles(subNum).ses(sessNum).name, meanALL_ICAFile(1).name};
        else
            varsToCheck = {'pcasig', 'icasig', 'compSet', subjectICAFiles(1).ses(1).name, subjectICAFiles(1).ses(1).name};
        end
        
    end
    
    tmpStep = steps(1, inds);
    
    if ((tmpStep == 3) && steps(2, inds) == 1)
        varsToCheck{1} = 'V';
    end
    
    status2 = checkFiles(file, tmpStep, varsToCheck);
    
    if (~status2)
        status = 0;
        break;
    end
    
end

if (~status)
    
    nStep = steps(1, inds);
    groupNo = steps(2, inds);
    dataSetNo = steps(3, inds);
    resume_info.stepsToRun = stepsToCheck(stepsToCheck >= nStep);
    resume_info.groupNo = groupNo;
    resume_info.dataSetNo = dataSetNo;
    
    return;
    
end


disp('Done checking. Analysis is complete.');
fprintf('\n\n');


function status = checkFiles(file, nStep, varsToCheck)
% Check if the files can be loaded or not

[outDir, fN, extn] = fileparts(file);

file = fullfile(outDir, [fN, extn]);

status = 0;

tmpRelDir = ['tmpZipDir_', fN];
tmpDir = fullfile(outDir, tmpRelDir);

try
    
    if (exist(deblank(file(end, :)), 'file'))
        
        if (strcmpi(extn, '.mat'))
            % MAT file
            %load(file, varsToCheck{nStep - 2});
            D = whos('-file', file);
            checkVar = strmatch(varsToCheck{nStep - 2}, char(D.name), 'exact');
            if (isempty(checkVar))
                error(['Variable ', varsToCheck{nStep - 2}, ' not found']);
            end
        elseif (strcmpi(extn, '.zip'))
            % Zip file
            if (~exist(tmpDir, 'dir'))
                mkdir(outDir, tmpRelDir);
            end
            icatb_unzip(file, tmpDir);
            tmpF = icatb_fullFile('directory', tmpDir, 'files', char(regexprep(cellstr(varsToCheck{nStep - 2}), ['.*\', filesep], '')));
            if (exist(deblank(tmpF(end, :)), 'file'))
                tmpF = icatb_rename_4d_file(tmpF);
                icatb_loadData(deblank(tmpF(end, :)));
            end
        else
            % IMG or NII files
            file = icatb_rename_4d_file(file);
            icatb_loadData(deblank(file(end, :)));
        end
        
        status = 1;
    end
    
catch
    
end


if (strcmpi(extn, '.zip'))
    try
        rmdir(tmpDir, 's') ;
    catch
    end
end

function isZip = checkZip(zipContents, outfile)

isZip = 0;
if (~isempty(zipContents.zipFiles))
    isZip = ~isempty(find(strcmpi(zipContents.zipFiles, outfile) == 1));
end


function [isPCAChanged, optionsChanged] = checkPCAOpts(oldSessInfo, sesInfo)
%% Detect change in PCA type and options
%

isPCAChanged = 0;
optionsChanged = 0;

oldPCAType = 'standard';
newPCAType = oldPCAType;

if (isfield(oldSessInfo, 'pcaType'))
    oldPCAType = lower(oldSessInfo.pcaType);
end

if (isfield(sesInfo, 'pcaType'))
    newPCAType = lower(sesInfo.pcaType);
end

oldSessInfo.pcaType = oldPCAType;
sesInfo.pcaType = newPCAType;

if (isfield(oldSessInfo, 'pca_opts'))
    if (iscell(oldSessInfo.pca_opts))
        oldSessInfo.pca_opts = oldSessInfo.pca_opts{1}.pca_opts;
    end
end

if (isfield(sesInfo, 'pca_opts'))
    if (iscell(sesInfo.pca_opts))
        sesInfo.pca_opts = sesInfo.pca_opts{1};
    end
end


oldSessInfo = icatb_check_pca_opts(oldSessInfo);
sesInfo = icatb_check_pca_opts(sesInfo);

if (~strcmpi(oldSessInfo.pcaType, sesInfo.pcaType))
    isPCAChanged = 1;
    return;
end

fN = fieldnames(sesInfo.pca_opts);

for i = 1:length(fN)
    a = getfield(oldSessInfo.pca_opts, fN{i});
    b = getfield(sesInfo.pca_opts, fN{i});
    if (ischar(a))
        optionsChanged = ~strcmpi(a, b);
    else
        optionsChanged = (a ~= b);
    end
    if (optionsChanged)
        break;
    end
end
