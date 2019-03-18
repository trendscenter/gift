function sesInfo = icatb_dataReduction(sesInfo, statusHandle)
%% Reduces each subject's data using pca. If more than one reduction
% step is specified by parameters then a concatenation step is done along
% with another pca step.

% Input: sesInfo - structure containing all parameters necessary for group
% ica analysis.

if (~exist('sesInfo', 'var'))
    P = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Parameter File', 'filter', '*param*.mat');
    if (isempty(P))
        error('Parameter file is not selected for analysis');
    end
    outputDir = fileparts(P);
    % Make sure parameter file exists
    load(P);
    if (~exist('sesInfo', 'var'))
        error(['The selected file ', P, ' does not contain the sesInfo variable']);
    end
else
    outputDir = sesInfo.outputDir;
end

sesInfo.outputDir = outputDir;

drawnow;

if (~exist('statusHandle', 'var'))
    statusHandle = [];
end

if (sesInfo.isInitialized == 0)
    icatb_error('Parameter file has not been initialized');
end

modalityType = icatb_get_modality;

appDataName = 'gica_waitbar_app_data';

algoList = icatb_icaAlgorithm;
algoName = deblank(algoList(sesInfo.algorithm, :));
if (strcmpi(algoName, 'moo-icar'))
    algoName = 'gig-ica';
end

if (strcmpi(algoName, 'gig-ica') || strcmpi(algoName, 'constrained ica (spatial)'))
    return;
end


if ~isempty(statusHandle)
    % get the status handles
    statusData = getappdata(statusHandle, appDataName);
    statusData.perCompleted = statusData.perCompleted + statusData.unitPerCompleted;
    setappdata(statusHandle, appDataName, statusData);
    set(statusHandle, 'name', [num2str(round(statusData.perCompleted*100)), '% analysis done']); waitbar(statusData.perCompleted, statusHandle);
end


%% Mask Indices
mask_ind = sesInfo.mask_ind;
numVoxels = length(mask_ind);

preproc_type = 'remove mean per timepoint';

if (isfield(sesInfo, 'preproc_type'))
    preproc_type = sesInfo.preproc_type;
end

groupPCAOpts = char('Subject Specific', 'Grand Mean');

group_pca_type = 'subject specific';
if isfield(sesInfo, 'group_pca_type')
    group_pca_type = sesInfo.group_pca_type;
end

if (isnumeric(group_pca_type))
    group_pca_type = lower(deblank(groupPCAOpts(group_pca_type, :)));
    sesInfo.group_pca_type = group_pca_type;
end

doGrandMeanPCA = (strcmpi(group_pca_type, 'grand mean') && (sesInfo.numOfSub*sesInfo.numOfSess > 1));

disp(' ');
disp('------------------------------------------------------------------------------------');
disp('STARTING DATA REDUCTION (PRINCIPAL COMPONENTS ANALYSIS)');
disp('------------------------------------------------------------------------------------');

if (~isfield(sesInfo, 'dataType'))
    sesInfo.dataType = 'real';
end

[sesInfo, complexInfo] = icatb_name_complex_images(sesInfo, 'read');
numReductionSteps = sesInfo.numReductionSteps;

%% Group PCA Options
gpca_opts = getPCAOpts(sesInfo);

%% Reduction Steps To Run
reductionStepsToRun = (1:numReductionSteps);
if (isfield(sesInfo, 'reductionStepsToRun'))
    reductionStepsToRun = sesInfo.reductionStepsToRun;
    sesInfo = rmfield(sesInfo, 'reductionStepsToRun');
end

if (max(reductionStepsToRun) > 2)
    error('GIFT no longer runs three data-reduction steps. STP is a variation of 3 step PCA. Use STP instead.');
end

reductionStepsToRun(reductionStepsToRun > numReductionSteps) = [];


try
    if (strcmpi(algoName, 'iva-gl') || strcmpi(algoName, 'iva-l'))
        reductionStepsToRun = 1;
    end
catch
end


useTemporalICA = 0;

if (strcmpi(modalityType, 'fmri'))
    
    try
        useTemporalICA = strcmpi(sesInfo.group_ica_type, 'temporal');
    catch
    end
    
    if (useTemporalICA)
        reductionStepsToRun = 1;
        if (strcmpi(algoName, 'iva-gl') || strcmpi(algoName, 'iva-l') || strcmpi(algoName, 'gig-ica') || strcmpi(algoName, 'constrained ica (spatial)'))
            error(['Temporal ICA cannot be run using the algorithm ', algoName]);
        end
    end
    
end


intermediatePCA = 1;

if (~strcmpi(algoName, 'iva-gl') && ~strcmpi(algoName, 'iva-l'))
    if ((numReductionSteps == 1) && (sesInfo.numOfSub*sesInfo.numOfSess > 1))
        intermediatePCA = 0;
    end
end


mask = zeros(sesInfo.HInfo.DIM);
mask(sesInfo.mask_ind) = 1;
mask = (mask == 1);

%% Parallel info
num_workers = 4;
parallelMode = 'serial';

try
    parallelMode = sesInfo.parallel_info.mode;
catch
end

try
    num_workers = sesInfo.parallel_info.num_workers;
catch
end

conserve_disk_space = getConserveDiskSpaceInfo(sesInfo, gpca_opts, algoName);

toolboxNames = ver;
parallelCluster = ~isempty(find(strcmpi(cellstr(char(toolboxNames.Name)), 'parallel computing toolbox') == 1));

runGroupICAParallel = 0;
if (strcmpi(parallelMode, 'parallel') && parallelCluster)
    runGroupICAParallel = 1;
end

if (intermediatePCA)
    %% Whiten data and do second data-reduction if specified
    
    %% Compute mean of data-sets and get the tranformation vector
    if (doGrandMeanPCA)
        VM = grandMeanPCA(sesInfo, preproc_type, gpca_opts{1}.pca_opts.precision);
    end
    
    if (strcmpi(parallelMode, 'parallel') && any(reductionStepsToRun == 1))
        sesInfo = parFirstPCA(sesInfo, conserve_disk_space, num_workers);
        reductionStepsToRun(reductionStepsToRun==1)=[];
    end
    
    for nR = reductionStepsToRun
        
        numOfPC = sesInfo.reduction(nR).numOfPCAfterReduction;
        disp(['--Extracting principal components for data reduction( time #', num2str(nR), ' )']);
        
        pcaType = gpca_opts{nR}.pcaType;
        pca_opts = gpca_opts{nR}.pca_opts;
        
        if (nR == 1)
            
            %% First data reduction
            numOfGroupsAfter = sesInfo.reduction(nR).numOfGroupsAfterCAT;
            dataSetNo = 1;
            if (isfield(sesInfo, 'dataSetNo'))
                dataSetNo = sesInfo.dataSetNo;
                sesInfo = rmfield(sesInfo, 'dataSetNo');
            end
            
            dataSetsToRun = (1:numOfGroupsAfter);
            dataSetsToRun(dataSetsToRun < dataSetNo) = [];
            
            for nF = dataSetsToRun
                
                fprintf('\n');
                subNum = ceil(nF/sesInfo.numOfSess);
                ses = mod(nF-1, sesInfo.numOfSess) + 1;
                if (~strcmpi(modalityType, 'smri'))
                    if (~doGrandMeanPCA)
                        msg_string = ['--Doing pca on Subject #', num2str(subNum), ' Session #', num2str(ses)];
                    else
                        msg_string = ['-- Projecting Subject #', num2str(subNum), ' Session #', num2str(ses), ' to the eigen space of the mean of data-sets'];
                    end
                else
                    msg_string = 'Doing pca';
                end
                disp(msg_string);
                
                tmpPreprocType = 'none';
                if strcmpi(modalityType, 'smri')
                    tmpPreprocType = preproc_type;
                end
                
                if ~strcmpi(modalityType, 'smri')
                    
                    %% Load data
                    data = preprocData(sesInfo.inputFiles(nF).name, mask_ind, preproc_type, pca_opts.precision);
                    
                    %% Project data on to the eigen space of the mean data
                    if (doGrandMeanPCA)
                        data = data*VM;
                    end
                    
                else
                    
                    data  = cellstr(sesInfo.inputFiles(nF).name);
                    
                end
                
                if (conserve_disk_space ~= 1)
                    %% Do PCA with whitening
                    [pcasig, dewhiteM, Lambda, V, whiteM] = icatb_calculate_pca(data, numOfPC, 'type', pcaType, 'mask', mask, 'whiten', 1, 'verbose', 1, 'preproc_type', tmpPreprocType, ...
                        'pca_options', pca_opts);
                    
                else
                    %% Do pca
                    [V, Lambda] = icatb_calculate_pca(data, numOfPC, 'type', pcaType, 'mask', mask, 'whiten', 0, 'verbose', 1, 'preproc_type', tmpPreprocType, ...
                        'pca_options', pca_opts);
                    
                end
                
                %% Apply the transformations to Eigen vectors, dewhitening
                % and whitening matrices
                if (doGrandMeanPCA)
                    V = VM*V;
                    if (exist('dewhiteM', 'var'))
                        dewhiteM = VM*dewhiteM;
                        whiteM = pinv(dewhiteM);
                    end
                end
                
                if (conserve_disk_space ~= 1)
                    pcaout = [sesInfo.data_reduction_mat_file, '1-', num2str(nF), '.mat'];
                else
                    pcaout = [sesInfo.data_reduction_mat_file, '1-', num2str(1), '.mat'];
                end
                
                pcaout = fullfile(outputDir, pcaout);
                
                if (conserve_disk_space == 1)
                    Lambda = diag(Lambda);
                    Lambda = Lambda(:)';
                    if (nF > 1)
                        info = load(pcaout, 'V', 'Lambda');
                        info.V = info.V(1:sum(sesInfo.diffTimePoints(1:nF-1)), :);
                        info.Lambda = info.Lambda(1:nF-1, :);
                        V = [info.V; V];
                        Lambda = [info.Lambda; Lambda];
                        clear info;
                    end
                end
                
                if (conserve_disk_space ~= 1)
                    icatb_save(pcaout, 'pcasig', 'dewhiteM', 'V', 'Lambda', 'whiteM');
                else
                    icatb_save(pcaout, 'V', 'Lambda');
                end
                
                %if (strcmpi(modalityType, 'smri'))
                 %   sesInfo.pca_variances = compute_var(data, V, %Lambda);
                %end
                
                fprintf('\n\n');
                
            end
        else
            %% Second data reduction
            files = cell(1, length(sesInfo.inputFiles));
            for nF = 1:length(files)
                files{nF} = fullfile(outputDir, [sesInfo.data_reduction_mat_file, '1-', num2str(nF), '.mat']);
            end
            
            if (conserve_disk_space == 1)
                files = projectEigVecOnData(sesInfo, preproc_type, pca_opts.precision);
            end
            
            if (~runGroupICAParallel)
                [pcasig, dewhiteM, Lambda, V, whiteM] = icatb_calculate_pca(files, numOfPC, 'type', pcaType, 'whiten', 1, 'verbose', 1, 'preproc_type', 'none', ...
                    'pca_options', pca_opts);
            else
                [pcasig, dewhiteM, Lambda, V, whiteM] = icatb_parCalculatePCA(files, numOfPC, 'type', pcaType, 'whiten', 1, 'verbose', 1, 'preproc_type', 'none', ...
                    'pca_options', pca_opts);
            end
            pcaout = [sesInfo.data_reduction_mat_file, num2str(nR), '-1.mat'];
            pcaout = fullfile(outputDir, pcaout);
            icatb_save(pcaout, 'pcasig', 'dewhiteM', 'V', 'Lambda', 'whiteM');
            fprintf('\n\n');
            
        end
        
        
        if ~isempty(statusHandle)
            % get the status handles
            statusData = getappdata(statusHandle, appDataName);
            statusData.perCompleted = statusData.perCompleted + (statusData.unitPerCompleted / length(reductionStepsToRun));
            setappdata(statusHandle, appDataName, statusData);
            set(statusHandle, 'name', [num2str(round(statusData.perCompleted*100)), '% analysis done']);
            waitbar(statusData.perCompleted, statusHandle);
            
        end
        
    end
    
    disp(['Done with data reduction( time # ', num2str(nR),')' ]);
    disp(' ');
    
else
    %% Do group pca on pre-processed data
    disp('Group pca is done by stacking data across subjects ...');
    pcaType = gpca_opts{1}.pcaType;
    pca_opts = gpca_opts{1}.pca_opts;
    numOfPC = sesInfo.reduction(1).numOfPCAfterReduction;
    files = cell(1, length(sesInfo.inputFiles));
    for nF = 1:length(files)
        files{nF} = sesInfo.inputFiles(nF).name;
    end
    
    varToLoad = '';
    dims = [length(sesInfo.mask_ind), sum(sesInfo.diffTimePoints)];
    
    pcaout = [sesInfo.data_reduction_mat_file, '1-1.mat'];
    pcaout = fullfile(outputDir, pcaout);
    if (~useTemporalICA)
        if (~runGroupICAParallel)
            [pcasig, dewhiteM, Lambda, V, whiteM] = icatb_calculate_pca(files, numOfPC, 'type', pcaType, 'mask', mask, 'whiten', 1, 'verbose', 1, 'preproc_type', preproc_type, ...
                'pca_options', pca_opts, 'varToLoad', varToLoad, 'dims', dims);
        else
            [pcasig, dewhiteM, Lambda, V, whiteM] = icatb_parCalculatePCA(files, numOfPC, 'type', pcaType, 'mask', mask, 'whiten', 1, 'verbose', 1, 'preproc_type', preproc_type, ...
                'pca_options', pca_opts, 'varToLoad', varToLoad, 'dims', dims);
        end
        icatb_save(pcaout, 'pcasig', 'dewhiteM', 'V', 'Lambda', 'whiteM');
    else
        if (~runGroupICAParallel)
            [V, Lambda] = icatb_calculate_pca(files, numOfPC, 'type', pcaType, 'mask', mask, 'whiten', 0, 'verbose', 1, 'preproc_type', preproc_type, ...
                'pca_options', pca_opts, 'varToLoad', varToLoad, 'dims', dims);
        else
            [V, Lambda] = icatb_parCalculatePCA(files, numOfPC, 'type', pcaType, 'mask', mask, 'whiten', 0, 'verbose', 1, 'preproc_type', preproc_type, ...
                'pca_options', pca_opts, 'varToLoad', varToLoad, 'dims', dims);
        end
        icatb_save(pcaout, 'V', 'Lambda');
    end
    fprintf('\n\n');
    
    if ~isempty(statusHandle)
        % get the status handles
        statusData = getappdata(statusHandle, appDataName);
        statusData.perCompleted = statusData.perCompleted + statusData.unitPerCompleted;
        setappdata(statusHandle, appDataName, statusData);
        set(statusHandle, 'name', [num2str(round(statusData.perCompleted*100)), '% analysis done']);
        waitbar(statusData.perCompleted, statusHandle);
        
    end
    
end


%% Save parameter file
[pp, fileName] = fileparts(sesInfo.userInput.param_file);
drawnow;
icatb_save(fullfile(outputDir, [fileName, '.mat']), 'sesInfo');

disp('------------------------------------------------------------------------------------');
disp('ENDING DATA REDUCTION (PRINCIPAL COMPONENTS ANALYSIS)');
disp('------------------------------------------------------------------------------------');
disp('');


if ~isempty(statusHandle)
    % get the status handles
    statusData = getappdata(statusHandle, appDataName);
    statusData.perCompleted = statusData.perCompleted + statusData.unitPerCompleted;
    setappdata(statusHandle, appDataName, statusData);
    set(statusHandle, 'name', [num2str(round(statusData.perCompleted*100)), '% analysis done']); waitbar(statusData.perCompleted, statusHandle);
    
end

function pca_opts = getPCAOpts(sesInfo)
%% Get PCA opts
%

isOptsCell = 0;
try
    isOptsCell = iscell(sesInfo.pca_opts);
catch
end

if (~isOptsCell)
    
    if (isfield(sesInfo, 'pca_opts'))
        tmp_pca_opts = sesInfo.pca_opts;
    end
    
    %% Group PCA Options
    pcaType = 'standard';
    if (isfield(sesInfo, 'pcaType'))
        pcaType = sesInfo.pcaType;
    end
    
    sesInfo.pcaType = pcaType;
    sesInfo = icatb_check_pca_opts(sesInfo);
    pcaType = sesInfo.pcaType;
    tmp_pca_opts = sesInfo.pca_opts;
    
    pca_opts{1} = struct('pcaType', pcaType, 'pca_opts', tmp_pca_opts);
    
else
    
    pca_opts = sesInfo.pca_opts;
    for nP = 1:length(pca_opts)
        pca_opts{nP} = icatb_check_pca_opts(pca_opts{nP});
    end
    
end

if (length(pca_opts) ~= sesInfo.numReductionSteps)
    pca_opts = repmat(pca_opts(1), 1, sesInfo.numReductionSteps);
end



function VM = grandMeanPCA(sesInfo, preProcType, precisionType)
%% Compute grand mean
%

disp('Computing mean of data-sets ...');

meanData = zeros(length(sesInfo.mask_ind), min(sesInfo.diffTimePoints));

for nD = 1:length(sesInfo.inputFiles)
    subNum = ceil(nD/sesInfo.numOfSess);
    ses = mod(nD-1, sesInfo.numOfSess) + 1;
    disp(['Loading Subject #', num2str(subNum), ' Session #', num2str(ses)]);
    tmp = preprocData(sesInfo.inputFiles(nD).name, sesInfo.mask_ind, preProcType, precisionType);
    meanData = meanData + tmp(:, 1:size(meanData, 2));
    clear tmp;
end

meanData = meanData/length(sesInfo.inputFiles);

disp('Done');

fprintf('\n');

disp('Calculating PCA on mean of data-sets ...');

[VM, LambdaM] = icatb_calculate_pca(meanData, sesInfo.reduction(1).numOfPCAfterReduction, 'preproc_type', 'none', 'whiten', 0, 'type', 'standard');

pcaFile = fullfile(sesInfo.outputDir, [sesInfo.data_reduction_mat_file, 'mean.mat']);

icatb_save(pcaFile, 'LambdaM', 'VM');

fprintf('Done\n');


function data = preprocData(fileN, mask_ind, preProcType, precisionType)
%% Preprocess data
%

% Load data
data = icatb_read_data(fileN, [], mask_ind, precisionType);

% Call pre-processing function
data = icatb_preproc_data(data, preProcType);

if (~strcmpi(preProcType, 'remove mean per timepoint'))
    % Remove mean per timepoint
    data = icatb_remove_mean(data, 1);
end


function sesInfo = parFirstPCA(sesInfo, conserve_disk_space, num_workers)
%% Run pca in parallel mode
%

cValue = clock;
% remove a minute
cValue(end-1) = max([0, cValue(end-1) - 1]);
currentTime = datenum(cValue);

[dd, paramFile, extn] = fileparts(sesInfo.userInput.param_file);
paramFile = fullfile(sesInfo.outputDir, [paramFile, extn]);
numOfGroupsAfter = sesInfo.reduction(1).numOfGroupsAfterCAT;
dataSetNo = 1;
if (isfield(sesInfo, 'dataSetNo'))
    dataSetNo = sesInfo.dataSetNo;
    sesInfo = rmfield(sesInfo, 'dataSetNo');
end

dataSetsToRun = (1:numOfGroupsAfter);
dataSetsToRun(dataSetsToRun < dataSetNo) = [];

toolboxNames = ver;
parallelCluster = ~isempty(find(strcmpi(cellstr(char(toolboxNames.Name)), 'parallel computing toolbox') == 1));

disp('--Extracting principal components for data reduction( time #1)');
fprintf('\n');

if (parallelCluster)
    %% Use PCT
    icatb_parFirstPCA_cluster(sesInfo, dataSetsToRun, conserve_disk_space);
else
    %% Run in background mode
    giftPath = fileparts(which('gift.m'));
    dummyScriptPath = fullfile(giftPath, 'icatb_parallel_files', 'icatb_dummyScript.m');
    totalSubjectsIn = length(dataSetsToRun);
    if (totalSubjectsIn < num_workers)
        num_workers = totalSubjectsIn;
    end
    pcaFiles = cellstr(strcat(fullfile(sesInfo.outputDir, [sesInfo.data_reduction_mat_file, '1-']), icatb_numberToString(dataSetsToRun(:)), '.mat'));
    increments = ceil(totalSubjectsIn/num_workers);
    eW = 0;
    for nF = 1:num_workers
        sW = eW + 1;
        eW = eW + increments;
        eW = min([totalSubjectsIn, eW]);
        tmpDataSetsToRun = dataSetsToRun(sW:eW);
        % Run separate matlab sessions in background mode
        % (Dummyscript i.e., no text required to run in linux OS)
        commandStr = ['!matlab -nodesktop -nosplash -r "addpath(genpath(''', giftPath, '''));icatb_parFirstPCA(''', paramFile,  ''',[', ...
            num2str(tmpDataSetsToRun), '],0);exit" < "', dummyScriptPath, '" &'];
        eval(commandStr);
    end
    
    icatb_waitForTaskCompletion(pcaFiles, currentTime);
    
    LambdaV = []; VV = [];
    pcaFiles = cellstr(strcat(fullfile(sesInfo.outputDir, [sesInfo.data_reduction_mat_file, '1-']), icatb_numberToString((1:numOfGroupsAfter)'), '.mat'));
    
    if (conserve_disk_space)
        for nO = 1:length(pcaFiles)
            load(pcaFiles{nO}, 'Lambda', 'V');
            Lambda = diag(Lambda);
            Lambda = Lambda(:)';
            LambdaV = [LambdaV;Lambda];
            VV = [VV;V];
        end
        delete(fullfile(sesInfo.outputDir, [sesInfo.data_reduction_mat_file, '1-*mat']));
        pcaout = [sesInfo.data_reduction_mat_file, '1-1.mat'];
        pcaout = fullfile(sesInfo.outputDir, pcaout);
        icatb_parSave(pcaout, {VV, LambdaV}, {'V', 'Lambda'});
    end
end

disp('Done extracting principal components for data reduction( time #1)');
fprintf('\n');


function conserve_disk_space = getConserveDiskSpaceInfo(sesInfo, gpca_opts, algoName)

conserve_disk_space = 0;

% Write eigen vectors for 2 data reductions with stack data set to on
if (sesInfo.numReductionSteps == 2)
    stackData = 'yes';
    try
        stackData = gpca_opts{2}.pca_opts.stack_data;
    catch
    end
    
    if (strcmpi(gpca_opts{2}.pcaType, 'svd'))
        stackData = 'yes';
    end
    
    if (isfield(sesInfo, 'conserve_disk_space') && strcmpi(stackData, 'yes'))
        conserve_disk_space = sesInfo.conserve_disk_space;
    end
end

if ((sesInfo.numOfSub*sesInfo.numOfSess == 1) && (strcmpi(algoName, 'iva-gl') || strcmpi(algoName, 'iva-l')))
    try
        conserve_disk_space = sesInfo.conserve_disk_space;
    catch
    end
end


function [whiteningMatrix, dewhiteningMatrix] = get_pca_info(V, Lambda)
%% Get Whitening and de-whitening matrix
%
% Inputs:
% 1. V - Eigen vectors
% 2. Lambda - Eigen values diagonal matrix
%
% Outputs:
% 1. whiteningMatrix - Whitening matrix
% 2. dewhiteningMatrix - Dewhitening matrix
%


whiteningMatrix = sqrtm(Lambda) \ V';
dewhiteningMatrix = V * sqrtm(Lambda);


function data = projectEigVecOnData(sesInfo, preProcType, precisionType)
%% Project eigen vectors on to data
%

outputDir = sesInfo.outputDir;
mask_ind = sesInfo.mask_ind;

pcain = [sesInfo.data_reduction_mat_file, '1-1.mat'];
load(fullfile(outputDir, pcain), 'V', 'Lambda');
VStacked = V;
LambdaStacked = Lambda;
clear V Lambda;

data = zeros(length(mask_ind), length(sesInfo.inputFiles)*size(VStacked, 2));

endTp = 0;
endPC = 0;
for nF = 1:length(sesInfo.inputFiles)
    subNum = ceil(nF/sesInfo.numOfSess);
    ses = mod(nF - 1, sesInfo.numOfSess) + 1;
    disp(['Loading subject #', num2str(subNum), ' session #', num2str(ses), ' ...']);
    
    startTp = endTp + 1;
    endTp = endTp + sesInfo.diffTimePoints(nF);
    wM = get_pca_info(VStacked(startTp:endTp, :), diag(LambdaStacked(nF, :)));
    dat = preprocData(sesInfo.inputFiles(nF).name, mask_ind, preProcType, precisionType);
    dat = dat*wM';
    
    startPC = endPC + 1;
    endPC = endPC + size(VStacked, 2);
    
    data(:, startPC:endPC) = dat;
    
    clear dat;
    
end


function pf = compute_var(data, V, Lambda)
% Compute percent variances
[whiteM, dewhiteM] = get_pca_info(V, Lambda);
pcasig = data*whiteM';
pf = icatb_compute_var_evd(data, dewhiteM, pcasig);
