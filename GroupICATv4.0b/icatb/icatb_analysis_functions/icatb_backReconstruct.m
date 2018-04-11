function varargout = icatb_backReconstruct(sesInfo, statusHandle)
%% Back reconstruction is done to reconstruct individual subject image maps and time courses
% Input: sesInfo - structure containing all parameters for
% analysis

if(~exist('sesInfo','var'))
    P = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Parameter File', 'filter', '*param*.mat');
    if isempty(P)
        error('Parameter file is not selected for analysis');
    end
    [pathstr, fileName]=fileparts(P);
    outputDir = pathstr;
    sesInfo.outputDir = outputDir;
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
    error('Parameter file has not been initialized. You must Initialize the parameters file');
end

% Get modality, data title and compset fields
[modalityType, dataTitle, compSetFields] = icatb_get_modality;

icaStr = icatb_icaAlgorithm;

algorithmName = deblank(icaStr(sesInfo.algorithm, :));

if strcmpi(algorithmName, 'moo-icar')
    algorithmName = 'gig-ica';
end

if (strcmpi(algorithmName, 'gig-ica') || strcmpi(algorithmName, 'constrained ica (spatial)') )
    return;
end

useTemporalICA = 0;
if (strcmpi(modalityType, 'fmri'))
    try
        useTemporalICA = strcmpi(sesInfo.group_ica_type, 'temporal');
    catch
    end
end

%% Parallel info
num_workers = 4;
parallelMode = 'serial';

try
    parallelMode = sesInfo.parallel_info.mode;
catch
end

toolboxNames = ver;
parallelCluster= ~isempty(find(strcmpi(cellstr(char(toolboxNames.Name)), 'parallel computing toolbox') == 1));

if (strcmpi(parallelMode, 'parallel'))
    statusHandle = [];
end


% if (useTemporalICA)
%     return;
% end


appDataName = 'gica_waitbar_app_data';
if (~isempty(statusHandle))
    
    % get the status handles
    statusData = getappdata(statusHandle, appDataName);
    statusData.perCompleted = statusData.perCompleted + statusData.unitPerCompleted;
    setappdata(statusHandle, appDataName, statusData);
    set(statusHandle, 'name', [num2str(round(statusData.perCompleted*100)), '% analysis done']); waitbar(statusData.perCompleted, statusHandle);
    
end

verbose = (nargout < 2);

if (verbose)
    disp(' ');
    disp('---------------------------------------------------------------------');
    disp('STARTING BACK RECONSTRUCTION STEP');
    disp('---------------------------------------------------------------------');
end

preproc_type = 'remove mean per timepoint';

if (isfield(sesInfo, 'preproc_type'))
    preproc_type = sesInfo.preproc_type;
end


icain =[sesInfo.ica_mat_file, '.mat'];
backReconType = 'regular';

if (~useTemporalICA)
    
    if (~strcmpi(algorithmName, 'iva-gl') && ~strcmpi(algorithmName, 'iva-l'))
        
        if isfield(sesInfo, 'backReconType')
            backReconType = sesInfo.backReconType;
        end
        
        if strcmpi(backReconType, 'moo-icar')
            backReconType = 'gig-ica';
        end
        
        %         if ((sesInfo.numOfSub*sesInfo.numOfSess > 1) && (sesInfo.numReductionSteps == 1))
        %             backReconType = 'spatial-temporal regression';
        %         end
        %
        %         if (strcmpi(backReconType, 'str'))
        %             backReconType = 'spatial-temporal regression';
        %         end
        
        backReconOptions = cellstr(icatb_backReconOptions);
        backReconInd = strmatch(backReconType, lower(backReconOptions), 'exact');
        
        if (verbose)
            fprintf('\n');
            disp(['Using ', backReconOptions{backReconInd}, ' Back Reconstruction Approach ...']);
            fprintf('\n');
        end
        
        
        if (~strcmpi(backReconType, 'spatial-temporal regression') && ~strcmpi(backReconType, 'gig-ica'))
            load(fullfile(outputDir, icain), 'W');
        else
            load(fullfile(outputDir, icain), 'icasig');
        end
        
    else
        
        load(fullfile(outputDir, icain), 'W');
        
    end
    
else
    % Temporal ica signal
    load(fullfile(outputDir, icain), 'temporal_icasig');
    
end

if ~strcmpi(modalityType, 'eeg')
    compSetFields = {'ic', 'tc'};
else
    compSetFields = {'timecourse', 'topography'};
end

numReductionSteps = sesInfo.numReductionSteps;

dataSetsToRun = (1:sesInfo.numOfSub*sesInfo.numOfSess);

if (isfield(sesInfo, 'dataSetNo'))
    if (nargout < 2)
        dataSetsToRun(dataSetsToRun < sesInfo.dataSetNo) = [];
    else
        dataSetsToRun = sesInfo.dataSetNo;
    end
    sesInfo = rmfield(sesInfo, 'dataSetNo');
end

dataSetsToRun(dataSetsToRun > sesInfo.numOfSub*sesInfo.numOfSess) = [];

if (isempty(dataSetsToRun))
    error('Check the dataSetNo field passed');
end


if (length(dataSetsToRun) == 1)
    parallelMode = 'serial';
end

if (strcmpi(backReconType, 'moo-icar'))
    backReconType = 'gig-ica';
end

%% ICA back-reconstruction
if (~strcmpi(backReconType, 'spatial-temporal regression') && ~strcmpi(backReconType, 'gig-ica') && ~useTemporalICA)
    
    if (~strcmpi(deblank(icaStr(sesInfo.algorithm, :)), 'iva-gl') && ~strcmpi(deblank(icaStr(sesInfo.algorithm, :)), 'iva-l'))
        if (~isfield(sesInfo, 'tcInfo') || ~isfield(sesInfo, 'icInfo'))
            [tcInfo, icInfo] = icatb_groupBackReconInfo(sesInfo, W);
        else
            tcInfo = sesInfo.tcInfo;
            icInfo = sesInfo.icInfo;
        end
    end
    
    fS = whos('-file', fullfile(outputDir, [sesInfo.data_reduction_mat_file, '1-1.mat']));
    fNames = cellstr(char(fS.name));
    chkPcasig = isempty(strmatch('pcasig', fNames, 'exact'));
    
    if (chkPcasig)
        info = load(fullfile(outputDir, [sesInfo.data_reduction_mat_file, '1-', num2str(1), '.mat']));
        VStacked = info.V;
        LambdaStacked = info.Lambda;
        clear info;
    end
    
    %% Load first level data reduction files
    for j = dataSetsToRun
        
        if (verbose)
            disp(['Back reconstructing  set ', num2str(j)]);
        end
        %% PCAInfo
        if (~chkPcasig)
            [dewhiteM, whiteM, pcasig] = load_vars(fullfile(outputDir, [sesInfo.data_reduction_mat_file, '1-', num2str(j), '.mat']));
        else
            if (j == 1)
                startT = 1;
            else
                startT = sum(sesInfo.diffTimePoints(1:j-1)) + 1;
            end
            endT = sum(sesInfo.diffTimePoints(1:j));
            [whiteM, dewhiteM] = get_pca_info(VStacked(startT:endT, :), diag(LambdaStacked(j, :)));
            if (verbose)
                disp(['Loading  dataset ', num2str(j)]);
            end
            dat = icatb_read_data(sesInfo.inputFiles(j).name, [], sesInfo.mask_ind);
            % Call pre-processing function
            dat = icatb_preproc_data(dat, preproc_type, verbose);
            % Remove mean per timepoint
            dat = icatb_remove_mean(dat, 0);
            pcasig = dat*whiteM';
            clear wM dat;
        end
        
        pcasig = pcasig';
        
        if (~strcmpi(deblank(icaStr(sesInfo.algorithm, :)), 'iva-gl') && ~strcmpi(deblank(icaStr(sesInfo.algorithm, :)), 'iva-l'))
            ic = W*icInfo{j}*pcasig;
            tc = dewhiteM*tcInfo{j}*pinv(W);
        else
            % IVA
            ic = squeeze(W(:, :, j))*pcasig;
            tc = dewhiteM*pinv(squeeze(W(:, :, j)));
        end
        
        if (sesInfo.numOfSub*sesInfo.numOfSess == 1)
            tc = tc';
        end
        
        compSet = struct(compSetFields{1}, ic, compSetFields{2}, tc);
        
        if (verbose)
            disp(['-done back reconstructing  set ', num2str(j)]);
        end
        
        if ~isempty(statusHandle)
            
            % get the status handles
            statusData = getappdata(statusHandle, appDataName);
            statusData.perCompleted = statusData.perCompleted + (statusData.unitPerCompleted / length(dataSetsToRun));
            setappdata(statusHandle, appDataName, statusData);
            set(statusHandle, 'name', [num2str(round(statusData.perCompleted*100)), '% analysis done']); waitbar(statusData.perCompleted, statusHandle);
            
        end
        
        drawnow;
        
        %% Save subject components
        subFile = [sesInfo.back_reconstruction_mat_file, num2str(j), '.mat'];
        
        if (verbose)
            msgString = ['-saving back reconstructed ica data for set ', num2str(j),' -> ',subFile];
            disp(msgString);
            drawnow;
            icatb_save(fullfile(outputDir, subFile), 'compSet');
            clear compSet;
        end
        
    end
    %% End for loading first level data reduction files
    
else
    
    %% Spatial-temporal regression,  GIG-ICA, temporal ica
    if (strcmpi(parallelMode, 'serial'))
        if (nargout == 0)
            icatb_parBackReconstruct(sesInfo, dataSetsToRun, verbose, statusHandle);
        else
            [ddd, compSet] = icatb_parBackReconstruct(sesInfo, dataSetsToRun, verbose, statusHandle);
        end
    else
        if (parallelCluster)
            icatb_parBackReconstruct_cluster(sesInfo, dataSetsToRun, verbose);
        else
            parBackRecon(sesInfo, dataSetsToRun, verbose);
        end
    end
    
end


[pp, fileName] = fileparts(sesInfo.userInput.param_file);
drawnow;

icatb_save(fullfile(outputDir, [fileName, '.mat']), 'sesInfo');


if (verbose)
    disp('---------------------------------------------------------------------');
    disp('DONE WITH BACK RECONSTRUCTION STEP');
    disp('---------------------------------------------------------------------');
    disp(' ');
end

varargout{1} = sesInfo;

if (nargout == 2)
    varargout{2} = compSet;
end

if (~isempty(statusHandle))
    
    % get the status handles
    statusData = getappdata(statusHandle, appDataName);
    statusData.perCompleted = statusData.perCompleted + statusData.unitPerCompleted;
    setappdata(statusHandle, appDataName, statusData);
    set(statusHandle, 'name', [num2str(round(statusData.perCompleted*100)), '% analysis done']); waitbar(statusData.perCompleted, statusHandle);
    
end


function varargout = load_vars(fileIn)
% Load dewhiteM and whiteM variables

load(fileIn);

if (~exist('dewhiteM', 'var'))
    if (numel(Lambda) == length(Lambda))
        Lambda = diag(Lambda);
    end
    [whiteM, dewhiteM] = get_pca_info(V, Lambda);
end

if (~exist('whiteM', 'var'))
    whiteM = pinv(dewhiteM);
end

varargout{1} = dewhiteM;
varargout{2} = whiteM;

if (exist('pcasig', 'var'))
    varargout{3} = pcasig;
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



function parBackRecon(sesInfo, dataSetsToRun, verbose)
%% Back-reconstruct components using independent sessions
%

cValue = clock;
% remove a minute
cValue(end-1) = max([0, cValue(end-1) - 1]);
currentTime = datenum(cValue);

outputDir = sesInfo.outputDir;

num_workers = 4;
try
    num_workers = sesInfo.parallel_info.num_workers;
catch
end

[ddd, paramFile, extn] = fileparts(sesInfo.userInput.param_file);
paramFile = fullfile(sesInfo.outputDir, [paramFile, extn]);

giftPath = fileparts(which('gift.m'));
dummyScriptPath = fullfile(giftPath, 'icatb_parallel_files', 'icatb_dummyScript.m');
totalSubjectsIn = length(dataSetsToRun);
if (num_workers > totalSubjectsIn)
    num_workers = totalSubjectsIn;
end

backReconFiles = cellstr(strcat(fullfile(outputDir, [sesInfo.back_reconstruction_mat_file, '1-']), icatb_numberToString(dataSetsToRun(:)), '.mat'));
increments = ceil(totalSubjectsIn/num_workers);
eW = 0;
%currentTime = datenum(clock);
for nF = 1:num_workers
    sW = eW + 1;
    eW = eW + increments;
    eW = min([totalSubjectsIn, eW]);
    tmpDataSetsToRun = dataSetsToRun(sW:eW);
    % Run separate matlab sessions in background mode
    % (Dummyscript i.e., no text required to run in linux OS)
    commandStr = ['!matlab -nodesktop -nosplash -r "addpath(genpath(''', giftPath, '''));icatb_parBackReconstruct(''', paramFile, ''',[', ...
        num2str(tmpDataSetsToRun), '],', num2str(verbose), ''');exit" < "', dummyScriptPath, '" &'];
    eval(commandStr);
end
icatb_waitForTaskCompletion(backReconFiles, currentTime);