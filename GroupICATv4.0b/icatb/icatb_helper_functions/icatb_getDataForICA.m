function data = icatb_getDataForICA(sesInfo, algorithmName, statusHandle)
%% Load data before applying ICA (reduced data from pca excluding spatially constrained ica algorithms)
%

if (ischar(sesInfo))
    load (sesInfo);
end


appDataName = 'gica_waitbar_app_data';

if (~exist('statusHandle', 'var'))
    statusHandle = [];
end

modalityType = 'fMRI';
try
    modalityType = sesInfo.modality;
catch
end

outputDir = sesInfo.outputDir;
mask_ind = sesInfo.mask_ind;

useTemporalICA = 0;
if (strcmpi(modalityType, 'fmri'))
    try
        useTemporalICA = strcmpi(sesInfo.group_ica_type, 'temporal');
    catch
    end
end

if (~strcmpi(algorithmName, 'iva-gl') && ~strcmpi(algorithmName, 'iva-l'))
    
    pcain = [sesInfo.data_reduction_mat_file,num2str(sesInfo.numReductionSteps),'-',num2str(1), '.mat'];
    
    if (~useTemporalICA)
        load(fullfile(outputDir, pcain));
        data = pcasig;
        if size(data, 1) == prod(sesInfo.HInfo.DIM(1:3))
            data = data(mask_ind, :);
        end
    else
        load(fullfile(outputDir, pcain), 'V');
        data = V;
        data = icatb_remove_mean(V);
    end
    
    % Changed this code to allow the users add their own ICA algorithm
    % transpose data to equal components by volume
    data = data';
    
else
    
    fS = whos('-file', fullfile(outputDir, [sesInfo.data_reduction_mat_file, '1-1.mat']));
    fNames = cellstr(char(fS.name));
    chkPcasig = isempty(strmatch('pcasig', fNames, 'exact'));
    
    disp('Stacking data across subjects ...');
    
    if (chkPcasig)
        
        pcain = [sesInfo.data_reduction_mat_file, '1-1.mat'];
        info = load(fullfile(outputDir, pcain));
        VStacked = info.V;
        LambdaStacked = info.Lambda;
        
        for nD = 1:sesInfo.numOfSub*sesInfo.numOfSess
            if (nD == 1)
                startT = 1;
            else
                startT = sum(sesInfo.diffTimePoints(1:nD-1)) + 1;
            end
            endT = sum(sesInfo.diffTimePoints(1:nD));
            [whiteM, dewhiteM] = get_pca_info(VStacked(startT:endT, :), diag(LambdaStacked(nD, :)));
            dat = icatb_read_data(sesInfo.inputFiles(nD).name, [], mask_ind);
            % Call pre-processing function
            dat = icatb_preproc_data(dat, sesInfo.preproc_type, 0);
            % Remove mean per timepoint
            dat = icatb_remove_mean(dat, 0);
            pcasig = dat*whiteM';
            
            if (nD == 1)
                data = zeros(sesInfo.numComp, length(mask_ind), sesInfo.numOfSub*sesInfo.numOfSess);
            end
            
            data(:, :, nD) = pcasig';
            
            clear wM dat pcasig;
            
        end
        
    else
        
        % stack data of all subjects for doing IVA
        for nD = 1:sesInfo.numOfSub*sesInfo.numOfSess
            
            pcain = [sesInfo.data_reduction_mat_file,num2str(sesInfo.numReductionSteps),'-',num2str(nD), '.mat'];
            load(fullfile(outputDir, pcain), 'pcasig');
            
            if (nD == 1)
                data = zeros(sesInfo.numComp, length(mask_ind), sesInfo.numOfSub*sesInfo.numOfSess);
            end
            
            data(:, :, nD) = pcasig';
            
            clear pcasig;
            
            if ~isempty(statusHandle)
                
                % get the status handles
                statusData = getappdata(statusHandle, appDataName);
                statusData.perCompleted = statusData.perCompleted + (statusData.unitPerCompleted / (sesInfo.numOfSub*sesInfo.numOfSess));
                setappdata(statusHandle, appDataName, statusData);
                set(statusHandle, 'name', [num2str(round(statusData.perCompleted*100)), '% analysis done']); waitbar(statusData.perCompleted, statusHandle);
                
            end
        end
        
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