function [tcInfo, icInfo] = icatb_groupBackReconInfo(sesInfo, W)

j = sesInfo.numReductionSteps;

% Use size of Weights matrix instead of numComp to handle Constrained
% ICA (Spatial)
icInfo = {eye(size(W, 2), size(W, 2))};
tcInfo = icInfo;

%% Load higher level data reduction files
while (j > 1)
    
    numGroupsAfterCAT = sesInfo.reduction(j).numOfGroupsAfterCAT;
    numOfPCBeforeCAT = sesInfo.reduction(j).numOfPCBeforeCAT;
    numOfPrevGroupsInEachNewGroupAfterCAT = sesInfo.reduction(j).numOfPrevGroupsInEachNewGroupAfterCAT;
    
    tmpICInfo = cell(1, sum(numOfPrevGroupsInEachNewGroupAfterCAT));
    tmpTCInfo = tmpICInfo;
    
    countN = 0;
    %% Loop over no. of groups after concatenation
    for groupIndex = 1:numGroupsAfterCAT
        
        %% Load jth pca step file
        pcain = [sesInfo.data_reduction_mat_file, num2str(j), '-', num2str(groupIndex), '.mat'];
        fileIn = fullfile(sesInfo.outputDir, pcain);
        
        if (~exist(fileIn, 'file'))
            error(['File ', fileIn, ' doesn''t exist. Please set CONSERVE_DISK_SPACE in icatb_defaults.m to 0 or 1 prior to running setup analysis tool.']);
        end
        
        [dewhiteM, whiteM] = load_vars(fileIn);
        
        endT = 0;
        %% Loop over no. of sets in the concatenation
        for nP = 1:numOfPrevGroupsInEachNewGroupAfterCAT(groupIndex)
            countN = countN + 1;
            startT = endT + 1;
            endT = endT + numOfPCBeforeCAT;
            tmpICInfo{countN} = icInfo{groupIndex}*whiteM(:, startT:endT);
            if (strcmpi(sesInfo.backReconType, 'regular') || strcmpi(sesInfo.backReconType, 'gica'))
                %% Regular
                tmpTCInfo{countN} = dewhiteM(startT:endT, :)*tcInfo{groupIndex};
            else
                %% GICA3
                tmpTCInfo{countN} = pinv(whiteM(:, startT:endT))*tcInfo{groupIndex};
            end
            
            if strcmpi(sesInfo.backReconType, 'gica')
                tmpICInfo{countN} = icInfo{groupIndex}*pinv(dewhiteM(startT:endT, :));
            elseif (strcmpi(sesInfo.backReconType, 'gica3') || strcmpi(sesInfo.backReconType, 'regular'))
                tmpICInfo{countN} = icInfo{groupIndex}*whiteM(:, startT:endT);
            end
            
        end
        %% End of loop over no. of sets in the concatenation
        
    end
    %% End of loop over no. of groups after concatenation
    
    icInfo = tmpICInfo;
    clear tmpICInfo;
    tcInfo = tmpTCInfo;
    clear tmpICInfo;
    
    j = j - 1;
    
end
%% End for loading higher level data reduction files

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

