function regressParam = icatb_parse_regression_mat(matFile, numSub, numSess, numCond, avgRuns, writeExcel)
%% Parse regression MAT file
%
% Inputs:
% 1. matFile - Regression parameters mat file name.
% 2. numSub - Number of subjects.
% 3. numSess - Number of sessions.
% 4. numCond - Number of conditions.
% 5. avgRuns - Average runs. Options are 0 and 1. 0 means don't average
% runs whereas 1 means average over runs.
% 6. writeExcel - Write to Excel file. Options are 0 and 1. A value of 1
% means write regression parameters to excel file. If you don't average
% across runs, the number of excel files written will be equal to the no.
% of runs. Note that xlswrite function is required to write to excel files.
% This option may not work for versions less than MATLAB 7.
%
% Outputs:
% regressParam - Regression parameters matrix of dimensions components by
% subjects by runs by conditions. If you averaged across runs, the
% dimension of regressParam will be components by subjects by conditions.
%

pathstr = [];

if ~isstruct(matFile)
    
    [pathstr, fName, extn] = fileparts(matFile);       
    
    load (matFile);      
    
else
    
    regressInfo = matFile;
    clear matFile;
    
end

if (~isfield(regressInfo, 'regressionParameters'))
    error('You need to select a valid regression parameters file');
end

if isempty(pathstr)
    pathstr = pwd;
end

if (~exist('avgRuns', 'var'))
    avgRuns = 0;
end

if (~exist('writeExcel', 'var'))
    writeExcel = 0;
end


regressParameters = regressInfo.regressionParameters;
compNum = regressInfo.componentNumbers;
[compNum, sortInd] = sort(compNum);
regressParameters = regressParameters(:, sortInd);
numComp = length(compNum);

countS = 0;
% Initialise regressParam
regressParam = zeros(numComp, numSub, numSess, numCond);
% Loop over subjects
for nSub = 1:numSub
    % Loop over sessions
    for nSess = 1:numSess
        % Loop over conditions
        for nCond = 1:numCond
            countS = countS + 1;
            regressParam(:, nSub, nSess, nCond) = regressParameters(countS, :);
        end
        % End loop over conditions
    end
    % End loop over sessions
end
% End loop over subjects

clear regressParameters;

if avgRuns
    % Mean of regression parameters
    regressParam = squeeze(mean(regressParam, 3));
    numSess = 1;
    excelFileName = fullfile(pathstr, [fName, '.xls']);
end

if (writeExcel)
    
    warning off MATLAB:xlswrite:AddSheet;
    
    for nSess = 1:numSess
        if ~avgRuns
            excelFileName = fullfile(pathstr, [fName, '_sess_', num2str(nSess), '.xls']);
            data = squeeze(regressParam(:, :, nSess, :));
        else
            data = squeeze(regressParam(:, :, :));
        end
        disp(['Writing excel file ', excelFileName, '...']);
        for nComp = 1:size(data, 1)
            currentData = squeeze(data(nComp, :, :));
            [stat, mmsg] = xlswrite(excelFileName, currentData, ['Sheet', num2str(nComp)]);
        end
        % End for writing
    end
end
disp('Done');
