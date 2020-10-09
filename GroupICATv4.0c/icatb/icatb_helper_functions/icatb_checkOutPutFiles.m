function [outputFiles, status] = icatb_checkOutPutFiles(outputFiles, outputDir, prefix)
% Replace subject ICA output files with calibrated output files
% Status variable contains whether MAT file can be used or not

icatb_defaults;

global CALIBRATE_MAT_FILE;
global SUBJECT_ICA_AN3_FILE;
global COMPONENT_NAMING;
global SESSION_POSTFIX;

% Initialise status
status = 0;

regexpPattern = [prefix, SUBJECT_ICA_AN3_FILE, '(\d+)', COMPONENT_NAMING, SESSION_POSTFIX, '(\d+)'];

calibrateFileName = [prefix, CALIBRATE_MAT_FILE];

if ~iscell(outputFiles)
    outputFiles = cellstr(outputFiles);
end

[startInd, endInd, tokenExtents] = regexp(outputFiles{1}, regexpPattern);

inds = icatb_good_cells(startInd);

if ~isempty(find(inds ~= 0))
    currentStr = outputFiles{1};
    temp = tokenExtents{1};
    subNum = str2num(currentStr(temp(1, 1):temp(1, 2)));
    sessNum = str2num(currentStr(temp(2, 1):temp(2, 2)));
    currentStr = [calibrateFileName, num2str(subNum), '-', num2str(sessNum), '.mat'];    
    if exist(fullfile(outputDir, currentStr), 'file')
        outputFiles = outputFiles(1);
        status = 1;
        outputFiles{1} = currentStr;
    end
end


outputFiles = str2mat(outputFiles);