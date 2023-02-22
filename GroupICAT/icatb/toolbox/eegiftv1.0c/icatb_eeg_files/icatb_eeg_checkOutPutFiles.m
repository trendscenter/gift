function outputFiles = icatb_eeg_checkOutPutFiles(outputFiles, prefix)
% Replace subject ICA output files with calibrated output files

icatb_defaults;

global CALIBRATE_MAT_FILE;
global SUBJECT_ICA_AN3_FILE;
global COMPONENT_NAMING;
global SESSION_POSTFIX;

subNaming = regexprep(SUBJECT_ICA_AN3_FILE, '^(_)', '');

regexpPattern = [prefix, '.*', subNaming, '(\d+)', COMPONENT_NAMING, SESSION_POSTFIX, '(\d+)'];

calibrateFileName = [prefix, CALIBRATE_MAT_FILE];

if ~iscell(outputFiles)
    outputFiles = cellstr(outputFiles);
end

outputFiles = outputFiles(1);

[startInd, endInd, tokenExtents] = regexp(outputFiles, regexpPattern);

inds = icatb_good_cells(startInd);

if ~isempty(find(inds ~= 0))
    currentStr = outputFiles{1};
    % Fixed in version GroupICATv2.0b
    temp = tokenExtents{1};
    if iscell(temp)
        temp = temp{1};
    end
    subNum = str2num(currentStr(temp(1, 1):temp(1, 2)));
    sessNum = str2num(currentStr(temp(2, 1):temp(2, 2)));
    currentStr = [calibrateFileName, num2str(subNum), '-', num2str(sessNum), '.mat'];
    outputFiles{1} = currentStr;
end

outputFiles = str2mat(outputFiles);