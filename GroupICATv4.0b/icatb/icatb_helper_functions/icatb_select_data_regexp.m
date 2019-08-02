function [allDirs, numOfSub, numOfSess] = icatb_select_data_regexp(inputDir, subDirRegExp, sessDirRegExp)
%% Function to select the data using regular expression pattern match.
%
% Inputs:
% 1. inputDir - Directory where the subjects and sessions are located.
% 2. subDirRegExp - Subject directories regular expression. If subDirRegExp is empty, GIFT will treat this as one subject multiple sessions
% provided sessDirRegExp is not empty. This can have highly nested paths like 'Sub\d+; Study\w+; Visit\w+' where semi-colon is used as a path
% separator.
% 3. sessDirRegExp - Session directories regular expression. If sessDirRegExp is empty, GIFT will treat this as single session analysis.
% This variable cannot have nested paths.
%
%
% Outputs:
% 1. allDirs - All directories obtained after applying regular expressions
% 2. numOfSub - Number of subjects
% 3. numOfSess - Number of sessions
%

if ~exist('inputDir', 'var')
    error('Input directory is missing');
end

if ~exist('subDirRegExp', 'var')
    subDirRegExp = '';
end

if ~exist('sessDirRegExp', 'var')
    sessDirRegExp = '';
end

%% Get full path of the input directory
oldDir = pwd;

cd(inputDir);

inputDir = pwd;

cd(oldDir);

%% Find subject directories
if (isempty(subDirRegExp))
    subDirs = inputDir;
else
    %% List dirs based on subject dir regular expression (This can have
    % highly nested paths)
    subDirs = listSubDirs(inputDir, subDirRegExp);
end

if (isempty(subDirs))
    error('Error:InputData', 'No subject directories found with regular expression (%s)\n', subDirRegExp);
end

%% Number of subjects
numOfSub = size(subDirs, 1);

if ~isempty(sessDirRegExp)

    %% Don't handle nested directories for sessions
    allDirs = cell(size(subDirs, 1), 1);
    % Loop over subject directories
    for nSub = 1:size(subDirs, 1)
        allDirs{nSub} = getDirsRegExp(deblank(subDirs(nSub, :)), sessDirRegExp);
    end
    % End loop over subject directories

    good_inds = icatb_good_cells(allDirs);

    allDirs = allDirs(good_inds);

    if isempty(allDirs)
        error('Error:InputData', 'No session directories found with regular expression (%s)\n', sessDirRegExp);
    end

    % Session numbers
    sessNum = cellfun('size', allDirs, 1);

    if (any(sessNum ~= sessNum(1)))
        error('Error:InputData', 'Session numbers must be the same over subjects.\nPlease check regular expression (%s)\n', sessDirRegExp);
    end

    % Number of sessions
    numOfSess = sessNum(1);
    % All directories
    allDirs = str2mat(allDirs);

else

    % Number of sessions
    numOfSess = 1;
    % All directories
    allDirs = subDirs;

end

function dirs = listSubDirs(inpDir, regExpToMatch)
%% Get all possible directories based on regular expression pattern match
%

regExpToMatch2 = regExpToMatch;

try
    regExpToMatch2 = strtrim(regExpToMatch2);
catch
end

%% Parse regular expression using semi colon as delimiter
regExChars = strread(regExpToMatch2, '%s', 'delimiter', ';');

dirs = inpDir;

count = 1;
while ((count <= length(regExChars)) && ~isempty(dirs))

    % Temp dirs
    tmpDirs = cell(size(dirs, 1), 1);
    for nDir = 1:size(dirs, 1)
        tmpDirs{nDir} = getDirsRegExp(deblank(dirs(nDir, :)), regExChars{count});
    end

    % Get good cells
    good_inds = icatb_good_cells(tmpDirs);

    % Remove bad indices
    dirs = tmpDirs(good_inds);

    clear tmpDirs;

    if ~isempty(dirs)
        dirs = str2mat(dirs);
    else
        dirs = [];
    end

    % Update count
    count = count + 1;

end


function dirs = getDirsRegExp(inpDir, regExpToMatch)
%% Get directories based on the regular expression to match
%

%% List sub folders
dirs = icatb_listSubFolders(deblank(inpDir), 'relative');

if (isempty(dirs))
    return;
end

%% Make a cell array from character array
dirs = cellstr(dirs);

%% Use regexp for OS other than Windows
if ~ispc
    start_ind = regexp(dirs, regExpToMatch);
else
    start_ind = regexpi(dirs, regExpToMatch);
end

good_inds = icatb_good_cells(start_ind);

%% Return matched dirs with full path
dirs = dirs(good_inds);

if ~isempty(dirs)
    dirs = str2mat(dirs);
    if (~strcmp(inpDir(end), filesep))
        inpDir = [inpDir, filesep];
    end
    %dirs = strcat(inpDir, dirs);
    dirs = [repmat(inpDir, size(dirs, 1), 1), dirs];
else
    dirs = [];
end