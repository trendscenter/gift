function parameters = icatb_parseRegressParamTextFile(fileName)
% Parse text file and return the results in a parameters structure


%%%%% Read the contents of text file %%%%%%%%%%%%%
icatb_defaults;
global PARAMETER_INFO_MAT_FILE;


readContents = textread(fileName, '%s', 'delimiter', '\n');
ind = good_cells(readContents);
readContents = readContents(ind);
clear ind;

%%%% End for reading the contents of text file %%%%%

% check the associated parameter file
[pathstr, fileN, extn] = fileparts(fileName);
regressFileSuffix = '_temporal_regression';


try
    checkParameterFile = icatb_findstr(fileN, regressFileSuffix);
    if ~isempty(checkParameterFile)
        outputPrefix = strrep(fileN, regressFileSuffix, '');
        parameterFile = fullfile(outputPrefix, [PARAMETER_INFO_MAT_FILE, '.mat']);
        if exist(parameterFile, 'file')
            load(parameterFile);
            spmMatFlag = sesInfo.userInput.spmMatFlag;
        end
    end
catch
end

% All viewing sets
allViewingSets = {'Subject', 'Session', 'All data-sets'};

titleStr = ['TEMPORAL SORTING OF COMPONENTS OF '];

inds = regexpi(readContents, titleStr);

inds = find(good_cells(inds));

if isempty(inds)
    error('Error:RegressionParameters', 'File %s \nis not a valid regression parameters text file', fileName);
end

% Truncate the text file
readContents = readContents(inds(end):end);
inds = good_cells(readContents);
readContents = readContents(inds);

viewingSet = 'single subject';
for nR = 1:length(allViewingSets)
    checkViewingSet = findstr(lower(readContents{1}), lower([titleStr, allViewingSets{nR}]));
    if ~isempty(checkViewingSet)
        viewingSet = allViewingSets{nR};
    end
end

if strcmpi(viewingSet, 'single subject')
    error('Statistical testing of time courses cannot be done for viewing set single subject single session');
end

inds = regexp(readContents, 'Component Number');
inds = find(good_cells(inds));

if isempty(inds)
    error('Check the regression text file as Component Number variable doesn''t exist');
end

if length(readContents) < 5
    error('Need atleast two data-sets to do statistical testing on time courses');
end

componentStr = readContents(inds(end));

componentStr = strread(componentStr{1}, '%s', 'delimiter', '\t');
inds = regexp(componentStr, '\d');
inds = find(good_cells(inds));

componentStr = componentStr(inds);
compNum = str2num(str2mat(componentStr))';

origComp = compNum;
[compNum, sortCompInd] = sort(compNum);

% Regression parameter string
regressionParamStr = readContents(4:end);

% Initialise regression info structure
regressionInfo = repmat(struct('name', '', 'values', []), length(regressionParamStr), 1);

% Loop over number of regression parameters
for nRegression = 1:length(regressionParamStr)
    strParts = strread(regressionParamStr{nRegression}, '%s', 'delimiter', '\t');
    regressionInfo(nRegression).name = deblank(strParts{1});
    regressionInfo(nRegression).values = str2num(str2mat(strParts(2:end)));
    %regressionInfo(nRegression).values = zeros(1, length(strParts) - 1);
    %     % Loop over string parts
    %     for nn = 1:length(strParts)
    %         if nn == 1
    %             regressionInfo(nRegression).name = deblank(strParts{1});
    %         else
    %             regressionInfo(nRegression).values(nn-1) = str2num(strParts{nn});
    %         end
    %     end
    %     % End loop over string parts

    regressionInfo(nRegression).values = regressionInfo(nRegression).values(sortCompInd);

end
% End loop over number of regression parameters


% Regressor names
regressorNames = cellstr(str2mat(regressionInfo.name));

if strcmpi(viewingSet, 'session')
    %regExpPattern = 'Subject 1 ';
    regExpPattern = 'Subject (\d+) ';
    [startInd, endInd, tokens] = regexp(regressorNames, regExpPattern);
    subNum = zeros(1, length(tokens));
    for nn = 1:length(tokens)
        temp = regressorNames{nn};
        myToken = tokens{nn}{1};
        subNum(nn) = str2num(temp(myToken(1):myToken(2)));
    end

    numOfSub = max(subNum);
    numOfSess = 1;
elseif strcmpi(viewingSet, 'subject')
    %regExpPattern = 'Session 1 ';
    regExpPattern = 'Session (\d+) ';
    [startInd, endInd, tokens] = regexp(regressorNames, regExpPattern);
    sessNum = zeros(1, length(tokens));
    for nn = 1:length(tokens)
        temp = regressorNames{nn};
        % Session number
        myToken = tokens{nn}{1};
        sessNum(nn) = str2num(temp(myToken(1):myToken(2)));
    end

    numOfSub = 1;
    numOfsess = max(sessNum);

else
    %     regExpPattern = 'Subject 1 Session 1 ';
    %     inds = regexp(regressorNames, 'Subject 1 Session 1');
    %     inds = find(good_cells(inds));
    %     if isempty(inds)
    %         regExpPattern = 'Subject 1 ';
    %     end

    regExpPattern = 'Subject (\d+) Session (\d+) ';
    [startInd, endInd, tokens] = regexp(regressorNames, regExpPattern);
    ind = find(good_cells(startInd));
    if isempty(ind)
        regExpPattern = 'Subject (\d+) ';
        [startInd, endInd, tokens] = regexp(regressorNames, regExpPattern);
        ind2 = find(good_cells(startInd));
        if ~isempty(ind2)
            subNum = zeros(1, length(tokens));
            sessNum = 1;
            for nn = 1:length(tokens)
                temp = regressorNames{nn};
                % Subject number
                myToken = tokens{nn}{1};

                subNum(nn) = str2num(temp(myToken(1, 1):myToken(1, 2)));
            end
        else
            subNum = 1;
            sessNum = 1;
        end
    else
        tokens = tokens(ind);
        subNum = zeros(1, length(tokens));
        sessNum = zeros(1, length(tokens));
        for nn = 1:length(tokens)
            temp = regressorNames{nn};
            % Subject number
            myToken = tokens{nn}{1};

            subNum(nn) = str2num(temp(myToken(1, 1):myToken(1, 2)));
            % Session number
            sessNum(nn) = str2num(temp(myToken(2, 1):myToken(2, 2)));
        end
    end

    % Number of subjects and sessions
    numOfSub = max(subNum);
    numOfSess = max(sessNum);

end


if numOfSub*numOfSess == 1
    error('Stats cannot be done on single subject single session analysis');
end

if (numOfSub*numOfSess > 1) && (~strcmpi(viewingSet, 'all data-sets'))
    error('Error:StatTestTc', ['Statistical testing of time courses (beta weights) is done only when ', ...
        '\nall data-sets are used for temporal sorting']);
end

% Use regular expr
conditionStrings = regexprep(regressorNames, regExpPattern, '');

% Number of conditions
numConditions = ceil(length(regressorNames) / (numOfSub*numOfSess));

% Condition string structure
condStringStruct = repmat(struct('cond', repmat(struct('name', ''), numConditions, 1)), numOfSub*numOfSess, 1);

% Loop over number of conditions
countC = 0;
for nD = 1:length(condStringStruct)
    for nC = 1:numConditions
        countC = countC + 1;
        condStringStruct(nD).cond(nC).name = conditionStrings{countC};
    end
end
% End loop over number of conditions


% % Initialise regressInfo variable
regressInfo = repmat(struct('cond', repmat(struct('name', '', 'values', []), numConditions, 1)), ...
    numOfSub*numOfSess, 1);

% Initialise count for data-set and condition number
countDataSet = 0;
condNum = 0;

% Loop over number of subjects
for nSub = 1:numOfSub
    % Loop over number of sessions
    for nSess = 1:numOfSess
        countDataSet = countDataSet + 1;
        % Loop over conditions
        for nCond = 1:numConditions
            condNum = condNum + 1;
            % Regress Info structure
            regressInfo(countDataSet).cond(nCond).name = regressionInfo(condNum).name;
            regressInfo(countDataSet).cond(nCond).values = regressionInfo(condNum).values;
            regressInfo(countDataSet).regressorName(nCond).name = condStringStruct(countDataSet).cond(nCond).name;
        end
        % End loop over number of conditions
    end
    % End loop over number of sessions
end
% End loop over number of subjects

clear regressionInfo;

% Return parameters structure
parameters.viewingSet = viewingSet;
parameters.allViewingSets = allViewingSets;
parameters.numOfSub = numOfSub;
parameters.numOfSess = numOfSess;
parameters.numConditions = numConditions;
parameters.componentNum = compNum;
parameters.regressInfo = regressInfo;

if exist('spmMatFlag', 'var')
    parameters.spmMatFlag = spmMatFlag;
end


function ind = good_cells(mycell)
% Find good cells

if ~iscell(mycell)
    mycell = {mycell};
end

for j = 1:length(mycell)
    ind(j) = ~isempty(mycell{j});
end


function [matchedStr, inds, startInd, endInd] = findStringMatch(myString, regExpPattern)
% Return matched string

[startInd, endInd] = regexp(myString, regExpPattern);

inds = find(good_cells(startInd));

% Matched string
matchedStr = myString(inds);

% Start and End indices
startInd = startInd(inds);
endInd = endInd(inds);