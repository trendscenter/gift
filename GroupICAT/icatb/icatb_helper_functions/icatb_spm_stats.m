function icatb_spm_stats(parameter_file, desVal, group1Name, group1, group2Name, group2, compNumbers, writeTal, avgRuns)
% SPM stats on subject component maps. Currently wrapper is provided only
% for one sample t-test and two sample t-test.
%

% Load defaults.
icatb_defaults;

global SPM_STATS_TTEST2_EXPLICIT_MASK;
global SPM_STATS_TTEST_THRESHOLD;

spmPath = which('spm.m');

if isempty(spmPath);
    error('SPM does not exist on MATLAB path');
end

verNum = str2num(strrep(lower(spm('ver')), 'spm', ''));

if (verNum < 5)
    error('This utility works with SPM5 and higher');
end

% Defaults
icatb_defaults;
global PARAMETER_INFO_MAT_FILE;

filterP = ['*', PARAMETER_INFO_MAT_FILE, '*.mat'];
if ~exist('parameter_file', 'var') || isempty(parameter_file)
    parameter_file = icatb_selectEntry('typeEntity', 'file', 'typeSelection', 'single', 'filter', ...
        filterP, 'title', 'Select a valid parameter file');
end

if isempty(parameter_file)
    error('Parameter file is not selected');
end

load(parameter_file);

if ~exist('sesInfo', 'var')
    error(['Selected file: ', parameter_file, ' is not a valid parameter file']);
end

% Results directory
outputDir = fileparts(parameter_file);

if isempty(outputDir)
    outputDir = pwd;
end

cd(outputDir);

if ~isfield(sesInfo, 'icaOutputFiles')
    error('You need to run ICA analysis in order to do SPM stats on subject component maps');
end

% No. of subjects, sessions and components
numOfSub = sesInfo.numOfSub;
numOfSess = sesInfo.numOfSess;
numComp = sesInfo.numComp;
prefix = sesInfo.userInput.prefix;

if (numOfSub*numOfSess == 1)
    error('Cannot do stats as there is only data set present');
end

% Subject ICA files
subjectICAFiles = icatb_parseOutputFiles('icaOutputFiles', sesInfo.icaOutputFiles, 'numOfSub', numOfSub, 'numOfSess', numOfSess, ...
    'flagTimePoints', sesInfo.flagTimePoints);

zipContents = sesInfo.zipContents;

choiceString = char('One sample t-test', 'Two sample t-test', 'Paired T-test');

if ~exist('desVal', 'var')
    InputHandle = icatb_getGraphics('Select Design', 'normal', 'Select Design'); % figure handle
    % set menubar none
    set(InputHandle, 'menubar', 'none');
    value = icatb_promptUI('popup', 'Select Design', choiceString, 'numeric', InputHandle);
    if ((numOfSub > 1) && (numOfSess > 1) && (value ~= 3))
        avgRuns = icatb_promptUI('popup', 'Do You Want To Use Average Of Runs?', char('Yes', 'No'), 'numeric', InputHandle);
        if (avgRuns == 2)
            avgRuns = 0;
        end
    end
    close(InputHandle);
else
    if size(choiceString, 1) < desVal
        desVal = 1;
    end
    value = desVal;
end


desType = deblank(choiceString(value, :));

listStr = num2str((1:numComp)');

if ~exist('compNumbers', 'var')
    title_fig = 'Select components';
    compNumbers = icatb_listdlg('PromptString', title_fig, 'SelectionMode', 'multiple', 'ListString', listStr, ...
        'movegui', 'center', 'windowStyle', 'modal', 'title_fig', title_fig);
end

if isempty(compNumbers)
    error('Components are not selected');
end

if (size(listStr, 1) < max(compNumbers))
    error('Error:ComponentNumbers', 'Maximum of component numbers (%d) exceeds the maximum number of components (%d)\n',  max(compNumbers), size(listStr, 1));
end

clear listStr;

if (strcmpi(desType, 'two sample t-test') || strcmpi(desType, 'paired t-test'))
    % Check if one sample t-test results are done or not
    if (SPM_STATS_TTEST2_EXPLICIT_MASK)
        % Loop over components
        for nComp = compNumbers
            [fileIndex] = icatb_returnFileIndex(nComp);
            comp_ttest_dir = [prefix, '_one_sample_ttest_results',  filesep, 'Comp_', fileIndex];
            cc = load(fullfile(outputDir, comp_ttest_dir, 'SPM.mat'));
            ttest_map = fullfile(outputDir, comp_ttest_dir, cc.SPM.xCon(1).Vspm.fname);
            %ttest_map = fullfile(outputDir, comp_ttest_dir, 'spmT_0001.img');
            if ~exist(ttest_map, 'file')
                error('Error:CheckTmap', ['Please run one sample t-test for component %s if you want to use explicit mask.\nSee icatb_defaults.m ', ...
                    'for more information.'], fileIndex);
            end
        end
        % End loop over components
    end
    % End for checking
end

if (~exist('avgRuns', 'var') || (numOfSess == 1) || strcmpi(desType, 'paired t-test'))
    avgRuns = 0;
end

if (avgRuns == 0)
    
    % Select subjects and sessions for groups
    count = 0;
    listStr = repmat({''}, numOfSub*numOfSess, 1);
    for nSub = 1:numOfSub
        for nSess = 1:numOfSess
            count = count + 1;
            listStr{count} = ['Subject ', num2str(nSub), ' Session ', num2str(nSess)];
        end
    end
    
else
    
    % Select subjects for groups
    listStr = repmat({''}, numOfSub, 1);
    for nSub = 1:numOfSub
        listStr{nSub} = ['Subject ', num2str(nSub)];
    end
    
end

if strcmpi(desType, 'one sample t-test')
    if ~exist('group1Name', 'var')
        [groupName, group] = selectGroups(listStr, 'group');
    else
        groupName = group1Name;
        group = group1;
    end
    if isempty(group)
        error('Data sets are not selected for group');
    end
else
    
    g1Str = 'group 1';
    g2Str = 'group 2';
    
    if (strcmpi(desType, 'paired t-test'))
        g1Str = 'Condition 1';
    end
    
    if (strcmpi(desType, 'paired t-test'))
        g2Str = 'Condition 2';
    end
    
    if ~exist('group1Name', 'var')
        [group1Name, group1] = selectGroups(listStr, g1Str);
    end
    
    if isempty(group1)
        error('Data sets are not selected for group 1');
    end
    
    if ~exist('group2Name', 'var')
        [group2Name, group2] = selectGroups(listStr, g2Str);
    end
    
    if isempty(group2)
        error('Data sets are not selected for group 2');
    end
    
    group = [group1, group2];
    
end

if (length(group) == 1)
    error('Cannot do stats as you have selected only one data set');
end

if (size(listStr, 1) < max(group))
    error('Error:Groups', 'Maximum of groups (%d) exceeds the maximum of data sets (%d)\n',  max(group), size(listStr, 1));
end

drawnow;

if (~exist('writeTal', 'var') || isempty(writeTal))
    writeTal = 0;
end

warning off all;

filesToDelete = {};

if (avgRuns == 0)
    
    dataSetNum = group;
    disp('Checking subject component maps ...');
    % If component files are not present unzip files
    for nDataSet = dataSetNum
        nSub = ceil(nDataSet/numOfSess);
        nSess = nDataSet - (nSub - 1)*numOfSess;
        currentFile = deblank(subjectICAFiles(nSub).ses(nSess).name(1, :));
        if ~exist(currentFile, 'file')
            % Check if zip file is present or not
            [zipFileName, files_in_zip] = icatb_getViewingSet_zip(currentFile, [], 'real', zipContents);
            icatb_unzip(regexprep(zipFileName, ['.*\', filesep], ''), fullfile(outputDir, fileparts(currentFile)));
            %icatb_unzip(zipFileName, outputDir);
            filesToDelete{length(filesToDelete) + 1} = files_in_zip;
            clear files_in_zip;
        end
    end
    % End for checking if zip files are present or not
    disp('Done checking subject component maps');
    
    disp('Reading headers ...');
    for nDataSet = dataSetNum
        nSub = ceil(nDataSet/numOfSess);
        nSess = nDataSet - (nSub - 1)*numOfSess;
        compFiles = icatb_rename_4d_file(subjectICAFiles(nSub).ses(nSess).name);
        subjectICAFiles(nSub).ses(nSess).name = compFiles;
        clear compFiles;
    end
    disp('Done reading headers');
    
else
    
    if (~isfield(sesInfo, 'spm_stats') || ~exist(fullfile(outputDir, sesInfo.spm_stats.avg_runs_dir), 'dir'))
        icatb_spm_avg_runs(parameter_file);
        load(parameter_file);
    end
    
    % Get the stats files and directory from spm_stats field
    subjectICAFiles = sesInfo.spm_stats.avg_runs_files;
    avg_runs_dir = sesInfo.spm_stats.avg_runs_dir;
end

clear sesInfo;

try
    
    disp('Running SPM stats ...');
    
    if strcmpi(desType, 'one sample t-test')
        analysisDir = [prefix, '_one_sample_ttest_results'];
    elseif strcmpi(desType, 'two sample t-test')
        analysisDir = [prefix, '_two_sample_ttest_results'];
    else
        analysisDir = [prefix, '_paired_ttest_results'];
    end
    
    % Loop over components
    for nComp = compNumbers
        [fileIndex] = icatb_returnFileIndex(nComp);
        disp(['Doing ', lower(desType), ' on spatial maps for component ', fileIndex]);
        
        compDir = [analysisDir, filesep, 'Comp_', fileIndex];
        if (exist(fullfile(outputDir, compDir), 'dir') ~= 7)
            mkdir(outputDir, compDir);
        end
        
        temp_dir = fullfile(outputDir, compDir);
        
        if (strcmpi(desType, 'one sample t-test'))
            % One sample t-test
            
            % Group 1 files
            if ~avgRuns
                group1Files = get_files(subjectICAFiles, group, numOfSess, nComp, outputDir);
            else
                group1Files = icatb_fullFile('directory', fullfile(outputDir, avg_runs_dir), 'files', char(subjectICAFiles{nComp, group}));
            end
            
            % Do one sample t-test
            icatb_ttest_maps(group1Files, temp_dir, groupName);
            
            clear group1Files;
            
        else
            % Two sample t-test or paired t-test
            
            
            comp_ttest_dir = [prefix, '_one_sample_ttest_results',  filesep, 'Comp_', fileIndex];
            cc = load(fullfile(outputDir, comp_ttest_dir, 'SPM.mat'));
            ttest_map = fullfile(outputDir, comp_ttest_dir, cc.SPM.xCon(1).Vspm.fname);
            %ttest_map = fullfile(outputDir, comp_ttest_dir, 'spmT_0001.img');
            explicitMask = '';
            
            if (SPM_STATS_TTEST2_EXPLICIT_MASK)
                % Apply threshold and write images
                VT = spm_vol(ttest_map);
                dataT = spm_read_vols(VT);
                dataT(isnan(dataT)) = 0;
                ind = (abs(dataT(:)) >= SPM_STATS_TTEST_THRESHOLD);
                if isempty(find(ind ~= 0))
                    error('Error:Voxels', 'No voxels found in file %s\n for SPM_STATS_TTEST_THRESHOLD = %s.', ...
                        VT.fname, num2str(SPM_STATS_TTEST_THRESHOLD));
                else
                    dataT((ind == 0)) = 0;
                    dataT((ind == 1)) = 1;
                    VT.fname = fullfile(outputDir, comp_ttest_dir, 'spmT_0001_mask.img');
                    spm_write_vol(VT, dataT);
                    clear dataT;
                    explicitMask = VT.fname;
                end
            end
            
            if ~avgRuns
                % Group 1 files
                group1Files = get_files(subjectICAFiles, group1, numOfSess, nComp, outputDir);
                % Group 2 files
                group2Files = get_files(subjectICAFiles, group2, numOfSess, nComp, outputDir);
            else
                % Group 1 files
                group1Files = icatb_fullFile('directory', fullfile(outputDir, avg_runs_dir), 'files', char(subjectICAFiles{nComp, group1}));
                % Group 2 files
                group2Files = icatb_fullFile('directory', fullfile(outputDir, avg_runs_dir), 'files', char(subjectICAFiles{nComp, group2}));
            end
            
            if (strcmpi(desType, 'two sample t-test'))
                % Do two sample t-test
                icatb_ttest2_maps(group1Files, group2Files, temp_dir, group1Name, group2Name, explicitMask);
            else
                % Paired t-test
                icatb_paired_t_test_maps(group1Files, group2Files, explicitMask, temp_dir, [group1Name, ' - ', group2Name]);
            end
            
            clear group1Files group2Files;
            
        end
        
        if (writeTal >= 2)
            % Write talariach coords
            cc = load(fullfile(temp_dir, 'SPM.mat'));
            tmapFile = fullfile(temp_dir, cc.SPM.xCon(1).Vspm.fname);
            icatb_talairach(tmapFile);
        end
        
        fprintf('\n');
        
    end
    % End loop over components
    
    cd(outputDir);
    
    if ~isempty(filesToDelete)
        for nF = 1:length(filesToDelete)
            files_in_zip = filesToDelete{nF};
            files_in_zip = str2mat(files_in_zip);
            icatb_delete_file_pattern(files_in_zip, outputDir);
        end
    end
    
    disp(['Results are stored in ', fullfile(outputDir, analysisDir)]);
    
    fprintf('\n');
    warning on;
    
catch
    warning on;
    icatb_displayErrorMsg;
end

%%%%%% Sub functions %%%%%

function group1Files = get_files(subjectICAFiles, group1, numOfSess, compNum, outputDir)
% get files based on group

group1Files = repmat({''}, length(group1), 1);
count = 0;
for nDataSet = group1
    nSub = ceil(nDataSet/numOfSess);
    nSess = nDataSet - (nSub - 1)*numOfSess;
    count = count + 1;
    group1Files{count} = icatb_fullFile('directory', outputDir, 'files', deblank(subjectICAFiles(nSub).ses(nSess).name(compNum, :)));
end


function [groupName, groupVal] = selectGroups(subjectString, titleString)
% Select datasets for group

groupName = '';
groupVal = [];

icatb_defaults;
global UI_FS;

figTag = 'select_subjects_spm_stats';

% Delete figures that have tag select_subjects_spm_stats
if ~isempty(findobj(0, 'tag', figTag))
    delete(findobj(0, 'tag', figTag));
end

figTitle = ['Select data sets for ', titleString];

% Plot Figure
[graphicsHandle] = icatb_getGraphics(figTitle, 'normal', figTag, 'off');
set(graphicsHandle, 'menubar', 'none');

if ispc
    set(graphicsHandle, 'windowstyle', 'modal');
end

% Offsets
xOffset = 0.05; yOffset = 0.05;
editTextHeight = 0.05; editTextWidth = 0.45; yPos = 0.95;

editTextPos = [xOffset, yPos - yOffset - editTextHeight, editTextWidth, editTextHeight];

% Plot Text
editTextH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', ...
    'position', editTextPos, 'String', ['Enter name for ', titleString], 'fontsize', UI_FS - 1, ...
    'horizontalalignment', 'center');

[editTextH] = icatb_wrapStaticText(editTextH);

editTextPos = get(editTextH, 'position');

%%% Edit Box
editWidth = 0.4;
editPos = editTextPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = editWidth;

% Plot edit control
editH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'edit', ...
    'position', editPos, 'String', titleString, 'fontsize', UI_FS - 1, 'tag', 'group_name');

% List text width and height
listTextWidth = 0.6; listTextHeight = 0.05;
listTextPos = [0.5 - 0.5*listTextWidth, editTextPos(2) - 2*yOffset - 0.5*listTextHeight, listTextWidth, ...
    listTextHeight];

% Plot listbox
listTextH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', ...
    'position', listTextPos, 'String', 'Select data sets', 'fontsize', UI_FS - 1);


[listTextH] = icatb_wrapStaticText(listTextH);

listTextPos = get(listTextH, 'position');

% Listbox position
listWidth = listTextWidth; listHeight = 0.4;
listPos = listTextPos;
listPos(2) = listPos(2) - yOffset - listHeight;
listPos(3) = listWidth;
listPos(4) = listHeight;

% Plot listbox
listH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'listbox', ...
    'position', listPos, 'String', subjectString, 'fontsize', UI_FS - 1, 'tag', 'subject_listbox', 'min', 0, ...
    'max', 2, 'value', []);

editWidth = 0.3;
editHeight = 0.05;
editPos = [listPos(1) + 0.5*listPos(3) - 0.5*editWidth, listPos(2) - yOffset - 0.5*editHeight, ...
    editWidth, editHeight];

% Plot edit control
editH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'edit', ...
    'position', editPos, 'String', num2str(get(listH, 'value')), 'fontsize', UI_FS - 1, 'tag', ...
    'subject_editbox');


set(listH, 'callback', {@listSubjectsCallback, editH});
set(editH, 'callback', {@editSubjectsCallback, listH});

% Ok button
okWidth = 0.2; okHeight = 0.05;
okPos = [0.5 - 0.5*okWidth, yOffset + 0.5*okHeight, okWidth, okHeight];
okH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', okPos, 'String', 'OK', 'fontsize', UI_FS - 1, 'tag', 'ok_select_subjects', 'callback', ...
    {@okSubjectsCallback, graphicsHandle,listH,  titleString});

try
    % Set graphics handle visibility on
    set(graphicsHandle, 'visible', 'on');
    waitfor(graphicsHandle);
catch
end

if isappdata(0, 'groupValData')
    groupData = getappdata(0, 'groupValData');
    groupName = groupData.groupName;
    groupVal = groupData.groupVal;
    rmappdata(0, 'groupValData');
end


function okSubjectsCallback(hObject, event_data, handles, listH, titleString)
% Ok button for selecting data sets

try
    
    groupNameH = findobj(handles, 'tag', 'group_name');
    
    groupName = get(groupNameH, 'string');
    
    groupVal = get(listH, 'value');
    
    if isempty(groupVal)
        error(['Data sets are not selected for ', titleString]);
    end
    
    groupData.groupName = groupName;
    groupData.groupVal = groupVal;
    
    setappdata(0, 'groupValData', groupData);
    
    delete(handles);
    
catch
    disp(lasterr);
    icatb_errorDialog(lasterr, 'Error Found');
end

function editSubjectsCallback(hObject, event_data, subjectListH)
% Set the listbox value

editString = deblank(get(hObject, 'string'));

editVal = str2num(editString);

numItems = size(get(subjectListH, 'string'), 1);

%%%% Do Error checking %%%%
if isempty(editVal)
    error('Error:EditBox', 'Check the edit box string (%s) \nas it must generate a valid numeric value', editString);
end

if ~isempty(find(editVal == 0))
    error('Error:EditBox', 'Edit box string (%s) cannot accept a value of 0', editString);
end

if numItems < length(editVal)
    error('Error:EditBox', 'Number of items in edit box (%d) is larger than the \nnumber of items (%d) in listbox ', ...
        length(editVal), numItems);
end

if numItems < max(editVal)
    error('Error:EditBox', 'Maximum value in editbox string (%s) is larger than the \nnumber of items (%d) in listbox ', ...
        editString, numItems);
end
%%%% End for doing error checking %%%%

checkInteger = [];
try
    checkInteger = strread(num2str(editVal), '%d');
catch
end
if isempty(checkInteger)
    error('Error:EditBox', 'Editbox string (%s) does not contain integer items', editString);
end

set(subjectListH, 'value', editVal);


function listSubjectsCallback(hObject, event_data, editH)
% Set the listbox value

editString = num2str(get(hObject, 'value'));

set(editH, 'string', editString);