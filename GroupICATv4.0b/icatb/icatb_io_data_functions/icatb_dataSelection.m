function [files, designMatrix, numOfSub, numOfSess, dataSelMethod, diffTimePoints, spmMatFlag]  = icatb_dataSelection(...
    inputFile, outputDir, output_prefix, read_complex_file_naming, read_complex_images)
%% Data selection step in GIFT toolbox
%
% Inputs:
% 1. inputFile - inputFile containing the directory and file pattern for
% subjects
% 2. outputDir - Output directory for the analysis.
% 3. output_prefix - File output prefix
% 4. read_complex_file_naming - File naming for reading complex images
% 5. read_complex_images - Complex data type (real&imaginary) or
% (magnitude&phase)
%
% Output:
% 1. files - data structure containing the subject images. Data will be stored
% in sessions.
% 2. designMatrix - design matrix structure
% 3. numOfSub - Number of subjects selected
% 4. numOfSess - Number of sessions selected
% 5. dataSelMethod - Data selection method type
% 6. diffTimePoints - Time points vector.
% 7. spmMatFlag - flag to distinguish how to select the design matrix.

if ~exist('inputFile', 'var')
    inputFile = [];
end

% Get modality type
[modalityType, dataTitle] = icatb_get_modality;

fileType = 'any';
if strcmpi(modalityType, 'eeg')
    filterTextStr = '*.mat';
    type_file_selection = 'single';
else
    filterTextStr = '*.img;*.nii';
    type_file_selection = 'multiple';
    fileType = 'image';
end

if ispc
    windowStyle = 'modal';
else
    windowStyle = 'normal';
end


% Initialise output vars
files.name = [];
designMatrix.name = [];
dataSelMethod = 2;
diffTimePoints = [];
spmMatFlag = 'no';
numOfSub = 0;
numOfSess = 0;

if ~exist('outputDir', 'var')
    outputDir = pwd;
end

if ~exist('output_prefix', 'var')
    output_prefix = '';
end

if ~exist('inputFile', 'var')
    inputFile = [];
end

if ~exist('dataType', 'var')
    dataType = 'real';
end

if ~exist('read_complex_file_naming', 'var')
    read_complex_file_naming = [];
end

if ~exist('read_complex_images', 'var')
    read_complex_images = [];
end

fileNumbers = [];

% get the subject matrix file naming
subjectFile = [output_prefix, 'Subject.mat'];
subjectFile = fullfile(outputDir, subjectFile);

if isempty(inputFile)
    % Use GUI to select the data selection method
    
    
    if (~strcmpi(modalityType, 'smri'))
        
        % open a popup window asking the user to select the data
        questionString = 'Is your data stored in one group folder?';
        choiceString = str2mat('Yes', 'No');
        datasel_Handle = icatb_getGraphics('Selecting data method', 'normal', 'data selection'); % figure handle
        set(datasel_Handle, 'menubar', 'none', 'windowStyle', windowStyle);
        
        
        % plot menu on this figure window
        % Help Menu
        helpMenu = uimenu('parent', datasel_Handle, 'label', 'GIFT-Help');
        htmlHelpMenu = uimenu(helpMenu, 'label', 'Data Selection', 'callback', ...
            'icatb_openHTMLHelpFile(''icatb_select_data.htm'');');
        
        popupAnswer = icatb_promptUI('popup', questionString, choiceString, 'numeric', datasel_Handle);
        delete(datasel_Handle);
        
        % Data selection method
        if popupAnswer == 1
            dataSelMethod = 1;
        else
            dataSelMethod = 2;
        end
        
    else
        
        dataSelMethod = 2;
        
    end
    
else
    
    if ischar(inputFile)
        
        %% Evaluate function only once as there may be lot of subjects and
        % sessions.
        inputData = icatb_eval_script(inputFile);
        
    else
        inputData = inputFile;
        inputFile = inputData.inputFile;
    end
    
    if strcmpi(modalityType, 'smri')
        dataSelMethod = 4;
    else
        try
            dataSelMethod = checkCriteria(inputData, 'dataSelectionMethod', 'integer', inputFile);
        catch
            dataSelMethod = 2;
        end
    end
    
    % Validate data selection method
    if (dataSelMethod ~= 1) && (dataSelMethod ~= 2) && (dataSelMethod ~= 3) && (dataSelMethod ~= 4)
        error('Please provide a valid option for data selection. Data selection options are 1, 2, 3 and 4');
    end
    
    try
        spmMatFlag = checkCriteria(inputData, 'keyword_designMatrix', 'character', inputFile);
        spmMatFlag = lower(spmMatFlag);
    catch
        spmMatFlag = 'no';
    end
    
    if ~strcmpi(modalityType, 'fmri')
        spmMatFlag = 'no';
    end
    
    % change the spm mat flag
    if strcmpi(spmMatFlag, 'one')
        spmMatFlag = 'same_sub_same_sess';
    end
    
    % change the spm mat flag
    if strcmpi(spmMatFlag, 'all')
        spmMatFlag = 'diff_sub_diff_sess';
    end
    
    % get the spm mat flag
    if ~strcmp(spmMatFlag, 'same_sub_same_sess') && ~strcmp(spmMatFlag, 'same_sub_diff_sess') && ...
            ~strcmp(spmMatFlag, 'diff_sub_diff_sess') && ~strcmp(spmMatFlag, 'no')
        error('Error:SPMDesign', 'Please check the keyword_designMatrix variable in input file %s\n', inputFile);
    end
    
    % For all subjects search for the suffix designMat
    if  strcmpi(spmMatFlag, 'same_sub_diff_sess') || strcmpi(spmMatFlag, 'same_sub_same_sess')
        designMatrix.name = checkCriteria(inputData, 'OnedesignMat', 'file', inputFile);
    end
    
end


% answer to question 1
if (dataSelMethod == 1)
    
    if isempty(inputFile)
        data_setDir = icatb_selectEntry('typeEntity', 'directory', 'title', ...
            'Select root folder for subjects and sessions');
        drawnow;
        if ~isempty(data_setDir)
            % dialog Title
            dlg_title = [dataTitle, ' data information.'];
            
            numParameters = 1;
            
            %                 inputText(numParameters).promptString = 'Select data type.';
            %                 inputText(numParameters).uiType = 'popup';
            %                 inputText(numParameters).answerString = {'Real', 'Complex'};
            %                 inputText(numParameters).dataType = 'string';
            %                 inputText(numParameters).tag = 'datatype';
            %                 inputText(numParameters).enable = 'on';
            %
            %                 numParameters = numParameters + 1;
            
            inputText(numParameters).promptString = 'Select file pattern for reading data.';
            inputText(numParameters).uiType = 'edit';
            inputText(numParameters).answerString = filterTextStr;
            inputText(numParameters).dataType = 'string';
            inputText(numParameters).tag = 'filepattern';
            inputText(numParameters).enable = 'on';
            
            numParameters = numParameters + 1;
            
            inputText(numParameters).promptString = 'Are session folders inside subject folders?';
            inputText(numParameters).uiType = 'popup';
            inputText(numParameters).answerString = str2mat('Yes', 'No');
            inputText(numParameters).dataType = 'string';
            inputText(numParameters).tag = 'data_folder';
            inputText(numParameters).enable = 'on';
            
            if strcmpi(modalityType, 'fmri')
                numParameters = numParameters + 1;
                
                inputText(numParameters).promptString = 'Enter file numbers to include. Leave empty if you want to select all.';
                inputText(numParameters).uiType = 'edit';
                inputText(numParameters).answerString = '';
                inputText(numParameters).dataType = 'string';
                inputText(numParameters).tag = 'file_numbers';
                inputText(numParameters).enable = 'on';
            end
            
            %numParameters = numParameters + 1;
            
            %                 inputText(numParameters).promptString = 'Select an option to read complex images.';
            %                 inputText(numParameters).uiType = 'popup';
            %                 inputText(numParameters).answerString = {'Real&Imaginary', 'Magnitude&Phase'};
            %                 inputText(numParameters).dataType = 'string';
            %                 inputText(numParameters).tag = 'read_complex_images';
            %                 inputText(numParameters).enable = 'off';
            %
            %                 numParameters = numParameters + 1;
            %
            %
            %                 inputText(numParameters).promptString = 'Select an option to write complex images.';
            %                 inputText(numParameters).uiType = 'popup';
            %                 inputText(numParameters).answerString = {'Real&Imaginary', 'Magnitude&Phase'};
            %                 inputText(numParameters).dataType = 'string';
            %                 inputText(numParameters).tag = 'write_complex_images';
            %                 inputText(numParameters).enable = 'off';
            
            numUIControls = length(inputText);
            
            % Input dialog box (get the necessary numbers)
            answer = icatb_inputDialog('inputtext', inputText, 'Title', dlg_title, 'windowStyle', windowStyle);
            
            getSPMMatrix = 'no';
            if ~isempty(answer)
                % data type
                %dataType = answer{1};
                % file pattern
                filePattern = answer{1};
                if strcmpi(answer{2}, 'yes')
                    % flag folder
                    flagFolder = 'data_in_subject_subfolder';
                else
                    flagFolder = 'data_in_subject_folder';
                end
                % read complex images
                %read_complex_images = answer{4};
                % write complex images
                %write_complex_images = answer{5};
            else
                error('information regarding data type and images should be specified');
            end
            
        else
            error('data-sets directory is not selected');
        end
        
        designMatrix.name = [];
        
        if strcmpi(modalityType, 'fmri')
            try
                if ~isempty(answer{3})
                    fileNumbers = str2num(answer{3});
                end
            catch
                disp('File numbers are not entered correctly');
            end
        end
        
    else
        
        %% read data-sets directory and filter file pattern
        sourceDir_fileP_flag = checkCriteria(inputData, 'sourceDir_filePattern_flagLocation', 'cell', inputFile);
        data_setDir = sourceDir_fileP_flag{1}; % Data-set directory
        filePattern = sourceDir_fileP_flag{2}; % File pattern
        flagFolder = sourceDir_fileP_flag{3}; % Flag for folder
        
        if strcmpi(modalityType, 'fmri')
            % Get the file numbers if possible
            if length(sourceDir_fileP_flag) == 4
                fileNumbers = sourceDir_fileP_flag{4};
            end
        end
        
    end
    % end for getting the information through input file or GUI
    
    drawnow;
    disp(['Reading data from source directory ', data_setDir, ' ...']);
    
    % get the data_sets automatically
    [files, numOfSub, numOfSess, selected_data_sets, data_folders] = icatb_get_sub_data(data_setDir, filePattern, flagFolder);
    
    % There are two ways to select the data get design matrix for each subject
    if strcmpi(spmMatFlag, 'diff_sub_diff_sess')
        
        designMatrix = repmat(struct('name', []), 1, numOfSub);
        
        spmDesignFilter = checkCriteria(inputData, 'spmDesignFilter', 'string', inputFile);
        
        %         keywd = 'spmDesignFilter';
        %         inputData = icatb_read_variables(inputFile, keywd, 'scalar', 'string');
        %         spmDesignFilter = getfield(inputData, keywd);
        %         clear inputData;
        
        % loop over subjects
        for ii = 1:numOfSub
            % current subject folder
            currentSub = deblank(data_folders((ii-1)*numOfSess + 1).name);
            % design matrix
            tempdesignMat = icatb_listFiles_inDir(currentSub, spmDesignFilter);
            if isempty(tempdesignMat)
                try
                    tempdesignMat = icatb_get_sub_data(currentSub, spmDesignFilter, 'data_in_subject_folder', 0);
                catch
                    error(['Cannot find SPM design matrix for subject ', num2str(ii), ...
                        ' in its folder or sub-folders where folder path is ', ...
                        currentSub]);
                end
                
                % get the full file path for SPM design matrix
                designMatrix(ii).name = deblank(tempdesignMat(1).name(1, :));
            else
                % get the full file path for SPM design matrix
                designMatrix(ii).name = fullfile(currentSub, deblank(tempdesignMat(1, :)));
            end
            clear tempdesignMat;
        end
        % end loop over subjects
    end
    % end for getting the design matrix
    
    selectedSubTxtFile = fullfile(outputDir, [output_prefix, 'SelectedDataFolders.txt']);
    % open the file for recording data-sets
    fid = fopen(selectedSubTxtFile, 'w');
    for nLines = 1:size(selected_data_sets, 1)
        tempStr = deblank(selected_data_sets(nLines, :));
        fprintf(fid, '%s\n', tempStr);
    end
    fclose(fid);
    fprintf('\n');
    % closing the file
    disp(['Please see the text file ', selectedSubTxtFile, ' for the selected data folders in order']);
    
    nF = 0;
    diffTimePoints = zeros(1, numOfSub*numOfSess);
    % Loop over subjects
    for nSub = 1:numOfSub
        % Loop over sessions
        for nSess = 1:numOfSess
            % Loop over files
            % for nF = 1:length(files)
            nF = nF + 1;
            tempFiles = icatb_rename_4d_file(files(nF).name);
            currentFileNum = fileNumbers;
            if ~isempty(currentFileNum)
                % Exclude the file number that exceed the number of files
                currentFileNum(currentFileNum > size(tempFiles, 1)) = [];
                if isempty(currentFileNum)
                    error(['Unable to find the images with the file numbers you have specified for subject ', num2str(nSub), ...
                        ' session ', num2str(nSess)]);
                end
            else
                % Use all files
                currentFileNum = (1:size(tempFiles, 1));
            end
            files(nF).name = tempFiles(currentFileNum, :);
            
            diffTimePoints(nF) = size(files(nF).name, 1);
            %end
            % End loop over files
        end
        % End loop over sessions
    end
    % End loop over subjects
    
elseif (dataSelMethod == 2)
    
    
    if isempty(inputFile)
        
        
        if (~strcmpi(modalityType, 'smri'))
            
            % dialog Title
            dlg_title = [dataTitle, ' data information.'];
            
            numParameters = 1;
            
            % define all the input parameters in a structure
            inputText(numParameters).promptString = 'Number of Subjects?';
            inputText(numParameters).uiType = 'edit';
            inputText(numParameters).answerString = '1';
            inputText(numParameters).dataType = 'numeric';
            inputText(numParameters).tag = 'num_subjects';
            inputText(numParameters).enable = 'on';
            
            numParameters = numParameters + 1;
            
            inputText(numParameters).promptString = 'Number of Sessions Per Subject?';
            inputText(numParameters).uiType = 'edit';
            inputText(numParameters).answerString = '1';
            inputText(numParameters).dataType = 'numeric';
            inputText(numParameters).tag = 'num_sessions';
            inputText(numParameters).enable = 'on';
            
            %             numParameters = numParameters + 1;
            %
            %             inputText(numParameters).promptString = 'Select data type.';
            %             inputText(numParameters).uiType = 'popup';
            %             inputText(numParameters).answerString = {'Real', 'Complex'};
            %             inputText(numParameters).dataType = 'string';
            %             inputText(numParameters).tag = 'datatype';
            %             inputText(numParameters).enable = 'on';
            %
            %             numParameters = numParameters + 1;
            %
            %
            %             inputText(numParameters).promptString = 'Select an option to read complex images.';
            %             inputText(numParameters).uiType = 'popup';
            %             inputText(numParameters).answerString = {'Real&Imaginary', 'Magnitude&Phase'};
            %             inputText(numParameters).dataType = 'string';
            %             inputText(numParameters).tag = 'read_complex_images';
            %             inputText(numParameters).enable = 'off';
            %
            %             numParameters = numParameters + 1;
            %
            %
            %             inputText(numParameters).promptString = 'Select an option to write complex images.';
            %             inputText(numParameters).uiType = 'popup';
            %             inputText(numParameters).answerString = {'Real&Imaginary', 'Magnitude&Phase'};
            %             inputText(numParameters).dataType = 'string';
            %             inputText(numParameters).tag = 'write_complex_images';
            %             inputText(numParameters).enable = 'off';
            
            numUIControls = length(inputText);
            
            % Input dialog box (get the necessary numbers)
            answer = icatb_inputDialog('inputtext', inputText, 'Title', dlg_title, 'windowStyle', windowStyle);
            
            if ~isempty(answer)
                % number of subjects
                numOfSub = answer{1};
                % number of sessions
                numOfSess = answer{2};
                getSPMMatrix = 'no'; %answer{3};
                %                 dataType = answer{3};
                %                 % read complex images
                %                 read_complex_images = answer{4};
                %                 % write complex images
                %                 write_complex_images = answer{5};
            else
                error('Number of data sets is not specified');
            end
            
        else
            
            numOfSub = 1;
            numOfSess = 1;
            
        end
        
        fileInfo.fileName = subjectFile;
        fileInfo.format = '';
        
        %             sesInfo.userInput.read_complex_images = lower(read_complex_images);
        %
        %             sesInfo.userInput.write_complex_images = lower(write_complex_images);
        
        
        % for one subject and one session
        if numOfSub*numOfSess == 1
            
            % select data
            if (~strcmpi(modalityType, 'smri'))
                dataWindowTitle = 'Select files for subject 1 session 1';
            else
                dataWindowTitle = 'Select subject images for SBM';
            end
            
            files(1).name = icatb_selectEntry('typeEntity', 'file', 'typeSelection', type_file_selection, ...
                'filter', filterTextStr, 'title', dataWindowTitle, 'fileType', fileType);
            drawnow;
            if isempty(files(1).name)
                error([dataTitle, ' data is not selected for the analysis']);
            else
                % get count for time points
                %[diffTimePoints] = icatb_get_countTimePoints(files);
                diffTimePoints = size(files(1).name, 1);
            end
            
        else
            
            if (strcmpi(modalityType, 'fmri'))
                %% select the data for fmri
                [files, diffTimePoints] = icatb_select_data('title', ['Select ', modalityType, ' data'], 'num_subjects', ...
                    numOfSub, 'num_sessions', numOfSess, 'files_specification', 'unequal', 'spm_check', 'no', ...
                    'filter_string', filterTextStr, 'type_file_selection', type_file_selection, 'fileInfo', fileInfo, 'figure_menu', 'data', ...
                    'datatype', dataType, 'complex_file_naming', read_complex_file_naming, ...,
                    'read_complex_images', read_complex_images, 'fileType', fileType, 'windowStyle', windowStyle);
                drawnow;
            else
                %% Select the data for eeg
                temporary_files = icatb_selectEntry('typeEntity', 'file', 'typeSelection', 'multiple', ...
                    'filter', filterTextStr, 'title', 'Select files for all data-sets');
                
                drawnow;
                if (isempty(temporary_files))
                    error('Error:InputFiles', 'Files are not selected');
                end
                
                if (size(temporary_files, 1) ~= numOfSub*numOfSess)
                    error('Error:NumDataSets', 'Number of files selected (%d) must match the number of datasets (%d)\n', ...
                        size(temporary_files, 1), numOfSub*numOfSess);
                end
                
                % Initialise files
                files = repmat(struct('name', ''), 1, numOfSub*numOfSess);
                count_files = 0;
                % Loop over subjects
                for nSub = 1:numOfSub
                    % Loop over sessions
                    for nSess = 1:numOfSess
                        count_files = count_files + 1;
                        files(count_files).name = deblank(temporary_files(count_files, :));
                    end
                    % End of loop over sessions
                end
                % End of loop over subjects
                
                clear temporary_files;
                
                diffTimePoints = ones(1, numOfSub*numOfSess);
                
            end
            
        end
        
    else
        %% Select data using the file paths for each data set
        
        % Selected subjectss
        selectedSubjects = checkCriteria(inputData, 'selectedSubjects', 'cell', inputFile);
        numOfSub = length(selectedSubjects);
        
        % Number of sessions
        numOfSess = checkCriteria(inputData, 'numOfSess', 'integer', inputFile);
        numOfSess = numOfSess(1);
        if (numOfSess == 0)
            error('Number of sessions cannot be zero');
        end
        
        disp('Reading parameters for getting input data files ...');
        
        if (strcmpi(modalityType, 'fmri') && strcmpi(spmMatFlag, 'diff_sub_diff_sess'))
            designMatrix = repmat(struct('name', []), 1, numOfSub);
        end
        
        % Initialise some variables
        allDirs = cell(numOfSub*numOfSess, 1);
        filePatterns = cell(numOfSub*numOfSess, 1);
        fileNums = cell(numOfSub*numOfSess, 1);
        
        countD = 0;
        %% Loop over number of subjects
        for i = 1 : numOfSub
            % Strings corresponding to the selected subjects
            subStr = selectedSubjects{i};
            
            % one design matrix for each subject
            if (strcmpi(modalityType, 'fmri') && strcmpi(spmMatFlag, 'diff_sub_diff_sess'))
                % For all subjects search for the suffix designMat
                keywd = [subStr, '_designMat'];
                designMatrix(i).name = checkCriteria(inputData, keywd, 'string', inputFile);
            end
            
            %% Loop over sessions
            for j = 1 : numOfSess
                countD = countD + 1;
                sessStr = ['_s', num2str(j)];
                keywd = [subStr, sessStr];
                Value_vector = checkCriteria(inputData, keywd, 'cell', inputFile);
                allDirs{countD} = Value_vector{1};
                filePatterns{countD} = Value_vector{2};
                if strcmpi(modalityType, 'fmri')
                    if (length(Value_vector) == 3)
                        fileNums{countD} = Value_vector{3};
                    else
                        fileNums{countD} = [];
                    end
                end
            end
            %% End for loop over sessions
        end
        %% End for loop over subjects
        
        [files, diffTimePoints] = listFilesWithFullPaths(allDirs, filePatterns, fileNums, numOfSub, numOfSess);
        
        disp('Done reading variables for getting input data files');
        
        fprintf('\n');
        
    end
    
elseif (dataSelMethod == 3)
    %% Select data using regular expressions
    
    % Get fields
    input_directory_name = inputData.input_directory_name;
    subject_dir_regexp = inputData.subject_dir_regexp;
    session_dir_regexp = inputData.session_dir_regexp;
    data_file_pattern = inputData.data_file_pattern;
    try
        file_numbers_to_include = inputData.file_numbers_to_include;
    catch
        file_numbers_to_include = [];
    end
    
    if strcmpi(modalityType, 'fmri') && strcmpi(spmMatFlag, 'diff_sub_diff_sess')
        try
            spm_stats_dir = inputData.spm_stats_dir;
        catch
            spm_stats_dir = '';
        end
    end
    
    disp('Listing files based on regular expressions ...');
    fprintf('\n');
    
    % Select data based on regular expressions
    [allDirs, numOfSub, numOfSess] = icatb_select_data_regexp(input_directory_name, subject_dir_regexp, session_dir_regexp);
    
    allDirs = cellstr(allDirs);
    
    % List files with full paths and get the number of time points
    % information
    [files, diffTimePoints] = listFilesWithFullPaths(allDirs, data_file_pattern, file_numbers_to_include, numOfSub, numOfSess);
    
    if strcmpi(modalityType, 'fmri') && strcmpi(spmMatFlag, 'diff_sub_diff_sess')
        
        if isempty(subject_dir_regexp) && isempty(session_dir_regexp)
            %% If the files are in the input directory itself
            currentDir = fullfile(allDirs{1}, spm_stats_dir);
            fileListWithDir = icatb_listFiles_inDir(currentDir, 'SPM.mat');
            if isempty(fileListWithDir)
                error('Error:SPMFile', 'SPM mat file doesn''t exist in directory (%s) \n', currentDir);
            end
            designMatrix(1).name = fullfile(currentDir, deblank(fileListWithDir(1, :)));
            clear currentDir fileListWithDir;
        elseif isempty(subject_dir_regexp) && ~isempty(session_dir_regexp)
            %% If the files are in session directory itself
            % Loop over sessions
            for nSess = 1:numOfSess
                currentDir = fullfile(allDirs{nSess}, spm_stats_dir);
                if (exist(currentDir, 'dir') == 7)
                    fileListWithDir = icatb_listFiles_inDir(currentDir, 'SPM.mat');
                    if ~isempty(fileListWithDir)
                        designMatrix(1).name = fullfile(currentDir, deblank(fileListWithDir(1, :)));
                        break;
                    end
                end
            end
            % End loop over sessions
            
            clear currentDir fileListWithDir;
            
            if isempty(designMatrix(1).name)
                error('Error:SPMFile', 'SPM mat file doesn''t exist in any of the session directories\n');
            end
            
        else
            
            designMatrix = repmat(struct('name', []), 1, numOfSub);
            %% Loop over subjects
            for nSub = 1:numOfSub
                fileListWithDir = [];
                if ~isempty(session_dir_regexp)
                    subDir = fileparts(allDirs{(nSub - 1)*numOfSess + 1});
                    currentDir = fullfile(subDir, spm_stats_dir);
                    if (exist(currentDir, 'dir') == 7)
                        fileListWithDir = icatb_listFiles_inDir(currentDir, 'SPM.mat');
                    end
                end
                if isempty(fileListWithDir)
                    %% Loop over sessions
                    for nSess = 1:numOfSess
                        currentDir = fullfile(allDirs{(nSub - 1)*numOfSess + nSess}, spm_stats_dir);
                        if (exist(currentDir, 'dir') == 7)
                            fileListWithDir = icatb_listFiles_inDir(currentDir, 'SPM.mat');
                            if ~isempty(fileListWithDir)
                                designMatrix(nSub).name = fullfile(currentDir, deblank(fileListWithDir(1, :)));
                                break;
                            end
                        end
                    end
                    %% End loop over sessions
                else
                    designMatrix(nSub).name = fullfile(currentDir, deblank(fileListWithDir(1, :)));
                end
                
                if isempty(designMatrix(nSub).name)
                    error('Error:SPMFile', 'Design matrix doesn''t exist for subject %s\n', allDirs{(nSub - 1)*numOfSess + 1});
                end
                
                clear fileListWithDir;
                
            end
            %% End loop over subjects
            
        end
        % End for handling design matrix
    end
    
    disp('Done');
    fprintf('\n');
    
    
else
    
    % Get fields
    input_data_file_patterns = inputData.input_data_file_patterns;
    
    if (~iscell(input_data_file_patterns))
        %         if (strcmpi(modalityType, 'smri'))
        %             input_data_file_patterns = {input_data_file_patterns};
        %         else
        input_data_file_patterns = cellstr(input_data_file_patterns);
        %end
    end
    
    good_inds = icatb_good_cells(input_data_file_patterns);
    if (~isempty(find(good_inds == 0)))
        error('Error:InputFilePatterns', 'Some of the filepatterns are empty. Please check variable input_data_file_patterns in file %s\n', ...
            inputFile);
    end
    
    if (strcmpi(modalityType, 'smri'))
        input_data_file_patterns = {char(input_data_file_patterns)};
        %input_data_file_patterns = input_data_file_patterns(1);
    end
    
    %% Number of subjects and sessions
    numOfSub = size(input_data_file_patterns, 1);
    numOfSess = size(input_data_file_patterns, 2);
    
    dummy_scans = 0;
    
    %% Check dummy scans and design matrices for fmri
    if (strcmpi(modalityType, 'fmri'))
        
        dummy_scans = inputData.dummy_scans;
        
        %% Different design between subjects
        if (strcmpi(spmMatFlag, 'diff_sub_diff_sess'))
            
            input_design_matrices = inputData.input_design_matrices;
            if (~iscell(input_design_matrices))
                input_design_matrices = cellstr(input_design_matrices);
            end
            good_inds = icatb_good_cells(input_design_matrices);
            if (~isempty(find(good_inds == 0)))
                error('Error:DesignMatrices', 'Some of the design matrices are empty. Please check variable input_design_matrices in file %s\n', ...
                    inputFile);
            end
            
            if (length(input_design_matrices) ~= numOfSub)
                error('Error:NumDesignMat', 'Length of design matrices must equal no. of subjects (%d)\n', numOfSub);
            end
            
            designMatrix = repmat(struct('name', []), 1, numOfSub);
            %% Loop over subjects
            for nSub = 1:numOfSub
                designMatrix(nSub).name = deblank(input_design_matrices{nSub});
            end
            
        end
    end
    
    % Initialise input file information
    files = repmat(struct('name', []), 1, numOfSub*numOfSess);
    diffTimePoints = zeros(1, numOfSub*numOfSess);
    counter = 0;
    %% Loop over subjects
    for nSub = 1:numOfSub
        %% Loop over sessions
        for nSess = 1:numOfSess
            counter = counter + 1;
            fileP = input_data_file_patterns{nSub, nSess};
            if (size(fileP, 1) == 1)
                [pathstr, fp, extn] = fileparts(fileP);
                fileContents = icatb_listFiles_inDir(pathstr, [fp, extn]);
                if (isempty(fileContents))
                    error('Error:FilePattern', 'Please check file pattern %s as there are no files found\n', fileP);
                end
                fileListWithDir = icatb_fullFile('directory', pathstr, 'files', fileContents);
            else
                fileListWithDir = fileP;
            end
            
            fileListWithDir = icatb_rename_4d_file(fileListWithDir);
            
            if (length(dummy_scans) > 1)
                fileNum = dummy_scans;
                fileNum(fileNum > size(fileListWithDir, 1)) = [];
                if isempty(fileNum)
                    error('Error:DummyScans', 'Cannot find the files specified with file numbers. Please check dummy_scans variable for \nfile %s\n', fileP);
                end
                fileListWithDir = fileListWithDir(fileNum, :);
            else
                if (dummy_scans > 0)
                    if (dummy_scans >= size(fileListWithDir, 1))
                        error('Error:DummyScans', 'Please check dummy_scans variable (%d) as it exceeds or equals the no. of time points (%d) for\nfile %s\n',  ...
                            dummy_scans, size(fileListWithDir, 1), fileP);
                    end
                    fileNum = (1:size(fileListWithDir, 1));
                    fileNum(1:dummy_scans) = [];
                    fileListWithDir = fileListWithDir(fileNum, :);
                    
                end
            end
            
            %% For EEG only one file is allowed per dataset
            if strcmpi(modalityType, 'eeg')
                fileListWithDir = deblank(fileListWithDir(1, :));
            end
            
            files(counter).name = fileListWithDir; % append the file list with directory
            diffTimePoints(counter) = size(fileListWithDir, 1);
            
        end
        %% End of loop over sessions
    end
    %% End of loop over subjects
    
    
end

drawnow;

if (strcmpi(modalityType, 'eeg'))
    fprintf('\n');
    disp('Checking data dimensions for eeg data ...');
    filesN = str2mat(files.name);
    if size(filesN, 1) ~= numOfSub*numOfSess
        error('You are allowed to enter only one file per subject. The no. of files is %s and number of datasets is %s', ...
            num2str(size(filesN, 1)), num2str(numOfSub*numOfSess));
    end
    diffTimePoints = icatb_get_num_electrodes(filesN);
    disp('End for checking data dimensions for eeg data');
    fprintf('\n');
else
    % Check spm design matrix
    icatb_check_spm_design_matrix(designMatrix, numOfSub, numOfSess, diffTimePoints, spmMatFlag, inputFile);
end

function varOut = checkCriteria(inputData, keywd, varCriteria, inputFile)
%% Check the variable criteria
%

if (~isfield(inputData, keywd))
    error('Error:InputFile', 'Variable %s doesn''t exist in file %s\n', keywd, inputFile);
end

varOut = getfield(inputData, keywd);

varCheck = varOut;
if isnumeric(varOut)
    varCheck = num2str(varCheck);
end

%% Do error check
[status, message] = icatb_errorCheck(varCheck, varCriteria, keywd);

if (status == 0)
    error('Error:InputFile', '%s. Please check the input file %s\n', message, inputFile);
end


function [files, diffTimePoints] = listFilesWithFullPaths(allDirs, filePatterns, fileNums, numOfSub, numOfSess)
%% List files with full paths
%

if ~iscell(filePatterns)
    filePatterns = cellstr(filePatterns);
end

if ~iscell(fileNums)
    fileNums = {fileNums};
end

if (length(filePatterns) == 1)
    filePatterns = repmat(filePatterns, numOfSub*numOfSess, 1);
end

if (length(fileNums) == 1)
    fileNums = repmat(fileNums, numOfSub*numOfSess, 1);
end

counter = 0;
diffTimePoints = zeros(1, numOfSub*numOfSess);

% Initialise input file information
files = repmat(struct('name', []), 1, numOfSub*numOfSess);
for j = 1 : numOfSub
    for k = 1 : numOfSess
        counter = counter + 1;
        % get scans and parse
        fileDir = allDirs{counter};
        filePattern = filePatterns{counter};
        fileNum = fileNums{counter};
        
        % list Files with the matching pattern
        fileList = icatb_listFiles_inDir(fileDir, filePattern);
        
        if (isempty(fileList))
            icatb_error('Could not find any files matching pattern', {fileDir, filePattern});
        end
        % get the full file path for the files
        fileListWithDir = icatb_fullFile('directory', fileDir, 'files', fileList);
        fileListWithDir = icatb_rename_4d_file(fileListWithDir);
        if ~isempty(fileNum)
            fileNum(fileNum > size(fileListWithDir, 1)) = [];
            if isempty(fileNum)
                error(['Unable to find the files with the file numbers you have specified for subject ', num2str(j), ...
                    ' session ', num2str(k)]);
            end
            fileListWithDir = fileListWithDir(fileNum, :);
        end
        files(counter).name = fileListWithDir; % append the file list with directory
        diffTimePoints(counter) = size(fileListWithDir, 1);
    end
end
