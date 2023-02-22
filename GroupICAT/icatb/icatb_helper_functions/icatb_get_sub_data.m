function [data_sets, numOfSub, numOfSess, selected_data_sets, data_folders] = icatb_get_sub_data(sourceDir, filePattern, ...
    flagFolder, verbose)
% purpose: gets the subject data based on the filepattern and flag for
% folder
% 
% Input:
% 1. sourceDir - root directory
% 2. filePattern - wild card 
% 3. flagFolder - 'data_in_subject_folder' or 'data_in_subject_subfolder'
%
%
% Ouput: 
% 1. data_sets - array of structures with field name
% 2. numOfSub - Number of subjects
% 3. numOfSess - Number of sessions
%
% Example 1: C:\MATLAB6p5\work\Example Subjects\new_example_subjects root
% folder contains three example subjects folders like sub01_vis, sub02_vis,
% sub03_vis and the data are in those folders.
%
% Usage:
% sourceDir = C:\MATLAB6p5\work\Example Subjects\new_example_subjects; 
% filePattern = nsrs*.img; flagFolder = 'data_in_subject_folder';
% [data_sets, numOfSub, numOfSess] = icatb_get_sub_data(sourceDir,
% filePattern, flagFolder);
% 
% Example 2: I:\drive pilot data (last two) root folder contains three folders and each
% folder contains two sub-folders with the data. The folders are drive_7364400, 
% drive_4014069 and stats. drive_7364400 contains two session folders like
% 2 and 3 and drive_4014069 contains folders like 3 and 4. In the example, stats doesn't
% contain the data specified by the file pattern and therefore there are
% 2 subjects and sessions.
%
% Usage:
% sourceDir = 'I:\drive pilot data (last two)'; filePattern = 'sw*.img';
% flagFolder = 'data_in_subject_subfolder';
% [data_sets, numOfSub, numOfSess] = icatb_get_sub_data(sourceDir,
% filePattern, flagFolder);
%

% check the existence
if ~exist('flagFolder', 'var')
    flagFolder = 'data_in_subject_subfolder';
    numOfSess = 1;
end

% keyword for printing data-sets
if ~exist('verbose', 'var')
    verbose = 1;
end
% end for keyword for printing data-sets

% get subject folders
subFolders = icatb_listSubFolders(sourceDir);

% loop over subject folders
data_count = 0; sessionCount = 0;
data_sets.name = []; numOfSub = 0; numOfSess = 0;
sessionVec = zeros(1, size(subFolders, 1));
data_folders.name = [];

if ~isempty(subFolders)
    
    for ii = 1:size(subFolders, 1)
        if strcmpi(flagFolder, 'data_in_subject_folder')
            tempSub = deblank(subFolders(ii, :));
        else
            try
                tempSub = icatb_listSubFolders(deblank(subFolders(ii, :)));
            catch
                tempSub = [];
            end
        end
        % end for checking the flag
        sessionCount = 0;
        if ~isempty(tempSub)
            for jj = 1:size(tempSub, 1)
                try
                    tempF = icatb_listFiles_inDir(deblank(tempSub(jj, :)), filePattern);
                catch
                    tempF = [];
                end
                if ~isempty(tempF)
                    sessionCount = sessionCount + 1; % update session count
                    %                     if jj == 1
                    %                         subCount = subCount + 1;
                    %                     end
                    data_count = data_count + 1;
                    data_sets(data_count).name = icatb_fullFile('files', tempF, 'directory', deblank(tempSub(jj, :)));
                    data_folders(data_count).name = deblank(subFolders(ii, :));
                end
            end
        end
        % end for checking
        sessionVec(ii) = sessionCount;
    end
    % end loop over subject folders
    
else
    
    % subject folder is current folder
    subFolder = sourceDir;
    
    % get the file names in the current directory itself
    currentFiles = icatb_listFiles_inDir(subFolder, filePattern);
    if ~isempty(currentFiles)
        % full file names
        data_sets(1).name = icatb_fullFile('files', currentFiles, 'directory', sourceDir);
        sessionVec = 1;
        data_folders(1).name = subFolder;
    end
    
end
% end for checking the current folder


checkSess = find(sessionVec ~= 0);

if ~isempty(checkSess)
    % return number of subjects
    numOfSub = length(checkSess);
    sessionVec = sessionVec(checkSess);
    checkDiffSess = find(sessionVec(1) ~= sessionVec(end));
    if isempty(checkDiffSess)
        % return number of sessions
        numOfSess = sessionVec(1);
    else
        error('Number of sessions should be the same over subjects.');
    end
    
else
    
    numOfSub = 0;
    numOfSess = 0;
    
end
% end for getting number of subjects and sessions

if numOfSub*numOfSess == 0
    error('data doesn''t exist');
else
    % display selected folders in order
    % selected subjects and sessions
    for ii = 1:length(data_sets)
        sel_sub(ii).name = fileparts(deblank(data_sets(ii).name(1, :)));
    end
    fprintf('\n');
    if verbose
        fprintf(['The selected data folders are in the following order:\n']);
    end
    selected_data_sets = str2mat(sel_sub.name);
    if verbose
        disp(selected_data_sets);
    end
end