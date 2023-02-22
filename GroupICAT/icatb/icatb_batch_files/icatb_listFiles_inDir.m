function [fileContents] = icatb_listFiles_inDir(input_directory, filesPattern)

% List files in a directory with the matching
% file pattern
fileContents = [];

oldDir = pwd;

cd(input_directory);

if ~exist('filesPattern', 'var')
    filesPattern = '*';
end

if isempty(filesPattern)
    filesPattern = '*';
end

if ~exist('fileNumber', 'var')
    fileNumber = [];
end

filesPattern = strread(filesPattern, '%s', 'delimiter', ';');

count = 0;
% loop over file pattern
for ii = 1:length(filesPattern)

    specifyDir = dir(filesPattern{ii}); % List directories and files

    % Return zeros and ones for the directories or not
    listDir = [specifyDir.isdir];

    % Find only files
    tempFiles = find(listDir == 0);

    % List all the folder or file names
    listContents = str2mat(specifyDir.name);

    if ~isempty(tempFiles)
        count = count + 1;
        fileC = listContents(tempFiles, :);
        fileContents(count).name = fileC;
        clear fileC;
    end

    clear listDir; clear tempFiles;
end
% end loop over filter pattern

if ~isempty(fileContents)
    temp = fileContents;
    clear fileContents;
    fileContents = str2mat(temp.name);
    % get unique strings only for more than one file pattern
    if length(filesPattern) > 1
        clear temp;
        [uniqueStr] = icatb_getUniqueStrings(fileContents);
        clear fileContents;
        fileContents = str2mat(uniqueStr);
    end
    % end for getting file contents

end


% change directory
cd(oldDir);