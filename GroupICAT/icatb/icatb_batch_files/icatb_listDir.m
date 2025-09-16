function [dirContents] = icatb_listDir(inputDir)

% List sub folders in a directory
oldDir = pwd;

if ~exist('inputDir', 'var')
    inputDir = pwd;
end

cd(inputDir);

specifyDir = dir; % List directories and files

% Return zeros and ones for the directories or not
listDir = [specifyDir.isdir];

% Find only directories
tempDir = find(listDir == 1);  

% List all the folders
listContents = str2mat(specifyDir.name);

dirContents = [];

if ~isempty(tempDir)
    dirContents = listContents(tempDir, :);
end   

% change directory
cd(oldDir); 