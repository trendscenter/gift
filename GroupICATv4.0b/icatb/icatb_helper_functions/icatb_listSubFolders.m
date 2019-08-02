function subFolders = icatb_listSubFolders(sourceDir, optional)
%% List sub folders
%
% Inputs:
% 1. sourceDir - Source Directory
% 2. optional - Options are 'fullpath' or 'relative'. By default full path
% will be listed.
%
% Outputs:
% subFolders - Sub folders in the source directory
%

%% List sub folders
subFolders = icatb_listDir(sourceDir);

index1 = strmatch('.', subFolders, 'exact');
index2 = strmatch('..', subFolders, 'exact');

index = [index1, index2];

CheckOne = ones(size(subFolders, 1), 1);

CheckOne(index) = 0;

tempVec = find(CheckOne ~= 0);

if ~exist('optional', 'var')
    optional = 'fullpath';
end

if ~isempty(tempVec)
    subFolders = subFolders(tempVec, :);

    %% Add full path
    if strcmpi(optional, 'fullpath')
        if (~strcmp(sourceDir(end), filesep))
            sourceDir = [sourceDir, filesep];
        end
        %subFolders = strcat(sourceDir, subFolders);
        subFolders = [repmat(sourceDir, size(subFolders, 1), 1), subFolders];
    end
else
    subFolders = [];
end



