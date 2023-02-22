function icatb_delete_file_pattern(fileN, outputDir)
% form component file patterns

if iscell(fileN)
    fileN = str2mat(fileN);
end

if ~exist('outputDir', 'var') || isempty(outputDir)
    outputDir = pwd;
end

% delete component files after compressing
temp = icatb_spm_get('fileSummary', fileN);
temp = icatb_fullFile('directory', outputDir, 'files', temp);
files_to_deleteP = cellstr(temp);
clear temp;

% delete the files
for ii = 1:length(files_to_deleteP)
    drawnow;
    delete(files_to_deleteP{ii});       
end   