function [zipFileName, files_in_zip] = icatb_getViewingSet_zip(compFiles, complexInfoWrite, dataType, zipContents)

if ~strcmpi(dataType, 'real')
    [checkPath, checkCompFile, extn] = fileparts(deblank(compFiles(1, :)));
    checkCompFile = fullfile(checkPath, [complexInfoWrite.complex_file_naming{1}, ...
        checkCompFile, extn]);
else
    checkCompFile = [deblank(compFiles(1, :))];
end


zipCount = 0;
zipFileName = [];
files_in_zip = {};
for ii = 1:length(zipContents.zipFiles)
    % get the first file
    %firstFile = zipContents.files_in_zip(ii).name{1};
    if ispc
        compareFiles = strmatch(lower(checkCompFile), lower(zipContents.files_in_zip(ii).name), 'exact');
    else
        compareFiles = strmatch(checkCompFile, zipContents.files_in_zip(ii).name, 'exact');
    end
    if ~isempty(compareFiles)
        zipCount = ii;
        break;
    end
end

% get the
if zipCount > 0
    zipFileName = zipContents.zipFiles{zipCount};
    files_in_zip = zipContents.files_in_zip(zipCount).name;
end