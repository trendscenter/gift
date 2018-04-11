function [zipContents] = icatb_getUniqueZipFiles(zipContents)
% remove repeated sequence of zip files

lenZipFile = length(zipContents.zipFiles);
if lenZipFile > 0
    lenZipFile = length(zipContents.zipFiles{1});
end


if lenZipFile ~= 0
    % ignore the repeated zip file names
    [zipContents.zipFiles, countVec] = icatb_getUniqueStrings(zipContents.zipFiles);
    tempZipFiles = zipContents.files_in_zip;
    zipContents = rmfield(zipContents, 'files_in_zip');
    zipContents.files_in_zip = repmat(struct('name', {}), 1, length(countVec));
    % loop over zip files
    for nZip = 1:length(countVec)
        zipContents.files_in_zip(nZip).name = tempZipFiles(nZip).name;
    end
    clear tempZipFiles;
end
