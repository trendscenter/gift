function icatb_cleanupFiles(filesToDelete, outputDir)
%% Cleanup files
%

for nF = 1:length(filesToDelete)
    files_in_zip = filesToDelete{nF};
    files_in_zip = str2mat(files_in_zip);
    icatb_delete_file_pattern(files_in_zip, outputDir);
end