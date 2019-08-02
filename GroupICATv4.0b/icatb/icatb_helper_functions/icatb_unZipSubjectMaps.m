function filesToDelete = icatb_unZipSubjectMaps(sesInfo, subjectICAFiles)
%% Uncompress subject files
%

if (~exist('subjectICAFiles', 'var'))
    subjectICAFiles = icatb_parseOutputFiles('icaOutputFiles', sesInfo.icaOutputFiles, 'numOfSub', sesInfo.numOfSub, 'numOfSess', sesInfo.numOfSess, 'flagTimePoints', sesInfo.flagTimePoints);
end

if ~isfield(sesInfo, 'zipContents')
    zipContents.zipFiles = {};
    zipContents.files_in_zip(1).name = {};
    sesInfo.zipContents = zipContents;
end


filesToDelete = {};
% If component files are not present unzip files
% Loop over subjects
for nSub = 1:sesInfo.numOfSub
    % Loop over sessions
    for nSess = 1:sesInfo.numOfSess
        currentFile = deblank(subjectICAFiles(nSub).ses(nSess).name(1, :));
        if ~exist(fullfile(sesInfo.outputDir, currentFile), 'file')
            % Check if zip file is present or not
            [zipFileName, files_in_zip] = icatb_getViewingSet_zip(currentFile, [], 'real', sesInfo.zipContents);
            if (~isempty(zipFileName))
                icatb_unzip(regexprep(zipFileName, ['.*\', filesep], ''), fullfile(sesInfo.outputDir, fileparts(currentFile)));
                filesToDelete{length(filesToDelete) + 1} = files_in_zip;
            end
            clear files_in_zip;
        end
    end
    % End loop over sessions
end
% End loop over subjects
% End for checking if zip files are present or not