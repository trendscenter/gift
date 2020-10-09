function sesInfo = icatb_parCalibrateComponents_cluster(param_file, dataSetsToRun, conserve_disk_space)
%% Scale components. Run code in parallel
%

icatb_defaults;
global FUNCTIONAL_DATA_FILTER;
global ZIP_IMAGE_FILES;

if (ischar(param_file))
    load(param_file);
else
    sesInfo = param_file;
end

if (~exist('conserve_disk_space', 'var'))
    conserve_disk_space = 0;
end


[dd1, dd2, imExtn] = fileparts(FUNCTIONAL_DATA_FILTER);

modalityType = sesInfo.modality;
componentSigns = sesInfo.scale_opts.componentSigns;
meanImOffsets = sesInfo.scale_opts.meanImOffsets;
compSetFields = sesInfo.scale_opts.compSetFields;
outputDir = sesInfo.outputDir;
mask_ind = sesInfo.mask_ind;
numOfSess = sesInfo.numOfSess;
numOfSub = sesInfo.numOfSub;
back_reconstruction_mat_file = sesInfo.back_reconstruction_mat_file;
calibrate_components_mat_file = sesInfo.calibrate_components_mat_file;

[subjectICAFiles, meanICAFiles, tmapICAFiles, meanALL_ICAFile] = icatb_parseOutputFiles('icaOutputFiles', sesInfo.icaOutputFiles, 'numOfSub', sesInfo.numOfSub, 'numOfSess', sesInfo.numOfSess, ...
    'flagTimePoints', sesInfo.flagTimePoints);

outputFileNames = getOutputFileNames(subjectICAFiles, dataSetsToRun, numOfSess);

zipFileName = cell(1, length(dataSetsToRun));
files_in_zip = cell(1, length(dataSetsToRun));
tmpSessInfo = sesInfo;
zipImages = ZIP_IMAGE_FILES;

parfor i = 1:length(dataSetsToRun)
    
    tmpComponentSigns = componentSigns;
    tmpOffsets = meanImOffsets;
    
    fieldName  = compSetFields;
    subjectNumberStr = num2str(ceil(dataSetsToRun(i) / numOfSess));
    sessNumberStr = num2str(mod(dataSetsToRun(i) - 1, numOfSess) + 1);
    
    disp(['--Subject ', subjectNumberStr,' Session ', sessNumberStr, '''s Component Set']);
    if (conserve_disk_space ~= 1)
        subFile = [back_reconstruction_mat_file, num2str(dataSetsToRun(i)), '.mat'];
        dd = load(fullfile(outputDir, subFile));
        compSet = dd.compSet;
    else
        compSet = getBackReconSet(tmpSessInfo, dataSetsToRun(i));
    end
    
    tempICs = compSet.(fieldName{1});
    tempTCs = compSet.(fieldName{2});
    
    % Since multi-subject BR comps are written as Time by components.
    % Reshape it as components by time
    if (tmpSessInfo.numOfSub*tmpSessInfo.numOfSess > 1)
        tempTCs = tempTCs';
    end
    
    
    tempICs = bsxfun(@minus, tempICs, tmpOffsets(:));
    tempICs = bsxfun(@times, tempICs, tmpComponentSigns(:));
    tempTCs = bsxfun(@times, tempTCs, tmpComponentSigns(:));
    
    
    if (strcmpi(modalityType, 'fmri'))
        
        if (~tmpSessInfo.center_images && (tmpSessInfo.scaleType == 2))
            % Remove mean of component images
            disp('Removing mean of component images ...');
            tempICs = icatb_remove_mean(tempICs')';
            fprintf('Done\n');
        end
    end
    
    
    [tempICs, tempTCs] = icatb_scaleICA(tempICs, tempTCs, tmpSessInfo.inputFiles(i).name, ...
        tmpSessInfo.scaleType, tmpSessInfo.dataType, [], mask_ind);
    
    
    
    if ~strcmpi(modalityType, 'eeg')
        
        if (conserve_disk_space ~= 1)
            %--- save data
            if strcmpi(imExtn, '.nii')
                dispStr = ['...saving scaled ica data for subject ', subjectNumberStr, ' session ', sessNumberStr,...
                    ' in nifti format and as matlab file'];
            else
                dispStr = ['...saving scaled ica data for subject ', subjectNumberStr, ' session ', sessNumberStr,...
                    ' in analyze format and as matlab file'];
            end
        else
            if strcmpi(imExtn, '.nii')
                dispStr = ['...saving scaled ica data for subject ', subjectNumberStr, ' session ', sessNumberStr,...
                    ' in nifti format'];
            else
                dispStr = ['...saving scaled ica data for subject ', subjectNumberStr, ' session ', sessNumberStr,...
                    ' in analyze format'];
            end
            
        end
        
    else
        
        dispStr = ['...saving scaled ica data for subject ', subjectNumberStr, ' session ', sessNumberStr, ...
            ' as matlab file'];
        
    end
    
    disp(dispStr);
    
    if ((conserve_disk_space ~= 1) || strcmpi(modalityType, 'eeg'))
        
        % save data in matlab file
        subFileOut = [calibrate_components_mat_file, subjectNumberStr, '-', sessNumberStr, '.mat'];
        subFileOut = fullfile(outputDir, subFileOut);
        drawnow;
        
        icatb_parSave(subFileOut, {tempICs, tempTCs}, compSetFields);
        
    end
    
    %     out_file = tmpSubjectICAFiles(subjectNumber).ses(sessNumber).name;
    %     out_file = out_file(1, :);
    
    if ~strcmpi(modalityType, 'eeg')
        if ((conserve_disk_space == 1) && tmpSessInfo.center_images)
            saveImgFiles(tmpSessInfo, outputFileNames{i}, tempICs, tempTCs, 'no', 1);
        else
            [zipFileName{i}, files_in_zip{i}] = saveImgFiles(tmpSessInfo, outputFileNames{i}, tempICs, tempTCs, zipImages, 1);
        end
    end
    
end


zipFileName(cellfun('isempty', zipFileName) == 1)  = [];
files_in_zip(cellfun('isempty', files_in_zip) == 1) = [];

sesInfo.zipContents.zipFiles = zipFileName;
for nF = 1:length(zipFileName)
    sesInfo.zipContents.files_in_zip(nF).name = files_in_zip;
end



function [zipFileName, files_in_zip] = saveImgFiles(sesInfo, outfile, tempIC, tempTC, zipFiles, deleteFiles)
%% Save images and compress images if necessary
%

icatb_defaults;
global ZIP_IMAGE_FILES;

if (~exist('zipFiles', 'var') || isempty(zipFiles))
    zipFiles = ZIP_IMAGE_FILES;
end

if (~exist('deleteFiles', 'var'))
    deleteFiles = 1;
end

[fileNames, zipFileName, files_in_zip] = icatb_saveICAData(outfile, tempIC, tempTC, sesInfo.mask_ind, sesInfo.numComp, ...
    sesInfo.HInfo, sesInfo.dataType, [], sesInfo.outputDir, zipFiles, deleteFiles);


function compSet = getBackReconSet(sesInfo, nDataSet)
%% Get back-reconstructed component of the subject
%

sesInfo.dataSetNo = nDataSet;
[sesInfo, compSet] = icatb_backReconstruct(sesInfo);


function outputFileNames = getOutputFileNames(subjectICAFiles, dataSetsToRun, numOfSess)

outputFileNames = cell(1, length(dataSetsToRun));

for i = 1:length(dataSetsToRun)
    subjectNumber = (ceil(dataSetsToRun(i) / numOfSess));
    sessNumber = (mod(dataSetsToRun(i) - 1, numOfSess) + 1);
    outputFileNames{i} = subjectICAFiles(subjectNumber).ses(sessNumber).name(1, :);
end