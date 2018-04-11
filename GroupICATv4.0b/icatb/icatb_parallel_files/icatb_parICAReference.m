function icatb_parICAReference(param_file, dataSetsToRun, algorithmName, statusHandle)
%% Run ica with reference in parallel mode
%

if (ischar(param_file))
    load(param_file);
else
    sesInfo = param_file;
end

appDataName = 'gica_waitbar_app_data';

if (~exist('sesInfo', 'var'))
    error('File selected is not a valid ICA parameter file');
end

ICA_Options = sesInfo.userInput.ICA_Options;

if (~exist('dataSetsToRun', 'var'))
    dataSetsToRun = 1:sesInfo.numOfSub*sesInfo.numOfSess;
end

files = sesInfo.inputFiles;
preproc_type = sesInfo.preproc_type;
mask_ind = sesInfo.mask_ind;
outputDir = sesInfo.outputDir;
numOfSub = sesInfo.numOfSub;
numOfSess = sesInfo.numOfSess;


try
    modalityType = sesInfo.modality;
    setappdata(0, 'group_ica_modality', modalityType);
catch
end

[modalityType, dataTitle, compSetFields] = icatb_get_modality;

if (~exist('statusHandle', 'var'))
    statusHandle = [];
end

%% ICA with reference (Write back-reconstructed files directly)
for nD = dataSetsToRun
    
    disp(['Calculating ICA on data-set ', num2str(nD)]);
    data = icatb_remove_mean(icatb_preproc_data(icatb_read_data(files(nD).name, [], mask_ind), preproc_type));
    
    if (strcmpi(algorithmName, 'constrained ica (spatial)'))
        num_comps = rank(data);
        [tmpW, tmpDw] = icatb_calculate_pca(data, num_comps, 'remove_mean', 0, 'whiten', 1);
        [dd, ddW, tmpA, tmpS] = icatb_icaAlgorithm(algorithmName, tmpW', ICA_Options);
        tmpA = tmpDw*tmpA;
    else
        [dd, ddW, tmpA, tmpS] = icatb_icaAlgorithm(algorithmName, data', ICA_Options);
    end
    
    if (numOfSub*numOfSess == 1)
        % to maintain consistency with the dimensions of timecourses
        % when number of subjects is 1.
        tmpA = tmpA';
    end
    
    compSet = struct(compSetFields{1}, tmpS, compSetFields{2}, tmpA);
    
    % Save subject components
    subFile = [sesInfo.back_reconstruction_mat_file, num2str(nD), '.mat'];
    msgString = ['-saving back reconstructed ica data for set ', num2str(nD),' -> ',subFile];
    disp(msgString);
    drawnow;
    icatb_parSave(fullfile(outputDir, subFile), {compSet}, {'compSet'});
    
    if ~isempty(statusHandle)
        
        % get the status handles
        statusData = getappdata(statusHandle, appDataName);
        statusData.perCompleted = statusData.perCompleted + (statusData.unitPerCompleted / sesInfo.numOfSub*sesInfo.numOfSess);
        setappdata(statusHandle, appDataName, statusData);
        set(statusHandle, 'name', [num2str(round(statusData.perCompleted*100)), '% analysis done']); waitbar(statusData.perCompleted, statusHandle);
        
    end
    
    fprintf('Done\n');
    
end
