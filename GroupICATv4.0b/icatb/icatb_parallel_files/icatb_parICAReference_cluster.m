function icatb_parICAReference_cluster(param_file, dataSetsToRun, algorithmName)
%% Run ica with reference in parallel mode
%

if (ischar(param_file))
    load(param_file);
else
    sesInfo = param_file;
end

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

back_reconstruction_mat_file = sesInfo.back_reconstruction_mat_file;

%% ICA with reference (Write back-reconstructed files directly)
parfor nDataSet = 1:length(dataSetsToRun)
    
    tmpDataSetsToRun = dataSetsToRun;
    nD = tmpDataSetsToRun(nDataSet);
    tmpFiles = files;
    
    tmpCompSetFields = compSetFields;
    
    disp(['Calculating ICA on data-set ', num2str(nD)]);
    data = icatb_remove_mean(icatb_preproc_data(icatb_read_data(tmpFiles(nD).name, [], mask_ind), preproc_type));
    
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
    
    compSet = struct(tmpCompSetFields{1}, tmpS, tmpCompSetFields{2}, tmpA);
    
    % Save subject components
    subFile = [back_reconstruction_mat_file, num2str(nD), '.mat'];
    msgString = ['-saving back reconstructed ica data for set ', num2str(nD),' -> ',subFile];
    disp(msgString);
    drawnow;
    icatb_parSave(fullfile(outputDir, subFile), {compSet}, {'compSet'});
    
    
    fprintf('Done\n');
    
end
