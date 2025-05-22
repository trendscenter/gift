function param_file = icatb_read_batch_file(file_name)
%% Read batch file.
%
% Inputs:
% file_name - Input file
%
% Outputs:
% param_file - ICA parameter file
%

%% Load defaults
icatb_defaults;

global SCALE_DEFAULT;
global PARAMETER_INFO_MAT_FILE;
global WRITE_ANALYSIS_STEPS_IN_DIRS;
global CONSERVE_DISK_SPACE;

scaleType = SCALE_DEFAULT;
write_analysis_steps_in_dirs = WRITE_ANALYSIS_STEPS_IN_DIRS;
conserve_disk_space = CONSERVE_DISK_SPACE;

if (~isstruct(file_name))
    
    pathstr = fileparts(file_name);
    if (isempty(pathstr))
        tempF = file_name;
        file_name = which(file_name);
        if (isempty(file_name))
            error('Error:InpFile', 'File %s doesn''t exist', tempF);
        end
    end
    
    disp('Reading parameters for group ICA ...');
    
    % Evaluate script
    inputData = icatb_eval_script(file_name);
    inputData.inputFile = file_name;
    
else
    inputData = file_name;
    clear file_name;
end

%% Modality type
modalityType = 'fmri';

if (isfield(inputData, 'modalityType'))
    modalityType = inputData.modalityType;
end

% Set application data
setappdata(0, 'group_ica_modality', lower(modalityType));

if (strcmpi(modalityType, 'fnc'))
    param_file = read_fnc(inputData);
    return;
end


%% Check for bids
if (strcmpi(modalityType, 'fmri'))
    if (isfield(inputData, 'bids_info') && ~isempty(inputData.bids_info.root_dir))
        input_data_file_patterns = icatb_parseBIDS(inputData.bids_info);
        inputData.input_data_file_patterns = input_data_file_patterns;
        inputData.dataSelectionMethod = 4;
        sesInfo.userInput.bids_info = inputData.bids_info;
    end
end


%% Group ICA type
if (strcmpi(modalityType, 'fmri'))
    group_ica_type = 'spatial';
    try
        group_ica_type = inputData.group_ica_type;
    catch
    end
    sesInfo.userInput.group_ica_type = group_ica_type;
end

%% Parallel info
analysisMode = 'serial';
num_workers = 4;
try
    analysisMode = lower(inputData.parallel_info.mode);
catch
end

try
    num_workers = inputData.parallel_info.num_workers;
catch
end

sesInfo.userInput.parallel_info.mode = analysisMode;
sesInfo.userInput.parallel_info.num_workers = num_workers;

%% Conserve disk space
try
    conserve_disk_space = inputData.conserve_disk_space;
catch
end

%% Analysis Type
which_analysis = 1;
if (isfield(inputData, 'which_analysis'))
    which_analysis = inputData.which_analysis;
end

sesInfo.userInput.which_analysis = which_analysis;

% ICASSO Options
if (which_analysis == 2)
    
    icasso_opts = inputData.icasso_opts;
    sel_mode = lower(icasso_opts.sel_mode);
    num_ica_runs = icasso_opts.num_ica_runs;
    
    if (~strcmpi(sel_mode, 'randinit') && ~strcmpi(sel_mode, 'bootstrap') && ~strcmpi(sel_mode, 'both'))
        error('Please provide a valid option for icasso_opts.sel_mode. Valid option must be in {randinit, bootstrap, both}');
    end
    
    if (num_ica_runs < 2)
        error('Error:ICASSO', ['You need to run ICA algorithm atleast two times inorder to use ICASSO.', ...
            '\nPlease check variable icasso_opts.num_ica_runs.']);
    end
    
    sesInfo.userInput.icasso_opts = inputData.icasso_opts;
    
end

% MST options
if (which_analysis == 3)
    sesInfo.userInput.mst_opts = inputData.mst_opts;
end

% Cross ISI options
if (which_analysis == 4)
    sesInfo.userInput.cross_isi_opts = inputData.cross_isi_opts;
end

%% Group PCA performance settings
perfOptions = icatb_get_analysis_settings;
perfType = 'user specified settings';
if (isfield(inputData, 'perfType'))
    perfType = inputData.perfType;
end

if (isnumeric(perfType))
    perfType = perfOptions{perfType};
end

perfType = lower(perfType);
sesInfo.userInput.perfType = perfType;

%% Output prefix
prefix = '';
if (isfield(inputData, 'prefix'))
    prefix = inputData.prefix;
end

[status, message] = icatb_errorCheck(prefix, 'output_prefix', 'prefix');
if (~status)
    error(message);
end

sesInfo.userInput.prefix = prefix;


%% Pre-processing type
preproc_options = icatb_preproc_data;
if ~(strcmpi(modalityType, 'smri'))
    preproc_type = 'remove mean per timepoint';
    alias_options = {'rt', 'rv', 'in', 'vn'};
else
    preproc_type = 'remove mean per subject';
    alias_options = {'rt', 'none'};
end

%if (~strcmpi(modalityType, 'smri'))
if (isfield(inputData, 'preproc_type'))
    preproc_type = inputData.preproc_type;
end
%end

try
    ind = getIndex(preproc_type, preproc_options, 'Data Pre-processing Type');
catch
    ind = getIndex(preproc_type, alias_options, 'Data Pre-processing Type');
end

if (isempty(ind))
    ind = 1;
end

sesInfo.userInput.preproc_type = lower(preproc_options{ind});

%% Output directory
outputDir = pwd;
if (isfield(inputData, 'outputDir'))
    outputDir = inputData.outputDir;
end

if (~exist(outputDir, 'dir'))
    sLastChar =  outputDir(length(outputDir));
    if (sLastChar == '/' || sLastChar == '\') %slash bugs out the folder creation
        [ParDir, ChildDir] = fileparts(outputDir(1:end-1)); %fix
    else
        [ParDir, ChildDir] = fileparts(outputDir);
    end
    mkdir(ParDir, ChildDir);
end

cd(outputDir);
outputDir = pwd;

sesInfo.userInput.pwd = outputDir;

% store these fields regarding dataType, complex naming
dataType = 'real'; read_complex_images = 'real&imaginary'; write_complex_images = 'real&imaginary';

sesInfo.userInput.param_file = [sesInfo.userInput.prefix, PARAMETER_INFO_MAT_FILE, '.mat'];
sesInfo.userInput.param_file = fullfile(outputDir, sesInfo.userInput.param_file);
param_file = sesInfo.userInput.param_file;
sesInfo.userInput.dataType = lower(dataType);
sesInfo.userInput.read_complex_images = lower(read_complex_images);
sesInfo.userInput.write_complex_images = lower(write_complex_images);

%% Read data
% Support to handle motion regression in batch as well as in the GUI
if ~isfield(inputData, 'input_data_hd_patterns')
    inputData.input_data_hd_patterns = '';
end

[sesInfo] = icatb_name_complex_images(sesInfo, 'read');

if (strcmpi(modalityType, 'conn'))
    inputData.input_data_file_patterns = cellstr(inputData.input_data_file_patterns);
    sessionsInvolved = size(inputData.input_data_file_patterns, 2);
    filePatterns = inputData.input_data_file_patterns;
    filePatterns = filePatterns';
    filePatterns = filePatterns(:);
    sesInfo.userInput.dataInfo.filesInfo.filesList = filePatterns;
    sesInfo.userInput.dataInfo.filesInfo.numOfSess = sessionsInvolved;
    
    dummy_scans = 0;
    try
        dummy_scans = inputData.dummy_scans;
    catch
    end
    sesInfo.userInput.dataInfo.filesInfo.file_numbers = dummy_scans;
    
    subsampling_depth = 1;
    try
        subsampling_depth = inputData.subsampling_depth;
    catch
    end
    sesInfo.userInput.dataInfo.subsampling_depth = subsampling_depth;
    maskFile = [];
    try
        maskFile = inputData.maskFile;
    catch
    end
    sesInfo.userInput.dataInfo.maskFile = maskFile;
    conn_type = 'ENLwFC';
    try
        conn_type = inputData.conn_type;
    catch
    end
    
    sesInfo.userInput.dataInfo.conn_type = conn_type;
    sesInfo = icatb_gen_data_conn_ica(sesInfo);
    inputData.maskFile = sesInfo.userInput.maskFile;
    files = sesInfo.userInput.files;
    numOfSub = sesInfo.userInput.numOfSub;
    numOfSess = sesInfo.userInput.numOfSess;
    numofDatasets = numOfSub*numOfSess;
    SPMFiles = [];
    diffTimePoints =  repmat(length(sesInfo.userInput.mask_ind), 1, numofDatasets);
    sesInfo.userInput.diffTimePoints = diffTimePoints;
else
    
    [files, designMatrix, numOfSub, numOfSess, dataSelMethod, diffTimePoints, spmMatFlag] = icatb_dataSelection(...
        inputData, sesInfo.userInput.pwd, sesInfo.userInput.prefix, ...
        sesInfo.userInput.read_complex_file_naming, sesInfo.userInput.read_complex_images);
    sesInfo.userInput.files = files;
    sesInfo.userInput.hd    = inputData.input_data_hd_patterns; %support for motion reg
    SPMFiles = designMatrix;
    sesInfo.userInput.dataSelMethod = dataSelMethod;
    sesInfo.userInput.designMatrix = designMatrix;
    sesInfo.userInput.spmMatFlag = spmMatFlag;
    sesInfo.userInput.diffTimePoints = diffTimePoints;
    sesInfo.userInput.numOfSub = numOfSub;
    sesInfo.userInput.numOfSess = numOfSess;
    
end

drawnow;

subjectFile = [sesInfo.userInput.prefix, 'Subject.mat'];
subjectFile = fullfile(outputDir, subjectFile);

icatb_save(subjectFile, 'files', 'numOfSub', 'numOfSess', 'SPMFiles', 'modalityType');



%% Mask file
maskFile = [];
if (isfield(inputData, 'maskFile'))
    maskFile = inputData.maskFile;
end

if (~strcmpi(modalityType, 'conn'))
    
    if (~isempty(maskFile) && (strcmpi(modalityType, 'fmri') || strcmpi(modalityType, 'smri')))
        
        if (strcmpi(maskFile, 'default') || strcmpi(maskFile, 'default mask'))
            maskFile = [];
        elseif (strcmpi(lower(maskFile), 'default&icv'))
            maskFile = ['default&icv'];
        else
            if (strcmpi(maskFile, 'average') || strcmpi(maskFile, 'average mask'))
                
                inpFiles = cell(1, length(files));
                for nF = 1:length(files)
                    inpFiles{nF} = files(nF).name;
                end
                
                icatb_generateMask(inpFiles, 'outputDir', outputDir, 'corr_threshold', 0.8, 'prefix', sesInfo.userInput.prefix);
                
                maskFile = fullfile(outputDir, [sesInfo.userInput.prefix, 'Mask.nii']);
                
            end
            
            [maskPath, maskF, extn] = fileparts(maskFile);
            if (isempty(maskPath))
                maskPath = pwd;
            end
            maskFile = fullfile(maskPath, [maskF, extn]);
            if (~exist(icatb_parseExtn(maskFile), 'file'))
                error([maskFile, ' doesn''t exist']);
            end
            
            % If mask resolution do not match the mask will be resliced
            [dT, extns, maskDim] = icatb_get_countTimePoints(icatb_parseExtn(maskFile));
            [dT, extns, dims] = icatb_get_countTimePoints(icatb_parseExtn(deblank(sesInfo.userInput.files(1).name(1, :))));
            if length(find(maskDim == dims)) ~= length(maskDim)
                fprintf('Mask dimensions ([%s]) are not equal to image dimensions ([%s]). Resizing mask image/images to match functional image\n\n', ...
                    num2str(maskDim), num2str(dims));
                
                firstFile = deblank(sesInfo.userInput.files(1).name(1, :));
                
                % Handle gz files
                firstFileTmp = deblank(icatb_parseExtn(firstFile));
                if (strcmpi(firstFileTmp(end-2:end), '.gz'))
                    gzfn = gunzip (firstFileTmp, tempdir);
                    gzfn = char(gzfn);
                    firstFile = icatb_rename_4d_file(gzfn);
                    firstFile = deblank(firstFile(1, :));
                end
                
                sTemp = noisecloud_spm_coregister(firstFile, deblank(maskFile(1, :)), maskFile, sesInfo.userInput.pwd);
                [sFPath,sFName,sFExt] = fileparts(sTemp);
                % Rename tmp output to the mask name
                maskFile = [sesInfo.userInput.pwd filesep sesInfo.userInput.prefix 'Mask' sFExt];
                movefile(sTemp, maskFile, 'f')
            end
            
        end
        
    else
        maskFile = [];
    end
    
    sesInfo.userInput.maskFile = maskFile;
    %% Create mask
    sesInfo = icatb_update_mask(sesInfo);
end



%% Algorithm
algoType = 1;

if (isfield(inputData, 'algoType'))
    algoType = inputData.algoType;
end

ica_algo = lower(cellstr(icatb_icaAlgorithm));

if (ischar(algoType))
    if strcmpi(algoType, 'gig-ica')
        algoType = 'moo-icar';
    end
end

algoType = getIndex(algoType, ica_algo, 'ICA Algorithm');
sesInfo.userInput.algorithm = algoType;

if strcmpi(ica_algo{algoType}, 'semi-blind infomax')
    %% Error check for SBICA
    if (sesInfo.userInput.numOfSub*sesInfo.userInput.numOfSess > 1)
        error('Error:SBICA', ['Presently SBICA works with only one data set.', ...
            '\nPlease select another algorithm if you want to run multiple data-sets']);
    end
end

useTemporalICA = 0;
if (strcmpi(modalityType, 'fmri'))
    useTemporalICA = strcmpi(sesInfo.userInput.group_ica_type, 'temporal');
end

%% Data Reduction steps
numReductionSteps = 2;
if (isfield(inputData, 'numReductionSteps'))
    numReductionSteps = inputData.numReductionSteps;
end

if (numOfSub*numOfSess == 1)
    numReductionSteps = 1;
end


if (useTemporalICA || ~isempty(icatb_findstr(lower(ica_algo{algoType}),'iva')) ...
        || strcmpi(ica_algo{algoType}, 'moo-icar') || icatb_string_compare(ica_algo{algoType}, 'constrained'))
    numReductionSteps = 1;
end

if (numReductionSteps > 2)
    disp('!!!You have specified more than 2 reduction steps. Number of reduction steps is set to 2.');
    numReductionSteps = 2;
end

sesInfo.userInput.numReductionSteps = numReductionSteps;

%% Do Estimation
doEstimation = 0;
if (isfield(inputData, 'doEstimation'))
    doEstimation = inputData.doEstimation;
end

if (strcmpi(ica_algo{algoType}, 'moo-icar') || icatb_string_compare(ica_algo{algoType}, 'constrained'))
    doEstimation = 0;
end

if (doEstimation > 0)
    
    icatb_save(param_file, 'sesInfo');
    
    sesInfo.userInput.dim_est_opts.method = doEstimation;
    
    estimation_opts = struct;
    
    if (isfield(inputData, 'estimation_opts'))
        estimation_opts = inputData.estimation_opts;
    end
    
    if (~isfield(estimation_opts, 'PC1'))
        estimation_opts.PC1 = 'mean';
    end
    
    if (~isfield(estimation_opts, 'PC2'))
        estimation_opts.PC2 = 'mean';
    end
    
    try
        fwhm = estimation_opts.fwhm;
        sesInfo.userInput.dim_est_opts.fwhm = fwhm;
    catch
    end
    
    disp('Doing dimensionality estimation ...');
    
    %% Do estimation
    [estimated_comps, c, m, a, sesInfo] = icatb_estimateCompCallback([], [], sesInfo);
    
    clear c m a;
    
    if any(estimated_comps == 1)
        error('Only one component is estimated for one or more subjects. Try setting components manually.');
    end
    
    %% Evalulate PC1
    numPC1 = round(eval([lower(estimation_opts.PC1), '(estimated_comps);']));
    
    %% Evalulate PC2
    if (numReductionSteps == 2)
        numPC2 = round(eval([lower(estimation_opts.PC2), '(estimated_comps);']));
    end
    
    % Set the appropriate data reduction numbers
    inputData.numOfPC1 = numPC1;
    
    if (numReductionSteps == 2)
        inputData.numOfPC2 = numPC2;
    end
    
    disp('Done with the dimesionality estimation ...');
    
    fprintf('\n');
    
end

%% Principal component numbers
numOfDataSets = numOfSub*numOfSess;
try
    sesInfo.userInput.numOfPC1 = inputData.numOfPC1;
catch
    sesInfo.userInput.numOfPC1 = min(diffTimePoints);
end

sesInfo.userInput.numComp = sesInfo.userInput.numOfPC1;
sesInfo.userInput.numOfPC2 = 0;
sesInfo.userInput.numOfPC3 = 0;

sesInfo.userInput.numOfGroups1 = (numOfSub*numOfSess);
sesInfo.userInput.numOfGroups2 = 0;
sesInfo.userInput.numOfGroups3 = 0;

if (numReductionSteps == 2)
    sesInfo.userInput.numOfPC2 = inputData.numOfPC2;
    sesInfo.userInput.numComp = inputData.numOfPC2;
    sesInfo.userInput.numOfGroups2 = 1;
end

if (sesInfo.userInput.numComp < 2)
    error('Select Number of IC to be more than or equal to 2.');
end

%% TR
if (strcmpi(modalityType, 'fmri'))
    if (isfield(inputData, 'TR'))
        TR = inputData.TR;
        if (length(TR) > 1)
            if (length(TR) ~= numOfSub)
                error('TR length must match the number of subjects');
            end
        end
        sesInfo.userInput.TR = TR;
    end
end


%% PCA Type
pcaType = 'standard';
pcaOptions = icatb_pca_options;
if (isfield(inputData, 'pcaType'))
    pcaType = inputData.pcaType;
end
ind = getIndex(pcaType, pcaOptions, 'PCA Type');
pcaType = pcaOptions{ind};
sesInfo.userInput.pcaType = lower(pcaType);

pca_opts = icatb_pca_options(pcaType);

if (isfield(inputData, 'pca_opts'))
    pca_opts = inputData.pca_opts;
else
    if (isfield(inputData, 'covariance_opts'))
        pca_opts = inputData.covariance_opts;
    end
end

sesInfo.userInput.pca_opts = pca_opts;

%% Group PCA Type
groupPCAOpts = {'Subject Specific', 'Grand Mean'};
group_pca_type = 'subject specific';
if (isfield(inputData, 'group_pca_type'))
    group_pca_type = inputData.group_pca_type;
end
ind = getIndex(group_pca_type, groupPCAOpts, 'Group PCA Type');
sesInfo.userInput.group_pca_type = lower(groupPCAOpts{ind});

if (strcmpi(sesInfo.userInput.group_pca_type, 'grand mean') && any(sesInfo.userInput.diffTimePoints ~= sesInfo.userInput.diffTimePoints(1)))
    if (strcmpi(icatb_get_modality, 'fmri'))
        error('Please select the same no. of timepoints if you want to select the grand mean group PCA');
    elseif(strcmpi(icatb_get_modality, 'eeg'))
        error('Please select the same no. of electrodes if you want to select the grand mean group PCA');
    end
end

%% Back reconstruction type
backReconType = 'regular';
backReconOptions = icatb_backReconOptions;
backReconOptions = cellstr(backReconOptions);
if (isfield(inputData, 'backReconType'))
    backReconType = inputData.backReconType;
end

if (ischar(backReconType) && strcmpi(backReconType, 'str'))
    backReconType = 'spatial-temporal regression';
end

if (ischar(backReconType) && strcmpi(backReconType, 'moo-icar'))
    backReconType = 'gig-ica';
end

ind = getIndex(backReconType, backReconOptions, 'Back Reconstruction');

sesInfo.userInput.backReconType = lower(backReconOptions{ind});


%% Scale Type
scaleOptions = icatb_scaleICA;
if (isfield(inputData, 'scaleType'))
    scaleType = inputData.scaleType;
end

if (isnumeric(scaleType))
    scaleType = scaleType + 1;
else
    if (strcmpi(scaleType, 'no'))
        scaleType = 'no scaling';
    elseif (findstr(scaleType, 'norm ics'))
        if (strcmpi(modalityType, 'fmri'))
            scaleType = 'scaling in timecourses';
        else
            scaleType = 'scaling in topographies';
        end
    elseif (findstr(scaleType, 'joint scaling'))
        if (strcmpi(modalityType, 'fmri'))
            scaleType = 'scaling in maps and timecourses';
        else
            scaleType = 'scaling in timecourses and topographies';
        end
    elseif (strcmpi(scaleType, 'calibrate'))
        if (strcmpi(modalityType, 'fmri'))
            scaleType = 'scale to original data(%)';
        else
            scaleType = 'scale to original data';
        end
    end
end

ind = getIndex(scaleType, scaleOptions, 'Scaling Type');
sesInfo.userInput.scaleType = (ind - 1);


%% ICA Algorithm callback

% Size of the reduced data
dataSize = [sesInfo.userInput.numComp, length(sesInfo.userInput.mask_ind)];

% ICA Options
ICA_Options = icatb_icaOptions(dataSize, lower(ica_algo{algoType}), 'off');

ICA_Options = chkICAOptions(ICA_Options, inputData);

sesInfo.userInput.ICA_Options = ICA_Options;

if strcmpi(ica_algo{algoType}, 'semi-blind infomax')
    % SBICA
    sesInfo = sbICACallback(sesInfo, inputData.refFunNames);
    
elseif (icatb_string_compare(ica_algo{algoType}, 'constrained') || strcmpi(ica_algo{algoType}, 'moo-icar') || ...
        strcmpi(ica_algo{algoType}, 'iva-l-sos-adaptive'))
    % Constrained ICA (Spatial)
    sesInfo = constrainedICACallback(sesInfo, inputData.refFiles);
end

if (isempty(write_analysis_steps_in_dirs))
    write_analysis_steps_in_dirs = 0;
end

if (isempty(conserve_disk_space))
    conserve_disk_space = 0;
end

sesInfo.userInput.write_analysis_steps_in_dirs = write_analysis_steps_in_dirs;
sesInfo.userInput.conserve_disk_space = conserve_disk_space;

sesInfo.isInitialized = 0;
display_results = 0;
try
    display_results = inputData.display_results;
catch
end

try
    network_summary_opts = inputData.network_summary_opts;
catch
end


if (isnumeric(display_results))
    tmp_display_results = display_results;
    clear display_results;
    if (tmp_display_results ~= 0)
        display_results.formatName = tmp_display_results;
    end
end

if (exist('network_summary_opts', 'var'))
    display_results.network_summary_opts = network_summary_opts;
end

if (exist('display_results', 'var'))
    sesInfo.userInput.display_results = display_results;
end

%% Save Parameter file
icatb_save(param_file, 'sesInfo');
disp('');
displayString = [' Parameters are saved in ', param_file];
% display the parameters
disp(displayString);
fprintf('\n');

%% Create global variable also
global GIFT_PARAM_FILE;
GIFT_PARAM_FILE = param_file;

function selOption = getIndex(selOption, options, titleStr)
%% Get selected index
%

if (ischar(selOption))
    ind = strmatch(lower(selOption), lower(cellstr(options)), 'exact');
    if (isempty(ind))
        error(['Unknown ', titleStr, '(', selOption, ') passed']);
    end
    selOption = ind(1);
end

function sesInfo = sbICACallback(sesInfo, refFunNames)
%% SBICA Callback
%

% check if the design matrix is specified for semi-blind ICA
% check whether the index is passed or not,
if isfield(sesInfo.userInput, 'designMatrix')
    % open the design matrix file
    checkFile = sesInfo.userInput.designMatrix.name;
    if isempty(checkFile)
        error('Design matrix should be specified to run semi-blind ICA');
    end
    
    % do error checking for SPM design matrix
    [spmData] = icatb_loadSPM_new('spmName', checkFile, 'countTimePoints', ...
        sesInfo.userInput.diffTimePoints(1), 'data_sessionNumber', 1, ...
        'check_spm_design_matrix', 'yes', 'check_regressors', 'yes');
    allNames = spmData.sel_refnames;
    allNames = deblank(allNames);
    %refFunNames = inputData.refFunNames;
    if length(refFunNames) > 2
        disp('-- By default using only two reference function names');
        refFunNames = refFunNames(1:2);
    end
    
    getIndex = zeros(1, length(refFunNames));
    % check the index
    for ii = 1:length(refFunNames)
        tempIndex = strmatch(lower(refFunNames{ii}), lower(allNames), 'exact');
        if isempty(tempIndex)
            error(['Specified reference function with name ', lower(refFunNames{ii}), ' doesn''t exist']);
        end
        % store the indices
        getIndex(ii) = tempIndex;
    end
    
    refFunNames = char(refFunNames);
    
    spmData.selectedRegressors = refFunNames; % store the selected regressors
    spmData.getIndex = getIndex; % store the selected information in spm data
    
    % get the selected time course information
    spmData = icatb_loadSPM_new('spmData', spmData, 'get_timecourses', 'yes', ...
        'check_spm_design_matrix', 'no', 'check_regressors', 'no', 'flag_selecting_regressors', 'no');
    
    timecourse = spmData.timecourse; % get the selected time courses information
    
    clear spmData;
    
    sesInfo.userInput.ICA_Options = [sesInfo.userInput.ICA_Options, {'TC', timecourse}];
    
else
    error('The field designMatrix doesn''t exist for sesInfo');
end
% end for checking design matrix


function sesInfo = constrainedICACallback(sesInfo, spatial_references)
%% Constrained ICA Callback
%

spatial_references = char(spatial_references);

if (size(spatial_references, 1) == 1)
    [pathstr, fp, extn] = fileparts(spatial_references);
    fileContents = icatb_listFiles_inDir(pathstr, [fp, extn]);
    if (isempty(fileContents))
        error('Error:FilePattern', 'Please check file pattern %s as there are no files found\n', spatial_references);
    end
    spatial_references = icatb_fullFile('directory', pathstr, 'files', fileContents);
end

% get the spatial reference data
[images, imHInfo] = icatb_loadData(spatial_references);
numSpatialFiles = size(images, 4);

imDims = imHInfo(1).DIM(1:3);
funcDims = sesInfo.userInput.HInfo.DIM(1:3);

if (length(find((imDims == funcDims) ~= 0)) ~= length(funcDims))
    fprintf('Spatial reference image dimensions ([%s]) are not equal to functional image dimensions ([%s]). Resizing template image/images to match functional image\n\n', ...
        num2str(imDims), num2str(funcDims));
    
    firstFile = deblank(sesInfo.userInput.files(1).name(1, :));
    
    firstFileTmp = deblank(icatb_parseExtn(firstFile));
    if (strcmpi(firstFileTmp(end-2:end), '.gz'))
        gzfn = gunzip (firstFileTmp, tempdir);
        gzfn = char(gzfn);
        firstFile = icatb_rename_4d_file(gzfn);
        firstFile = deblank(firstFile(1, :));
    end
    
    spatial_references = noisecloud_spm_coregister(firstFile, deblank(spatial_references(1, :)), spatial_references, sesInfo.userInput.pwd);
end

%[images, imHInfo] = icatb_loadData(spatial_references);


% if (length(find((imDims == funcDims) ~= 0)) ~= length(funcDims))
%     error('Error:Dimensions', 'Spatial reference image dimensions ([%s]) are not equal to functional image dimensions ([%s])', ...
%         num2str(imDims), num2str(funcDims));
% end

%images = reshape(images, prod(funcDims), numSpatialFiles);
%images = (images(sesInfo.userInput.mask_ind, :))';

ICAOptions = sesInfo.userInput.ICA_Options;

% Update ICA Options
sesInfo.userInput.ICA_Options = [{'ref_data', {spatial_references, sesInfo.userInput.mask_ind}}, ICAOptions];
sesInfo.userInput.numComp = numSpatialFiles;
sesInfo.userInput.numOfPC1 = sesInfo.userInput.numComp;

function ICA_Options = chkICAOptions(ICA_Options, inputData)
%% Check ICA Options
%

if (~isempty(ICA_Options))
    if (isfield(inputData, 'icaOptions'))
        tmp = inputData.icaOptions;
        [matched_fields, ia, ib] = intersect(ICA_Options(1:2:end), tmp(1:2:end));
        if (~isempty(matched_fields))
            ICA_Options(2*ia) = tmp(2*ib);
        end
    end
end


function param_file = read_fnc(inputData)
%% Read FNC file
%

contrast_vector = [];
try
    contrast_vector = inputData.contrast_vector;
catch
end

if (~isnumeric(contrast_vector))
    contrast_vector = str2num(contrast_vector);
end

dataInfo = [];

if ~iscell(inputData.input_data_file_patterns)
    inputData.input_data_file_patterns = cellstr(inputData.input_data_file_patterns);
end

count = 0;
for nSess = 1:size(inputData.input_data_file_patterns, 2)
    count = count + 1;
    dataInfo(count).name = ['Session ', num2str(nSess)];
    dataInfo(count).files = inputData.input_data_file_patterns(:, nSess);
end

fnc_matrix = icatb_fnc_input(inputData);

sesInfo.isInitialized = 0;

sesInfo.userInput.prefix = inputData.prefix;
sesInfo.userInput.modality = inputData.modalityType;
sesInfo.userInput.outputDir = inputData.outputDir;
sesInfo.userInput.pwd = inputData.outputDir;
sesInfo.userInput.dataInfo = dataInfo;

scaleType = 'no';
try
    scaleType  = inputData.scaleType;
catch
end

sesInfo.userInput.scaleType = scaleType;

cd(sesInfo.userInput.outputDir);

doEstimation = 0;
try
    doEstimation = inputData.doEstimation;
catch
end

sesInfo.userInput.doEstimation = doEstimation;

ICA_Options = {};
try
    ICA_Options = sesInfo.userInput.ICA_Options;
catch
end

sesInfo.userInput.ICA_Options = ICA_Options;

fnc_file_name = fullfile(sesInfo.userInput.outputDir, [sesInfo.userInput.prefix, '_fnc_data.mat']);

save(fnc_file_name, 'fnc_matrix');

fnc_matrix = icatb_mat2vec(fnc_matrix);

if (doEstimation == 1)
    [numOfPC1, mdl] = order_selection(fnc_matrix);
    sesInfo.userInput.est_opts.mdl = mdl;
    sesInfo.userInput.est_opts.numOfPC1 = numOfPC1;
else
    numOfPC1 = inputData.numOfPC1;
end


which_analysis = 1;
try
    which_analysis = inputData.which_analysis;
catch
end

sesInfo.userInput.which_analysis = which_analysis;

% ICASSO Options
if (which_analysis == 2)
    
    icasso_opts = inputData.icasso_opts;
    sel_mode = lower(icasso_opts.sel_mode);
    num_ica_runs = icasso_opts.num_ica_runs;
    
    if (~strcmpi(sel_mode, 'randinit') && ~strcmpi(sel_mode, 'bootstrap') && ~strcmpi(sel_mode, 'both'))
        error('Please provide a valid option for icasso_opts.sel_mode. Valid option must be in {randinit, bootstrap, both}');
    end
    
    if (num_ica_runs < 2)
        error('Error:ICASSO', ['You need to run ICA algorithm atleast two times inorder to use ICASSO.', ...
            '\nPlease check variable icasso_opts.num_ica_runs.']);
    end
    
    sesInfo.userInput.icasso_opts = inputData.icasso_opts;
    
end

% MST options
if (which_analysis == 3)
    sesInfo.userInput.mst_opts = inputData.mst_opts;
end

% Cross ISI options
if (which_analysis == 4)
    sesInfo.userInput.cross_isi_opts = inputData.cross_isi_opts;
end

sesInfo.userInput.doEstimation = doEstimation;
sesInfo.userInput.numComp = numOfPC1;
sesInfo.userInput.numOfPC1 = numOfPC1;
sesInfo.userInput.contrast_vector = contrast_vector;


%% Algorithm
algoType = 1;

if (isfield(inputData, 'algoType'))
    algoType = inputData.algoType;
end

ica_algo = lower(cellstr(icatb_icaAlgorithm));

algoType = getIndex(algoType, ica_algo, 'ICA Algorithm');
sesInfo.userInput.algorithm = algoType;

if isfield(inputData, 'reference_file')
    reference_file = inputData.reference_file;
    if ~ischar(reference_file)
        labels_file_name = reference_file{2};
        reference_file_name = reference_file{1};
    else
        reference_file_name = reference_file;
        if ~(exist(reference_file_name, 'file'))
            error([reference_file_name, ' doesn''t exist']);
        end
        [pathstr, fn, extn] = fileparts(reference_file_name);
        labels_file_name = fullfile(pathstr, [fn, '.txt']);
        if ~(exist(labels_file_name, 'file'))
            error([labels_file_name, ' doesn''t exist']);
        end
    end
    
    sesInfo.userInput.reference_file = reference_file_name;
    compLabels = getCompLabels(labels_file_name);
    networkOpts = cell(length(compLabels), 2);
    for n = 1:length(compLabels)
        networkOpts{n, 1} = compLabels(n).name;
        networkOpts{n, 2} = compLabels(n).value;
    end
    sesInfo.userInput.network_summary_opts = networkOpts;
end

display_results = 0;
try
    display_results = inputData.display_results;
catch
end
sesInfo.userInput.display_results = display_results;

param_file = fullfile(sesInfo.userInput.outputDir, [sesInfo.userInput.prefix, '_ica_parameter_info.mat']);
save(param_file, 'sesInfo');

disp('');
displayString = [' Parameters are saved in ', param_file];
% display the parameters
disp(displayString);
fprintf('\n');


function [comp_est, mdl, aic, kic] = order_selection(data, fwhm)

disp('Estimating dimension ...');

%% Remove mean
data = detrend(data, 0);

%% Arrange data based on correlation threshold
[V1, D1] = icatb_svd(data);

lam = diag(D1);

lam = sort(lam);

lam = lam(end:-1:1);

lam = lam(:)';

N = (size(data, 2));

if (exist('fwhm', 'var') && ~isempty(fwhm))
    N = N  / prod(fwhm);
    N = ceil(N);
end

%% Make eigen spectrum adjustment
tol = max(size(lam)) * eps(max(lam));

if (lam(end) < tol)
    lam(end) = [];
end

%% Correction on the ill-conditioned results (when tdim is large, some
% least significant eigenvalues become small negative numbers)
lam(real(lam) <= tol) = tol;

p = length(lam);
aic = zeros(1, p - 1);
kic = zeros(1, p - 1);
mdl = zeros(1, p - 1);
for k = 1:p-1
    LH = log(prod(lam(k+1:end).^(1/(p-k)) )/mean(lam(k+1:end)));
    mlh = 0.5*N*(p-k)*LH;
    df = 1 + 0.5*k*(2*p-k+1);
    aic(k) =  -2*mlh + 2*df;
    kic(k) =  -2*mlh + 3*df;
    mdl(k) =  -mlh + 0.5*df*log(N);
end

% Find the first local minimum of each ITC
itc = zeros(3, length(mdl));
itc(1,:) = aic;
itc(2,:) = kic;
itc(3,:) = mdl;

%% Use only mdl
dlap = squeeze(itc(end, 2:end)-itc(end, 1:end-1));
a = find(dlap > 0);
if isempty(a)
    comp_est = length(squeeze(itc(end, :)));
else
    comp_est = a(1);
end

disp(['Estimated components is found to be ', num2str(comp_est)]);


if ((comp_est < 2) || all(diff(mdl) < 0))
    comp_est = min([6, min(size(data))]);
    warning('FNC:DimensionalityEstimation', 'Estimated components is 1 or MDL function is monotonically decreasing. Using %d instead', comp_est);
end

disp('Done');


function compLabels = getCompLabels(txtFile)

fid = fopen(txtFile, 'r');
if (fid == -1)
    error(['File ', txtFile, ' cannot be opened for reading']);
end
try
    dd = textscan(fid, '%s', 'delimiter', '\t\n,', 'multipleDelimsAsOne', 1, 'whitespace', ' ');
    val = dd{1};
catch
    val = [];
end
fclose(fid);
val = val(icatb_good_cells(val));
chk = cellfun('isempty', regexp(val, '^\d+$'));

inds = find(chk == 1);

compLabels = repmat(struct('name', '', 'value', []), 1, length(inds));
for nI = 1:length(inds)
    compLabels(nI).name = val{inds(nI)};
    if (nI == length(inds))
        endT = length(val);
    else
        endT = inds(nI + 1) - 1;
    end
    dd = str2num(char(val{inds(nI) + 1:endT}));
    compLabels(nI).value = dd(:)';
end