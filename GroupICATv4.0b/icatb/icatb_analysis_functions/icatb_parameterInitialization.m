function [sesInfo] = icatb_parameterInitialization(sesInfo, statusHandle)
%% Initializes the variables from the user input and perfoms error checking
% Inputs: sesInfo - structure containing all parameters for analysis

if ~exist('sesInfo','var')
    [P] = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Parameter File', 'filter', '*param*.mat');
    if isempty(P)
        error('Parameter file is not selected for analysis');
    end
    [pathstr, fileName] = fileparts(P);
    outputDir = pathstr;
    % Make sure parameter file exists
    load(P);
    if ~exist('sesInfo', 'var')
        %         infoCell{1} = P;
        %         icatb_error('The selected file does not contain sesInfo variable', infoCell);
        error(['The selected file ', P, ' does not contain the sesInfo variable']);
    end
else
    outputDir = sesInfo.outputDir; % get the output directory information
end


sesInfo.outputDir = outputDir;

if ~exist('statusHandle', 'var')
    statusHandle = [];
end

[modalityType, dataTitle] = icatb_get_modality;

% Load defaults
icatb_defaults;

global PARAMETER_INFO_MAT_FILE; % parameter file
global DATA_REDUCTION_MAT_FILE; %A file for each subject's principal components
global ICA_MAT_FILE;
global BACK_RECONSTRUCTION_MAT_FILE;
global CALIBRATE_MAT_FILE;
global AGGREGATE_AN3_FILE;
global FUNCTIONAL_DATA_FILTER;
global FLIP_ANALYZE_IMAGES;
global NUM_RUNS_GICA;

[pp, bb, functional_data_extn] = fileparts(FUNCTIONAL_DATA_FILTER);

if ~strcmpi(modalityType, 'eeg')
    if ~strcmpi(functional_data_extn, '.nii') & ~strcmpi(functional_data_extn, '.img')
        error(['Functional data filter used for writing images is neither nifti nor analyze format.' ...
            ' Check FUNCTIONAL_DATA_FILTER variable in icatb_defaults.m']);
    end
end


sesInfo.flip_analyze_images = FLIP_ANALYZE_IMAGES;


%% Parallel info
num_workers = 4;
parallelMode = 'serial';

try
    parallelMode = sesInfo.userInput.parallel_info.mode;
catch
end

try
    num_workers = sesInfo.userInput.parallel_info.num_workers;
catch
end

sesInfo.parallel_info.mode = parallelMode;
sesInfo.parallel_info.num_workers = num_workers;



icaStr = icatb_icaAlgorithm;
algorithmName = deblank(icaStr(sesInfo.userInput.algorithm, :));

if strcmpi(algorithmName, 'moo-icar')
    algorithmName = 'gig-ica';
end

if (strcmpi(modalityType, 'fmri'))
    
    try
        sesInfo.group_ica_type = sesInfo.userInput.group_ica_type;
    catch
        sesInfo.group_ica_type = 'spatial';
    end
    
    useTemporalICA = strcmpi(sesInfo.group_ica_type, 'temporal');
    if (useTemporalICA)
        if (strcmpi(algorithmName, 'iva-gl') || strcmpi(algorithmName, 'gig-ica') || strcmpi(algorithmName, 'constrained ica (spatial)') || ...
                strcmpi(algorithmName, 'semi-blind infomax'))
            error(['Temporal ica cannot be run using algorithm ', algorithmName]);
        end
    end
    
end

if (~strcmpi(algorithmName, 'iva-gl') && ~strcmpi(algorithmName, 'gig-ica') && ~strcmpi(algorithmName, 'constrained ica (spatial)'))
    
    if (sesInfo.userInput.numOfSub*sesInfo.userInput.numOfSess == 1)
        sesInfo.userInput.numReductionSteps = 1;
    elseif (sesInfo.userInput.numOfSub*sesInfo.userInput.numOfSess > 1)
        if (sesInfo.userInput.numReductionSteps > 2)
            sesInfo.userInput.numReductionSteps = 2;
        end
    end
else
    sesInfo.userInput.numReductionSteps = 1;
end


%setup rest of reduction step values
sesInfo.numReductionSteps = sesInfo.userInput.numReductionSteps;

PCs = [sesInfo.userInput.numOfPC1, sesInfo.userInput.numOfPC2];

sesInfo.userInput.numComp = PCs(sesInfo.numReductionSteps);

sesInfo.numComp = sesInfo.userInput.numComp;

%get number of subjects and sessions
sesInfo.numOfSub = sesInfo.userInput.numOfSub;
sesInfo.numOfSess = sesInfo.userInput.numOfSess;
sesInfo.numOfDataSets = sesInfo.userInput.numOfSub*sesInfo.userInput.numOfSess;
if isfield(sesInfo.userInput, 'dataType')
    sesInfo.dataType = sesInfo.userInput.dataType;
else
    sesInfo.dataType = 'real';
end

write_analysis_steps_in_dirs = 0;
if (isfield(sesInfo.userInput, 'write_analysis_steps_in_dirs'))
    write_analysis_steps_in_dirs = sesInfo.userInput.write_analysis_steps_in_dirs;
end

conserve_disk_space = 0;
if (isfield(sesInfo.userInput, 'conserve_disk_space'))
    conserve_disk_space = sesInfo.userInput.conserve_disk_space;
end

sesInfo.conserve_disk_space = conserve_disk_space;
sesInfo.write_analysis_steps_in_dirs = write_analysis_steps_in_dirs;


%naming convention for matlab output files
sesInfo.param_file = [sesInfo.userInput.prefix, PARAMETER_INFO_MAT_FILE];
sesInfo.data_reduction_mat_file = [sesInfo.userInput.prefix,DATA_REDUCTION_MAT_FILE];
sesInfo.ica_mat_file = [sesInfo.userInput.prefix,ICA_MAT_FILE];
sesInfo.back_reconstruction_mat_file = [sesInfo.userInput.prefix,BACK_RECONSTRUCTION_MAT_FILE];
sesInfo.calibrate_components_mat_file = [sesInfo.userInput.prefix,CALIBRATE_MAT_FILE];
sesInfo.aggregate_components_an3_file = [sesInfo.userInput.prefix,AGGREGATE_AN3_FILE];


% get the count for time points
if isfield(sesInfo.userInput, 'diffTimePoints')
    countFiles = sesInfo.userInput.diffTimePoints; % different time points
else
    % get the count of the files
    [countFiles] = icatb_get_countTimePoints(sesInfo.userInput.files);
end
% end for checking the time points

checkTimePoints = find(countFiles ~= countFiles(1));

if ~isempty(checkTimePoints)
    flagTimePoints = 'different_time_points';
else
    flagTimePoints = 'same_time_points';
end

% store these two parameters
sesInfo.flagTimePoints = flagTimePoints;
sesInfo.diffTimePoints = countFiles;

sesInfo.icaOutputFiles = icatb_getOutputFileNames(sesInfo.userInput.prefix, sesInfo.numComp, sesInfo.numOfSub, sesInfo.numOfSess, sesInfo.write_analysis_steps_in_dirs);


if (sesInfo.write_analysis_steps_in_dirs)
    
    sesInfo.data_reduction_mat_file = fullfile([sesInfo.userInput.prefix, '_data_reduction_files'], regexprep(DATA_REDUCTION_MAT_FILE, '^(_)', ''));
    sesInfo.ica_mat_file = fullfile([sesInfo.userInput.prefix, '_ica_files'], regexprep(ICA_MAT_FILE, '^(_)', ''));
    sesInfo.aggregate_components_an3_file = fullfile([sesInfo.userInput.prefix, '_ica_files'], regexprep(AGGREGATE_AN3_FILE, '^(_)', ''));
    sesInfo.back_reconstruction_mat_file = fullfile([sesInfo.userInput.prefix, '_back_reconstruction_files'], regexprep(BACK_RECONSTRUCTION_MAT_FILE, '^(_)', ''));
    sesInfo.calibrate_components_mat_file = fullfile([sesInfo.userInput.prefix, '_scaling_components_files'], regexprep(CALIBRATE_MAT_FILE, '^(_)', ''));
    
end

appDataName = 'gica_waitbar_app_data';
if~isempty(statusHandle)
    % Display Waitbar
    
    % get the status handles
    statusData = getappdata(statusHandle, appDataName);
    statusData.perCompleted = statusData.perCompleted + statusData.unitPerCompleted;
    set(statusHandle, 'name', [num2str(round(statusData.perCompleted*100)), '% analysis done']);
    setappdata(statusHandle, appDataName, statusData);
    waitbar(statusData.perCompleted, statusHandle);
    
end

%input files
sesInfo.inputFiles = sesInfo.userInput.files;
sesInfo.numOfScans = sesInfo.diffTimePoints(1);

% get the naming of the complex images
[sesInfo, complexInfo] = icatb_name_complex_images(sesInfo, 'read');

if ~strcmpi(modalityType, 'eeg')
    % return HInfo
    [V, sesInfo.HInfo] = icatb_returnHInfo(sesInfo.userInput.files(1).name(1, :));
    clear V;
else
    sesInfo.HInfo = sesInfo.userInput.HInfo;
end

if ~isempty(statusHandle)
    
    % get the status handles
    statusData = getappdata(statusHandle, appDataName);
    statusData.perCompleted = statusData.perCompleted + statusData.unitPerCompleted;
    setappdata(statusHandle, appDataName, statusData);
    set(statusHandle, 'name', [num2str(round(statusData.perCompleted*100)), '% analysis done']); waitbar(statusData.perCompleted, statusHandle);
    
end


% Analysis Type
which_analysis = 1;
if (isfield(sesInfo.userInput, 'which_analysis'))
    which_analysis = sesInfo.userInput.which_analysis;
end
sesInfo.which_analysis = which_analysis;

if (which_analysis == 2)
    if isfield(sesInfo.userInput, 'icasso_opts')
        sesInfo.icasso_opts = sesInfo.userInput.icasso_opts;
    else
        sesInfo.icasso_opts = struct('sel_mode', 'randinit', 'num_ica_runs', max([2, NUM_RUNS_GICA]));
    end
elseif (which_analysis == 3)
    if isfield(sesInfo.userInput, 'mst_opts')
        sesInfo.mst_opts = sesInfo.userInput.mst_opts;
    else
        sesInfo.mst_opts.num_ica_runs = max([2, NUM_RUNS_GICA]);
    end
end

%--algorithm
sesInfo.algorithm = sesInfo.userInput.algorithm;

% TR in seconds
try
    sesInfo.TR = sesInfo.userInput.TR;
catch
end

% Data pre-processing option
preprocType = 'remove mean per timepoint';
if (isfield(sesInfo.userInput, 'preproc_type'))
    preprocType = lower(sesInfo.userInput.preproc_type);
end

sesInfo.preproc_type = preprocType;

%-- Back reconstruction type
backReconType = 'regular';
if (isfield(sesInfo.userInput, 'backReconType'))
    backReconType = sesInfo.userInput.backReconType;
end

if ((sesInfo.numOfSub*sesInfo.numOfSess > 1) && (sesInfo.numReductionSteps == 1))
    backReconType = 'spatial-temporal regression';
end

if (strcmpi(backReconType, 'str'))
    backReconType = 'spatial-temporal regression';
end

sesInfo.backReconType = backReconType;

if isfield(sesInfo, 'covariance_opts')
    sesInfo.covariance_opts = rmfield(sesInfo, 'covariance_opts');
end

% %-- Covariance Options
if ~isfield(sesInfo.userInput, 'pcaType')
    sesInfo.userInput.pcaType = 'standard';
end

sesInfo.pcaType = sesInfo.userInput.pcaType;

if (isfield(sesInfo, 'covariance_opts'))
    sesInfo = rmfield(sesInfo, 'covariance_opts');
end

sesInfo.userInput = icatb_check_pca_opts(sesInfo.userInput);
sesInfo.pca_opts = sesInfo.userInput.pca_opts;


%-- Group PCA settings
if ~isfield(sesInfo.userInput, 'group_pca_type')
    sesInfo.userInput.group_pca_type = 'subject specific';
end

sesInfo.group_pca_type = sesInfo.userInput.group_pca_type;

perfOptions = icatb_get_analysis_settings;
perfType = 'user specified settings';
if (isfield(sesInfo.userInput, 'perfType'))
    perfType = sesInfo.userInput.perfType;
end

if (isnumeric(perfType))
    perfType = perfOptions{perfType};
end

perfType = lower(perfType);

sesInfo.perfType = perfType;

% Get ica options
ICA_Options = {};
if (isfield(sesInfo.userInput, 'ICA_Options'))
    ICA_Options = sesInfo.userInput.ICA_Options;
end

sesInfo.ICA_Options = ICA_Options;

%--scaleType
sesInfo.scaleType = sesInfo.userInput.scaleType;

%--Convert To Z
if isfield(sesInfo.userInput, 'convertToZ')
    convertToZ = sesInfo.userInput.convertToZ;
    if convertToZ
        sesInfo.scaleType = 2;
    end
    sesInfo.userInput = rmfield(sesInfo.userInput, 'convertToZ');
    if isfield(sesInfo, 'convertToZ')
        sesInfo = rmfield(sesInfo, 'convertToZ');
    end
end

% get the mask indices
sesInfo.mask_ind = sesInfo.userInput.mask_ind;

sesInfo = icatb_dataReductionSetup(sesInfo);

%check to make sure valid parameters
icatb_parameterErrorCheck(sesInfo);


% Store flip parameter
%sesInfo.flip_analyze_images = FLIP_ANALYZE_IMAGES;

%parameters are now initialized and ready for analysis
sesInfo.isInitialized = 1;

saveParamFile = 1;
if (isfield(sesInfo, 'saveParamFile'))
    saveParamFile = sesInfo.saveParamFile;
end

if (saveParamFile)
    %save group information in matlab file
    [pp, fileName] = fileparts(sesInfo.userInput.param_file);
    icatb_save(fullfile(outputDir, [fileName, '.mat']), 'sesInfo');
end

if ~isempty(statusHandle)
    % get the status handles
    statusData = getappdata(statusHandle, appDataName);
    statusData.perCompleted = statusData.perCompleted + statusData.unitPerCompleted;
    setappdata(statusHandle, appDataName, statusData);
    set(statusHandle, 'name', [num2str(round(statusData.perCompleted*100)), '% analysis done']); waitbar(statusData.perCompleted, statusHandle);
end
