function schronnInfo = icatb_run_spatial_chronnectome(param_file)
%% Run dfnc

icatb_defaults;
global DETRENDNUMBER;
global PREPROC_DEFAULT;

if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('title', 'Select spatial chronnectome Parameter File', 'typeEntity', 'file', 'typeSelection', ...
        'single', 'filter', '*schronn.mat');
    drawnow;
end

if (isempty(param_file))
    error('Please select spatial chronnectome parameter file');
end

if (ischar(param_file))
    load(param_file);
    outputDir = fileparts(param_file);
    if (isempty(outputDir))
        outputDir = pwd;
    end
    if (~exist('schronnInfo', 'var'))
        error(['Selected file ', param_file, ' is not a valid spatial chronnectome parameter file']);
    end
else
    schronnInfo = param_file;
    outputDir = schronnInfo.userInput.outputDir;
    clear param_file;
end

cd (outputDir);

%% Initialise
outputDir = pwd;
detrend_no = DETRENDNUMBER;
doDespike = 'no';
tc_filter = 0;
method = 'none';
wsize = 30;
window_alpha = 3;
window_type = 'gaussian';
num_repetitions = 10;
covariates = 'none';


try
    detrend_no = schronnInfo.userInput.feature_params.final.tc_detrend;
    doDespike = schronnInfo.userInput.feature_params.final.tc_despike;
    tc_filter = schronnInfo.userInput.feature_params.final.tc_filter;
    wsize = schronnInfo.userInput.feature_params.final.wsize;
    window_alpha = schronnInfo.userInput.feature_params.final.window_alpha;
    method = schronnInfo.userInput.feature_params.final.method;
    covariates = schronnInfo.userInput.feature_params.final.tc_covariates;
    num_repetitions = schronnInfo.userInput.feature_params.final.num_repetitions;
catch
end


TR = schronnInfo.userInput.TR;

load(schronnInfo.userInput.ica_param_file);

comps = schronnInfo.userInput.comp;

if (length(comps) ~= length(unique(comps)))
    error('One or more components are replicated across different component networks');
end

if ((length(TR) > 1) && (length(TR) ~= sesInfo.numOfSub))
    error(['You have specified multiple TRs. TRs should match the length of number of subjects (', num2str(sesInfo.numOfSub), ')']);
end


schronnInfo.HInfo = sesInfo.HInfo;
schronnInfo.mask_ind = sesInfo.mask_ind;
schronnInfo.prefix = schronnInfo.userInput.prefix;
schronnInfo.comps = comps;
schronnInfo.outputDir = outputDir;
schronnInfo.TR = TR;
schronnInfo.detrend_no = detrend_no;
schronnInfo.doDespike = doDespike;
schronnInfo.tc_filter = tc_filter;
schronnInfo.method = method;
schronnInfo.wsize = wsize;
schronnInfo.window_alpha = window_alpha;
schronnInfo.window_type = window_type;
schronnInfo.num_repetitions = num_repetitions;

outDir = fileparts(schronnInfo.userInput.ica_param_file);
sesInfo.outputDir = outDir;


tic;
logFile = fullfile(outputDir, [schronnInfo.prefix, '_results.log']);
diary(logFile);

disp('----------------------------------------------------');
disp('-------------- STARTING SPATIAL CHRONNECTOME  --------------');
disp('----------------------------------------------------');

disp('Loading timecourses ....');

covariate_files = '';
scansToInclude = [];
if (~strcmpi(covariates, 'none'))
    try
        covariate_files = cellstr(schronnInfo.userInput.feature_params.final.tc_covariates_userdata.filesList);
        scansToInclude = schronnInfo.userInput.feature_params.final.tc_covariates_userdata.file_numbers;
        disp('Variance associated with the covariates will be removed from the timecourses ... ');
    catch
    end
end

fprintf('\n');

%nT = min(sesInfo.diffTimePoints);

if (any(sesInfo.diffTimePoints ~= sesInfo.diffTimePoints(1)))
    error('Timepoints must be the same across data-sets');
end

minTR = min(TR);

outFiles = cell(sesInfo.numOfSub*sesInfo.numOfSess, length(comps));

numOfSub = sesInfo.numOfSub;
numOfSess = sesInfo.numOfSess;
schronnPrefix = schronnInfo.userInput.prefix;

if (length(TR) == 1)
    TR = repmat(TR, 1, sesInfo.numOfSub);
end

time_points = sesInfo.diffTimePoints;
TRN = repmat(TR(:)', sesInfo.numOfSess, 1);
TRN = TRN(:)';
interpFactor = TRN./minTR;
[num, denom] = rat(interpFactor);
minTP = min(ceil((time_points.*num) ./ denom));

preprocess_params.detrend_no = detrend_no;
preprocess_params.doDespike = doDespike;
preprocess_params.tc_filter = tc_filter;

windowing_params.method = schronnInfo.method;
windowing_params.wsize = schronnInfo.wsize;
windowing_params.window_alpha = schronnInfo.window_alpha;
windowing_params.num_L1_repetitions = schronnInfo.num_repetitions;
windowing_params.window_type = schronnInfo.window_type;


for nc = 1:length(comps)
    dirname = [schronnPrefix, '_comp_', icatb_returnFileIndex(comps(nc))];
    if (exist(dirname, 'dir') ~= 7)
        mkdir(outputDir, dirname);
    end
end

clear dirname;


mask_ind = sesInfo.mask_ind;
inputFiles = sesInfo.inputFiles;
preproc_default = PREPROC_DEFAULT;

tmpSesInfo = sesInfo;

if (~isempty(covariate_files))
    covariate_files = cellstr(covariate_files);
end

%% Loop over subjects
parfor dataSetCount = 1:numOfSub*numOfSess
    
    nSub = ceil(dataSetCount/numOfSess);
    nSess = mod(dataSetCount - 1, numOfSess) + 1;
    
    disp(['Computing spatial chronnectome on subject ', num2str(nSub), ' session ', num2str(nSess)]);
    
    if (strcmpi(doDespike, 'yes') && (tc_filter > 0))
        disp('Despiking and filtering timecourses ...');
    elseif (strcmpi(doDespike, 'yes') && (tc_filter == 0))
        disp('Despiking timecourses ...');
    elseif (strcmpi(doDespike, 'no') && (tc_filter > 0))
        disp('Filtering timecourses ...');
    end
    
    
    covariateInfo = [];
    if (~isempty(covariate_files))
        covariateInfo.file = covariate_files{dataSetCount};
        covariateInfo.file_numbers = scansToInclude;
    end
    
    tc = icatb_loadComp(tmpSesInfo, comps, 'vars_to_load', 'tc', 'detrend_no', [], 'subjects', nSub, 'sessions', nSess);
    
    % Load BOLD data
    bold_data = icatb_remove_mean(icatb_preproc_data(icatb_read_data(inputFiles(dataSetCount).name, [], mask_ind, 'double'), preproc_default));
    bold_data = bold_data';
    
    corrInfo = icatb_compute_dfnc(tc, TR(nSub), 'Y', bold_data, 'mintr', minTR, 'minTP', minTP, 'preprocess_params', ...
        preprocess_params, 'windowing_params', windowing_params, 'covariateInfo', covariateInfo);
    
    dc_maps_all = corrInfo.FNCdyn;
    
    
    %% save results
    oFiles = cell(1, length(comps));
    for ncomp = 1:length(comps)
        dirname = [schronnPrefix, '_comp_', icatb_returnFileIndex(comps(ncomp))];
        dynamic_coupling_maps = squeeze(dc_maps_all(:, ncomp, :));
        compResultsFile = [schronnPrefix, '_sub_', icatb_returnFileIndex(nSub), '_s', num2str(nSess), '_comp_', icatb_returnFileIndex(comps(ncomp)), '_results.mat'];
        results_file = fullfile(dirname, compResultsFile);
        oFiles{ncomp} = results_file;
        disp(['.... saving file ', results_file]);
        icatb_parSave(fullfile(outputDir, results_file), {dynamic_coupling_maps}, {'dynamic_coupling_maps'});
    end
    
    outFiles(dataSetCount, :) = oFiles;
    
    disp('Done');
    fprintf('\n');
    
    %end
    %% End of loop over sessions
end
%% End of loop over subjects

schronnInfo.outputFiles = outFiles;
fileN = fullfile(outputDir, [schronnInfo.prefix, '.mat']);
disp(['Saving parameter file ', fileN, ' ...']);
save(fileN, 'schronnInfo');

totalTime = toc;

fprintf('\n');

disp('Analysis Complete');

disp(['Total time taken to complete the analysis is ', num2str(totalTime/60), ' minutes']);

diary('off');

fprintf('\n');

