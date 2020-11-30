function icatb_run_dfc_roi(param_file)
%% Setup dynamic FC for ROI-ROI or ROI-voxel
%

icatb_defaults;
global DETRENDNUMBER;
global DFNC_DEFAULTS;
if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('title', 'Select dFC ROI Parameter File', 'typeEntity', 'file', 'typeSelection', ...
        'single', 'filter', '*dfc*roi*.mat');
    drawnow;
end

if (isempty(param_file))
    error('Please select dFC ROI parameter file');
end

if (ischar(param_file))
    load(param_file);
    outputDir = fileparts(param_file);
    if (isempty(outputDir))
        outputDir = pwd;
    end
    if (~exist('dfcRoiInfo', 'var'))
        error(['Selected file ', param_file, ' is not a valid dFNC parameter file']);
    end
else
    dfcRoiInfo = param_file;
    outputDir = dfcRoiInfo.userInput.outputDir;
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

window_step_size = 1;

try
    window_step_size = DFNC_DEFAULTS.step_size;
catch
end

if (~window_step_size)
    window_step_size = 1;
end

%% load mask
atlas_mask = icatb_loadData(dfcRoiInfo.userInput.atlas_mask_file);
if isempty(dfcRoiInfo.userInput.maskFile)
    mask = atlas_mask;
else
    mask =  icatb_loadData(dfcRoiInfo.userInput.maskFile);
end

mask = (abs(mask) > eps);

%% Load timecourses
roi_labels = dfcRoiInfo.userInput.roi_info{1};
roi_indices = dfcRoiInfo.userInput.roi_info{2};
analysisType = dfcRoiInfo.userInput.analysisType;


try
    detrend_no = dfcRoiInfo.userInput.feature_params.final.tc_detrend;
    doDespike = dfcRoiInfo.userInput.feature_params.final.tc_despike;
    tc_filter = dfcRoiInfo.userInput.feature_params.final.tc_filter;
    wsize = dfcRoiInfo.userInput.feature_params.final.wsize;
    window_alpha = dfcRoiInfo.userInput.feature_params.final.window_alpha;
    num_repetitions = dfcRoiInfo.userInput.feature_params.final.num_repetitions;
    method = dfcRoiInfo.userInput.feature_params.final.method;
    covariates = dfcRoiInfo.userInput.feature_params.final.tc_covariates;
    window_type = dfcRoiInfo.userInput.feature_params.final.window_type;
catch
end


TR = dfcRoiInfo.userInput.TR;
numOfDataSets = dfcRoiInfo.userInput.numOfDataSets;

if (length(TR) == 1)
    TR = repmat(TR, 1, numOfDataSets);
end

if ((length(TR) > 1) && (length(TR) ~=  numOfDataSets))
    error(['You have specified multiple TRs. TRs should match the length of number of subjects (', num2str(numOfDataSets), ')']);
end

if (strcmpi(analysisType, 'roi-voxel'))
    method = 'none';
end

time_points = dfcRoiInfo.userInput.time_points;
dfcRoiInfo.prefix = dfcRoiInfo.userInput.prefix;
dfcRoiInfo.outputDir = outputDir;
dfcRoiInfo.TR = TR;
dfcRoiInfo.detrend_no = detrend_no;
dfcRoiInfo.doDespike = doDespike;
dfcRoiInfo.tc_filter = tc_filter;
dfcRoiInfo.method = method;
dfcRoiInfo.wsize = wsize;
dfcRoiInfo.window_alpha = window_alpha;
dfcRoiInfo.window_type = window_type;
dfcRoiInfo.num_repetitions = num_repetitions;


dfcRoiInfo.param_file = dfcRoiInfo.userInput.param_file;
dfc_prefix = [dfcRoiInfo.userInput.prefix, '_dfc_roi'];

tic;
logFile = fullfile(outputDir, [dfc_prefix, '_results.log']);
diary(logFile);

disp('----------------------------------------------------');
disp('-------------- STARTING DYNAMIC FC ROI --------------');
disp('----------------------------------------------------');

covariate_files = '';
scansToInclude = [];
if (~strcmpi(covariates, 'none'))
    try
        covariate_files = dfcRoiInfo.userInput.feature_params.final.tc_covariates_userdata.filesList;
        scansToInclude = dfcRoiInfo.userInput.feature_params.final.tc_covariates_userdata.file_numbers;
        disp('Variance associated with the covariates will be removed from the timecourses ... ');
    catch
    end
end

fprintf('\n');

outFiles = cell(1, numOfDataSets);


if (strcmpi(doDespike, 'yes') && (tc_filter > 0))
    disp('Despiking and filtering timecourses ...');
elseif (strcmpi(doDespike, 'yes') && (tc_filter == 0))
    disp('Despiking timecourses ...');
elseif (strcmpi(doDespike, 'no') && (tc_filter > 0))
    disp('Filtering timecourses ...');
end

minTR =  min(TR);

preprocess_params.detrend_no = detrend_no;
preprocess_params.doDespike = doDespike;
preprocess_params.tc_filter = tc_filter;

windowing_params.method = dfcRoiInfo.method;
windowing_params.wsize = dfcRoiInfo.wsize;
windowing_params.window_alpha = dfcRoiInfo.window_alpha;
windowing_params.num_L1_repetitions = dfcRoiInfo.num_repetitions;
windowing_params.window_type = dfcRoiInfo.window_type;


interpFactor = TR./minTR;
[num, denom] = rat(interpFactor);
minTP = min(ceil((time_points.*num) ./ denom));

varsToSave = {'FNCdyn'};

if (strcmpi(dfcRoiInfo.method, 'L1'))
    varsToSave(end + 1:end + 2) = {'Pdyn', 'Lambdas'};
end

masks_y = cell(1, length(roi_indices));
for nR = 1:length(roi_indices)
    tmp_msk = (atlas_mask == roi_indices(nR));
    if (strcmpi(analysisType, 'roi-voxel'))
        mask2 = (((tmp_msk==0).*mask)==1);
        masks_y{nR} = mask2;
    end
end

best_lambda = zeros(1, numOfDataSets);
dims = size(atlas_mask);

filesList = dfcRoiInfo.userInput.filesInfo.filesList;

if (~isempty(covariate_files))
    covariate_files = cellstr(covariate_files);
end

%% Loop over subjects
parfor dataSetCount = 1:numOfDataSets
    
    vars = cell(1, length(varsToSave));
    
    results_file = [dfc_prefix, '_dataset_', icatb_returnFileIndex(dataSetCount), '_results.mat'];
    
    covariateInfo = [];
    if (~isempty(covariate_files))
        covariateInfo.file = covariate_files{dataSetCount};
        covariateInfo.file_numbers = scansToInclude;
    end
    
    data_sub = icatb_read_data(filesList{dataSetCount}, [], (1:prod(dims)));
    
    X = zeros(size(data_sub, 2), length(roi_indices));
    Y = [];
    for nR = 1:length(roi_indices)
        tmp_msk = (atlas_mask == roi_indices(nR));
        X(:, nR) = mean(data_sub(tmp_msk==1, :))';
    end
    
    disp(['Computing dynamic FC on dataset ', num2str(dataSetCount)]);
    corrInfo = [];
    
    if (strcmpi(analysisType, 'roi-voxel'))
        FNCdyn = cell(length(masks_y), 1);
        for nY = 1:length(masks_y)
            Y = data_sub(masks_y{nY}, :)';
            disp(['Using ROI-Voxel (', roi_labels{roi_indices(nY)}, ') ...']);
            corrInfo = icatb_compute_dfnc(X(:, nY), TR(dataSetCount), 'Y', Y, 'mintr', minTR, 'minTP', minTP, ...
                'preprocess_params', preprocess_params, 'windowing_params', windowing_params, 'covariateInfo', ...
                covariateInfo);
            tmp_fnc = zeros(size(corrInfo.FNCdyn, 1), prod(dims));
            tmp_fnc(:,  masks_y{nY}) = squeeze(corrInfo.FNCdyn);
            FNCdyn{nY} = tmp_fnc;
        end
        
    else
        disp('Using ROI-ROI...');
        corrInfo = icatb_compute_dfnc(X, TR(dataSetCount), 'mintr', minTR, 'minTP', minTP, ...
            'preprocess_params', preprocess_params, 'windowing_params', windowing_params, 'covariateInfo', ...
            covariateInfo);
        FNCdyn = corrInfo.FNCdyn;
    end
    
    vars{1} = FNCdyn;
    
    
    varCount = 1;
    if (strcmpi(windowing_params.method, 'L1'))
        best_lambda(dataSetCount) = corrInfo.best_lambda;
        vars(varCount+1:varCount+2) = {corrInfo.Pdyn, corrInfo.Lambdas};
    end
    
    disp(['.... saving file ', results_file]);
    icatb_parSave(fullfile(outputDir, results_file), vars, varsToSave);
    
    outFiles{dataSetCount} = results_file;
    disp('Done');
    fprintf('\n');
    
    %end
    %% End of loop over sessions
end
%% End of loop over subjects

dfcRoiInfo.best_lambda = best_lambda;
dfcRoiInfo.outputFiles = outFiles;
dfcRoiInfo.masks_y = masks_y;
fileN = fullfile(outputDir, dfcRoiInfo.param_file);
disp(['Saving parameter file ', fileN, ' ...']);
save(fileN, 'dfcRoiInfo');

totalTime = toc;

fprintf('\n');

disp('Analysis Complete');

disp(['Total time taken to complete the analysis is ', num2str(totalTime/60), ' minutes']);

diary('off');

fprintf('\n');




