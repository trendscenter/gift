function dfncInfo = icatb_run_dfnc(param_file)
%% Run dfnc

icatb_defaults;
global DETRENDNUMBER;
global DFNC_DEFAULTS;
if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('title', 'Select dFNC Parameter File', 'typeEntity', 'file', 'typeSelection', ...
        'single', 'filter', '*dfnc.mat');
    drawnow;
end

if (isempty(param_file))
    error('Please select dFNC parameter file');
end

if (ischar(param_file))
    load(param_file);
    outputDir = fileparts(param_file);
    if (isempty(outputDir))
        outputDir = pwd;
    end
    if (~exist('dfncInfo', 'var'))
        error(['Selected file ', param_file, ' is not a valid dFNC parameter file']);
    end
else
    dfncInfo = param_file;
    outputDir = dfncInfo.userInput.outputDir;
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
% average sliding window correlation
aswc = 0;

try
    aswc = DFNC_DEFAULTS.aswc;
catch
end

try
    aswc = dfncInfo.userInput.aswc;
catch
end


window_step_size = 1;

try
    window_step_size = DFNC_DEFAULTS.step_size;
catch
end

if (~window_step_size)
    window_step_size = 1;
end

% run tvdfnc
tvdfnc = 0;
try
    tvdfnc = DFNC_DEFAULTS.tvdfnc;
catch
end

try
    tvdfnc = dfncInfo.userInput.tv_dfnc;
catch
end

try
    detrend_no = dfncInfo.userInput.feature_params.final.tc_detrend;
    doDespike = dfncInfo.userInput.feature_params.final.tc_despike;
    tc_filter = dfncInfo.userInput.feature_params.final.tc_filter;
    wsize = dfncInfo.userInput.feature_params.final.wsize;
    window_alpha = dfncInfo.userInput.feature_params.final.window_alpha;
    num_repetitions = dfncInfo.userInput.feature_params.final.num_repetitions;
    method = dfncInfo.userInput.feature_params.final.method;
    covariates = dfncInfo.userInput.feature_params.final.tc_covariates;
    window_type = dfncInfo.userInput.feature_params.final.window_type;
catch
end

%% Tukey window option is disabled for now.
if (strcmpi(window_type, 'tukey'))
    %if (window_alpha <= 1)
    window_alpha = 3;
    %end
    disp(['Tukey window option is currently disabled. Using gaussian window with alpha = ', ...
        num2str(window_alpha), ' ...']);
end


TR = dfncInfo.userInput.TR;

load(dfncInfo.userInput.ica_param_file);

comps = [dfncInfo.userInput.comp.value];

if (length(comps) ~= length(unique(comps)))
    error('One or more components are replicated across different component networks');
end

if (length(TR) == 1)
    TR = repmat(TR, 1, sesInfo.numOfSub);
end

if ((length(TR) > 1) && (length(TR) ~= sesInfo.numOfSub))
    error(['You have specified multiple TRs. TRs should match the length of number of subjects (', num2str(sesInfo.numOfSub), ')']);
end


dfncInfo.prefix = dfncInfo.userInput.prefix;
dfncInfo.comps = comps;
dfncInfo.outputDir = outputDir;
dfncInfo.TR = TR;
dfncInfo.detrend_no = detrend_no;
dfncInfo.doDespike = doDespike;
dfncInfo.tc_filter = tc_filter;
dfncInfo.method = method;
dfncInfo.wsize = wsize;
dfncInfo.window_alpha = window_alpha;
dfncInfo.window_type = window_type;
dfncInfo.num_repetitions = num_repetitions;

outDir = fileparts(dfncInfo.userInput.ica_param_file);
sesInfo.outputDir = outDir;


tic;
logFile = fullfile(outputDir, [dfncInfo.prefix, '_results.log']);
diary(logFile);

disp('----------------------------------------------------');
disp('-------------- STARTING DYNAMIC FNC --------------');
disp('----------------------------------------------------');

covariate_files = '';
scansToInclude = [];
if (~strcmpi(covariates, 'none'))
    try
        covariate_files = dfncInfo.userInput.feature_params.final.tc_covariates_userdata.filesList;
        scansToInclude = dfncInfo.userInput.feature_params.final.tc_covariates_userdata.file_numbers;
        disp('Variance associated with the covariates will be removed from the timecourses ... ');
    catch
    end
end

fprintf('\n');

modelTC = {};
if (isfield(dfncInfo.userInput, 'selectedRegressors'))
    selectedRegressors = dfncInfo.userInput.selectedRegressors;
    spm_files  = cellstr(char(sesInfo.userInput.designMatrix.name));
    
    countTp = 0;
    modelTC = cell(sesInfo.numOfSub, sesInfo.numOfSess);
    if (length(spm_files) == 1)
        spmInfo = load (deblank(spm_files{1}));
    end
    
    for nSub = 1:sesInfo.numOfSub
        
        if (length(spm_files) > 1)
            spmInfo = load (deblank(spm_files{nSub}));
        end
        
        for nSess = 1:sesInfo.numOfSess
            countTp = countTp + 1;
            disp(['Loading model timecourse for subject ', num2str(nSub), ' Session ', num2str(nSess), ' ...']);
            modelX = icatb_getModelTimeCourse(sesInfo.diffTimePoints(countTp), spmInfo, selectedRegressors, nSess);
            modelTC{nSub, nSess} = modelX;
        end
    end
    
    
end

outFiles = cell(1, sesInfo.numOfSub);

numOfSub = sesInfo.numOfSub;
numOfSess = sesInfo.numOfSess;
dfncPrefix = dfncInfo.userInput.prefix;


if (strcmpi(doDespike, 'yes') && (tc_filter > 0))
    disp('Despiking and filtering timecourses ...');
elseif (strcmpi(doDespike, 'yes') && (tc_filter == 0))
    disp('Despiking timecourses ...');
elseif (strcmpi(doDespike, 'no') && (tc_filter > 0))
    disp('Filtering timecourses ...');
end

minTR =  min(TR);
averageLengthTR =  round(aswc /minTR) ;
time_points = sesInfo.diffTimePoints;

TRN = repmat(TR(:)', sesInfo.numOfSess, 1);
TRN = TRN(:)';
interpFactor = TRN./minTR;
[num, denom] = rat(interpFactor);
minTP = min(ceil((time_points.*num) ./ denom));

preprocess_params.detrend_no = detrend_no;
preprocess_params.doDespike = doDespike;
preprocess_params.tc_filter = tc_filter;

windowing_params.method = dfncInfo.method;
windowing_params.wsize = dfncInfo.wsize;
windowing_params.window_alpha = dfncInfo.window_alpha;
windowing_params.num_L1_repetitions = dfncInfo.num_repetitions;
windowing_params.window_type = dfncInfo.window_type;

varsToSave = {'FNCdyn', 'tc'};

if (strcmpi(dfncInfo.method, 'L1'))
    varsToSave(end + 1:end + 2) = {'Pdyn', 'Lambdas'};
end

if (tvdfnc)
    varsToSave(end + 1) = {'FNC_tvdfnc'};
end

if (~isempty(modelTC))
    varsToSave(end + 1) = {'task_connectivity'};
end

if (~isempty(covariate_files))
    covariate_files = cellstr(covariate_files);
end


best_lambda = zeros(1, numOfSub*numOfSess);

%% Loop over subjects
parfor dataSetCount = 1:numOfSub*numOfSess
    
    vars = cell(1, length(varsToSave));
    
    nSub = ceil(dataSetCount/numOfSess);
    nSess = mod(dataSetCount - 1, numOfSess) + 1;
    
    results_file = [dfncPrefix, '_sub_', icatb_returnFileIndex(nSub), '_sess_', icatb_returnFileIndex(nSess), '_results.mat'];
    
    modelX = [];
    try
        modelX = modelTC{nSub, nSess};
    catch
    end
    
    covariateInfo = [];
    if (~isempty(covariate_files))
        covariateInfo.file = covariate_files{dataSetCount};
        covariateInfo.file_numbers = scansToInclude;
    end
    
    disp(['Computing dynamic FNC on subject ', num2str(nSub), ' session ', num2str(nSess)]);
    
    tc = icatb_loadComp(sesInfo, comps, 'subjects', nSub, 'sessions', nSess, 'vars_to_load', 'tc', 'detrend_no', detrend_no, 'covariates', ...
        covariate_files, 'scansToInclude', scansToInclude);
    
    %preprocess_params.detrend_no = [];
    corrInfo = icatb_compute_dfnc(tc, TR(nSub), 'mintr', minTR, 'minTP', minTP, 'preprocess_params', preprocess_params, ...
        'windowing_params', windowing_params, 'modelTC', modelX, 'covariateInfo', covariateInfo);
    
    tc = corrInfo.X;
    FNCdyn = corrInfo.FNCdyn;
    
    if (averageLengthTR ~= 0)
        disp('Computing average sliding window correlation ...');
        FNCdyn = f_aswc(FNCdyn, averageLengthTR);
    end
    
    vars(1:2) = {FNCdyn, tc};
    
    varCount = 2;
    if (strcmpi(windowing_params.method, 'L1'))
        best_lambda(dataSetCount) = corrInfo.best_lambda;
        vars(varCount+1:varCount+2) = {corrInfo.Pdyn, corrInfo.Lambdas};
        varCount = varCount + 2;
    end
    
    if (tvdfnc)
        disp('Computing temporal variation dfnc and stacking the matrices FNC and tvdfnc ...');
        FNC_tvdfnc = icatb_compute_tvdfnc(FNCdyn);
        varCount = varCount + 1;
        vars(varCount) = {FNC_tvdfnc};
    end
    
    if (~isempty(modelTC))
        varCount = varCount + 1;
        vars(varCount) = {corrInfo.task_connectivity};
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

dfncInfo.best_lambda = best_lambda;
dfncInfo.outputFiles = outFiles;
fileN = fullfile(outputDir, [dfncInfo.prefix, '.mat']);
disp(['Saving parameter file ', fileN, ' ...']);
save(fileN, 'dfncInfo');

totalTime = toc;

fprintf('\n');

disp('Analysis Complete');

disp(['Total time taken to complete the analysis is ', num2str(totalTime/60), ' minutes']);

diary('off');

fprintf('\n');


function outData = f_aswc( data, aSz )
%
% Victor M. Vergara, PhD
% Created: 2018-08-01
%
% This function takes SWC data [ Time X nodes (FNC, etc..)] and performs an
% sliding averaging of the data for each column
%
N = size(data,1)-aSz;
outData = zeros(N,size(data,2));
for kk=1:N
    outData(kk,:) = mean(data(kk:kk+aSz-1,:));
end

%end

% function [modelX, model_inds] = getModelTimeCourse(tp, spmInfo, selectedRegressors, nSess)
%
% modelX = repmat(NaN, tp, length(selectedRegressors));
%
%
% XX = spmInfo.SPM.xX.X;
% names = cellstr(char(spmInfo.SPM.xX.name));
%
% try
%     spmInfo.SPM.nscan(nSess);
% catch
%     nSess = 1;
% end
%
% if (nSess == 1)
%     startTp = 1;
% else
%     startTp = sum(spmInfo.SPM.nscan(1:nSess-1)) + 1;
% end
%
% endTp = sum(spmInfo.SPM.nscan(1:nSess));
% model_inds = [];
% for nR = 1:length(selectedRegressors)
%     regressName = ['Sn(', num2str(nSess), ') ', selectedRegressors{nR}];
%     inds = strmatch(lower(regressName), lower(spmInfo.SPM.xX.name), 'exact');
%     if (~isempty(inds))
%         model_inds = [model_inds, nR];
%         modelX(:, nR) = XX(startTp:endTp, inds(1));
%     end
% end

