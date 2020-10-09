function icatb_run_dynamic_coherence(param_file)
%% Run dynamic coherence
%

icatb_defaults;
global DETRENDNUMBER;

%% Select dynamic coherence param file
if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('title', 'Select dynamic coherence parameter file', 'typeEntity', 'file', 'typeSelection', 'single', 'filter', '*dyn_coh.mat');
    drawnow;
    if (isempty(param_file))
        error('Dynamic Coherence parameter file is not selected');
    end
end

drawnow;

if (ischar(param_file))
    load(param_file);
    outDir = fileparts(param_file);
else
    cohInfo = param_file;
    outDir = cohInfo.userInput.outputDir;
end

if (~exist('cohInfo', 'var'))
    error('Selected file is not a valid dynamic coherence parameter file');
end


if (isempty(outDir))
    outDir = pwd;
end


cd(outDir);


cohInfo.outputDir = outDir;
load(cohInfo.userInput.ica_param_file);

icaDir = fileparts(cohInfo.userInput.ica_param_file);
if (isempty(icaDir))
    icaDir = pwd;
end

sesInfo.outputDir = icaDir;
sesInfo.userInput.pwd = icaDir;
numOfSub = sesInfo.numOfSub;
numOfSess = sesInfo.numOfSess;


numClusters = cohInfo.userInput.kmeans.num_clusters;
numRuns = cohInfo.userInput.kmeans.num_runs;
regress_tc = '';
try
    regress_tc = cohInfo.regress_tc;
catch
end

despike_tc = 'no';
try
    despike_tc = cohInfo.userInput.feature_params.final.tc_despike;
catch
end

tc_detrend = DETRENDNUMBER;

try
    tc_detrend = cohInfo.userInput.feature_params.final.tc_detrend;
catch
end


covariates = 'none';
try
    covariates = cohInfo.userInput.feature_params.final.tc_covariates;
catch
end

covariate_files = '';
scansToInclude = [];
if (~strcmpi(covariates, 'none'))
    try
        covariate_files = cohInfo.userInput.feature_params.final.tc_covariates_userdata.filesList;
        scansToInclude = cohInfo.userInput.feature_params.final.tc_covariates_userdata.file_numbers;
        disp('Variance associated with the covariates will be removed from the timecourses ... ');
    catch
    end
end

TR = cohInfo.userInput.TR;

cohInfo.TR = TR;
cohInfo.kmeans.num_clusters = numClusters;
cohInfo.kmeans.num_runs = numRuns;

comps = [cohInfo.userInput.comp.value];
cohInfo.comps = comps;
cohInfo.prefix = cohInfo.userInput.prefix;

tic;
logFile = fullfile(outDir, [cohInfo.prefix, '_results.log']);
diary(logFile);

disp('----------------------------------------------------');
disp('-------------- STARTING DYNAMIC COHERENCE --------------');
disp('----------------------------------------------------');


[tc, TR] = icatb_loadAndInterpTC(sesInfo, comps, TR, 'detrend_no', tc_detrend, 'covariates', covariate_files, 'scansToInclude', scansToInclude);
cohInfo.TRN = TR;

numHistBins = 5;
fscales = linspace(0.01, 1/TR/2, numHistBins);
cohInfo.fscales = fscales;

%tc = icatb_loadComp(sesInfo, cohInfo.comps, 'vars_to_load', 'tc', 'truncate_tp', 1, 'covariates', regress_tc);

%tc = reshape(tc, numOfSub, numOfSess, min(sesInfo.diffTimePoints), length(cohInfo.comps));

size_tc = size(tc);
%tc = reshape(tc, numOfSub*numOfSess,  min(sesInfo.diffTimePoints), length(cohInfo.comps));
tc = reshape(tc, numOfSub*numOfSess,  size_tc(end-1), size_tc(end));


if (strcmpi(despike_tc, 'yes'))
    fprintf('Despiking timecourses\n');
    for nD = 1:numOfSub*numOfSess
        tmp = squeeze(tc(nD, :, :));
        tmp = icatb_despike_tc(tmp, TR);
        tc(nD, :, :) = tmp;
    end
end


settings_file = which('WaveletInfo_TR_2.mat');

disp('Computing dynamic coherence ...');

[Corrected_Out_Distribution, Null_Distribution, mask_v_inds, coin] = icatb_compute_dynamic_coherence(tc, settings_file);

[IDX, centroids, sumD] = icatb_compute_complex_kmeans(Corrected_Out_Distribution, coin, numClusters, numRuns);

clusterInfo.IDX = IDX;
clusterInfo.centroids = centroids;
clusterInfo.sumD = sumD;

outFile = fullfile(outDir, [cohInfo.prefix, '_results.mat']);

save(outFile, 'clusterInfo', 'Corrected_Out_Distribution', 'Null_Distribution', 'coin', 'mask_v_inds');
icatb_save(fullfile(outDir, [cohInfo.prefix, '.mat']), 'cohInfo');

totalTime = toc;

fprintf('\n');

disp('Analysis Complete');

disp(['Total time taken to complete the analysis is ', num2str(totalTime/60), ' minutes']);

diary('off');

fprintf('\n');


