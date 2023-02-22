function icatb_run_sdh(sdhInfo)
%% Run spatial dynamics hierarchy
%

%% Initialize variables
icatb_defaults;
global EXPERIMENTAL_TR;
global PREPROC_DEFAULT;

TR = EXPERIMENTAL_TR;

resultsDir = sdhInfo.outputDir;
load(sdhInfo.ica_parameter_file);

doDespike = 'yes';
detrend_no = 3;
tc_filter = 0.15;
covariates = 'none';

try
    detrend_no = sdhInfo.feature_params.final.tc_detrend;
    doDespike = sdhInfo.feature_params.final.tc_despike;
    tc_filter = sdhInfo.feature_params.final.tc_filter;
    covariates = sdhInfo.feature_params.final.tc_covariates;
catch
end


covariate_files = '';
scansToInclude = [];
if (~strcmpi(covariates, 'none'))
    try
        covariate_files = cellstr(sdhInfo.feature_params.final.tc_covariates_userdata.filesList);
        scansToInclude = sdhInfo.feature_params.final.tc_covariates_userdata.file_numbers;
    catch
    end
end

cd(resultsDir);

output_LogFile = fullfile(resultsDir, [sesInfo.userInput.prefix, '_sdh_results.log']);

tic;

% Print output to a file
diary(output_LogFile);

load(fullfile(sesInfo.outputDir, [sesInfo.ica_mat_file, '.mat']), 'icasig');
icasig = icasig';

numOfSub = sesInfo.numOfSub;
numOfSess = sesInfo.numOfSess;
subjectFiles = cell(numOfSub*numOfSess, 1);
inputFiles = sesInfo.inputFiles;
prefix = sesInfo.userInput.prefix;
mask_ind = sesInfo.mask_ind;
sdhInfo.prefix = [prefix, '_sdh'];
sdhInfo.param_file = [prefix, '_sdh_info.mat'];
sdhInfo.V = sesInfo.HInfo.V(1);
sdhInfo.mask_ind = mask_ind;
sdhInfo.time_points = sesInfo.diffTimePoints;

if (~isempty(covariate_files))
    covariate_files = cellstr(covariate_files);
end

%% Save components
parfor countD = 1:numOfSub*numOfSess
    
    nSub = ceil(countD/numOfSess);
    nSess = mod(countD - 1, numOfSess) + 1;
    disp(['Loading subject ', num2str(nSub), ' session ', num2str(nSess), ' ...']);
    tmp_dat = icatb_remove_mean(icatb_preproc_data(icatb_read_data(inputFiles(countD).name, [], mask_ind), PREPROC_DEFAULT));
    [tc, ic] = icatb_dual_regress(tmp_dat, icasig);
    
    tc = icatb_detrend(tc, 1, [], detrend_no);
    
    
    if (~isempty(covariate_files))
        disp('Variance associated with the covariates will be removed from the timecourses ... ');
        tc = regress_cov(tc, covariate_files{countD}, scansToInclude);
    end
    
    % Despiking timecourses
    if (strcmpi(doDespike, 'yes'))
        disp('Despiking timeseries ...');
        tc = icatb_despike_tc(tc, TR);
    end
    
    % Filter timecourses
    if (tc_filter > 0)
        disp('Filtering timeseries ...');
        tc = icatb_filt_data(tc, TR, tc_filter);
    end
    
    subdir = [prefix, '_sdh_subject_comps'];
    
    if (exist(fullfile(resultsDir, subdir), 'dir') ~= 7)
        mkdir (resultsDir, subdir);
    end
    
    sub_comps_file = fullfile(subdir, [prefix, '_sdh_sub_', icatb_returnFileIndex(nSub), '_sess_', ...
        icatb_returnFileIndex(nSess), '.mat']);
    
    icatb_parSave(fullfile(resultsDir, sub_comps_file), {tc, ic}, {'tc', 'ic'});
    
    subjectFiles{countD} = sub_comps_file;
end

sdhInfo.subjectFiles = subjectFiles;

fname = fullfile(resultsDir, [sesInfo.userInput.prefix, '_sdh_info.mat']);
disp(['Saving parameters info in ', fname]);
save(fname, 'sdhInfo');
fprintf('\n\n');

diary('off');


totalTime = toc;

fprintf('\n');

disp('Analysis Complete');

disp(['Total time taken to complete the analysis is ', num2str(totalTime/60), ' minutes']);

diary('off');

fprintf('\n');