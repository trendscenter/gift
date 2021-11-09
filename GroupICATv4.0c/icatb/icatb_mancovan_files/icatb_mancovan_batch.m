function icatb_mancovan_batch(file_name)
%% Batch for mancovan analysis
%

file_name = deblank(file_name);

addpath(fileparts(which(file_name)));

%% Read input parameters
mancovanInfo = readFile(file_name);

skipCalc = 0;
if (isfield(mancovanInfo, 'skipCalc'))
    skipCalc = mancovanInfo.skipCalc;
end

if (~skipCalc)
    
    %% Create design matrix
    mancovanInfo = create_design(mancovanInfo);
    
    %% Run mancovan
    mancovanInfo = icatb_run_mancovan(mancovanInfo, 1);
    
    if (mancovanInfo.write_stats_info)
        icatb_mancovan_agg(mancovanInfo);
    end
    
else
    
    display_results = [];
    if (isfield(mancovanInfo, 'display'))
        display_results = mancovanInfo.display;
    end
    
    outDir = pwd;
    try
        outDir = mancovanInfo.outputDir;
    catch
        
    end
    
    mancovanInfo = icatb_mancovan_reduce(mancovanInfo.files, outDir, display_results);
    
    
end

%% Display results
display_results = 0;
try
    mancovanInfo.display.t_threshold;
    display_results = 1;
catch
end

if (display_results)
    mancovan_results_summary(mancovanInfo);
end




function mancovanInfo = create_design(mancovanInfo)
%% Create Design

fprintf('Creating design matrix ...\n');

mancovanInfo = icatb_mancovan_full_design(mancovanInfo, mancovanInfo.userInput.interactions);

if (isfield(mancovanInfo.userInput, 'time'))
    mancovanInfo.time = mancovanInfo.userInput.time;
    mancovanInfo.X = mancovanInfo.X(mancovanInfo.time.subjects{1}, :);
    mancovanInfo.userInput = rmfield(mancovanInfo.userInput, 'time');
end

fileN = fullfile(mancovanInfo.outputDir, [mancovanInfo.prefix, '.mat']);
icatb_save(fileN, 'mancovanInfo');

fprintf('Done\n\n');

function mancovanInfo = readFile(file_name)
%% Read file
%

icatb_defaults;
global MANCOVA_DEFAULTS;

fprintf('Reading input parameters for mancovan analysis \n');

inputData = icatb_eval_script(file_name);

%% distributed mancova
write_stats_info = 0;
try
    write_stats_info = MANCOVA_DEFAULTS.write_stats_info;
catch
end

if (isfield(inputData, 'write_stats_info'))
    write_stats_info = inputData.write_stats_info;
end

mancovanInfo.write_stats_info = write_stats_info;

ica_param_file = inputData.ica_param_file;

if (iscell(ica_param_file))
    variablesIn = whos('-file', ica_param_file{1});
    if (~isempty(strmatch('uni_results_info', cellstr(char(variablesIn.name)), 'exact')))
        mancovanInfo.skipCalc = 1;
        mancovanInfo.files = inputData.ica_param_file;
        
        try
            mancovanInfo.display = inputData.display;
        catch
        end
        
        try
            mancovanInfo.userInput.outputDir = inputData.outputDir;
            mancovanInfo.outputDir = mancovanInfo.userInput.outputDir;
        catch
        end
        
        return;
    end
end




%% Store some information
features = {'spatial maps', 'timecourses spectra', 'fnc correlations', 'fnc correlations (lag)'};
inputData.features = lower(cellstr(inputData.features));
[dd, ia] = intersect(inputData.features, lower(features));
if (isempty(dd))
    error('Please check features variable');
end

ia = sort(ia);
inputData.features = inputData.features(ia);

mancovanInfo.userInput.features = inputData.features;

mancovanInfo.userInput.outputDir = inputData.outputDir;

if (exist(mancovanInfo.userInput.outputDir, 'dir') ~= 7)
    mkdir(mancovanInfo.userInput.outputDir);
end

comp_network_names = inputData.comp_network_names;

if (iscell(ica_param_file))
    if (~isfield(inputData, 'merge_type'))
        inputData.merge_type = 'stack_subjects';
    end
    [ica_param_file, comp_network_names] = icatb_merge_analyses(ica_param_file, 'merge_type', inputData.merge_type, 'outputDir', ...
        mancovanInfo.userInput.outputDir, 'comp_network_names', comp_network_names);
end

mancovanInfo.userInput.ica_param_file = ica_param_file;
load( mancovanInfo.userInput.ica_param_file);


if (~isfield(inputData, 'avg_runs_info') || isempty(inputData.avg_runs_info))
    avg_runs_info = cellfun(@(subj, sess)(subj - 1).*sesInfo.numOfSess + 1:subj.*sesInfo.numOfSess, num2cell(1:sesInfo.numOfSub), ...
        repmat({sesInfo.numOfSess}, 1, sesInfo.numOfSub), 'UniformOutput', false);
else
    avg_runs_info = inputData.avg_runs_info;
end
mancovanInfo.userInput.avg_runs_info = avg_runs_info;

mancovanInfo.userInput.numOfSub = length(avg_runs_info);
mancovanInfo.userInput.numOfSess = sesInfo.numOfSess;
mancovanInfo.userInput.prefix = [sesInfo.userInput.prefix, '_mancovan'];
compFiles = sesInfo.icaOutputFiles(1).ses(1).name;
mancovanInfo.userInput.compFiles = icatb_rename_4d_file(icatb_fullFile('directory', fileparts(mancovanInfo.userInput.ica_param_file), 'files', compFiles));
mancovanInfo.userInput.numICs = sesInfo.numComp;
mancovanInfo.userInput.HInfo = sesInfo.HInfo.V(1);


if (isfield(inputData, 'univariate_tests') && ~isempty(inputData.univariate_tests))
    
    mancovanInfo.userInput.univariate_tests = inputData.univariate_tests;
    mancovanInfo.userInput.numOfPCs = [2,2,2];
else
    
    %% Number of PCs
    numOfPCs = inputData.numOfPCs;
    
    if (length(numOfPCs) == 1)
        numOfPCs = ones(1, length(mancovanInfo.userInput.features))*numOfPCs(1);
    end
    
    if (length(numOfPCs) > length(mancovanInfo.userInput.features))
        numOfPCs = numOfPCs(1:length(mancovanInfo.userInput.features));
    end
    
    if (length(numOfPCs) ~= length(mancovanInfo.userInput.features))
        error('Please check variable numOfPCs. The length of numOfPCs must match the no. of features');
    end
    
    mancovanInfo.userInput.numOfPCs = numOfPCs;
    
    if (any(mancovanInfo.userInput.numOfPCs > mancovanInfo.userInput.numOfSub))
        error(['One/more PCs exceed the no. of subjects (', num2str(mancovanInfo.userInput.numOfSub), ')']);
    end
    
    
    % Set estimation to none
    doEst = 0;
    try
        doEst = inputData.doEstimation;
    catch
    end
    mancovanInfo.userInput.doEstimation = doEst;
    
    
end

univ_test_name = '';

try
    univ_test_name = mancovanInfo.userInput.univariate_tests{1, 1};
catch
end

%% Covariates
if ~(strcmpi(univ_test_name, 'ttest') || strcmpi(univ_test_name, 'ttest2'))
    
    covariates = inputData.covariates;
    covM = repmat(struct('name', '', 'value', [], 'type', 'continuous', 'transformation', ''), 1, size(covariates, 1));
    
    for n = 1:size(covariates, 1)
        covM(n).name = covariates{n, 1};
        if (~strcmpi(covariates{n, 2}, 'continuous'))
            covM(n).type = 'categorical';
        end
        
        val = covariates{n, 3};
        
        if (isnumeric(val))
            val = num2str(val(:));
        else
            if (ischar(val) && (size(val, 1) == 1))
                val = icatb_mancovan_load_covariates(val, covM(n).type);
            end
        end
        
        val = strtrim(cellstr(val));
        
        if (length(val) ~=  mancovanInfo.userInput.numOfSub)
            error(['Covariate vector must match the no. of subjects in the analysis. Please check ', covM(n).name, ' covariates']);
        end
        
        covM(n).value = val(:)';
        clear val;
        
        try
            covM(n).transformation = covariates{n, 4};
        catch
        end
        
        if (strcmpi(covM(n).type, 'categorical'))
            covM(n).transformation = '';
        end
    end
    
else
    
    covM = [];
    
end

mancovanInfo.userInput.cov = covM;


interactions = [];
if (isempty(univ_test_name))
    try
        interactions = inputData.interactions;
    catch
    end
else
    if ~(strcmpi(univ_test_name, 'ttest') || strcmpi(univ_test_name, 'ttest2'))
        covariateNames = cellstr(char(mancovanInfo.userInput.cov.name));
        temp_cell_a = mancovanInfo.userInput.univariate_tests(:, 1);
        all_test_names = temp_cell_a;
        for nTest = 1:size(mancovanInfo.userInput.univariate_tests, 1)
            if ~isempty(mancovanInfo.userInput.univariate_tests{nTest, 2})
                tmp_cell_b = mancovanInfo.userInput.univariate_tests{nTest, 2};
                if ~iscell(tmp_cell_b)
                    tmp_cell_b = cellstr(tmp_cell_b);
                end
                all_test_names = cat(1, all_test_names, tmp_cell_b');
            end
        end
        
        all_test_names = unique(all_test_names);
        good_test_inds = icatb_good_cells(regexpi(all_test_names, '\_X\_'));
        
        all_test_names = all_test_names(good_test_inds);
        all_interactions = cell(1, length(all_test_names));
        for nT = 1:length(all_test_names)
            
            dd = regexpi(all_test_names{nT}, '\_X\_','split');
            inda = find(strcmpi(covariateNames, dd{1}) == 1);
            indb = find(strcmpi(covariateNames, dd{2}) == 1);
            
            tmp = [inda, indb];
            
            all_interactions{nT} = tmp;
            
        end
        
        interactions = cat(1, all_interactions{:});
        
        
    end
    
end

mancovanInfo.userInput.interactions = interactions;

%% Component network names
%comp_network_names = inputData.comp_network_names;
comp = repmat(struct('name', '', 'value', []), 1, size(comp_network_names, 1));

for n = 1:size(comp_network_names, 1)
    comp(n).name = comp_network_names{n, 1};
    val = comp_network_names{n, 2};
    if (ischar(val))
        val = load(val, '-ascii');
    end
    val = val(:)';
    comp(n).value = val;
end

value = [comp.value];

if (length(value) ~= length(unique(value)))
    error('There are duplicate entries of component/components. Please check variable comp_network_names');
end

if (max(value) > sesInfo.numComp)
    error('Max value of components is greater than the no. of components present');
end

mancovanInfo.userInput.comp = comp;

%% P threshold
p_threshold = 0.01;
try
    p_threshold = inputData.p_threshold;
catch
end
mancovanInfo.userInput.p_threshold = p_threshold;


%% TR
mancovanInfo.userInput.TR = inputData.TR;

default_params = icatb_mancovan_feature_options('tr', mancovanInfo.userInput.TR, 'mask_dims', mancovanInfo.userInput.HInfo(1).dim(1:3));
feature_params = default_params;

try
    fparams = inputData.feature_params;
catch
    fparams = [];
end

if (~isempty(fparams))
    
    %% Spatial maps
    fieldsToCheck = {'sm_center', 'stat_threshold_maps', 'z_threshold_maps'};
    tags = cellstr( char(feature_params.inputParameters(1).options.tag));
    
    chk = strmatch('sm_mask', tags, 'exact');
    
    try
        if (~isempty(chk))
            
            if (~isempty(fparams.sm_mask))
                
                feature_params.inputParameters(1).options(chk).value = 2;
                feature_params.inputParameters(1).options(chk).userdata = fparams.sm_mask;
                
            end
            
        end
    catch
    end
    
    
    tags = cellstr( char(feature_params.inputParameters(1).options.tag));
    feature_params.inputParameters(1) = fillValuesControls(feature_params.inputParameters(1), fieldsToCheck, fparams, tags);
    
    %% Spectra
    tags = cellstr( char(feature_params.inputParameters(2).options.tag));
    fieldsToCheck = {'spectra_detrend', 'spectra_tapers', 'spectra_sampling_freq', 'spectra_freq_band', ...
        'spectra_normalize_subs', 'spectra_transform'};
    
    feature_params.inputParameters(2) = fillValuesControls(feature_params.inputParameters(2), fieldsToCheck, fparams, tags);
    
    
    %% FNC
    tags = cellstr( char(feature_params.inputParameters(3).options.tag));
    fieldsToCheck = {'fnc_tc_detrend', 'fnc_tc_despike', 'fnc_tc_filter'};
    feature_params.inputParameters(3) = fillValuesControls(feature_params.inputParameters(3), fieldsToCheck, fparams, tags);
    
end

out = icatb_OptionsWindow(feature_params.inputParameters, feature_params.defaults, 'off', 'title', 'Feature Options');
feature_params.inputParameters = out.inputParameters;
feature_params.final = out.results;
mancovanInfo.userInput.feature_params = feature_params;

if (isfield(inputData, 'display'))
    fieldNames = cellstr(fieldnames(inputData.display));
    for nF = 1:length(fieldNames)
        mancovanInfo.display.(fieldNames{nF}) = inputData.display.(fieldNames{nF});
    end
end

if (isfield(inputData, 'time'))
    try
        if (~isempty(inputData.time) && ~isempty(inputData.time.subjects))
            mancovanInfo.userInput.time = inputData.time;
        end
    catch
    end
end

fprintf('Done\n\n');


function inputParameters = fillValuesControls(inputParameters, fieldsToCheck, fparams, tags)
%% Fill values in feature defaults window
%

for nF = 1:length(fieldsToCheck)
    
    chk = strmatch(fieldsToCheck{nF}, tags, 'exact');
    
    try
        if (~isempty(chk))
            
            val = fparams.(fieldsToCheck{nF});
            
            if (isnumeric(val))
                val = num2str(val);
            end
            
            if (strcmpi(inputParameters.options(chk).uiType, 'popup') || strcmpi(inputParameters.options(chk).uiType, 'popupmenu'))
                chkInd = strmatch(lower(val), lower(cellstr(inputParameters.options(chk).answerString)), 'exact');
                if (isempty(chkInd))
                    chkInd = 1;
                end
                inputParameters.options(chk).value = chkInd;
            else
                inputParameters.options(chk).answerString = val;
            end
        end
    catch
    end
    
end

function mancovan_results_summary(mancovanInfo)
% display results

if (~isdeployed)
    
    outDir = fullfile(mancovanInfo.outputDir, [mancovanInfo.prefix, '_results_summary']);
    opts.outputDir = outDir;
    opts.showCode = false;
    opts.useNewFigure = false;
    opts.format = 'html';
    opts.createThumbnail = true;
    if (strcmpi(opts.format, 'pdf'))
        opts.useNewFigure = false;
    end
    assignin('base', 'mancovanInfo', mancovanInfo);
    opts.codeToEvaluate = 'icatb_mancovan_results_summary(mancovanInfo);';
    disp('Generating reults summary. Please wait ....');
    drawnow;
    publish('icatb_mancovan_results_summary', opts);
    close all;
    
    if (strcmpi(opts.format, 'html'))
        icatb_openHTMLHelpFile(fullfile(outDir, 'icatb_mancovan_results_summary.html'));
    else
        open(fullfile(outDir, 'icatb_mancovan_results_summary.pdf'));
    end
    
    disp('Done');
    
else
    
    icatb_mancovan_results_summary(mancovanInfo);
    
end