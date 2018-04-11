function icatb_mancovan_batch(file_name)
%% Batch for mancovan analysis
%

file_name = deblank(file_name);

addpath(fileparts(which(file_name)));

%% Read input parameters
mancovanInfo = readFile(file_name);

%% Create design matrix
mancovanInfo = create_design(mancovanInfo);

%% Run mancovan
icatb_run_mancovan(mancovanInfo, 1);


function mancovanInfo = create_design(mancovanInfo)
%% Create Design

fprintf('Creating design matrix ...\n');

mancovanInfo = icatb_mancovan_full_design(mancovanInfo, mancovanInfo.userInput.interactions);
fileN = fullfile(mancovanInfo.outputDir, [mancovanInfo.prefix, '.mat']);
icatb_save(fileN, 'mancovanInfo');

fprintf('Done\n\n');

function mancovanInfo = readFile(file_name)
%% Read file
%

fprintf('Reading input parameters for mancovan analysis \n');

inputData = icatb_eval_script(file_name);

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

mancovanInfo.userInput.ica_param_file = inputData.ica_param_file;
load( mancovanInfo.userInput.ica_param_file);

mancovanInfo.userInput.numOfSub = sesInfo.numOfSub;
mancovanInfo.userInput.numOfSess = sesInfo.numOfSess;
mancovanInfo.userInput.prefix = [sesInfo.userInput.prefix, '_mancovan'];
compFiles = sesInfo.icaOutputFiles(1).ses(1).name;
mancovanInfo.userInput.compFiles = icatb_rename_4d_file(icatb_fullFile('directory', fileparts(mancovanInfo.userInput.ica_param_file), 'files', compFiles));
mancovanInfo.userInput.numICs = sesInfo.numComp;
mancovanInfo.userInput.HInfo = sesInfo.HInfo.V(1);

%% Covariates
covariates = inputData.covariates;
cov = repmat(struct('name', '', 'value', [], 'type', 'continuous', 'transformation', ''), 1, size(covariates, 1));

for n = 1:size(covariates, 1)
    cov(n).name = covariates{n, 1};
    if (~strcmpi(covariates{n, 2}, 'continuous'))
        cov(n).type = 'categorical';
    end
    
    val = covariates{n, 3};
    
    if (isnumeric(val))
        val = num2str(val(:));
    else
        if (ischar(val) && (size(val, 1) == 1))
            val = icatb_mancovan_load_covariates(val, cov(n).type);
        end
    end
    
    val = strtrim(cellstr(val));
    
    if (length(val) ~=  mancovanInfo.userInput.numOfSub)
        error(['Covariate vector must match the no. of subjects in the analysis. Please check ', cov(n).name, ' covariates']);
    end
    
    cov(n).value = val(:)';
    clear val;
    
    try
        cov(n).transformation = covariates{n, 4};
    catch
    end
    
    if (strcmpi(cov(n).type, 'categorical'))
        cov(n).transformation = '';
    end
end

mancovanInfo.userInput.cov = cov;

interactions = [];
try
    interactions = inputData.interactions;
catch
end

mancovanInfo.userInput.interactions = interactions;

%% Component network names
comp_network_names = inputData.comp_network_names;
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