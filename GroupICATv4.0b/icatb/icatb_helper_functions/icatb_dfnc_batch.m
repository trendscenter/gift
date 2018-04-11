function icatb_dfnc_batch(file_name)
%% Batch for dfnc analysis
%

file_name = deblank(file_name);

addpath(fileparts(which(file_name)));

%% Read input parameters
dfncInfo = readFile(file_name);

%% Run dFNC
dfncInfo = icatb_run_dfnc(dfncInfo);

%% Postprocess and display
if (isfield(dfncInfo, 'postprocess'))
    dfncInfo = icatb_post_process_dfnc(dfncInfo);
    %% Display
    if (dfncInfo.postprocess.display_results)
        icatb_dfnc_results_html(dfncInfo);
    end
    
end


function dfncInfo = readFile(file_name)
%% Read file
%

fprintf('Reading input parameters for dFNC analysis \n');

inputData = icatb_eval_script(file_name);

dfncInfo.userInput.outputDir = inputData.outputDir;

if (exist(dfncInfo.userInput.outputDir, 'dir') ~= 7)
    mkdir(dfncInfo.userInput.outputDir);
end

dfncInfo.userInput.ica_param_file = inputData.ica_param_file;
load( dfncInfo.userInput.ica_param_file);

dfncInfo.userInput.numOfSub = sesInfo.numOfSub;
dfncInfo.userInput.numOfSess = sesInfo.numOfSess;
dfncInfo.userInput.prefix = [sesInfo.userInput.prefix, '_dfnc'];
compFiles = sesInfo.icaOutputFiles(1).ses(1).name;
dfncInfo.userInput.compFiles = icatb_rename_4d_file(icatb_fullFile('directory', fileparts(dfncInfo.userInput.ica_param_file), 'files', compFiles));
dfncInfo.userInput.numICs = sesInfo.numComp;
dfncInfo.userInput.HInfo = sesInfo.HInfo.V(1);

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

dfncInfo.userInput.comp = comp;


%% TR
if (isfield(sesInfo, 'TR'))
    dfncInfo.userInput.TR = sesInfo.TR;
else
    dfncInfo.userInput.TR = inputData.TR;
end

covInfo.numOfDataSets = dfncInfo.userInput.numOfSub*dfncInfo.userInput.numOfSess;

default_params = icatb_dfnc_options('covInfo', covInfo);

feature_params = default_params;

try
    fparams = inputData.dfnc_params;
catch
    fparams = [];
end

if (~isempty(fparams))
    
    %% Preprocessing
    tags = cellstr( char(feature_params.inputParameters(1).options.tag));
    fieldsToCheck = {'tc_detrend', 'tc_despike', 'tc_filter'};
    feature_params.inputParameters(1) = fillValuesControls(feature_params.inputParameters(1), fieldsToCheck, fparams, tags);
    
    chk = strmatch('tc_covariates', tags, 'exact');
    
    try
        if (~isempty(chk))
            if (~isempty(fparams.tc_covariates.filesList))
                feature_params.inputParameters(1).options(chk).value = 2;
                
                covInfo.filesList = fparams.tc_covariates.filesList;
                covInfo.file_numbers = [];
                try
                    covInfo.file_numbers = fparams.tc_covariates.file_numbers;
                catch
                end
                
                feature_params.inputParameters(1).options(chk).userdata = covInfo;
                
            end
        end
    catch
    end
    
    %% dfnc
    tags = cellstr( char(feature_params.inputParameters(2).options.tag));
    fieldsToCheck = {'method', 'wsize', 'window_alpha', 'num_repetitions'};
    
    feature_params.inputParameters(2) = fillValuesControls(feature_params.inputParameters(2), fieldsToCheck, fparams, tags);
    
end

out = icatb_OptionsWindow(feature_params.inputParameters, feature_params.defaults, 'off', 'title', 'Feature Options');
feature_params.inputParameters = out.inputParameters;
feature_params.final = out.results;
dfncInfo.userInput.feature_params = feature_params;

% read post-processing
try
    dfncInfo.postprocess = inputData.postprocess;
catch
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