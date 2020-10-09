function icatb_icasso(param_file, selMode, numRuns)
%% ICASSO plugin for GIFT. ICASSO centrotype results are used and stored in
% ica.mat which will be used later in back reconstruction, scaling
% components steps.
%
% Inputs:
% 1. param_file - ICA parameter file
% 2. selMode - Selection mode. Options are 'randinit', 'bootstrap' and
% 'both'.
% 3. numRuns - Number of ICA runs.
%

% Defaults
icatb_defaults;

global PARAMETER_INFO_MAT_FILE;
global NUM_RUNS_GICA;

giftPath = fileparts(which('gift.m'));
if ~isempty(giftPath)
    addpath(fullfile(giftPath, 'toolbox', 'icasso122'));
end

% load parameters file
filterP = ['*', PARAMETER_INFO_MAT_FILE, '*.mat'];

if ~exist('param_file', 'var') || isempty(param_file)
    param_file = icatb_selectEntry('typeEntity', 'file', 'title', 'Select a valid parameter file', 'filter', filterP);
end

drawnow;

if (isempty(param_file))
    error('Parameter file is not selected');
end

outputDir = fileparts(param_file);

if (isempty(outputDir))
    outputDir = pwd;
end

cd(outputDir);

load(param_file);

if ~exist('sesInfo', 'var')
    error('Not a valid parameter file');
else
    disp('Parameters file succesfully loaded');
end

icaAlgo = icatb_icaAlgorithm; % available ICA algorithms

algoVal = sesInfo.userInput.algorithm; % algorithm index

% selected ICA algorithm
algorithmName = deblank(icaAlgo(algoVal, :));

if strcmpi(algorithmName, 'moo-icar')
    algorithmName = 'gig-ica';
end

if (strcmpi(algorithmName, 'gig-ica') || strcmpi(algorithmName, 'constrained ica (spatial)'))
    error(['ICASSO option is not available for algorithm ', algorithmName]);
end

%% output directory
sesInfo.outputDir = outputDir; % set the output directory
sesInfo.userInput.pwd = outputDir;

% Check selection mode and number of ICA runs
if (nargin < 2)
    
    if (isfield(sesInfo.userInput, 'icasso_opts'))
        icasso_opts = icatb_get_icasso_opts(sesInfo.userInput.icasso_opts);
    else
        icasso_opts = icatb_get_icasso_opts;
    end
    
else
    
    if (~strcmpi(selMode, 'randinit') && ~strcmpi(selMode, 'bootstrap') && ~strcmpi(selMode, 'both'))
        error('Please provide a valid option for selMode. Valid option must be in {randinit, bootstrap, both}');
    end
    
    if (~exist('numRuns', 'var'))
        numRuns = NUM_RUNS_GICA;
    end
    
    if (numRuns < 2)
        error('Error:ICASSO', ['You need to run ICA algorithm atleast two times inorder to use ICASSO.', ...
            '\nPlease check variable numRuns.']);
    end
    
    icasso_opts.sel_mode = lower(selMode);
    icasso_opts.num_ica_runs = numRuns;
    
end
% End for checking selection mode and number of ICA runs

%% Set analysis type as ICASSO
sesInfo.userInput.which_analysis = 2;
sesInfo.userInput.icasso_opts = icasso_opts;
sesInfo.which_analysis = sesInfo.userInput.which_analysis;
sesInfo.icasso_opts = sesInfo.userInput.icasso_opts;


%% Run parameter initialization if needed
if (~sesInfo.isInitialized)
    sesInfo = icatb_runAnalysis(sesInfo, 2);
end

%% Run data reduction if needed

% Number of components
numOfIC = sesInfo.numComp;

% Load last PCA reduction step file
pcain = [sesInfo.data_reduction_mat_file, num2str(sesInfo.numReductionSteps), '-', num2str(1), '.mat'];

if ~exist(fullfile(outputDir, pcain), 'file')
    sesInfo = icatb_runAnalysis(sesInfo, 3);
end

load(fullfile(outputDir, pcain));

if (size(pcasig, 2) ~= numOfIC)
    sesInfo = icatb_runAnalysis(sesInfo, 3);
    load(fullfile(outputDir, pcain));
end

%% Run ICA
sesInfo = icatb_runAnalysis(sesInfo, 4);

%% Run back reconstruction
sesInfo = icatb_runAnalysis(sesInfo, 5);

%% Run scaling components
sesInfo = icatb_runAnalysis(sesInfo, 6);

%% Run group stats
icatb_runAnalysis(sesInfo, 7);

