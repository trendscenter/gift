function icatb_parIcassoEst(param_file, sel_mode, num_ica_runs, outName)
%% Run ICASSO in parallel (Sessions)
%

%% Initialise variables

load(param_file);

ICA_Options = {};
% get ica options
if (isfield(sesInfo, 'ICA_Options'))
    ICA_Options = sesInfo.ICA_Options;
else
    if (isfield(sesInfo.userInput, 'ICA_Options'))
        ICA_Options = sesInfo.userInput.ICA_Options;
    end
end

% convert to cell
if isempty(ICA_Options)
    if ~iscell(ICA_Options)
        ICA_Options = {};
    end
end

icaAlgo = icatb_icaAlgorithm;
algoVal = sesInfo.algorithm; % algorithm index
% selected ICA algorithm
algorithmName = deblank(icaAlgo(algoVal, :));

%% ICASSO
data = icatb_getDataForICA(sesInfo, algorithmName);
numOfIC = size(data, 1);

% PCA
[V, Lambda] = icatb_v_pca(data, 1, numOfIC, 0, 'transpose', 'yes');
% Whiten matrix
[w, White, deWhite] = icatb_v_whiten(data, V, Lambda, 'transpose');
clear V Lambda;

sR = icatb_icassoEst(sel_mode, data, num_ica_runs, 'numOfPC', numOfIC, 'algoIndex', sesInfo.algorithm, ...
    'dewhiteM', deWhite, 'whiteM', White, 'whitesig', w, 'icaOptions', ICA_Options);

%% Save info
icatb_save(outName, 'sR');