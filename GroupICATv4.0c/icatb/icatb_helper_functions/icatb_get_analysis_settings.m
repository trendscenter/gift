function varargout = icatb_get_analysis_settings(sesInfo, perfType)
%% Get group PCA options based on the performance type
%

opts = {'Maximize Performance', 'Less Memory Usage', 'User Specified Settings'};

if (nargin == 0)
    varargout{1} = opts;
    return;
end

opts = lower(opts);

if (nargin == 1)
    perfType = 1;
end

if (isnumeric(perfType))
    perfType = opts{perfType};
end

perfType = lower(perfType);

perfType = deblank(perfType);

pcaType = 'standard';
if (isfield(sesInfo.userInput, 'pcaType'))
    pcaType = lower(sesInfo.userInput.pcaType);
end

if (strcmpi(pcaType, 'svd'))
    pcaType = 'standard';
end

if (~isfield(sesInfo.userInput, 'pca_opts'))
    if (isfield(sesInfo.userInput, 'covariance_opts'))
        pca_opts_user_specified = sesInfo.userInput.covariance_opts;
    else
        pca_opts_user_specified = icatb_pca_options(pcaType);
    end
else
    pca_opts_user_specified = sesInfo.userInput.pca_opts;
end

% Fill any missing fields
pca_opts_user_specified = icatb_pca_options(pcaType, pca_opts_user_specified, 'off');

sesInfo.userInput.pca_opts = pca_opts_user_specified;

if ((strcmpi(perfType, 'maximize performance')))
    sesInfo.userInput.pca_opts.precision = 'single';
end

% Get pca options for each setting
pca_opts = icatb_mem_ica(sesInfo);

algoNames = cellstr(icatb_icaAlgorithm);
algoName = deblank(algoNames{sesInfo.userInput.algorithm});

minDimSize = 10000;
dims = [length(sesInfo.userInput.mask_ind), sum(sesInfo.userInput.diffTimePoints)];
if (sesInfo.userInput.numReductionSteps == 2)
    dims(2) = sesInfo.userInput.numOfPC1*sesInfo.userInput.numOfSub*sesInfo.userInput.numOfSess;
end

pcaTypes = cellstr(char(pca_opts.pca_type));

if (strcmpi(perfType, 'maximize performance'))
    %% Standard PCA or MPOWIT is used.
    
    
    chk = (strcmpi('standard', pcaTypes) | strcmpi('svd', pcaTypes));
    includeEVD = ones(length(pcaTypes), 1);
    if (min(dims) <= minDimSize)
        includeEVD(chk) = 1;
    else
        includeEVD(chk) = 0;
    end
    
    excludeEMPCA = ones(length(pcaTypes), 1);
    chk = (strcmpi('expectation maximization', pcaTypes) | strcmpi('empca', pcaTypes));
    excludeEMPCA(chk) = 0;
    
    includePCA = includeEVD .* excludeEMPCA;
    
    tmp_pca_opts = getBestMatch(pca_opts(includePCA == 1), {'stack_data', 'yes'});
    
elseif (strcmpi(perfType, 'less memory usage'))
    %% Mostly standard PCA without stacking data-sets is used unless the covariance matrix size
    % exceeds the stacked data size
    
    if (~strcmpi(algoName, 'iva-gl'))
        chk = (strcmpi(cellstr(char(pca_opts.pca_type)), 'mpowit') & strcmpi(cellstr(char(pca_opts.stack_data)), 'no'));
        tmp_pca_opts = pca_opts(chk);
        if (isempty(tmp_pca_opts))
            tmp_pca_opts = [getBestMatch(pca_opts, {'stack_data', 'yes'}), getBestMatch(pca_opts, {'stack_data', 'no'})];
        end
    else
        tmp_pca_opts = getBestMatch(pca_opts, {'stack_data', 'yes'});
    end
    
    [dd, inds] = min([tmp_pca_opts.max_mem]);
    tmp_pca_opts = tmp_pca_opts(inds);
    
else
    %% User settings
    
    fieldN = {'pca_type', pcaType, 'stack_data', 'yes'};
    
    for ii = 1:2:length(fieldN)
        if (isfield(pca_opts_user_specified, fieldN{ii}))
            fieldN{ii + 1} = getfield(pca_opts_user_specified, fieldN{ii});
        end
    end
    
    tmp_pca_opts = getBestMatch(pca_opts, fieldN);
    
end

try
    if (strcmpi(tmp_pca_opts.precision, 'single'))
        tmp_pca_opts.tolerance = 1e-4;
    end
catch
end


fieldN = fieldnames(pca_opts_user_specified);

for ii = 1:length(fieldN)
    if (isfield(tmp_pca_opts, fieldN{ii}))
        pca_opts_user_specified = setfield(pca_opts_user_specified, fieldN{ii}, getfield(tmp_pca_opts, fieldN{ii}));
    end
end

varargout{1} = tmp_pca_opts.max_mem;
varargout{2} = lower(tmp_pca_opts.pca_type);
varargout{3} = pca_opts_user_specified;


function tmp_pca_opts = getBestMatch(pca_opts, fieldN)
%% Get best match in terms of RAM usage

icatb_defaults;
global MAX_AVAILABLE_RAM;

% By default use 1 GB RAM
if (isempty(MAX_AVAILABLE_RAM))
    MAX_AVAILABLE_RAM = 1;
end

for ii = 1:2:length(fieldN)
    tmpS = double(strcmpi(cellstr(str2mat(pca_opts.(fieldN{ii}))), fieldN{ii + 1}));
    if (ii == 1)
        inds = tmpS;
    else
        inds = inds.*tmpS;
    end
end

tmp_pca_opts = pca_opts(inds == 1);
if (isempty(tmp_pca_opts))
    tmp_pca_opts = pca_opts(end);
end

[mems, inds] = sort([tmp_pca_opts.max_mem]);
inds = inds(end:-1:1);
tmp_pca_opts = tmp_pca_opts(inds);

check = find([tmp_pca_opts.max_mem] < MAX_AVAILABLE_RAM);

if (isempty(check))
    tmp_pca_opts = tmp_pca_opts(end);
else
    tmp_pca_opts = tmp_pca_opts(check(1));
end



