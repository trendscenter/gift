function varargout = icatb_parMST(param_file, algorithmName, num_ica_runs, outName, data)
%% Run MST in parallel (Sessions)
%

if (ischar(param_file))
    load(param_file);
else
    sesInfo = param_file;
end

if (~exist('sesInfo', 'var'))
    error('File selected is not a valid parameter file');
end

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


verbose = (nargout == 1);

%% MST
if (~exist('data', 'var'))
    data = icatb_getDataForICA(sesInfo, algorithmName);
end

if (verbose)
    fprintf('\n');
    disp(['Number of times ', algorithmName, ' will run is ', num2str(num_ica_runs)]);
end

icasigR = cell(1, num_ica_runs);
for nR = 1:length(icasigR)
    if (verbose)
        if (nR == 1)
            fprintf('\n');
            disp(['Run ', num2str(nR), ' / ', num2str(num_ica_runs)]);
        end
    end
    [dd1, dd2, dd3, icasigR{nR}]  = icatb_icaAlgorithm(algorithmName, data, ICA_Options);
    fprintf('\n');
end

if (nargout == 0)
    icatb_save(outName, 'icasigR');
    varargout = {};
else
    varargout{1} = icasigR;
end