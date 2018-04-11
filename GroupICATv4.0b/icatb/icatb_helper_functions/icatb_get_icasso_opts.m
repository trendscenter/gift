function icasso_opts = icatb_get_icasso_opts(icasso_opts)
%% Get ICASSO Options
%

icatb_defaults;
global NUM_RUNS_GICA;

modes = {'RandInit', 'Bootstrap', 'Both'}; % All modes

val = 1;

if (exist('icasso_opts', 'var'))

    tempRuns = icasso_opts.num_ica_runs;
    selMode = icasso_opts.sel_mode;

    try
        min_cluster_size = icasso_opts.min_cluster_size;
    catch
    end

    try
        max_cluster_size = icasso_opts.max_cluster_size;
    catch
    end

    val = strmatch(lower(selMode), lower(modes), 'exact');
    if (isempty(val))
        val = 1;
    end

else

    tempRuns = NUM_RUNS_GICA;
    if (tempRuns < 2)
        tempRuns = 2;
    end

end

if (~exist('min_cluster_size', 'var'))
    min_cluster_size = ceil(0.8*tempRuns);
end

if (~exist('max_cluster_size', 'var'))
    max_cluster_size = tempRuns;
end

numParameters = 1;
inputText(numParameters).promptString = 'Select Mode';
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerString = modes;
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'mode';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = val;

numParameters = numParameters + 1;

inputText(numParameters).promptString = 'Enter no. of times ( > 1) you want ICA to be run.';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = num2str(tempRuns);
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'numRuns';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;

numParameters = numParameters + 1;

inputText(numParameters).promptString = 'Enter min cluster size to get the most stable run. Recommended option is ceil(0.8*number of runs).';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = num2str(min_cluster_size);
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'min_cluster_size';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;


numParameters = numParameters + 1;

inputText(numParameters).promptString = 'Enter max cluster size to get the most stable run. Recommended option is number of runs.';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = num2str(max_cluster_size);
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'max_cluster_size';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;

% Input dialog box
answer = icatb_inputDialog('inputtext', inputText, 'Title', 'Input parameters to determine ICA algorithm reliability', ...
    'handle_visibility', 'on', 'windowStyle', 'normal');

drawnow;

if isempty(answer)
    error('Input parameters to determine ICA algorithm reliability are not selected');
end

% Selected mode
selMode = answer{1};

if (isempty(answer{2}) && ~isnumeric(answer{2}))
    error('Number of ICA runs must be numeric');
end

% Number of ICA runs
numRuns = answer{2};

if (numRuns < 0)
    numRuns = 1;
end

min_cluster_size = answer{3};
max_cluster_size = answer{4};

icasso_opts.sel_mode = lower(selMode);
icasso_opts.num_ica_runs = numRuns;
icasso_opts.min_cluster_size = min_cluster_size;
icasso_opts.max_cluster_size = max_cluster_size;
