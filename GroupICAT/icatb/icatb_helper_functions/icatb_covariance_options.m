function covOpts = icatb_covariance_options(covOpts, handle_visibility)
%% Covariance matrix options
%
% Inputs:
% 1. covOpts - Covariance matrix options
% 2. handle_visibility - 'on' or 'off'
%
% Outputs:
% 1. covOpts - Covariance matrix options
%

if (~exist('covOpts', 'var') || isempty(covOpts))
    covOpts = struct('stack_data', 'yes', 'storage', 'full', 'precision', 'double', 'eig_solver', 'selective');
    return;
end

if (~isstruct(covOpts))
    error('Covariance matrix options must be in a data structure');
end

if (~exist('handle_visibility', 'var'))
    handle_visibility = 'on';
end

%% Stack Data
numParameters = 1;

inputText(numParameters).promptString = 'Do You Want To Stack Datasets?';
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerString = {'Yes', 'No'};
inputText(numParameters).value = 1;
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'stack_data';
inputText(numParameters).enable = 'on';

%% Storage Type
numParameters = numParameters + 1;

inputText(numParameters).promptString = 'Select Matrix Storage Type';
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerString = {'Full', 'Packed'};
inputText(numParameters).value = 1;
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'storage';
inputText(numParameters).enable = 'on';

%% Precision
numParameters = numParameters + 1;

inputText(numParameters).promptString = 'Select Precision';
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerString = {'Double', 'Single'};
inputText(numParameters).value = 1;
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'precision';
inputText(numParameters).enable = 'on';

%% Eigen Solver Type
numParameters = numParameters + 1;

inputText(numParameters).promptString = 'Select Eigen Solver Type';
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerString = {'Selective', 'All'};
inputText(numParameters).value = 1;
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'eig_solver';
inputText(numParameters).enable = 'on';

allTags = lower(cellstr(str2mat(inputText.tag)));
field_names = fieldnames(covOpts);

%% Loop over field names
for nF = 1:length(field_names)
    cF = lower(field_names{nF});
    ind = strmatch(cF, allTags, 'exact');
    ind = ind(1);
    cVal = getfield(covOpts, cF);
    if (~isnumeric(cVal))
        cVal = lower(cVal);
        match_index = strmatch(cVal, lower(cellstr(str2mat(inputText(ind).answerString))), 'exact');
    else
        match_index = cVal;
    end
    inputText(ind(1)).value = match_index;
end
%% End of loop over field names

%% Input dialog box
answers = icatb_inputDialog('inputtext', inputText, 'Title', 'Select Options For Covariance Matrix in Data Reduction Stage', 'handle_visibility',  handle_visibility, ...
    'windowStyle', 'modal');

drawnow;

if ~isempty(answers)
    fN = cell(1, 2*length(allTags));
    fN(1:2:end) = allTags;
    fN(2:2:end) = answers;
    fN = lower(fN);
    covOpts = struct(fN{:});
end