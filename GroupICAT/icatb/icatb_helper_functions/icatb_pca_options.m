function varargout = icatb_pca_options(pcaType, pcaOpts, handle_visibility)
%% PCA Options
%
% Inputs:
% 1. pcaType - Options are 'standard' or 'expectation maximization'
% 2. pcaOpts - PCA Options Data structure
% 3. handle_visibility - Options are 'on' or 'off'
%
% Outputs:
% pcaOpts - PCA Options
%

pcaTypes = {'Standard', 'Expectation Maximization', 'SVD', 'MPOWIT', 'STP'};
if (nargin == 0)
    varargout{1} = pcaTypes;
    return;
end

if (isnumeric(pcaType))
    pcaType = pcaTypes{pcaType};
end

switch (lower(pcaType))
    
    case {'standard', 'evd'}
        %% Standard PCA
        
        % Stack Data
        numParameters = 1;
        
        inputText(numParameters).promptString = 'Do You Want To Stack Datasets?';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'Yes', 'No'};
        inputText(numParameters).value = 1;
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'stack_data';
        inputText(numParameters).enable = 'on';
        
        % Storage Type
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select Matrix Storage Type';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'Full', 'Packed'};
        inputText(numParameters).value = 1;
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'storage';
        inputText(numParameters).enable = 'on';
        
        % Precision
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select Precision';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'Double', 'Single'};
        inputText(numParameters).value = 1;
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'precision';
        inputText(numParameters).enable = 'on';
        
        % Eigen Solver Type
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select Eigen Solver Type';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'Selective', 'All'};
        inputText(numParameters).value = 1;
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'eig_solver';
        inputText(numParameters).enable = 'on';
        
    case 'svd'
        %% SVD
        
        % Precision
        numParameters = 1;
        
        inputText(numParameters).promptString = 'Select Precision';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'Double', 'Single'};
        inputText(numParameters).value = 1;
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'precision';
        inputText(numParameters).enable = 'on';
        
        % SVD Solver Type
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select SVD Solver Type';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'Selective', 'All'};
        inputText(numParameters).value = 1;
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'solver';
        inputText(numParameters).enable = 'on';
        
        
    case {'empca', 'expectation maximization'}
        %% Expectation Maximization
        
        % Stack Data
        numParameters = 1;
        
        inputText(numParameters).promptString = 'Do You Want To Stack Datasets?';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'Yes', 'No'};
        inputText(numParameters).value = 1;
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'stack_data';
        inputText(numParameters).enable = 'on';
        
        % Precision
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select Precision';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'Double', 'Single'};
        inputText(numParameters).value = 1;
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'precision';
        inputText(numParameters).enable = 'on';
        
        % Select stopping tolerance
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Enter Stopping Tolerance';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(1e-4);
        inputText(numParameters).value = 1;
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'tolerance';
        inputText(numParameters).enable = 'on';
        
        % Enter max iterations
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Enter Max No. Of Iterations';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(1000);
        inputText(numParameters).value = 1;
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'max_iter';
        inputText(numParameters).enable = 'on';
        
    case 'mpowit'
        %% MPOWIT
        
        % Stack Data
        numParameters = 1;
        
        inputText(numParameters).promptString = 'Do You Want To Stack Datasets?';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'Yes', 'No'};
        inputText(numParameters).value = 1;
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'stack_data';
        inputText(numParameters).enable = 'on';
        
        % Precision
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select Precision';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'Double', 'Single'};
        inputText(numParameters).value = 1;
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'precision';
        inputText(numParameters).enable = 'on';
        
        % Select stopping tolerance
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Enter Stopping Tolerance';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(1e-5);
        inputText(numParameters).value = 1;
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'tolerance';
        inputText(numParameters).enable = 'on';
        
        % Enter max iterations
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Enter Max No. Of Iterations';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(1000);
        inputText(numParameters).value = 1;
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'max_iter';
        inputText(numParameters).enable = 'on';
        
        % Enter max iterations
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Enter block multiplier';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(10);
        inputText(numParameters).value = 1;
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'block_multiplier';
        inputText(numParameters).enable = 'on';
        
        
    case 'stp'
        %% Subsampled time pca
        
        % Precision
        numParameters = 1;
        
        inputText(numParameters).promptString = 'Select Precision';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'Double', 'Single'};
        inputText(numParameters).value = 1;
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'precision';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        inputText(numParameters).promptString = 'Enter intermediate components to retain';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(500);
        inputText(numParameters).value = 1;
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'num_comp';
        inputText(numParameters).enable = 'on';
        
        % number of subjects in each group
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select number of subjects in each group';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(10);
        inputText(numParameters).value = 1;
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'numGroups';
        inputText(numParameters).enable = 'on';
        
end

%% Default PCA Options
fieldN = cell(1, 2*length(inputText));
for i = 1:length(inputText)
    answerString = inputText(i).answerString;
    
    if (strcmpi(inputText(i).uiType, 'edit'))
        val = answerString;
    else
        tmp = inputText(i).value;
        val = answerString{tmp};
    end
    
    if (strcmpi(inputText(i).dataType, 'numeric'))
        val = str2num(val);
    else
        val = lower(val);
    end
    
    fieldN{2*i-1} = inputText(i).tag;
    fieldN{2*i} = val;
    
end

DefaultPCAOpts = struct(fieldN{:});

if (nargin < 2)
    varargout{1} = DefaultPCAOpts;
    return;
end

if (~isstruct(pcaOpts))
    error('PCA options must be in a data structure');
end

if (~exist('handle_visibility', 'var'))
    handle_visibility = 'on';
end

allTags = cellstr(str2mat(inputText.tag));
field_names = fieldnames(pcaOpts);

[common_fields, ia] = intersect(allTags, field_names);

if (isempty(common_fields))
    error('Please check the PCA Options data structure you have passed');
end

ia = ia(:)';

%% Update uicontrols with the selected values
for i = ia
    tmp = getfield(pcaOpts, inputText(i).tag);
    if (strcmpi(inputText(i).uiType, 'popup'))
        try
            val = strmatch(lower(tmp), lower(inputText(i).answerString), 'exact');
        catch
            val = 1;
        end
        inputText(i).value = val;
    else
        if (isnumeric(tmp))
            tmp = num2str(tmp);
        end
        inputText(i).answerString = tmp;
    end
end

%% Get Answers
answers = getAnswers(inputText, pcaType, handle_visibility);

if (isempty(answers))
    answers = getAnswers(inputText, pcaType, 'off');
end

drawnow;

fN = cell(1, 2*length(allTags));
char_inds  = strmatch('string', lower(cellstr(str2mat(inputText.dataType))), 'exact');
answers(char_inds) = lower(answers(char_inds));
fN(1:2:end) = allTags;
fN(2:2:end) = answers;
varargout{1} = struct(fN{:});


function answers = getAnswers(inputText, pcaType, handle_visibility)
%% Get Answers
%

if strcmpi(handle_visibility, 'on')
    answers = icatb_inputDialog('inputtext', inputText, 'Title', ['Select options for ', upper(pcaType), ' PCA in the Data Reduction Stage'], 'windowStyle', 'modal');
else
    answers = cell(1, length(inputText));
    for i = 1:length(inputText)
        answerString = lower(inputText(i).answerString);
        if (strcmpi(inputText(i).uiType, 'popup'))
            answerString = answerString{inputText(i).value};
        end
        if (strcmpi(inputText(i).dataType, 'numeric'))
            answerString = str2num(answerString);
        end
        answers{i} = answerString;
    end
end