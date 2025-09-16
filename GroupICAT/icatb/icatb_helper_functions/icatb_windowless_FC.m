function icatb_windowless_FC(param_file, num_dict, varargin)
%% Use windowless approach to compute the number of functional connectivity states
%
%
% Inputs:
% 1. param_file - ICA parameter file
% 2. nnum_dict - NUmber of dictionary elements
%
% Outputs:
% Dictionary - Functional connectivity states
%


icatb_defaults;
global PARAMETER_INFO_MAT_FILE;

filterP = ['*', PARAMETER_INFO_MAT_FILE, '*.mat'];

verbose = 1;
max_iter = 75;
for n = 1:2:length(varargin)
    if (strcmpi(varargin{n}, 'verbose'))
        verbose = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'max_iter'))
        max_iter = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'outputdir'))
        outputDir = varargin{n + 1};
    end
end

if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('typeEntity', 'file', 'typeSelection', 'single', 'filter', filterP, 'title', 'Select a valid parameter file');
end

if (~exist('outputDir', 'var'))
    outputDir = icatb_selectEntry('typeEntity', 'directory', 'title', 'Select Analysis Output Directory');
end

if (isempty(outputDir))
    outputDir = pwd;   
end


if (~exist('num_dict', 'var'))
    
    % open input dialog box
    numParameters = 1;
    
    inputText(numParameters).promptString = 'Enter number of dictionary elements';
    inputText(numParameters).answerString = '3';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).dataType = 'numeric';
    inputText(numParameters).tag = 'num_dict';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString = 'Display statements?';
    inputText(numParameters).answerString = {'Yes', 'No'};
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).dataType = 'string';
    inputText(numParameters).tag = 'verbose';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString = 'Enter max no of iterations';
    inputText(numParameters).answerString = '75';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).dataType = 'numeric';
    inputText(numParameters).tag = 'max_iter';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
    
    answer = icatb_inputDialog('inputtext', inputText, 'Title', 'KSVD params', 'handle_visibility',  'on');
    
    num_dict = answer{1};
    verbose = strcmpi(answer{2}, 'yes');
    max_iter = answer{3};
    
end


drawnow;

disp('Computing KSVD on the subject timecourses ...');

%outDir = fileparts(param_file);

load(param_file, 'sesInfo');
sesInfo.outputDir = fileparts(param_file);
sesInfo.userInput.pwd = sesInfo.outputDir;

TC = icatb_loadComp(sesInfo, (1:sesInfo.numComp), 'vars_to_load', 'tc');

TC_ICN = zscore(TC, [], 2);

D = reshape(TC_ICN, [], size(TC_ICN,3))';

param.L = 1;   % number of elements in each linear combination.
param.K = num_dict; % number of dictionary elements
param.errorFlag = 0; % decompose signals until a certain error is reached. do not use fix number of coefficients.
%param.errorGoal = sigma;
param.preserveDCAtom = 0;

%%%%%%%% initial dictionary: Dictionary elements %%%%%%%%
param.InitializationMethod =  'DataElements';

param.displayProgress = verbose;
param.numIteration = max_iter;


max_iter_num = 32;

%dic_sz = 4;

%for dic_sz=4:4
param.K = num_dict; % number of dictionary elements
Dictionary_iter = cell(max_iter_num,1);
output_iter = cell(max_iter_num,1);

parfor iter_num = 1:max_iter_num
    [Dictionary_iter{iter_num},output_iter{iter_num}]  = icatb_ksvd(D, param);
end

best_err = Inf;

for alaki_iter=1:size(output_iter,1)
    
    if(sum(output_iter{alaki_iter}.totalerr) < best_err)
        
        best_err = sum(output_iter{alaki_iter}.totalerr);
        Dictionary = Dictionary_iter{alaki_iter};
        output = output_iter{alaki_iter};
    end
end

drawnow;

wfcInfo.ksvd.Dictionary = Dictionary;
wfcInfo.ksvd.totalerr = output.totalerr;
wfcInfo.ksvd.coeff = output.CoefMatrix;
wfcInfo.param_file = param_file;

fileName = fullfile(outputDir, [sesInfo.userInput.prefix, '_windowless_fc.mat']);
icatb_save(fileName, 'wfcInfo');

disp(['Information is saved in file ', fileName]);
