% Enter the values for the variables required for the ICA analysis.
% Variables are on the left and the values are on the right.
% Characters must be enterd in single quotes


%% Modality
modalityType = 'fnc';

%% Type of analysis
% Options are 1 and 2.
% 1 - Regular Group ICA
% 2 - Group ICA using icasso
which_analysis = 2;


%% ICASSO options.
% This variable will be used only when which_analysis variable is set to 2.
icasso_opts.sel_mode = 'randinit';  % Options are 'randinit', 'bootstrap' and 'both'
icasso_opts.num_ica_runs = 5; % Number of times ICA will be run


%% Data file pattern. Enter file patterns for each cell session in %a cell array of size number of subjects by sessions or pass the 
%gift parameter file name 
input_data_file_patterns = { 'C:\Users\srrac\gsu\srinivas\data\olin_fnc_data\Session1\fnc_Sub_001.txt', 'C:\Users\srrac\gsu\srinivas\data\olin_fnc_data\Session2\fnc_Sub_001.txt';
    'C:\Users\srrac\gsu\srinivas\data\olin_fnc_data\Session1\fnc_Sub_002.txt', 'C:\Users\srrac\gsu\srinivas\data\olin_fnc_data\Session2\fnc_Sub_002.txt';
    'C:\Users\srrac\gsu\srinivas\data\olin_fnc_data\Session1\fnc_Sub_003.txt', 'C:\Users\srrac\gsu\srinivas\data\olin_fnc_data\Session2\fnc_Sub_003.txt';
    'C:\Users\srrac\gsu\srinivas\data\olin_fnc_data\Session1\fnc_Sub_004.txt', 'C:\Users\srrac\gsu\srinivas\data\olin_fnc_data\Session2\fnc_Sub_004.txt';
    'C:\Users\srrac\gsu\srinivas\data\olin_fnc_data\Session1\fnc_Sub_005.txt', 'C:\Users\srrac\gsu\srinivas\data\olin_fnc_data\Session2\fnc_Sub_005.txt';
    'C:\Users\srrac\gsu\srinivas\data\olin_fnc_data\Session1\fnc_Sub_006.txt', 'C:\Users\srrac\gsu\srinivas\data\olin_fnc_data\Session2\fnc_Sub_006.txt';
    'C:\Users\srrac\gsu\srinivas\data\olin_fnc_data\Session1\fnc_Sub_007.txt', 'C:\Users\srrac\gsu\srinivas\data\olin_fnc_data\Session2\fnc_Sub_007.txt';
    'C:\Users\srrac\gsu\srinivas\data\olin_fnc_data\Session1\fnc_Sub_008.txt', 'C:\Users\srrac\gsu\srinivas\data\olin_fnc_data\Session2\fnc_Sub_008.txt';
    'C:\Users\srrac\gsu\srinivas\data\olin_fnc_data\Session1\fnc_Sub_009.txt', 'C:\Users\srrac\gsu\srinivas\data\olin_fnc_data\Session2\fnc_Sub_009.txt';
    'C:\Users\srrac\gsu\srinivas\data\olin_fnc_data\Session1\fnc_Sub_010.txt', 'C:\Users\srrac\gsu\srinivas\data\olin_fnc_data\Session2\fnc_Sub_010.txt';};

%input_data_file_patterns = 'C:\Users\srrac\gsu\srinivas\results\olin_rest\new\olin_rest_ica_parameter_info.mat'; 


%% If Mat file, provide variable of mat file
fnc_variable_mat_file = '';

%% Enter directory to put results of analysis
outputDir = 'C:\Users\srrac\gsu\srinivas\results\test_fnc\test_fnc2';

%% Enter Name (Prefix) Of Output Files
prefix = 'fnc_analysis';

%% Contrast vector (if multiple sessions enter a contrast vector like 1 -1 to subtract timepoints
contrast_vector = [1, -1];

%% Batch Estimation. If 1 is specified then estimation of 
% the components takes place and the corresponding PC numbers are associated
% Options are 1 or 0
doEstimation = 1; 

%% Number of pc to reduce each subject down to at each reduction step
% The number of independent components the will be extracted is the same as 
% the number of principal components after the final data reduction step.  
numOfPC1 = 16;

%% Scale Type (options are no and 'z-scores')
scaleType = 'z-scores';

%% 'Which ICA Algorithm Do You Want To Use';
% see icatb_icaAlgorithm for details or type icatb_icaAlgorithm at the
% command prompt.
% Note: Use only one subject and one session for Semi-blind ICA. Also specify atmost two reference function names

% 1 means infomax, 2 means fastICA, etc.
algoType = 1;

%% Network nifti file (Labels text file.txt extension with the same name will be automatically searched based on the nifti file name). Otherwise provide
% network file namd and labels as cell array like {file_name, labels}
reference_file = 'C:\Users\srrac\Desktop\gift_latest\gift-master\GroupICAT\icatb\icatb_templates\RSN_28.nii';

%% Display results
display_results = 1;
