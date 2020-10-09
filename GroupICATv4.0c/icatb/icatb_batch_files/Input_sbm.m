% Enter the values for the variables required for the ICA analysis.
% Variables are on the left and the values are on the right.
% Characters must be enterd in single quotes

%% Modality
modalityType = 'smri';

%% Type of analysis
% Options are 1 and 2.
% 1 - Regular Group ICA
% 2 - Group ICA using icasso
which_analysis = 1;

%% ICASSO options.
% This variable will be used only when which_analysis variable is set to 2.
icasso_opts.sel_mode = 'randinit';  % Options are 'randinit', 'bootstrap' and 'both'
icasso_opts.num_ica_runs = 5; % Number of times ICA will be run


%% Data file pattern. Include all subjects (3D images) in the file pattern or use char 
% array if you want to use a particular order like char('D:\myfile.img',
% D:\sub002\fil2.img');
input_data_file_patterns = 'D:\sbm_data\*.img';

%% Enter directory to put results of analysis
outputDir = 'D:\test_sbm\batch_results';

%% Enter Name (Prefix) Of Output Files
prefix = 'smri_test';

%% Enter location (full file path) of the image file to use as mask
% or use Default mask which is []
maskFile = [];

%% Batch Estimation. If 1 is specified then estimation of 
% the components takes place and the corresponding PC numbers are associated
% Options are 1 or 0
doEstimation = 0; 

%% Number of pc to reduce each subject down to at each reduction step
% The number of independent components the will be extracted is the same as 
% the number of principal components after the final data reduction step.  
numOfPC1 = 45;

%% Scale the components. Options are 0 and 2
% 0 - Don't scale
% 2 - Scale to Z scores
scaleType = 0;

%% 'Which ICA Algorithm Do You Want To Use';
% see icatb_icaAlgorithm for details or type icatb_icaAlgorithm at the
% command prompt.
% Note: Use only one subject and one session for Semi-blind ICA. Also specify atmost two reference function names

% 1 means infomax, 2 means fastICA, etc.
algoType = 1;