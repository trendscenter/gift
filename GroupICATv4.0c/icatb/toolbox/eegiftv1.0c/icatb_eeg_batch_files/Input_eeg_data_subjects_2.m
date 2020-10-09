% Modality type:
modalityType = 'EEG';

% Type of analysis
% Options are 1 and 2.
% 1 - Regular Group ICA
% 2 - Group ICA using icasso
which_analysis = 1;

% ICASSO options.
% This variable will be used only when which_analysis variable is set to 2.
icasso_opts.sel_mode = 'randinit';  % Options are 'randinit', 'bootstrap' and 'both'
icasso_opts.num_ica_runs = 5; % Number of times ICA will be run
% Most stable run estimate is based on these settings. 
icasso_opts.min_cluster_size = 2; % Minimum cluster size
icasso_opts.max_cluster_size = 15; % Max cluster size. Max is the no. of components


%% Group PCA performance settings. Best setting for each option will be selected based on variable MAX_AVAILABLE_RAM in icatb_defaults.m. 
% If you have selected option 3 (user specified settings) you need to manually set the PCA options. 
%
% Options are:
% 1 - Maximize Performance
% 2 - Less Memory Usage
% 3 - User Specified Settings
perfType = 1;

% There are two ways to enter the subject data
% Options are 1, 2, 3 and 4
dataSelectionMethod = 1;

%%%%%%%%%%%%% Method 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If you have all subjects in one directory and their sessions in a separate folder or in the subject folder then specify 
% root directory, filePattern and flag.
% Options for flag are: data_in_subject_folder, data_in_subject_subfolder
%
% 1. data_in_subject_subfolder - Data is selected from the subject sub
% folders. Number of sessions is equal to the number of sub-folders
% containing the specified file pattern.
%
% 2. data_in_subject_folder - Data is selected from the subject
% folders. Number of sessions is 1 and number of subjects is equal to the number of subject folders
% containing the specified file pattern.

% Note: Make sure the sessions are the same over subjects.

sourceDir_filePattern_flagLocation = {'D:\test_GIFT\new_version\GIFT_EEG_Results\oddball\Data', 'Oddball.mat', 'data_in_subject_folder'};

%%%%%%%%%%%%%%%%%%%%%%%%% end for Method 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%% Method 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If you have different filePatterns and location for subjects not in one
% root folder then enter the data here.
% Number of subjects is determined getting the length of the selected subjects. Specify the data set or data sets needed for 
% the analysis here.

selectedSubjects = {'s1', 's2', 's3'};  % naming for subjects s1 refers to subject 1, s2 means subject 2. Use cell array convention even in case of one subject one session

% Number of Sessions
numOfSess = 1;

% Data folder and file pattern 
s1_s1 = {'D:\test_GIFT\new_version\GIFT_EEG_Results\oddball\Data\oddball_0001', 'Oddball.mat'}; % subject 1 session 1

s2_s1 = {'D:\test_GIFT\new_version\GIFT_EEG_Results\oddball\Data\oddball_0002', 'Oddball.mat'}; % subject 2 session 1

s3_s1 = {'D:\test_GIFT\new_version\GIFT_EEG_Results\oddball\Data\oddball_0003', 'Oddball.mat'}; % subject 3 session 1

%%%%%%%%%%%%%%%%%%%%%%% end for Method 2 %%%%%%%%%%%%%%


%%%%%%%%%%%%%% Method 3 (Uses Regular expressions) %%%%%%%%%%%%%%%%%%%%%

% Input data directory name
input_directory_name = 'D:\test_GIFT\new_version\GIFT_EEG_Results\oddball\Data\';

% Subject directory regular expression. This variable can have nested paths
% like Sub01_vis\Study1. To match this Sub\w+; Study\w+ regular expression can be used where semi-colon
% is used as a path separator. If there are no subject directories inside the input directory, leave it as empty like ''
subject_dir_regexp = 'odd\w+';

% Session directory regular expression. This variable cannot have nested
% paths. If there are no session directories inside subject directories, leave it as empty.
session_dir_regexp = '';

% Data file pattern. Use wild card for this and not regular expression.
data_file_pattern = 'Oddball.mat';

%%%%%%%% End for Method 3 %%%%%%%%%%%%

% Method 4
% Input data file pattern for data-sets must be in a cell array. The no. of rows of cell array correspond to no. of subjects
% and columns correspond to sessions. In the below example, there are 3
% subjects and 1 session. If you have multiple sessions, please see
% Input_data_subjects_2.m file.
input_data_file_patterns = {'D:\test_GIFT\new_version\GIFT_EEG_Results\oddball\Data\oddball_0001\Oddball.mat';
    'D:\test_GIFT\new_version\GIFT_EEG_Results\oddball\Data\oddball_0002\Oddball.mat';
    'D:\test_GIFT\new_version\GIFT_EEG_Results\oddball\Data\oddball_0003\Oddball.mat'};

%%%%%%%% End for Method 4 %%%%%%%%%%%%

% Enter directory to put results of analysis
outputDir = 'D:\test_GIFT\new_version\GIFT_EEG_Results\oddball\Results\batch_results';

% Enter Name (Prefix) Of Output Files
prefix = 'Oddball';


%% Group PCA Type. Used for analysis on multiple subjects and sessions.
% Options are 'subject specific' and 'grand mean'. 
%   a. Subject specific - Individual PCA is done on each data-set before group
%   PCA is done.
%   b. Grand Mean - PCA is done on the mean over all data-sets. Each data-set is
%   projected on to the eigen space of the mean before doing group PCA.
%
% NOTE: Grand mean implemented is from FSL Melodic. Make sure that there are
% equal no. of electrodes between data-sets.
%
group_pca_type = 'subject specific';

%% Back reconstruction type. Options are str and gica
backReconType = 'gica';

%% Data Pre-processing options
% 1 - Remove mean per time point
% 2 - Remove mean per voxel
% 3 - Intensity normalization
% 4 - Variance normalization
preproc_type = 3;


%% PCA Type. Also see options associated with the selected pca option.
%% Standard PCA options are commented.
% Options are 1 and 2
% 1 - Standard 
% 2 - Expectation Maximization
pcaType = 2;

%% PCA options (Standard)

% % a. Options are yes or no
% % 1a. yes - Datasets are stacked. This option uses lot of memory depending
% % on datasets, voxels and components.
% % 2a. no - A pair of datasets are loaded at a time. This option uses least
% % amount of memory and can run very slower if you have very large datasets.
% pca_opts.stack_data = 'yes';
% 
% % b. Options are full or packed.
% % 1b. full - Full storage of covariance matrix is stored in memory.
% % 2b. packed - Lower triangular portion of covariance matrix is only stored in memory.
% pca_opts.storage = 'full';
% 
% % c. Options are double or single.
% % 1c. double - Double precision is used
% % 2c. single - Floating point precision is used.
% pca_opts.precision = 'double';
% 
% % d. Type of eigen solver. Options are selective or all
% % 1d. selective - Selective eigen solver is used. If there are convergence
% % issues, use option all.
% % 2d. all - All eigen values are computed. This might run very slow if you
% % are using packed storage. Use this only when selective option doesn't
% % converge.
% 
% pca_opts.eig_solver = 'selective';


%% PCA Options (Expectation Maximization)
% a. Options are yes or no
% 1a. yes - Datasets are stacked. This option uses lot of memory depending
% on datasets, voxels and components.
% 2a. no - A pair of datasets are loaded at a time. This option uses least
% amount of memory and can run very slower if you have very large datasets.
pca_opts.stack_data = 'yes';

% b. Options are double or single.
% 1b. double - Double precision is used
% 2b. single - Floating point precision is used.
pca_opts.precision = 'single';

% c. Stopping tolerance 
pca_opts.tolerance = 1e-4;

% d. Maximum no. of iterations
pca_opts.max_iter = 1000;

% Maximum reduction steps you can select is 3
% Note: This number will be changed depending upon the number of subjects
% and sessions
numReductionSteps = 2;

% number of pc to reduce each subject down to at each reduction step
% The number of independent components the will be extracted is the same as 
% the number of principal components after the final data reduction step.  
numOfPC1 = 20;
numOfPC2 = 16;
numOfPC3 = 16;


%% Scale the Results. Options are 0, 1, 2, 3 and 4
% 0 - Don't scale
% 1 - Scale to data units
% 2 - Scale to Z scores
% 3 - Normalize timecourses using the maximum value and multiply topographies using the maximum value
% 4 - Scale topographies using the maximum value and timecourses using the standard deviation of topographies
scaleType = 0;

% 'Which ICA Algorithm Do You Want To Use';
% see icatb_icaAlgorithm for details or type icatb_icaAlgorithm at the
% command prompt.
% Note: You cannot use constrained ICA algorithms for EEGIFT

% 1 means infomax, 2 means fastICA, etc.
algoType = 2;

%% ICA Options - Name by value pairs in a cell array. Options will vary depending on the algorithm. See icatb_icaOptions for more details. Some options are shown below.
% Infomax -  {'posact', 'off', 'sphering', 'on', 'bias', 'on', 'extended', 0}
% FastICA - {'approach', 'symm', 'g', 'tanh', 'stabilization', 'on'}

icaOptions =  {'approach', 'symm', 'g', 'tanh', 'stabilization', 'on'};