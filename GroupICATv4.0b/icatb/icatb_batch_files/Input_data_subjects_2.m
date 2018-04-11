% Enter the values for the variables required for the ICA analysis.
% Variables are on the left and the values are on the right.
% Characters must be enterd in single quotes
%
% After entering the parameters, use icatb_batch_file_run(inputFile); 


%% Modality. Options are fMRI and EEG
modalityType = 'fMRI';

%% Type of analysis
% Options are 1, 2 and 3.
% 1 - Regular Group ICA
% 2 - Group ICA using icasso
% 3 - Group ICA using MST
which_analysis = 1;

%% ICASSO options.
% This variable will be used only when which_analysis variable is set to 2.
icasso_opts.sel_mode = 'randinit';  % Options are 'randinit', 'bootstrap' and 'both'
icasso_opts.num_ica_runs = 5; % Number of times ICA will be run
% Most stable run estimate is based on these settings. 
icasso_opts.min_cluster_size = 2; % Minimum cluster size
icasso_opts.max_cluster_size = 15; % Max cluster size. Max is the no. of components


%% Enter TR in seconds. If TRs vary across subjects, TR must be a row vector of length equal to the number of subjects.
TR = 1;


%% Group ica type
% Options are spatial or temporal for fMRI modality. By default, spatial
% ica is run if not specified.
group_ica_type = 'spatial';


%% Parallel info
% enter mode serial or parallel. If parallel, enter number of
% sessions/workers to do job in parallel
parallel_info.mode = 'parallel';
parallel_info.num_workers = 4;


%% Group PCA performance settings. Best setting for each option will be selected based on variable MAX_AVAILABLE_RAM in icatb_defaults.m. 
% If you have selected option 3 (user specified settings) you need to manually set the PCA options. 
%
% Options are:
% 1 - Maximize Performance
% 2 - Less Memory Usage
% 3 - User Specified Settings
perfType = 1;

%% Design matrix selection 
% Design matrix (SPM.mat) is used for sorting the components
% temporally (time courses) during display. Design matrix will not be used during the
% analysis stage except for SEMI-BLIND ICA.
% options are ('no', 'same_sub_same_sess', 'same_sub_diff_sess', 'diff_sub_diff_sess')
% 1. 'no' - means no design matrix.
% 2. 'same_sub_same_sess' - same design over subjects and sessions
% 3. 'same_sub_diff_sess' - same design matrix for subjects but different
% over sessions
% 4. 'diff_sub_diff_sess' - means one design matrix per subject.

keyword_designMatrix = 'same_sub_same_sess';

%% Specify location of design matrix here if you have selected 'same_sub_same_sess' or
% 'same_sub_diff_sess' option for keyword_designMatrix variable
OnedesignMat = 'E:\drive pilot data (last three)\stats\SPM.mat';


%% There are three ways to enter the subject data
% options are 1, 2, 3 or 4
dataSelectionMethod = 3;

%% Method 1 

% If you have all subjects in one directory and their sessions in a separate folder or in the subject folder then specify 
% root directory, filePattern, flag and file numbers to include.
% Options for flag are: data_in_subject_folder, data_in_subject_subfolder
%
% 1. data_in_subject_subfolder - Data is selected from the subject sub
% folders. Number of sessions is equal to the number of sub-folders
% containing the specified file pattern.
%
% 2. data_in_subject_folder - Data is selected from the subject
% folders. Number of sessions is 1 and number of subjects is equal to the number of subject folders
% containing the specified file pattern.

% You can provide the file numbers ([1:220]) to include as a vector. If you want to
% select all the files then leave empty.

% Note: Make sure the sessions are the same over subjects.

sourceDir_filePattern_flagLocation = {'E:\drive pilot data (last three)\', 'sw*.img', ...
        'data_in_subject_subfolder'};

%% Specify design matrix filter pattern here if you have selected
% 'diff_sub_diff_sess' option for keyword_designMatrix variable for Method 1. It looks
% for the design matrix in the respective subject folder or session folders.
spmDesignFilter = 'SPM*.mat';    

%%%%%%%%%%%%%%%%%%%%%%%%% end for Method 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Method 2 

% If you have different filePatterns and location for subjects not in one
% root folder then enter the data here.
% Number of subjects is determined getting the length of the selected subjects. specify the data set or data sets needed for 
% the analysis here.

selectedSubjects = {'s1', 's2', 's3'};  % naming for subjects s1 refers to subject 1, s2 means subject 2. Use cell array convention even in case of one subject one session

% Number of Sessions
numOfSess = 2;

% functional data folder, file pattern and file numbers to include
% You can provide the file numbers ([1:220]) to include as a vector. If you want to
% select all the files then leave empty.

s1_s1 = {'E:\drive pilot data (last three)\drive_4014069\3', 'sw*.img'}; % subject 1 session 1
s1_s2 = {'E:\drive pilot data (last three)\drive_4014069\4', 'sw*.img'}; % subject 1 session 2

s2_s1 = {'E:\drive pilot data (last three)\drive_6367240\2', 'sw*.img'}; % subject 2 session 1
s2_s2 = {'E:\drive pilot data (last three)\drive_6367240\3', 'sw*.img'}; % subject 2 session 2

s3_s1 = {'E:\drive pilot data (last three)\drive_7364400\2', 'sw*.img'}; % subject 3 session 1
s3_s2 = {'E:\drive pilot data (last three)\drive_7364400\3', 'sw*.img'}; % subject 3 session 2

% specify design matrix for each subject here if you selected 'diff_sub_diff_sess' under 
% keyword_designMatrix variable for Method 2
s1_designMat = 'E:\drive pilot data (last three)\stats\SPM.mat'; % subject 1 design matrix
s2_designMat = 'E:\drive pilot data (last three)\stats\SPM.mat'; % subject 2 design matrix
s3_designMat = 'E:\drive pilot data (last three)\stats\SPM.mat'; % subject 3 design matrix

%%%%%%%%%%%%%%%%%%%%%%% end for Method 2 %%%%%%%%%%%%%%


%% Method 3 (Uses Regular expressions) 

% Input data directory name
input_directory_name = 'E:\drive pilot data (last three)\';

% Subject directory regular expression. This variable can have nested paths
% like Sub01_vis\Study1. To match this Sub\w+; Study\w+ regular expression can be used where semi-colon
% is used as a path separator. If there are no subject directories inside the input directory, leave it as empty like ''
subject_dir_regexp = 'drive\w+';

% Session directory regular expression. This variable cannot have nested
% paths. If there are no session directories inside subject directories, leave it as empty.
session_dir_regexp = '\<\d\>';

% Data file pattern. Use wild card for this and not regular expression.
data_file_pattern = 'sw*.img';

% File numbers to include. Leave it as empty if you want to include all of
% them.
file_numbers_to_include = [];

% SPM stats directory name relative to subject or session directories. Use this only when you specify
% 'diff_sub_diff_sess' as the value for keyword_designMatrix variable. GIFT
% will first search in the subject directories and later session
% directories to find SPM.mat files
spm_stats_dir = '';

%%%%%%%% End for Method 3 %%%%%%%%%%%%

%% Method 4
% Input data file pattern for data-sets must be in a cell array. The no. of rows of cell array correspond to no. of subjects
% and columns correspond to sessions. In the below example, there are 3
% subjects and 2 sessions.
input_data_file_patterns = {'E:\drive pilot data (last three)\drive_4014069\3\sw*.img', 'E:\drive pilot data (last three)\drive_4014069\4\sw*.img';
    'E:\drive pilot data (last three)\drive_6367240\2\sw*.img', 'E:\drive pilot data (last three)\drive_6367240\3\sw*.img';
    'E:\drive pilot data (last three)\drive_7364400\2\sw*.img', 'E:\drive pilot data (last three)\drive_7364400\3\sw*.img'};

% Input for design matrices will be used only if you have a design matrix
% for each subject i.e., if you have selected 'diff_sub_diff_sess' for
% variable keyword_designMatrix.
input_design_matrices = {'E:\drive pilot data (last three)\stats\SPM.mat';
    'E:\drive pilot data (last three)\stats\SPM.mat';
    'E:\drive pilot data (last three)\stats\SPM.mat'};

% Enter no. of dummy scans to exclude from the group ICA analysis. If you have no dummy scans leave it as 0.
dummy_scans = 0;

%%%%%%%% End for Method 4 %%%%%%%%%%%%


%% Enter directory to put results of analysis
outputDir = 'D:\test_GIFT\batch_results';

%% Enter Name (Prefix) Of Output Files
prefix = 'Driving';

%% Enter location (full file path) of the image file to use as mask
% or use Default mask which is []
maskFile = [];

%% Group PCA Type. Used for analysis on multiple subjects and sessions.
% Options are 'subject specific' and 'grand mean'. 
%   a. Subject specific - Individual PCA is done on each data-set before group
%   PCA is done.
%   b. Grand Mean - PCA is done on the mean over all data-sets. Each data-set is
%   projected on to the eigen space of the mean before doing group PCA.
%
% NOTE: Grand mean implemented is from FSL Melodic. Make sure that there are
% equal no. of timepoints between data-sets.
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
% Standard PCA and SVD PCA options are commented.
% PCA options are commented.
% Options are 1, 2, 3, 4 and 5.
% 1 - Standard 
% 2 - Expectation Maximization
% 3 - SVD
% 4 - MPOWIT
% 5 - STP
pcaType = 1;

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

% %% PCA Options (SVD)
% % a. Options are double or single.
% % 1a. double - Double precision is used
% % 2a. single - Floating point precision is used.
% pca_opts.precision = 'single';

% % b. Type of eigen solver. Options are selective or all
% % 1b. selective - svds function is used.
% % 2b. all - Economy size decomposition is used.
% pca_opts.solver = 'selective';


%% Maximum reduction steps you can select is 2. Options are 1 and 2. For temporal ica, only one data reduction step is
% used.
numReductionSteps = 2;

%% Batch Estimation. If 1 is specified then estimation of 
% the components takes place and the corresponding PC numbers are associated
% Options are 1 or 0
doEstimation = 0; 

%% MDL Estimation options. This variable will be used only if doEstimation is set to 1.
% Options are 'mean', 'median' and 'max' for each reduction step. The length of cell is equal to
% the no. of data reductions used.
estimation_opts.PC1 = 'mean';
estimation_opts.PC2 = 'mean';

%% Number of pc to reduce each subject down to at each reduction step
% The number of independent components the will be extracted is the same as 
% the number of principal components after the final data reduction step.  
numOfPC1 = 30;
numOfPC2 = 25;

%% Scale the Results. Options are 0, 1, 2, 3 and 4
% 0 - Don't scale
% 1 - Scale to Percent signal change
% 2 - Scale to Z scores
% 3 - Normalize spatial maps using the maximum intensity value and multiply timecourses using the maximum intensity value
% 4 - Scale timecourses using the maximum intensity value and spatial maps using the standard deviation of timecourses
scaleType = 0;


%% 'Which ICA Algorithm Do You Want To Use';
% see icatb_icaAlgorithm for details or type icatb_icaAlgorithm at the
% command prompt.
% Note: Use only one subject and one session for Semi-blind ICA. Also specify atmost two reference function names

% 1 means infomax, 2 means fastICA, etc.
algoType = 2;

%% Specify atmost two reference function names if you select Semi-blind ICA algorithm.
% Reference function names can be acessed by loading SPM.mat in MATLAB and accessing 
% structure SPM.xX.name.
refFunNames = {'Sn(1) right*bf(1)', 'Sn(1) left*bf(1)'};


%% Specify spatial reference files for constrained ICA (spatial) or gig-ica
refFiles = {which('ref_default_mode.nii'), which('ref_left_visuomotor.nii'), which('ref_right_visuomotor.nii')};

%% ICA Options - Name by value pairs in a cell array. Options will vary depending on the algorithm. See icatb_icaOptions for more details. Some options are shown below.
% Infomax -  {'posact', 'off', 'sphering', 'on', 'bias', 'on', 'extended', 0}
% FastICA - {'approach', 'symm', 'g', 'tanh', 'stabilization', 'on'}

icaOptions =  {'approach', 'symm', 'g', 'tanh', 'stabilization', 'on'};