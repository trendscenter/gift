%% GIFT Batch Template (generic)
% Fill in the USER SETTINGS section, then run:
%   icatb_batch_file_run('input_this_file.m');
% Date 2/10/2026

%% -----------------------------
% USER SETTINGS
% -----------------------------

[s_root, ~]=fileparts(which('gift'));
[s_root, ~]=fileparts(s_root);
s_root = [s_root filesep];
s_dir_tmp = tempdir;
disp(s_dir_tmp)

% Modality: 'fMRI', 'sMRI', or 'EEG'
modalityType = 'fMRI';

% TR in seconds (scalar or 1 x nSubjects vector)
TR = 2;

% Output
outputDir = [s_dir_tmp filesep 'out_test_icatb_batch_infomax'];
if isfolder(outputDir)
    rmdir(outputDir, 's');
end

prefix    = 'infomaxrsn';

%% PCA Type. Also see options associated with the selected pca option. EM
% PCA options and SVD PCA are commented.
% Options are 1, 2, 3, 4 and 5.
% 1 - Standard 
% 2 - Expectation Maximization
% 3 - SVD
% 4 - MPOWIT
% 5 - STP
pcaType = 1;

%% PCA options (Standard)

% a. Options are yes or no
% 1a. yes - Datasets are stacked. This option uses lot of memory depending
% on datasets, voxels and components.
% 2a. no - A pair of datasets are loaded at a time. This option uses least
% amount of memory and can run very slower if you have very large datasets.
pca_opts.stack_data = 'yes';

% b. Options are full or packed.
% 1b. full - Full storage of covariance matrix is stored in memory.
% 2b. packed - Lower triangular portion of covariance matrix is only stored in memory.
pca_opts.storage = 'full';

% c. Options are double or single.
% 1c. double - Double precision is used
% 2c. single - Floating point precision is used.
pca_opts.precision = 'double';

% d. Type of eigen solver. Options are selective or all
% 1d. selective - Selective eigen solver is used. If there are convergence
% issues, use option all.
% 2d. all - All eigen values are computed. This might run very slow if you
% are using packed storage. Use this only when selective option doesn't
% converge.
pca_opts.eig_solver = 'selective';

%% Maximum reduction steps you can select is 2. Options are 1 and 2. For temporal ica, only one data reduction step is
% used.
numReductionSteps = 2;

%% Batch Estimation. If 1 is specified then estimation of 
% the components takes place and the corresponding PC numbers are associated
% Options are 1 or 0
doEstimation = 0; 

%% Number of pc to reduce each subject down to at each reduction step
% The number of independent components the will be extracted is the same as 
% the number of principal components after the final data reduction step. 
numOfPC1 = 30;
numOfPC2 = 20;

% Data selection method (1/2/3/4). This template uses Method 4.
dataSelectionMethod = 4;

% Method 4: list subject files (rows = subjects, cols = sessions)
% Example: 2 subjects, 1 session
input_data_file_patterns = {[s_root 'tests/data/subjects/sub-007_lowres.nii']};

% Optional: per-subject design matrices (only used for certain keyword_designMatrix settings)
% for each subject i.e., if you have selected 'diff_sub_diff_sess' for variable keyword_designMatrix.
input_design_matrices = {};

% Dummy scans to drop
dummy_scans = 0;

% Mask: [] for default, or full path, or special strings (if your lab uses them)
maskFile = [s_root 'tests/data/test_icatb_nmark_batch/sub-007mask.nii'];  % or [] / 'C:\path\mask.nii'

% Preprocessing:
% 1 Remove mean per time point
% 2 Remove mean per voxel
% 3 Intensity normalization
% 4 Variance normalization
preproc_type = 1;

% Scaling:
% 0 none, 1 percent signal change, 2 Z-scores
scaleType = 2;

% ICA algorithm (string name or numeric, depending on your GIFT version)
% Examples: 'infomax', 'fastica', 'moo-icar', ...
algoType = 'infomax';

%% -----------------------------
% PERFORMANCE / PARALLEL SETTINGS
% -----------------------------

% Performance type:
% 1 Maximize performance
% 2 Less memory usage
% 3 User specified settings
perfType = 1;

% Parallel execution
% mode: 'serial' or 'parallel'
parallel_info.mode        = 'serial';
parallel_info.num_workers = 4;

%% -----------------------------
% REPORT / DISPLAY SETTINGS  (fmri and smri only)
% -----------------------------
display_results = 0;

