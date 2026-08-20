%% NeuroMark PET SBM Batch Input File Example
% TReNDS 8/19/26 Cyrus Eierud
%
% Before running GIFT/SBM using PET input, complete the following steps:
%
% 1) Make sure you are using GIFT v4.0.6.44 or later.
%
% 2) Copy this input file, input_neuromark_pet.m, to:
%
%       /Users/myname/output
%
% 3) Edit the following variables in input_neuromark_pet.m so that they
%    match your file system:
%
%       dir_subjects
%       input_data_file_patterns
%       outputDir
%
% 4) If you want to retain the mean in the resulting loading values,
%    set the following variable in icatb_defaults.m:
%
%       SBM_Z_SCORES_RMMEAN = 0;
%
% 5) In the MATLAB Command Window, change to the output directory:
%
%       cd /Users/myname/output
%
% 6) Run the NeuroMark PET SBM analysis:
%
%       icatb_batch_file_run('input_neuromark_pet.m');
%
% -------------------------------------------------------------------------
% This batch file specifies the parameters required for the NeuroMark PET
% spatially constrained SBM analysis.
%
% Variables are specified on the left and their values on the right.
% Character values must be enclosed in single quotes.
% -------------------------------------------------------------------------


%% Modality
% GIFT modality used for structural/SBM analysis.
modalityType = 'sMRI';


%% Group PCA performance
% The optimal internal settings are selected based on MAX_AVAILABLE_RAM
% in icatb_defaults.m.
%
% Options:
%   1 - Maximize performance
%   2 - Reduce memory usage
%   3 - User-specified settings
perfType = 2;


%% Data selection method
dataSelectionMethod = 4;


%% Input data
% Directory containing the PET images.
dir_subjects = '/Users/myname/mysubjects/';

% Hint: If you have a lot of subjects you may create a for and next statement to 
%          build up your input_data_file_patterns variable.

% Specify one 3D PET image per subject for 20 subjects.
input_data_file_patterns = {
      [dir_subjects 'subject01.nii'] ; [dir_subjects 'subject02.nii'] ;
      [dir_subjects 'subject03.nii'] ; [dir_subjects 'subject04.nii'] ;
      [dir_subjects 'subject05.nii'] ; [dir_subjects 'subject06.nii'] ;
      [dir_subjects 'subject07.nii'] ; [dir_subjects 'subject08.nii'] ;
      [dir_subjects 'subject09.nii'] ; [dir_subjects 'subject10.nii'] ;
      [dir_subjects 'subject11.nii'] ; [dir_subjects 'subject12.nii'] ;
      [dir_subjects 'subject13.nii'] ; [dir_subjects 'subject14.nii'] ;
      [dir_subjects 'subject15.nii'] ; [dir_subjects 'subject16.nii'] ;
      [dir_subjects 'subject17.nii'] ; [dir_subjects 'subject18.nii'] ;
      [dir_subjects 'subject19.nii'] ; [dir_subjects 'subject20.nii'] ;
};


%% Output
% Directory where analysis results will be written.
outputDir = '/Users/myname/output';

% Prefix used for output files.
prefix = 'neuromarkpet';


%% Mask
% Use the default intracranial-volume mask.
% Verify that the resulting mask is appropriate for the PET data.
maskFile = 'default&icv';


%% Data preprocessing
% For PET/SBM, no mean removal may be performed so that the original
% intensity information can be retained in resulting loading values.
preproc_type = 'none';


%% PCA settings
% PCA type:
%   1 - Standard PCA
%   2 - Expectation Maximization PCA
%   3 - SVD PCA
%   4 - MPOWIT
%   5 - STP
pcaType = 1;

% Disable automatic estimation of the number of principal components.
doEstimation = 0;

% Number of group data-reduction steps.
numReductionSteps = 1;

% Number of principal components retained at each reduction step.
%% Number of principal components retained at the first PCA step may be 1.5*number of components in template
% e.g., if template has 93 components, numOfPC1 = round(93*1.5) = 140
numOfPC1 = 140;


%% Component scaling
% Options:
%   0 - Do not scale
%   2 - Scale components to Z scores
scaleType = 2;


%% ICA algorithm
% Algorithm identifiers are defined by icatb_icaAlgorithm.
%
% For SBM, the indexes are e.g.:
%   1  - Infomax
%   2  - FastICA
%   15 - MOO-ICAR
%
% NeuroMark PET uses MOO-ICAR for spatially constrained SBM.
algoType = 15;


%% NeuroMark PET spatial reference
% Locate the NeuroMark PET FBP WM/GM reference distributed with GIFT.
refFiles = which('Neuromark_PET-FBP_WMGM_2.0_modelorder-100_3x3x3.nii');


%% Back reconstruction
backReconType = 'gig-ica';


%% Report generation
% Set to 1 to generate/display analysis results.
display_results = 0;


%% Dummy scans
% Not applicable to PET/SBM.
dummy_scans = 0;


%% Parallel processing
% Uncomment these lines to enable parallel processing.
%
% parallel_info.mode = 'parallel';
% parallel_info.num_workers = 40;


%% Analysis type
% Options:
%   1 - Regular group ICA
%   2 - Group ICA using ICASSO
which_analysis = 1;


%% ICASSO options
% These settings are used only when which_analysis = 2.
%
% sel_mode options:
%   'randinit'  - Random initialization
%   'bootstrap' - Bootstrap
%   'both'      - Random initialization and bootstrap
icasso_opts.sel_mode = 'randinit';
icasso_opts.num_ica_runs = 1;