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
outputDir = [s_dir_tmp filesep 'out_test_icatb_nmark_batch'];
if isfolder(outputDir)
    rmdir(outputDir, 's');
end

prefix    = 'nmarkrsn';

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
scaleType = 0;

% ICA algorithm (string name or numeric, depending on your GIFT version)
% Examples: 'infomax', 'fastica', 'moo-icar', ...
algoType = 'moo-icar';

% Spatial reference template for constrained / Neuromark-style ICA
% (only used by certain algorithms like 'moo-icar' / constrained spatial ICA)
% fMRI templates: Neuromark_fMRI_1.0.nii, Neuromark_fMRI_2.0_modelorder-175.nii,
%   Neuromark_fMRI_2.0_modelorder-25.nii, Neuromark_fMRI_2.1_modelorder-multi.nii,
%   Neuromark_fMRI_2.2_modelorder-multi.nii, Neuromark_fMRI_3.0_aging_modelorder-100.nii
%   Neuromark_fMRI_3.0_development_modelorder-100.nii, Neuromark_fMRI_3.0_infant_modelorder-100.nii
%   Neuromark_fMRI_WM_2.2_modelorder-multi.nii
% Templates for sMRI Neuromark: Neuromark_sMRI_1.0_modelorder-30_2x2x2.nii
%   Neuromark_sMRI_3.0_modelorder-100_3x3x3.nii, Neuromark_dMRI_3.0_modelorder-100_3x3x3.nii
%   Neuromark_PET-FBP_1.0_modelorder-40_2x2x2.nii
refFiles = [s_root 'tests/data/misc/RSN_28_lowres.nii'];

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

