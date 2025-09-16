
% Enter the values for the variables required for the ICA analysis.
% Variables are on the left and the values are on the right.
% Characters must be enterd in single quotes
%
% After entering the parameters, use icatb_batch_file_run(inputFile); 

%% Modality. Options are fMRI and EEG
modalityType = 'CONN';

%% Type of stability analysis
% Options are 1 and 2.
% 1 - Regular Group ICA
% 2 - Group ICA using icasso
% 3 - Group ICA using Minimum spanning tree (MST)
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
% If you have selected option 3 (user specified settings) you need to manually set the PCA options. See manual or other
% templates (icatb/icatb_batch_files/Input_data_subjects_1.m) for more information to set PCA options 
%
% Options are:
% 1 - Maximize Performance
% 2 - Less Memory Usage
% 3 - User Specified Settings
perfType = 1;


%% Data selection
% Input data file pattern for data-sets must be in a cell array. The no. of rows of cell array correspond to no. of subjects
% and columns correspond to sessions. In the below example, there are 3
% subjects and 1 session.
input_data_file_patterns = {'D:\Srinivas\gsu\data\visuomotor\sub01_vis\singleSlice\ns*.img';
    'D:\Srinivas\gsu\data\visuomotor\sub02_vis\singleSlice\ns*.img';
    'D:\Srinivas\gsu\data\visuomotor\sub03_vis\singleSlice\ns*.img'};

% Enter no. of dummy scans to exclude from the group ICA analysis. If you have no dummy scans leave it as 0.
dummy_scans = 0;


%% Enter directory to put results of analysis
outputDir = 'D:\Srinivas\gsu\results\vis_conn_ica';

%% Enter Name (Prefix) Of Output Files
prefix = 'Visuomotor';

%% Enter location (full file path) of the image file to use as mask
% or use Default mask which is []
maskFile = [];

%% Enter subsampling depth for reducing the voxel space
subsampling_depth = 1;

%% Back reconstruction type. Options are 1 and 2
% 1 - Regular
% 2 - Spatial-temporal Regression 
% 3 - GICA3
% 4 - GICA
% 5 - GIG-ICA
backReconType = 4;

%% Number of pc to reduce each subject down to at each reduction step
% The number of independent components the will be extracted is the same as 
% the number of principal components after the final data reduction step.  
numOfPC1 = 45;
numOfPC2 = 30;

%% Scale the Results. Options are 0, 1, 2
% 0 - Don't scale
% 2 - Scale to Z scores
scaleType = 2;


%% 'Which ICA Algorithm Do You Want To Use';
% see icatb_icaAlgorithm for details or type icatb_icaAlgorithm at the
% command prompt.
% Note: Use only one subject and one session for Semi-blind ICA. Also specify atmost two reference function names

% 1 means infomax, 2 means fastICA, etc.
algoType = 1;

%% Report generator (fmri and smri only)
display_results.formatName = 'html'; 
display_results.slices_in_mm = (-40:4:72);
display_results.convert_to_zscores = 'yes';
display_results.threshold = 1.0;
display_results.image_values = 'positive and negative';
display_results.slice_plane = 'axial';
display_results.anatomical_file = which('ch2bet_3x3x3.nii');


%% ICA Options - Name by value pairs in a cell array. Options will vary depending on the algorithm. See icatb_icaOptions for more details. Some options are shown below.
% Infomax -  {'posact', 'off', 'sphering', 'on', 'bias', 'on', 'extended', 0}
% FastICA - {'approach', 'symm', 'g', 'tanh', 'stabilization', 'on'}

icaOptions = {'posact', 'off', 'sphering', 'on', 'bias', 'on', 'extended', 0};