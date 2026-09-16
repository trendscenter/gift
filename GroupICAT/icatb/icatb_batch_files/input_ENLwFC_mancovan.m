%% ENL-wFC / Linear wFC MANCOVAN batch example
%
% Use full paths for directories and files wherever possible.
%
% After editing the parameters below, run:
%
%   icatb_mancovan_batch('input_file');
%
% where input_file is the name of this batch script without the ".m"
% extension.
%
% TReNDS
% 09/08/2026
% Cyrus Eierud
%
% This script can be used after Explicitly Nonlinear Windowed Functional
% Connectivity (ENL-wFC) has been calculated in GIFT.
%
% It can be used to perform statistical testing with MANCOVAN after the
% primary ENL-wFC analysis, following a workflow similar to that used in
% Spencer et al. (2026).
%
%
%% Suggested ENL-wFC workflow
%
% Before running this batch script, process the data as follows:
%
% 1) In MATLAB, enter:
%
%       groupica
%
% 2) Click the [Connectivity Domain] button. The Connectivity window will
%    appear.
%
% 3) Click [Import Data]. The data-selection window will appear.
%
% 4) Select an output directory and click [OK]. The Import Data window
%    will appear.
%
% 5) Enter a prefix (for example, "demo"), select the fMRI NIfTI files,
%    and set "Add Linear Comparison" to [Yes].
%
%    Select the remaining analysis settings as appropriate.
%
%    NOTE: The number of principal components should be approximately
%    1.5 times the number of independent components (ICs) that you select
%    in the next step.
%
%    This procedure creates two GIFT projects:
%
%       prefix-enl   Explicitly nonlinear connectivity analysis
%       prefix-lin   Linear connectivity analysis
%
% 6) Click [Done]. Wait while the subject data are converted to
%    PCA-compressed connectivity matrices. This may take several minutes.
%
% 7) Click [Setup Analysis], select the parameter file, and click [OK].
%    The Setup Analysis window will appear.
%
% 8) Select the desired analysis settings and click [Done].
%
% 9) Click [Run Analysis].
%
%    After the analysis finishes, a message box reports the number of
%    components matched between the linear and ENL analyses and indicates
%    the file in which the component-matching information was saved.
%
% 10) Click [Display], select the:
%
%       demo-enl_ica_parameter_info.mat
%
%     file, choose the desired report settings, and click [OK].
%     The ENL report will then be generated.
%
% 11) Click [Display] again, this time selecting:
%
%       demo-lin_ica_parameter_info.mat
%
%     Select analogous report settings and click [OK].
%     The linear-connectivity report will then be generated.
%
%     Using the component-matching file produced in step 9, you can create
%     a list of matched linear and ENL components similar to Fig. 3 in
%     Spencer et al. (2026).
%
% 12) Make two copies of this script (input_ENLwFC_mancovan.m) to:
%
%       input_my_ENLwFC_mancovan_enl.m
%       input_my_ENLwFC_mancovan_lin.m
%
% 13) In input_my_ENLwFC_mancovan_enl.m:
%
%       - Set outputDir to the desired MANCOVAN output directory.
%
%       - Set:
%
%           ica_param_file = '/your/workdir/demo-enl_ica_parameter_info.mat';
%
%       - Set comp_network_names to the components/networks that you want
%         to test.
%
%       - To perform a two-sample t-test between two groups, similar to
%         Spencer et al. (2026), set:
%
%           univariate_tests = ...
%               {'Ttest2', {(1:XX), (XX+1:end)}, {'Group 1', 'Group 2'}};
%
%         Here:
%
%           1:XX       = subjects belonging to Group 1
%           XX+1:end   = subjects belonging to Group 2
%
%       - Set the ENL spatial mask:
%
%           feature_params.sm_mask = '/your/workdir/demo-enlMask.nii';
%
% 14) Run:
%
%       icatb_mancovan_batch('input_my_ENLwFC_mancovan_enl');
%
% 15) Configure input_my_ENLwFC_mancovan_lin.m analogously:
%
%       - Set outputDir to the desired MANCOVAN output directory.
%
%       - Set:
%
%           ica_param_file = '/your/workdir/demo-lin_ica_parameter_info.mat';
%
%       - Set comp_network_names to the components/networks that you want
%         to test.
%
%       - For a two-sample t-test:
%
%           univariate_tests = ...
%               {'Ttest2', {(1:XX), (XX+1:end)}, {'Group 1', 'Group 2'}};
%
%         Here:
%
%           1:XX       = subjects belonging to Group 1
%           XX+1:end   = subjects belonging to Group 2
%
%       - Set the linear spatial mask:
%
%           feature_params.sm_mask = '/your/workdir/demo-linMask.nii';
%
% 16) Run:
%
%       icatb_mancovan_batch('input_my_ENLwFC_mancovan_lin');
%
%
%% Comparing linear and ENL results
%
% After completing these steps, use the component-matching information
% from step 9 to compare the linear and ENL statistical results, similar
% to Fig. 5 in Spencer et al. (2026).
%
% The MANCOVAN p-value maps may look different from the t-score maps shown
% in Spencer et al. (2026). The corresponding t-statistics are also saved
% within the MANCOVAN output directories and can be used to generate
% t-score maps if desired.
%
%
%% Reference
%
% Title:
%   Networks extracted from nonlinear fMRI connectivity exhibit unique
%   spatial variation and enhanced sensitivity to differences between
%   individuals with schizophrenia and controls
%
% Year:
%   2026
%
% Authors:
%   Spencer Kinsey, Katarzyna Kazimierczak, Pablo Andres Camazon,
%   Jiayu Chen, Tulay Adali, Peter Kochunov, Bhim M. Adhikari,
%   Judith Ford, Theo G. M. van Erp, Mukesh Dhamala,
%   Vince D. Calhoun, and Armin Iraji
%
% The instructions above are intended to reproduce a statistical workflow
% similar to that used in the paper.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Output directory
% Directory in which the MANCOVAN results will be saved.

outputDir = '/home/cyrus/ext4max/fromSsd/Documents/trends/work/2025/misc/ENL-wFC/out082626/09/mancout_lin2b';


%% ICA parameter file
% GIFT parameter file corresponding to the analysis being tested.
%
% For example:
%   ENL analysis:    demo-enl_ica_parameter_info.mat
%   Linear analysis: demo-lin_ica_parameter_info.mat

ica_param_file = '/home/cyrus/ext4max/fromSsd/Documents/trends/work/2025/misc/ENL-wFC/out082626/09/demo-lin_ica_parameter_info.mat';


%% Features
% Features to include in the statistical analysis.
%
% Available feature types include:
%   - spatial maps
%   - timecourses spectra
%   - FNC correlations

features = {'spatial maps'};


%% Component/network definitions
% Cell array with one row per network.
%
% Column 1: Network name
% Column 2: Component indices belonging to that network
%
% Do not assign the same component to multiple network groups.

comp_network_names = {'all', [1:20]};


%% Univariate statistical tests
% If univariate_tests is specified, the multivariate tests are skipped.
%
% Format:
%
%   {'TestType', {data-set indices}, {'Group/condition names'}}
%
% TestType can be:
%
%   'Ttest'   - one-sample or paired t-test
%   'Ttest2'  - independent two-sample t-test
%
% Examples:
%
% One-sample t-test:
%
%   univariate_tests = ...
%       {'Ttest', {(1:50)}, {'Group'}};
%
% Paired t-test:
%
%   univariate_tests = ...
%       {'Ttest', {(2:25), (27:50)}, ...
%       {'Condition 1', 'Condition 2'}};
%
% Independent two-sample t-test:
%
%   Group 1 = subjects 1-3
%   Group 2 = subjects 4-6

univariate_tests = ...
    {'Ttest2', {(1:3), (4:6)}, {'Group 1', 'Group 2'}};


%% Significance threshold
% Statistical significance threshold used by the analysis.

p_threshold = 0.05;


%% Repetition time (TR)
% fMRI repetition time in seconds.

TR = 2;


%% Spatial-map feature settings
%
% sm_center
%   Center the distribution of subject-specific component maps.
%   Options:
%       'yes'
%       'no'
%
% sm_mask
%   External spatial mask used for the analysis.
%   Leave empty to allow MANCOVAN to generate/use its default mask.
%
% stat_threshold_maps
%   Statistic used to threshold spatial maps.
%   Options:
%       'T'
%       'Z'
%
%   When 'T' is selected, the threshold is determined from the
%   t-distribution.
%
% z_threshold_maps
%   Z-score threshold used when stat_threshold_maps = 'Z'.

feature_params.sm_center = 'yes';

feature_params.sm_mask = ...
    '/home/cyrus/ext4max/fromSsd/Documents/trends/work/2025/misc/ENL-wFC/out082626/03/prefix-enlMask.nii';

feature_params.stat_threshold_maps = 'T';

feature_params.z_threshold_maps = 1;


%% Display settings

% Frequency range used when displaying spectral results.
% This setting is not relevant when only spatial maps are analyzed.
display.freq_limits = [0.09, 0.15];

% Structural image used as the anatomical background for spatial-map
% displays.
display.structFile = ...
    fullfile(fileparts(which('gift.m')), 'icatb_templates', 'ch2bet.nii');

% T-statistic threshold used when displaying spatial maps.
display.t_threshold = 2;

% P-value threshold used when displaying univariate maps, spectra, etc.
display.p_threshold = 0.05;

% Display both positive and negative image values.
display.image_values = 'Positive and Negative';

% Multiple-comparison correction for displayed univariate results.
% Common options include:
%   'fdr'
%   'none'
display.threshdesc = 'fdr';

% Display FNC connectogram:
%   1 = display
%   0 = do not display
display.display_connectogram = 0;
