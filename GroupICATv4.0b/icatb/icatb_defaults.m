function icatb_defaults
% All the GIFT toolbox defaults are stored in this file

if (isdeployed)
    
    gicaInstallDir = getenv('GICA_INSTALL_DIR');
    defsFile = fullfile(gicaInstallDir, 'icatb_defaults.m');
    
    if (~isempty(dir(defsFile)))
        icatb_eval_script(defsFile);
        return;
    end
    
end

%% Colormap info
global COLORMAP_FILE;
global COLORLIST;

%% Conventions for naming components, timecourses and structural
global COMPONENT_NAMING;
global TIMECOURSE_NAMING;

%% Matlab files
global PARAMETER_INFO_MAT_FILE; %Holds information for group session parameters
global DATA_REDUCTION_MAT_FILE; %A file for each subject's principal components
global ICA_MAT_FILE;
global BACK_RECONSTRUCTION_MAT_FILE;
global CALIBRATE_MAT_FILE;

%% Analyze files
global SUBJECT_ICA_AN3_FILE;
global MEAN_AN3_FILE;
global TMAP_AN3_FILE;
global STD_AN3_FILE;
global MEAN_ALL_AN3_FILE; % Mean for different sessions
global SESSION_POSTFIX;
global AGGREGATE_AN3_FILE;

%% File Indices
global SUBJECT_ICA_INDEX;
global MEAN_INDEX;
global TMAP_INDEX;
global STD_INDEX;
global MEAN_ALL_INDEX;  % Index for mean over different sessions

%% Screen Color Defaults
global BG_COLOR; % Figure background color
global BG2_COLOR; % User Interface background color except pushbutton
global FG_COLOR;
global AXES_COLOR;
global FONT_COLOR; % User Interface foreground color except pushbutton
global LEGENDCOLOR;
global BUTTON_COLOR; % BUTTON background color
global BUTTON_FONT_COLOR; % BUTTON foreground color
global HELP_FONT_COLOR;

%% Pictures
global COMPOSITE_VIEWER_PIC;
global COMPONENT_EXPLORER_PIC;
global ORTHOGONAL_VIEWER_PIC;

%% Display defaults
global SORT_COMPONENTS;
global IMAGE_VALUES;
global CONVERT_Z;
global THRESHOLD_VALUE;
global IMAGES_PER_FIGURE;
global COMPLEX_IMAGES_PER_FIGURE; % for complex data
global ANATOMICAL_PLANE;

%% Slice Range defaults
global USE_DEFAULT_SLICE_RANGE;
global SLICE_RANGE;

global ASPECT_RATIO_FOR_SQUARE_FIGURE;
global WS;
global MIN_SCREEN_DIM_IN_PIXELS;
global S0;

%% FONT DEFAULTS
global UI_FONTNAME;
global UI_FONTUNITS;
global UI_FS;

%% FUNTIONAL FILE DEFAULTS
% Functional data filter is used as filter in figure windows while reading
% images. This is also the keyword used in writing images (.nii or .img)
global FUNCTIONAL_DATA_FILTER;


%% Timecourse defaults
global DETRENDNUMBER;
global SMOOTHINGVALUE;
global SMOOTHPARA;

% Flag for displaying acknowledgements of ICA authors
global FLAG_ACKNOWLEDGE_CREATORS;

% ICA options window display default (displays ICA Options)
global ICAOPTIONS_WINDOW_DISPLAY;

% declare vars for storing directory information
global STORE_DIRECTORY_INFORMATION;
global DIRS_TO_BE_STORED;

% method for entering regressors
global METHOD_ENTERING_REGRESSORS;

% text file for entering regressors
global TXTFILE_REGRESSORS;

% Including a global variable to flip images
global FLIP_ANALYZE_IMAGES;

% variable for displaying text for slices
global TEXT_DISPLAY_SLICES_IN_MM;

%% Variables to detect the reading and writing of complex images
global READ_NAMING_COMPLEX_IMAGES;
global WRITE_NAMING_COMPLEX_IMAGES;

%% Event average defaults
global EVENTAVG_WINDOW_SIZE;
global EVENTAVG_INTERP_FACTOR;

% Defaults for printing regressors to text file
global PRINTTYPE_REGRESSORS;

% Zip image files default
global ZIP_IMAGE_FILES;

% default for scaling
global SCALE_DEFAULT;

% Default mask option
global DEFAULT_MASK_OPTION;

% Remove constant voxels
global REMOVE_CONSTANT_VOXELS;

% Default mask multiplier in SBM
global DEFAULT_MASK_SBM_MULTIPLIER;

% Number of times ICA must Run
global NUM_RUNS_GICA;

% Open display GUI
global OPEN_DISPLAY_GUI;

global fMRI_STATS_T;


% EEG defaults
global EEG_RMBASE;

%% EEG display_defaults
global EEG_IMAGE_VALUES;
global EEG_CONVERT_Z;
global EEG_THRESHOLD_VALUE;
global EEG_IMAGES_PER_FIGURE;
global EEG_TOPOPLOT_COLORMAP;

%% Talairach defaults
global TALAIRACH_DIST;
global TALAIRACH_THRESHOLD;

%% SPMS STATS defaults
% Compute one sample t-test and write talairach coords
global SPM_STATS_WRITE_TAL;

% Two sample t-test defaults
global SPM_STATS_TTEST2_EXPLICIT_MASK;
global SPM_STATS_TTEST_THRESHOLD;

% Average runs for subjects for all components using SPM5
global SPM_STATS_AVG_RUNS;

%% Experimental TR (in seconds for fmri data)
global EXPERIMENTAL_TR;

%% Timecourses Postprocess options (FNC, spectra)
global TIMECOURSE_POSTPROCESS;

%% FFT defaults for group comparison
global NPOINT_FFT_GROUP_COMPARISON;
global NUM_BINS_GROUP_COMPARISON;
global DEFAULT_TR_SPECTRAL_GROUP_COMPARE;

%% Enforce MAT file version
global ENFORCE_MAT_FILE_VER;

%% Center Image Distribution
global CENTER_IMAGES;

%% Data pre-processing default
global PREPROC_DEFAULT;

%% PCA Type default
global PCA_DEFAULT;

%% Back Reconstruction default
global BACKRECON_DEFAULT;

%% MAX Available RAM
global MAX_AVAILABLE_RAM;

%% Write analysis steps in dirs
global WRITE_ANALYSIS_STEPS_IN_DIRS;

%% Conserve disk space
global CONSERVE_DISK_SPACE;

%% Mancova defaults
global MANCOVA_DEFAULTS;

%% SBM defaults
global SBM_Z_SCORES_RMMEAN;

%% COMPOSITE COLOR defaults
global USE_UNIFORM_COLOR_COMPOSITE;

%% Connectogram options
global CONNECTOGRAM_OPTIONS;

%% Minimization function used in despike
global DESPIKE_SOLVER;

%% NIFTI GZ files
global NIFTI_GZ;
global GZIPINFO;

%% Dimensionality estimation options
global DIM_ESTIMATION_OPTS;

%% Shuffle random numbers
global RAND_SHUFFLE;

modalityType = icatb_get_modality;

%% Naming Convention Output Files( Analyze Format)
COMPONENT_NAMING = '_component_ica_';

if (~strcmpi(modalityType, 'smri'))
    TIMECOURSE_NAMING = '_timecourses_ica_';
else
    TIMECOURSE_NAMING = '_loading_coeff_';
end

%% Naming Convention Output Files( Matlab Format)
PARAMETER_INFO_MAT_FILE = '_ica_parameter_info';
DATA_REDUCTION_MAT_FILE = '_pca_r';
ICA_MAT_FILE = '_ica';
BACK_RECONSTRUCTION_MAT_FILE = '_ica_br';
CALIBRATE_MAT_FILE = '_ica_c';

% Specific Naming Convention for Component Files( Analyze Format)
SESSION_POSTFIX = 's';

% aggregate file naming
AGGREGATE_AN3_FILE = ['_agg_', COMPONENT_NAMING];

% subject component naming
if (~strcmpi(modalityType, 'smri'))
    SUBJECT_ICA_AN3_FILE = '_sub';
else
    SUBJECT_ICA_AN3_FILE = '_group';
end

%% Mean file naming
MEAN_AN3_FILE = '_mean'; % mean for a session
MEAN_ALL_AN3_FILE = '_mean'; % mean for all data sets

% tmap file naming
TMAP_AN3_FILE = '_tmap'; % tmap for a session

% std file naming
STD_AN3_FILE = '_std'; % std for a session

%% Indicies For Component Files
MEAN_ALL_INDEX = 1; % Mean for all sessions over subjects
MEAN_INDEX = 2;
TMAP_INDEX = 3;
STD_INDEX = 4;
SUBJECT_ICA_INDEX = 5; %5 to (5 + number of subjects)


%% Images of the composite viewer, component explorer, Orthoviewer
COMPONENT_EXPLORER_PIC = which('component_explorer_pic.tif');
ORTHOGONAL_VIEWER_PIC = which('orthogonal_viewer_pic.tif');
COMPOSITE_VIEWER_PIC =  which('composite_viewer_pic.tif');

% Name of file with colormaps
COLORMAP_FILE = 'icatb_colors';

%% Filter for selecting Functional Data.
% This is also the image extension used for writing component images
% options are ('*.img' or '*.nii')
% *.img - analyze
% *.nii - Nifti
FUNCTIONAL_DATA_FILTER = '*.nii';

% List of colors for plotting
COLORLIST = ['g','r','c','m','y'];

isMacOS = ~isempty(strfind(lower(mexext), 'mac'));

%% Screen Color Defaults
BG_COLOR = [0 0 0]; % Figure background color
BG2_COLOR = [.2 .2 .2]; % User Interface controls background color except push button
FG_COLOR = [.562 .562 .562];
FONT_COLOR = [1 1 1]; % User Interface controls foreground color except push button
AXES_COLOR = [0.25 0.25 0.25];
LEGENDCOLOR = [1 1 1];

%% Button background and font colors
BUTTON_COLOR = [.2 .2 .2]; %(default is black)
BUTTON_FONT_COLOR = [1 1 1]; %(default is white)
HELP_FONT_COLOR = [1 1 0]; % default is yellow

if (isMacOS)
    %% White background and black foreground for MAC OS
    BG_COLOR = [1, 1, 1];
    BG2_COLOR = BG_COLOR;
    BUTTON_COLOR = BG_COLOR;
    FG_COLOR = [0, 0, 0];
    FONT_COLOR = FG_COLOR;
    BUTTON_FONT_COLOR = FG_COLOR;
end

%% Display defaults
SORT_COMPONENTS = 'No'; % options are 'No' and 'Yes'
IMAGE_VALUES = 'Positive'; % options are 'Positive and Negative', 'Positive', 'Absolute Value' and 'Negative'
CONVERT_Z = 'Yes'; % options are 'No' and 'Yes'
THRESHOLD_VALUE = '1.0';
IMAGES_PER_FIGURE = '4'; % options are '1', '4', '9', '16' and '25'
COMPLEX_IMAGES_PER_FIGURE = '2';
ANATOMICAL_PLANE = 'Axial'; % options are 'Axial', 'Sagittal' and 'Coronal'

%% Slice Range defaults:
% Options are 1 or 0
% 1 means default slice range mentioned in variable SLICE_RANGE will be used.
% 0 means slice range will be calculated based on the data
USE_DEFAULT_SLICE_RANGE = 0;
SLICE_RANGE = '0:4:72';


% ASPECT RATIO TO MAKE SQUARE DIMENSIONS REALLY SQUARE
% WS = ASPECT RATIO OF THE COMPUTER YOUR LOOKING AT NOW COMPARED TO THE
% COMPUTER THAT GIFT WAS DELELOPED ON
% dimensions for square figure on development computer(in pixels)
dimForSquareFigure = [50 50 853 800];
S0   = get(0,'ScreenSize');
WS = [S0(3)/1280 (S0(4))/960 S0(3)/1280 (S0(4))/960];

xDiff = (dimForSquareFigure(3)*WS(3));
yDiff = (dimForSquareFigure(4)*WS(4));
xRatio = xDiff/yDiff;
yRatio = yDiff/xDiff;
if(xRatio >1)
    xRatio =1;
else
    yRatio =1;
end
%ASPECT_RATIO_FOR_SQUARE_FIGURE = [xRatio yRatio];
ASPECT_RATIO_FOR_SQUARE_FIGURE = [1, 1];
WS = WS;
MIN_SCREEN_DIM_IN_PIXELS=min([S0(3) S0(4)]);
PERCENT_SCREEN_OCUPIED = .9;

% Font Size
UI_FONTNAME = 'times';
UI_FONTUNITS = 'points';
UI_FS = 12;

% Sorting detrend defaults:
% DETRENDNUMBER - case 0 - Removes the mean
% DETRENDNUMBER - case 1 - Removes the mean and linear trend
% DETRENDNUMBER - case 2 - Uses sine and cosine one cycle, removes mean and
% linear trend
% DETRENDNUMBER - case 3 - Uses sine two cycles, cosine two cycles, sine one cycle, cosine one cycle,
% removes mean and inear trend

if (strcmpi(modalityType, 'fmri'))
    DETRENDNUMBER = 3;
else
    DETRENDNUMBER = 0;
end

%% Smoothing Defaults
% If you want to smooth time courses - Replace smooth para by yes(lower
% case)
SMOOTHPARA = 'No';
SMOOTHINGVALUE = 1.1;

% Acknowledge creators display ('off' turns off the acknowledgement display)
FLAG_ACKNOWLEDGE_CREATORS = 'off';

% ICA options window display ('off' turns off the options display)
ICAOPTIONS_WINDOW_DISPLAY = 'on';

% Options for storing subject directories in Directory History popup
% control in file selection window.
STORE_DIRECTORY_INFORMATION = 'No'; % flag for storing directories
% Enter the directories to be stored between curly brackets (use comma to
% separate the directories)
DIRS_TO_BE_STORED = {'C:\MATLAB6p5\work\Visuomotor_Data\visomot', 'C:\MATLAB6p5\work\Example Subjects'};

%% Method for entering regressors
% options are GUI, BATCH or AUTOMATIC
% 1. GUI - Means GUI will be used to enter regressors.
% 2. BATCH - Means text file specified in the variable TXTFILE_REGRESSORS (see below) will be used.
% 3. AUTOMATIC - Means the regressors will be used automatically and if the
% design matrix is session specific then session related regressors will be
% used. For correlation first regressor for the specified design matrix
% will be used.

% if BATCH is specified (text file can be entered by typing full file path
% for TXTFILE_REGRESSORS variable or by selecting 'Load file for temporal
% sorting' under DISPLAY GUI Options menu).
METHOD_ENTERING_REGRESSORS = 'GUI';

% text file for entering regressors. Default will not be used if Load file for temporal sorting is
% used in DISPLAY GUI OPTIONS MENU.
TXTFILE_REGRESSORS = which('Input_data_regressors_1.m');

%% Flip parameter for analyze images
FLIP_ANALYZE_IMAGES = 0;

%% Number of times ICA will be run
NUM_RUNS_GICA = 1;

% Text for showing the slices on component explorer and composite viewer
% options are 'on' and 'off'
TEXT_DISPLAY_SLICES_IN_MM = 'on';

% event average window size in seconds
EVENTAVG_WINDOW_SIZE = 30;

% event average interpolation factor
EVENTAVG_INTERP_FACTOR = 5;

% default Print Type for regression parameters
% options are: row_wise, column_wise
PRINTTYPE_REGRESSORS = 'row_wise';

%% Zip the image files default
% Options are No, yes
% Zips the image files by viewing set like subject 1 session 1 components,
% Mean for session 1 components, etc.

if (strcmpi(modalityType, 'fmri'))
    ZIP_IMAGE_FILES = 'No';
else
    ZIP_IMAGE_FILES = 'No';
end

%% Options for scaling are in icatb_scaleICA
% Options are 0, 1, 2, 3 and 4
% 0 - Don't scale components
% 1 - Convert components to percent signal change
% 2 - Convert components to z-scores
% 3 - Normalize independent component maps and multiply the maximum value of the component maps to timecourses.
% 4 - Scale component maps using standard deviation of timecourses and
% timecourses using maximum value of the component maps.
SCALE_DEFAULT = 2;


%% Mask Options (!!!! Use these settings before selecting the data)
% Mask is calculated by doing a Boolean AND of voxels that surpass or equal the mean.
% Options are 'all_files' and 'first_file'.
% 1. 'all_files' - Uses all files.
% 2. 'first_file' - Uses first file of each subject session.

DEFAULT_MASK_OPTION = 'first_file';

% Remove constant voxels from fmri. A value of 1 will remove constant
% voxels using all scans of all subjects.
REMOVE_CONSTANT_VOXELS = 0;

% Default mask multiplier in source based morphometry (SBM). Defaults is 1% of
% mean. Voxels >= 1% of mean will be used.
DEFAULT_MASK_SBM_MULTIPLIER = 0.01;

% Variable for reading complex images
% Complex images can be of real&imaginary or magnitude&phase type
% For real&imaginary first field of variable is read and for
% magnitude&phase second field of variable is read.
% Naming should contain _ term to distinguish images (R_ means the real
% images start with R_ and _R means images end with R)
READ_NAMING_COMPLEX_IMAGES.real_imag = {'R_', 'I_'};
READ_NAMING_COMPLEX_IMAGES.mag_phase = {'Mag_', 'Phase_'};
%
% % writing complex data
% % if the data is of type real_imag then first set is used
% % else second set is used
WRITE_NAMING_COMPLEX_IMAGES.real_imag = {'R_', 'I_'};
WRITE_NAMING_COMPLEX_IMAGES.mag_phase = {'Mag_', 'Phase_'};

% Same as defaults.stats.fmri.t used while calculating HRF
fMRI_STATS_T = 30;

%% Open display GUI after analysis
% Options are 1 and 0.
% 1 means open display GUI
OPEN_DISPLAY_GUI = 1;


%% EEG defaults

% Remove baseline default (Used in smooth trials input dialog box)
EEG_RMBASE = [-100 0];

%% EEG display defaults
EEG_IMAGE_VALUES = 'Positive and Negative'; % options are 'Positive and Negative', 'Positive', 'Absolute Value' and 'Negative'
EEG_CONVERT_Z = 'Yes'; % options are 'No' and 'Yes'
EEG_THRESHOLD_VALUE = '1.0';
EEG_IMAGES_PER_FIGURE = '4'; % options are '1', '4', '9'
EEG_TOPOPLOT_COLORMAP = jet(64);

%% Talairach defaults %
% Distance between contiguous voxels in mm
TALAIRACH_DIST = 4;
% Threshold. This pulls out positive and negative areas by applying data
% > abs(threshold) and -data > abs(threshold)
TALAIRACH_THRESHOLD = 3.5;

%% SPM5 stats defaults

% Compute one sample t-test on subject component maps using SPM5 and write
% talariach coords for the t-map.
%
% NOTE:
% a. You need SPM5 for this utility.
% b. Set ZIP_IMAGE_FILES as 'no' if you want this script to run faster.

% Options are:
% 0 - Don't compute one sample t-test and don't write talairach coords.
% 1 - Compute one sample t-test and don't write talairach coords.
% 2 - Compute one sample t-test and write talairach coordinates for the
% t-map.
SPM_STATS_WRITE_TAL = 0;


% Explicit mask when two sample t-test is used.
%
% 1 - A value of 1 means one sample t-test maps are
% thresholded based on SPM_STATS_TTEST_THRESHOLD and the values lying greater than
% threshold are used as a mask.
%
% 0 - A value of 0 means no explicit mask is used.

SPM_STATS_TTEST2_EXPLICIT_MASK = 1;

% Threshold on one sample t-test maps that will be used as a explicit mask
% when two sample t-test is done.
%
% NOTE: threshold criteria here used is abs(data) >= threshold.
SPM_STATS_TTEST_THRESHOLD = 1.5;


% Average runs for subjects for all components using spm_imcalc
% This will be computed automatically when you run group stats step and if you set SPM_STATS_WRITE_TAL to 1
% or 2.
SPM_STATS_AVG_RUNS = 1;

%% Experimental TR in seconds
EXPERIMENTAL_TR = 1;

%% Timecourse postprocess (FNC and spectra info)
TIMECOURSE_POSTPROCESS.write = 1; % 1 - Write FNC correlations and spectra info.
% spectra parameters
TIMECOURSE_POSTPROCESS.spectra.tapers = [3, 5];
TIMECOURSE_POSTPROCESS.spectra.sampling_frequency = 1/min(EXPERIMENTAL_TR);
TIMECOURSE_POSTPROCESS.spectra.frequency_band = [0, 1/(2*min(EXPERIMENTAL_TR))];
TIMECOURSE_POSTPROCESS.spectra.freq_limits = [0.1, 0.15];
% FNC parameters
TIMECOURSE_POSTPROCESS.fnc.despike_tc = 1; % 1 - Despike timecourses
TIMECOURSE_POSTPROCESS.fnc.cutoff_frequency = 0.15; % High frequency cutoff in Hz

%% FFT For Group Comparison
% N point FFT used for group comparison
NPOINT_FFT_GROUP_COMPARISON = 300;

% No. of frequency bins used for group comparison
NUM_BINS_GROUP_COMPARISON = 6;

% Default TR for spectral group comparision in seconds
DEFAULT_TR_SPECTRAL_GROUP_COMPARE = EXPERIMENTAL_TR;

%% Enforce MAT files versioning for MATLAB versions greater than 6.5. Use
% the correct option. For more help on version compatibility, please check
% MATLAB save command options for MAT files.
ENFORCE_MAT_FILE_VER = '-v6';

%% Center Image distribution
% Options are 0 and 1
% 1 means image distribution is centered
CENTER_IMAGES = 1;

%% Data Pre-processing options
% 1 - Remove mean per time point
% 2 - Remove mean per voxel
% 3 - Intensity normalization
% 4 - Variance normalization
PREPROC_DEFAULT = 1;

%% PCA Type
% 1 - Standard
% 2 - Expectation Maximization
% 3 - SVD
% 4 - MPOWIT
% 5 - STP
PCA_DEFAULT = 1;

%% Backreconstruction options.
% 2 - Spatial-temporal Regression
% 4 - GICA
BACKRECON_DEFAULT = 4;

%% Maximum available memory (RAM) in gigabytes.
% This variable will be used during data reduction or Group PCA step
% Default is set to 4 GB.
MAX_AVAILABLE_RAM = 4;


%% Set both WRITE_ANALYSIS_STEPS_IN_DIRS and CONSERVE_DISK_SPACE before setting up analysis.

% Organize analysis results in directories. Options are 0 and 1. The analysis directories are:
% a. Data reduction files - prefix_data_reduction_files
% b. ICA files - prefix_ica_files
% c. Back-reconstruction files - prefix_back_reconstruction_files
% d. Scaling components files - prefix_scaling_components_files
% e. Group stats files - prefix_group_stats_files
WRITE_ANALYSIS_STEPS_IN_DIRS = 0;

% Conserve disk space. Options are:
% 0 - Write all analysis files including intermediate files (PCA, Backreconstruction, scaled component MAT files)
% 1 - Write only necessary files required to resume the analysis. The files written are as follows:
%   a. Data reduction files - Only eigen vectors and eigen values are written in the first data reduction step. PCA components are written at the last reduction stage.
%   b. Back-reconstruction files - Back-reconstruction files are not written to the disk. The information is computed while doing scaling components step
%   c. Scaling components files - Scaling components MAT files are not written when using GIFT or SBM.
%   d. Group stats files - Only mean of all data-sets is written.
% 2 - Write all files till the group stats. Cleanup intermediate files at the end of the group stats (PCA, Back-reconstruct, Scaled component MAT files in GIFT). Analysis cannot be
% resumed if there are any changes to the setup parameters. Utilities that work with PCA and Backreconstruction files like Remove components, Percent Variance, etc won't work with this option .
CONSERVE_DISK_SPACE = 0;

%% Mancova defaults
%

% Default mask options:
% Standard deviation used when computing voxels using t-map statistic
MANCOVA_DEFAULTS.sm.t_std = 4;
% Image values to include after computing t-cutoff.
% 1 - Both positive and negative voxels
% 2  - Positive voxels
% 3 - Negative voxels
MANCOVA_DEFAULTS.sm.t_image_values = 2;

% FNC lag info used when computing FNC correlations (Same as in the FNC toolbox).
MANCOVA_DEFAULTS.fnc.lag = 3;
% shift resolution
MANCOVA_DEFAULTS.fnc.shift_resolution = 25;

% domain average (static fnc only)
MANCOVA_DEFAULTS.fnc.domain_average = 1;

%% SBM Defaults
% Remove mean of ICA loading coefficients
% Options are 0 and 1
SBM_Z_SCORES_RMMEAN = 1;


%% Colormap defaults
% 1 - Uniform color is used when displaying composite plots in network
% summary
% 0 - Intensity map is used.
USE_UNIFORM_COLOR_COMPOSITE = 0;

%% CONNECTOGRAM options
% CONNECTOGRAM SPATIAL MAPS WIDTH. If the spatial maps are smaller in size, you could increase it to a value like 0.08. By default, size is automatically determined.
CONNECTOGRAM_OPTIONS.sm_width = [];
% CONNECTOGRAM figure Position in pixels. Enter [xorigin, yorigin, figure_width, figure_height]. If this variable is set to empty, figure size is set to almost the screen size.
CONNECTOGRAM_OPTIONS.fig_pos = [];
% Radius at which spatial maps are plotted in a circle. By default, spatial maps are plotted in a circle with radius = 1.7
CONNECTOGRAM_OPTIONS.radius_sm = [];

%% Minimization function used in despike. Options are 'lsqcurvefit', 'fminsearch' and 'fminunc'
% lsqcurvefit - Faster compared to the other two methods
% fminsearch - Slower and requires general matlab version
% fminunc - Requires optimization toolbox
DESPIKE_SOLVER = 'lsqcurvefit';


%% NIFTI GZ file. Options are 0 and 1. By default, .nii.gz files are read using gzipinputstream. If set to 1, files are un-archived to temp directory and read using SPM.
NIFTI_GZ = 0;

%% GZIPInfo - Additional info 
% Option 0:
% isLargeFile - Options are 0 and 1. Option 0 assumes that enough heap space is
% there. For large data, use a value of 1.
% Option 1:
% tempdir - Temporary directory where gzip files can be unzipped if you set
% NIFTI_GZ = 1

% buffer size - 2^14 bytes
GZIPINFO.isLargeFile = 1;
GZIPINFO.buffer_size = 2^14;
GZIPINFO.tempdir = tempdir;

%% MDL Dimensionality estimation
% iid_sampling - Option are 1 and 0. By default, i.i.d sampling is run on the data. If set to 0, please provide data
% smoothness kernel to determine no. of i.i.d samples. 
DIM_ESTIMATION_OPTS.iid_sampling = 1;
DIM_ESTIMATION_OPTS.fwhm = [5, 5, 5];

%% shuffle random numbers
RAND_SHUFFLE = 1;
