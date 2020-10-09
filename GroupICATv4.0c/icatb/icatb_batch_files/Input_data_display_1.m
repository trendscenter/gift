% Input data for displaying components using display methods like
% component explorer, composite viewer or orthogonal viewer. The function
% that accepts the input file is icatb_batch_display.

% Parameters must be valid characters or integers


%%%%%%% fMRI data and component data %%%%%%%%

% This variable requires fMRI images directory and file pattern.
% and is used only for orthogonal viewer display method
sourceDir = 'C:\MATLAB6p5p2\work\Example Subjects\Visuomotor_data\sub01_vis\';

% fMRI file pattern
sourceFilePattern = 'ns*.img';


% Ouput directory where the component images are:
outputDir = 'D:\test_GIFT\Multiple_sub_Multiple_sess\';

% File pattern 
compFilePattern = '*mean*comp*s1*.img';

% Component Numbers
% [] means all components are selected for component explorer or
% specify the components in a string like [1 3 4] or [1:10]
%compNumbers = [];
compNumbers = [];

% anatomical file used for displaying component images:
structFile = which('nsingle_subj_T1_2_2_5.nii');

%%%%%%% fMRI data and component data %%%%%%%%


%%% Display Parameters %%%%

% Image Values
% Options are 1, 2, 3, 4
% 1 - 'Positive and Negative'
% 2 - 'Positive'
% 3 - 'Absolute Value'
% 4 - 'Negative'
returnValue = 2;

% Convert to z-scores

% Options are 1 or 0
% 1 means convert images to z-scores
% 0 means don't convert images to z-scores
convertToZ = 1;

% Threshold
thresholdValue = 1.0;

% Images per figure: (variable used only for component explorer)
% Options are 1, 4, 9, 16 or 25
imagesPerFigure = 1;

% Slice or anatomical plane (used for component explorer or composite
% viewer)
% Options are 'axial', 'sagittal' or 'coronal'
anatomicalPlane = 'axial';

%%%%%% End for display parameters %%%%%%

% display type:
% Options are 'component explorer', 'composite viewer', 'orthogonal viewer'
displayType = 'component explorer';