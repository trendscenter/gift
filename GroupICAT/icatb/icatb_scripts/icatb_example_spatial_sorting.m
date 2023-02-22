clear templateFiles; clear dispParameters;

% Specify ICA parameter file
param_file = 'E:\test_GIFT\Multiple_sub_Multiple_sess\Driving_ica_parameter_info.mat';

% Options for selectedStr are:
% 1.'Same set of spatial templates for all data-sets'
% 2. 'Different set of spatial templates for sessions'
% 3. 'Different set of spatial templates for subjects and sessions'
selectedStr = 'Different set of spatial templates for subjects and sessions';

% Count is the current data-set number
% Enter template files for subject 1 sessions followed by subject 2
% sessions, etc.

%%%% Note: Number of spatial templates must be the same between data-sets

count = 1; % Enter subject 1 session 1 templates
templateFiles(count).name = str2mat('E:\test_GIFT\Multiple_sub_Multiple_sess\Driving_mean_component_ica_s1_004.img', ...
                            'E:\test_GIFT\Multiple_sub_Multiple_sess\Driving_mean_component_ica_s1_007.img');                        
count = count + 1; % Enter subject 1 session 2 templates
templateFiles(count).name = str2mat('E:\test_GIFT\Multiple_sub_Multiple_sess\Driving_mean_component_ica_s2_004.img', ...
                            'E:\test_GIFT\Multiple_sub_Multiple_sess\Driving_mean_component_ica_s2_007.img');
                        
count = count + 1; % Enter subject 2 session 1 templates
templateFiles(count).name = str2mat('E:\test_GIFT\Multiple_sub_Multiple_sess\Driving_mean_component_ica_s1_003.img', ...
                            'E:\test_GIFT\Multiple_sub_Multiple_sess\Driving_mean_component_ica_s1_008.img');
count = count + 1; % Enter subject 2 session 2 templates
templateFiles(count).name = str2mat('E:\test_GIFT\Multiple_sub_Multiple_sess\Driving_mean_component_ica_s2_003.img', ...
                            'E:\test_GIFT\Multiple_sub_Multiple_sess\Driving_mean_component_ica_s2_008.img');


count = count + 1; % Enter subject 3 session 1 templates
templateFiles(count).name = str2mat('E:\test_GIFT\Multiple_sub_Multiple_sess\Driving_mean_component_ica_s1_002.img', ...
                            'E:\test_GIFT\Multiple_sub_Multiple_sess\Driving_mean_component_ica_s1_009.img');
count = count + 1; % Enter subject 3 session 2 templates
templateFiles(count).name = str2mat('E:\test_GIFT\Multiple_sub_Multiple_sess\Driving_mean_component_ica_s2_005.img', ...
                            'E:\test_GIFT\Multiple_sub_Multiple_sess\Driving_mean_component_ica_s2_010.img');                        

                        
%%%%%%% Specify Display parameters %%%%%%%%%

% Options for image values are
% 1 means positive and negative
% 2 means positive
% 3 means Absolute
% 4 means Negative
dispParameters.imagevalues = 2; 

% Anatomical plane options are axial, sagital, coronal
dispParameters.anatomicalplane = 'axial'; 

% slices in mm (vector of real numbers)
dispParameters.slicerange = [0:4:72]; 

% Number of images per figure
% Options are 1, 4, 9, 16, 25
dispParameters.imagesperfigure = 4;

% Convert to z scores:
% 1 means convert, 0 means don't convert to z-scores
dispParameters.convertToZ = 1;     

% Z Threshold: 
dispParameters.thresholdvalue = 1; 

% Anatomical file for overlaying component images:
% Image used from icatb/icatb_templates
dispParameters.structFile = which('nsingle_subj_T1_2_2_5.nii'); 

%%%%%%%% End for specifying display parameters %%%%%%%%%%

% Call spatial sorting function
icatb_spatialSorting(param_file, selectedStr, templateFiles, dispParameters);
