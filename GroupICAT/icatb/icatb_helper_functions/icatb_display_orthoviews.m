function icatb_display_orthoviews(anatomical_file, component_files, displayParameters)
%% Display Orthoviews of components
%
% Inputs:
% 1. anatomical_file - Anatomical file in .img or .nii format
% 2. component_files - Component file/files in .img or .nii format
% 3. displayParameters - Display parameters structure. Fields are as
% follows
%   a. image_values - Image values. Options are 'Positive and Negative',
%   'Positive', 'Absolute value' and 'Negative'
%   b. convert_to_zscores - COnvert to z-scores. Options are 'Yes' and No'
%   c. threshold - Threshold value
%

icatb_defaults;
global SUBJECT_ICA_AN3_FILE;
global COMPONENT_NAMING;
global IMAGE_VALUES
global CONVERT_Z
global THRESHOLD_VALUE

%% Select anatomical file
if (~exist('anatomical_file', 'var'))
    anatomical_file = icatb_selectEntry('typeEntity', 'file', 'title', 'Select anatomical file ...', 'filter', '*.img;*.nii', ...
        'typeselection', 'single', 'fileType', 'image', 'filenumbers', 1);
end

if (isempty(anatomical_file))
    error('Anatomical file is not selected');
end

[startPath, f_name, extn] = fileparts(icatb_parseExtn(anatomical_file));

% check the image extension
if ~strcmpi(extn, '.nii') && ~strcmpi(extn, '.img')
    error('Structural image should be in NIFTI or Analyze format');
end

compNaming = ['*', SUBJECT_ICA_AN3_FILE, COMPONENT_NAMING, '*.img;*', SUBJECT_ICA_AN3_FILE, COMPONENT_NAMING, '*.nii'];

%% Select component files
if (~exist('component_files', 'var'))
    component_files = icatb_selectEntry('typeEntity', 'file', 'title', 'Select component files ...', 'filter', compNaming, ...
        'typeselection', 'multiple', 'fileType', 'image', 'filenumbers', '');
end

if (isempty(component_files))
    error('Component files are not selected');
end

component_files = icatb_rename_4d_file(component_files);

%% Select component parameters
if (~exist('displayParameters', 'var') || isempty(displayParameters))
    
    numParameters = 1;
    
    if (~exist('image_values', 'var'))
        opts = char('Positive', 'Positive and Negative', 'Absolute Value', 'Negative');
        inputText(numParameters).promptString = 'Select image values';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = opts;
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'image_values';
        inputText(numParameters).enable = 'on';
    end
    
    if (~exist('convert_to_zscores', 'var'))
        numParameters = numParameters + 1;
        inputText(numParameters).promptString = 'Do you want to convert to z-scores?';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = char('Yes', 'No');
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'convert_to_zscores';
        inputText(numParameters).enable = 'on';
    end
    
    if (~exist('threshold', 'var'))
        numParameters = numParameters + 1;
        inputText(numParameters).promptString = 'Enter threshold';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = '1';
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'threshold';
        inputText(numParameters).enable = 'on';
    end
    
    
    
    if (exist('inputText', 'var'))
        answers = icatb_inputDialog('inputtext', inputText, 'Title', 'Component Parameters', 'handle_visibility',  'on', 'windowStyle', 'modal');
        if (isempty(answers))
            error('Component parameters are not selected');
        end
        for nA = 1:length(answers)
            val = answers{nA};
            eval([inputText(nA).tag, '=val;']);
        end
    end
    
else
    
    image_values = IMAGE_VALUES;
    convert_to_zscores = CONVERT_Z;
    threshold = str2num(THRESHOLD_VALUE);
    
    try
        image_values = displayParameters.image_values;
    catch
    end
    
    try
        convert_to_zscores = displayParameters.convert_to_zscores;
    catch
    end
    
    try
        threshold = displayParameters.threshold;
    catch
    end
    
end


helpH = helpdlg('Displaying Orthviews', 'Orthviews');
drawnow;


%% Display component files
component_files = cellstr(component_files);

for nF = 1:length(component_files)
    [dd, fN, extn] = fileparts(component_files{nF});
    [extn, num] = icatb_parseExtn(extn);
    
    % check the image extension
    if ~strcmpi(extn, '.nii') && ~strcmpi(extn, '.img')
        error(['File ', component_files{nF}, ' should be in NIFTI or Analyze format']);
    end
    gH(nF).H = icatb_orth_views(component_files{nF}, 'structfile', anatomical_file, 'threshold', threshold, 'image_values', image_values, ...
        'convert_to_zscores', convert_to_zscores, 'fig_title', [fN, '_comp_', icatb_returnFileIndex(num)], 'set_to_max_voxel', 1);
    drawnow;
end

icatb_plotNextPreviousExitButtons(gH);

try
    delete(helpH);
catch
end