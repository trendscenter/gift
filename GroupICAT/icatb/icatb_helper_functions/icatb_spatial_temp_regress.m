function icatb_spatial_temp_regress(compFiles, inputFiles, varargin)
%% Backreconstruct subject components using spatial-temporal regression
%
%
% Inputs:
% 1. compFiles - Group or aggregate component files (full path) in character array format.
% 2. inputFiles - Original data files (full path) in a cell array where each cell can contain a character array of images like
% 3D analyze or 3D Nifti.
% 3. varargin - Optional parameters are passed in pairs. Here is the
% description of optional parameters
%   a. outputDir - Output directory
%   b. format - Output format. Options are .nii, .img or .mat
%   c. outputPrefix - Output prefix.
%
%
% Outputs:
% Individual subject component images are written in the output directory in 4D Nifti format or 3D
% analyze format.
%
% Example:
% compFiles = {'E:\Multiple_sub_Multiple_sess\sens\Sensorimotor_agg__component_ica_018.img', 'E:\Multiple_sub_Multiple_sess\sens\Sensorimotor_agg__component_ica_019.img'};
% inputFiles = {'F:\smr1\swSm_nifti.nii', 'F:\smr2\swSm_nifti.nii'};
% icatb_spatial_temp_regress(compFiles, inputFiles, 'outputDir', 'E:\Multiple_sub_Multiple_sess\results', 'format', '.nii', 'outputPrefix', 'test');


try
    warning('off', 'MATLAB:dispatcher:ShadowedMEXExtension');
catch
end

%% Load defaults
icatb_defaults;
global SUBJECT_ICA_AN3_FILE;
global SESSION_POSTFIX;
global COMPONENT_NAMING;
global TIMECOURSE_NAMING;

outputPrefix = 'STR';
outputFormat = '.nii';
dummy_scans = 0;

%% Loop over variables
for n = 1:2:length(varargin)
    if (strcmpi(varargin{n}, 'outputdir'))
        outputDir = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'format'))
        outputFormat =  varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'outputprefix'))
        outputPrefix = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'preproc_type'))
        preproc_type = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'dummy_scans'))
        dummy_scans = varargin{n + 1};
    end
end
%% End of loop over variables

%% Output directory
if (~exist('outputDir', 'var'))
    outputDir = icatb_selectEntry('typeEntity', 'directory', 'title', 'Select output directory to save results ...');
end

if (isempty(outputDir))
    error('Output directory is not selected');
end

%% Component files selection
if (~exist('compFiles', 'var'))
    compFiles = icatb_selectEntry('typeEntity', 'file', 'typeSelection', 'multiple', 'filter', '*.img;*.nii', 'title', 'Select group or aggregate component files ...');
end

if (isempty(compFiles))
    error('Component files are not selected');
end

compFiles = deblank(char(compFiles));

if (size(compFiles, 1) == 1)
    pathstrComp = fileparts(compFiles);
    if (isempty(pathstrComp))
        pathstrComp = pwd;
    end
    tmp = dir(compFiles);
    if (isempty(tmp))
        error(['File ', compFiles, ' doesn''t exit']);
    end
    compFiles = icatb_fullFile('directory', pathstrComp, 'files', char(tmp.name));
    clear tmp;
end

%% Input files selection
if (~exist('inputFiles', 'var'))
    inputFiles = icatb_selectEntry('typeEntity', 'file', 'typeSelection', 'multiple', 'filter', '*.nii', 'title', 'Select original data files of subjects ...');
end

if (isempty(inputFiles))
    error('Original data files are not selected');
end

modalityType = icatb_get_modality;

if (~strcmpi(modalityType, 'smri'))
    if (~exist('preproc_type', 'var'))
        
        numParameters = 1;
        
        inputText(numParameters).promptString = 'Do you want to do Intensity Normalization on the data?';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'Yes', 'No'};
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'preproc_type';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Enter the file numbers to include';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = '';
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'dummy_scans';
        inputText(numParameters).enable = 'on';
        
        answers = icatb_inputDialog('inputtext', inputText, 'Title', 'Input Parameters', 'handle_visibility', 'on', 'windowStyle', 'modal');
        
        if (isempty(answers))
            error('Input dialog box is quit');
        end
        
        if (strcmpi(answers{1}, 'yes'))
            preproc_type = 'intensity normalization';
        else
            preproc_type = 'remove mean per timepoint';
        end
        
        dummy_scans = answers{2};
    end
else
    preproc_type = 'remove mean per timepoint';
    dummy_scans = [];
end


if (isnumeric(preproc_type))
    preprocOpts = icatb_preproc_data;
    preproc_type = lower(preprocOpts{preproc_type});
end

if (~strcmpi(modalityType, 'smri'))
    inputFiles = cellstr(inputFiles);
else
    if (~iscell(inputFiles))
        inputFiles = {inputFiles};
    end
end

if (~exist(outputDir, 'dir'))
    mkdir(outputDir);
end

disp('Components selected are: ');
disp(compFiles);
fprintf('\n');
disp(['Components will be prepended with ', outputPrefix]);
disp(['Components will be saved in ', outputFormat, ' format']);


fprintf('\nBackreconstructing subject components using Spatial-temporal regression ...\n');

fprintf('\n');
disp('Loading group components ...');
[group_icasig, HInfo] = icatb_loadData(compFiles);
V = HInfo.V(1);
clear HInfo;

comp_dims = V.dim(1:3);
mask = find(prod(group_icasig, 4) ~= 0);

if (isempty(mask))
    error('No in brain voxels found');
end

cd(outputDir);

outputDir = pwd;

varInfo.compFiles = compFiles;
varInfo.inputFiles = inputFiles;
varInfo.outputFormat = outputFormat;
varInfo.outputPrefix = outputPrefix;

icatb_save(fullfile(outputDir, [outputPrefix, '_input_parameters.mat']), 'varInfo');

%% Write mask
tmpDat = zeros(comp_dims);
tmpDat(mask) = 1;
V.fname = fullfile(outputDir, [outputPrefix, '_mask.nii']);
V.dt(1) = 4;
V.n(1) = 1;
icatb_write_vol(V, tmpDat);
clear tmpDat;

group_icasig = reshape(group_icasig, [prod(comp_dims), size(group_icasig, 4)]);
group_icasig = group_icasig(mask, :);
numOfIC = size(group_icasig, 2);

%% Loop over subjects
for row = 1:length(inputFiles)
    
    currentFiles = deblank(inputFiles{row});
    
    if (size(currentFiles, 1) == 1)
        [pathstr, fN, extn] = fileparts(currentFiles);
        if (isempty(pathstr))
            pathstr = pwd;
        end
        tmp = dir(currentFiles);
        if (isempty(tmp))
            error(['File ', currentFiles, ' doesn''t exit']);
        end
        currentFiles = icatb_fullFile('directory', pathstr, 'files', char(tmp.name));
    end
    
    firstFile = deblank(currentFiles(1, :));
    
    disp(['Loading file ', firstFile, ' ...']);
    
    [firstFile, number] = icatb_parseExtn(firstFile);
    
    firstV = icatb_spm_vol_nifti(firstFile, number);
    
    check = find((firstV.dim(1:3) == comp_dims) == 1);
    
    if (length(check) ~= 3)
        disp(['Dimensions of file ', deblank(currentFiles(1, :)), ' are different from component dimensions']);
        continue;
    end
    
    currentFiles = icatb_rename_4d_file(currentFiles);
    
    if (length(dummy_scans) > 1)
        fileNum = dummy_scans;
        fileNum(fileNum > size(currentFiles, 1)) = [];
        if isempty(fileNum)
            error('Error:DummyScans', 'Cannot find the files specified with file numbers. Please check dummy_scans variable for \nfile %s\n', firstFile);
        end
        currentFiles = currentFiles(fileNum, :);
    else
        if (dummy_scans > 0)
            if (dummy_scans >= size(currentFiles, 1))
                error('Error:DummyScans', 'Please check dummy_scans variable (%d) as it exceeds or equals the no. of time points (%d) for\nfile %s\n',  ...
                    dummy_scans, size(currentFiles, 1), firstFile);
            end
            fileNum = (1:size(currentFiles, 1));
            fileNum(1:dummy_scans) = [];
            currentFiles = currentFiles(fileNum, :);
        end
    end
    
    tmpMask = mask;
    
    %% Load data
    data = icatb_read_data(currentFiles, [], tmpMask);
    
    tmpM = any(diff(data, 1, 2), 2);
    
    data = data(tmpM, :);
    
    tmpMask = tmpMask(tmpM);
    
    if (~strcmpi(preproc_type, 'remove mean per timepoint'))
        data = icatb_preproc_data(data, preproc_type);
    end
    
    disp('Backreconstructing subject components ...');
    
    %% Backreconstruct using STR
    [tc, ic] = icatb_dual_regress(data, group_icasig(tmpM, :));
    
    clear data;
    
    %% Files naming
    if  (~strcmpi(modalityType, 'smri'))
        component_name = [outputPrefix, SUBJECT_ICA_AN3_FILE, icatb_returnFileIndex(row), COMPONENT_NAMING, SESSION_POSTFIX, '1_'];
    else
        component_name = [outputPrefix, SUBJECT_ICA_AN3_FILE, COMPONENT_NAMING];
    end
    timecourse_name = strrep(component_name, COMPONENT_NAMING, TIMECOURSE_NAMING);
    
    tmpOutDir = outputDir;
    outFile = fullfile(tmpOutDir, [component_name, outputFormat]);
    
    if (~exist(tmpOutDir, 'dir'))
        mkdir(tmpOutDir);
    end
    
    disp(['Saving file ', outFile, ' ...']);
    
    if (~strcmpi(outputFormat, '.mat'))
        
        %% Write timecourse
        V.dim(1:3) = [size(tc), 1];
        V.n(1) = 1;
        V.fname = fullfile(tmpOutDir, [timecourse_name, outputFormat]);
        icatb_write_vol(V, tc);
        
        %% Write spatial maps
        V.dim(1:3) = comp_dims;
        V.fname = outFile;
        
        for nComp = 1:numOfIC
            if strcmpi(outputFormat, '.nii')
                V.n(1) = nComp;
            else
                V.n(1) = 1;
                V.fname = fullfile(tmpOutDir, [component_name, icatb_returnFileIndex(nComp), outputFormat]);
            end
            tmpDat = zeros(comp_dims);
            tmpDat(tmpMask) = ic(nComp, :);
            icatb_write_vol(V, tmpDat);
            clear tmpDat;
        end
        
    else
        
        %% MAT format
        icatb_save(outFile, 'tc', 'ic');
        
    end
    
    clear tc ic;
    
    fprintf('\n');
    
end
%% End of loop over subjects

fprintf('\nDone\n\n');