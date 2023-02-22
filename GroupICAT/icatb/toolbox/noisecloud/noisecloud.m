function [features_norm, feature_labels] = noisecloud(TR_value, network_paths, timeseries_paths, varargin)
%% Noise cloud feature extraction. Adapted from original noise cloud function
%at https://github.com/vsoch/noisecloud
%
% Inputs:
% 1. TR_value - Experimental TR in seconds
% 2. network_paths - File names of spatial maps in a cell array
% 3. timeseries_paths - File names of timecourses in a cell array
% 4. varargin - Variable number of arguments passed
%   a. convert_to_z - Convert to z-scores. Options are 'no' and 'yes'
%   b. outDir - Output directory to place results and temporary
%   coregistered images
%   c. coregister - Options are 1 and 0. Set this option to 0 only if you coregistered the template
%   images to functional image space. Templates are stored in directory
%   mr/raw and tweak the appropriate file names in noisecloud_setup.
%

% -------------------------------------------------------------------------
% Add path for noisecloud
fullpath = fileparts(which('noisecloud'));
addpath(genpath(fullpath));
global TR;
global NOISE_ATLAS_DIRS;

%convert_to_z = 'no';
%outDir = pwd;
coregister_im = 1;
for nV = 1:length(varargin)
    if (strcmpi(varargin{nV}, 'convert_to_z'))
        convert_to_z = varargin{nV + 1};
    elseif (strcmpi(varargin{nV}, 'outDir'))
        outDir = varargin{nV + 1};
    elseif (strcmpi(varargin{nV}, 'coregister'))
        coregister_im = varargin{nV + 1};
    end
end

% Check for SPM
if isempty(which('spm'))
    error('Please download SPM and add it to your path!');
end

% Ask the user to input images TR
%TR = input('Enter TR value:');
if (~exist('TR_value', 'var'))
    answer_tr = inputdlg({'Enter TR value in seconds: '}, 'Enter TR', 1, {'1'});
    if (isempty(answer_tr))
        error('TR is not selected');
    end
    TR_value = str2num(answer_tr{1});
end

TR = TR_value;

drawnow;

if (~exist('network_paths', 'var'))
    network_paths = selectFiles('Select component maps ...', '*sub*comp*nii', 'on');
end

if (~exist('timeseries_paths', 'var'))
    timeseries_paths = selectFiles('Select timecourses ...', '*sub*time*nii', 'off');
end

if (~exist('convert_to_z', 'var'))
    convert_to_z = questdlg('Do You Want To Convert Components To Z-scores?', 'Z-scores?', 'Yes', 'No', 'Yes');
end


if (~exist('outDir', 'var'))
    outDir = selectDir('Select output directory ...');
end

addpath (outDir);

if (~isempty(network_paths))
    VV = spm_vol(char(network_paths));
    file_names = cell(length(VV), 1);
    for nV = 1:length(VV)
        file_names{nV} = [VV(nV).fname, ',', num2str(VV(nV).n(1))];
    end
    network_paths = file_names;
    network_row_names = network_paths;
end

numComp2 = [];
if (~isempty(timeseries_paths))
    TC = readTimecourses(timeseries_paths);
    numComp2 = size(TC, 2);
end

numComp1 = [];
if (~isempty(network_paths))
    
    numComp1 = length(network_paths);
    
    if (coregister_im)
        
        NOISE_ATLAS_DIRS = outDir;
        
        sourceDir = fullfile(fullpath, 'mr', 'raw');
        atlas_source_images = {'aal2mni152.nii', 'grey.nii', 'white.nii', 'csf.nii', 'MNI152_T1_2mm_edges.nii', 'MNI152_T1_2mm_strucseg.nii', ...
            'MNI152_T1_2mm_eye_mask.nii', 'MNI152_T1_2mm_skull.nii'};
        for nF = 1:length(atlas_source_images)
            disp(['Coregistering ', atlas_source_images{nF}, ' to functional image ...']);
            noisecloud_spm_coregister(network_paths{1}, fullfile(sourceDir, atlas_source_images{nF}), '', outDir);
        end
        
   end
    
end

numComp = [numComp1, numComp2];

if (~all(numComp == numComp(1)))
    error('Component numbers must match between features');
end

numComp = numComp(1);

% This will grab paths to atlases and tissue maps.  You must prepare these
% in advance to be registered to your data!  See INSTRUCTIONS.txt for details
noisecloud_setup;


%% Get feature labels
feature_labels = getFeatureLabels(network_paths, timeseries_paths);

features = zeros(numComp, length(feature_labels));

%% Step 2: Extract features
for s = 1:numComp
    
    sm = []; tc = [];
    if (~isempty(network_paths))
        current_image = network_paths{s};
        current_image = spm_read_vols(spm_vol(current_image));
        current_image(isfinite(current_image) == 0) = 0;
        sm = nc_spatial_features(current_image, convert_to_z);
    end
    
    if (~isempty(timeseries_paths))
        tc = nc_temporal_features(TC(:, s), convert_to_z);
    end
    
    tmp = [sm, tc];
    
    features(s, :) = tmp;
    
end


%% Step 3: Normalization

features_norm = zeros(size(features));
for i = 1:size(features_norm,2)
    single_feature = features(:,i);
    features_norm(:,i) = (single_feature - mean(single_feature)) / std(single_feature);
end

% Replace NaN with zero
features_norm(isnan(features_norm)) = 0;

fprintf('%s\n','Feature extraction complete.');
fprintf('%s%s%s\n','     ',num2str(length(feature_labels)),' total features');
fprintf('%s%s%s\n','     ',num2str(length(network_row_names)),' total networks');

%end

function network_paths = selectFiles(titleStr, filterPattern, multiSelect)
%% Select files
%

[network_paths, pathname] = uigetfile(filterPattern, titleStr, 'MultiSelect', multiSelect);
if (isnumeric(network_paths))
    if (~network_paths)
        error('Files are not selected');
    end
end
network_paths = char(network_paths);
if (~strcmp(pathname(end), filesep))
    pathname = [pathname, filesep];
end
network_paths = [repmat(pathname, size(network_paths, 1), 1), network_paths];

%network_paths = cellstr(network_paths);

function outDir = selectDir(titleStr)
%% Select directory
%

outDir = uigetdir(pwd, titleStr);

if (isempty(outDir))
    error('Output directory is not selected');
end

if (~outDir)
    error('Output directory is not selected');
end

function data = readTimecourses(time_text)
%% Read timecourses
%
[dd, fN, extn] = fileparts(deblank(time_text(1, :)));

commaPos = find(extn == ',');
if (~isempty(commaPos))
    extn = extn(1:commaPos(1)-1);
end

if (~strcmpi(extn, '.nii') && ~strcmpi(extn, '.img'))
    
    data = [];
    
    for nT = 1:size(time_text, 1)
        
        temp = load(deblank(time_text(nT, :)));
        
        if isstruct(temp)
            % MAT file
            fieldNames = fieldnames(temp);
            dat = getfield(temp, fieldNames{1});
        else
            % Ascii file
            dat = temp;
        end
        
        data = [data,dat];
    end
    
    
else
    data = [];
    
    for nT = 1:size(time_text, 1)
        dat = spm_read_vols(spm_vol(deblank(time_text(nT, :))));
        dat(isfinite(dat) == 0) = 0;
        data = [data,dat];
    end
end


function feature_labels = getFeatureLabels(network_paths, timeseries_paths)

load('noise_cloud_params.mat');

feature_labels = [];
idx = 1;

if (~isempty(network_paths))
    
    % % Write temporary scripts for spatial and temporal features
    for s=1:size(spatial,2)
        nfeatures = str2num(spatial{s}.n);
        labels = regexp(spatial{s}.label,',','split');
        if nfeatures ~= size(labels,2) % If there is one label to describe features
            for f=1:nfeatures
                feature_labels{idx} = [ spatial{s}.label '.' num2str(f) ];
                idx = idx +1;
            end
        else % If there is a unique label per feature
            for f=1:nfeatures
                feature_labels{idx} = labels{f};
                idx = idx +1;
            end
        end
    end
end


if (~isempty(timeseries_paths))
    for t=1:size(temporal,2)
        nfeatures = str2num(temporal{t}.n);
        labels = regexp(temporal{t}.label,',','split');
        if nfeatures ~= size(labels,2) % If there is one label to describe features
            for f=1:nfeatures
                feature_labels{idx} = [ temporal{t}.label '.' num2str(f) ];
                idx = idx + 1;
            end
        else % If there is a unique label per feature
            for f=1:nfeatures
                feature_labels{idx} = labels{f};
                idx = idx + 1;
            end
        end
    end
    
end
