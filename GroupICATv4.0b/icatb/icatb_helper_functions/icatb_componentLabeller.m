function icatb_componentLabeller(files, templateFiles, labelsFile, dispDefaults)
%% Label the components using the templates information
%
% Inputs:
% 1. files - Image files
% 2. templateFiles - Template files in a 4D nifti format. Please see
% icatb/icatb_templates/RSN.zip for an example
% 3. labelsFile - Template labels in a text file. Please see
% icatb/icatb_templates/RSN.txt
% 4. dispDefaults - Display defaults in a structure. Field names are
% 'structFile', 'threshold', 'image_values' and 'convert_to_zscores'
%
%

icatb_defaults;

global IMAGE_VALUES;
global CONVERT_Z;
global THRESHOLD_VALUE;

if (icatb_get_matlab_version < 14)
    error('Component labeller requires R14 and higher');
end

oldDir=pwd;

%% Component files
if (~exist('files', 'var') || isempty(files))
    files = icatb_selectEntry('typeSelection', 'multiple', 'filter', '*.img; *.nii', 'fileType', 'image', 'title', 'Select image files ...');
end

if (isempty(files))
    error('Image files are not selected');
end

files = deblank(files);

tmpDir = fileparts(deblank(files(1, :)));
if (isempty(tmpDir))
    tmpDir = pwd;
end

startPath = fullfile(fileparts(which('gift.m')), 'icatb_templates');

%% Template files
if (~exist('templateFiles', 'var') || isempty(templateFiles))
    templateFiles = icatb_selectEntry('typeSelection', 'multiple', 'filter', '*.img; *.nii; *.zip', 'title', 'Select template files ...', 'startPath', ...
        startPath);
end

if (isempty(templateFiles))
    error('Template files are not selected');
end

templateFiles = deblank(templateFiles);

chkZip = find(icatb_good_cells(regexp(cellstr(templateFiles), '\.zip$')) ~= 0);

if (~isempty(chkZip))
    if (size(templateFiles, 1) > 1)
        error('Select atmost one zip file when selecting templates');
    end
end

isZip = strcmpi(templateFiles(1, end-2:end), 'zip');
if (isZip)
    tmpOutDir = fullfile(tmpDir, 'tmp_gift_rsn_templates');
    if (~exist(tmpOutDir, 'dir'))
        mkdir(tmpOutDir);
    end
    removeFiles(tmpOutDir);
    templateFiles = deblank(templateFiles(1, :));
    templateFiles = listFiles(templateFiles, tmpOutDir);
end

files = icatb_rename_4d_file(files);
templateFiles = icatb_rename_4d_file(templateFiles);

promptString = 'Select component labels text file';

if (~exist('labelsFile', 'var') || isempty(labelsFile))
    labelsFile = icatb_selectEntry('typeEntity', 'file', 'typeselection', 'single', ...
        'title', promptString, 'filter', '*.txt', 'startPath', startPath);
end

if (isempty(labelsFile))
    error('Component labels file is not selected');
end

cd(oldDir);

%if (~exist('thresh', 'var') || isempty(thresh))
thresh = THRESHOLD_VALUE;
%end

if (ischar(thresh))
    thresh = str2num(thresh);
end

DISP_DEFS = struct('structFile', fullfile(fileparts(which('gift.m')), 'icatb_templates', 'nsingle_subj_T1_2_2_5.nii'), 'threshold', thresh, 'image_values', IMAGE_VALUES, ...
    'convert_to_zscores', CONVERT_Z);
if (~exist('dispDefaults', 'var'))
    dispDefaults = DISP_DEFS;
end

fnames = fieldnames(DISP_DEFS);

for i = 1:length(fnames)
    if (~isfield(dispDefaults, fnames{i}))
        dispDefaults.(fnames{i}) = DISP_DEFS.(fnames{i});
    end
end

runLabeller(files, templateFiles, labelsFile, dispDefaults);

function runLabeller(files, templateFiles, compLabelsTxtFile, dispDefaults)
%% Run labeller
%

%% Get component labels
comp_labels = getCompLabels(compLabelsTxtFile);
if (length(unique([comp_labels.value])) ~= size(templateFiles, 1))
    error(['No. of labels specified in file ', compLabelsTxtFile, ' doesn''t match the no. of spatial templates specified']);
end

%% Compute r-value
allLabels = [];
for nL = 1:length(comp_labels)
    allLabels = [allLabels, repmat({comp_labels(nL).name}, 1, length(comp_labels(nL).value))];
end

vals = [comp_labels.value];

thresh = eps; %dispDefaults.threshold;
labelsR = repmat({''}, 1, size(files, 1));

fileContents = repmat({''}, size(files, 1) + 2, 3);

fileContents(1, :) = {'Name', 'Label', 'Max correlation w.r.t template'};

for n = 1:size(files, 1)
    
    disp(['Determining label for file ', deblank(files(n, :)), ' ...']);
    
    tmp = zeros(1, length(comp_labels));
    cF = deblank(files(n, :));
    [pathstr, fN, extn] = fileparts(cF);
    structVol = icatb_spm_vol(cF);
    data = icatb_spm_read_vols(structVol);
    data(isfinite(data) == 0) = 0;
    maskM = abs(data(:)) > thresh;
    for nL = 1:size(templateFiles, 1)
        ic = icatb_resizeData(structVol, deblank(templateFiles(nL, :)), 1);
        if (max(abs(ic(maskM))) > eps)
            tmp(nL) = icatb_corr2(data(maskM), ic(maskM));
        end
    end
    [R, inds] = max(abs(tmp));
    name = sprintf(repmat('%s, ', 1, length(inds)), allLabels{vals == vals(inds)});
    rval = num2str(round(tmp(inds)*10000)/10000);
    labelsR{n} = [name(1:end-2), ', R = ', rval];
    fileContents(2 + n, :) = {[fN, extn], name(1:end-2), rval};
    
end

[p, fN] = fileparts(fN);
labels_file_name = fullfile(pathstr, [fN, '_labels.txt']);
icatb_print_table(fileContents, labels_file_name);
disp(['Component labels are stored in file ', labels_file_name]);
disp('');
clear data ic;


disp('Displaying files ...');

%% Display
graphicsH = [];
for nComp = 1:size(files, 1)
    file = deblank(files(nComp, :));
    [d, fN] = fileparts(file);
    H = icatb_orth_views(file, 'structfile', dispDefaults.structFile, 'threshold', dispDefaults.threshold, 'convert_to_zscores', dispDefaults.convert_to_zscores, ...
        'image_values', dispDefaults.image_values, 'labels', [fN, ' (', labelsR{nComp}, ')'], 'fig_title', ['Component ', icatb_returnFileIndex(nComp)], 'set_to_max_voxel', 1);
    graphicsH(length(graphicsH) + 1).H = H;
end

icatb_plotNextPreviousExitButtons(graphicsH);

fprintf('Done\n');

function compLabels = getCompLabels(txtFile)

fid = fopen(txtFile, 'r');
if (fid == -1)
    error(['File ', txtFile, ' cannot be opened for reading']);
end
try
    dd = textscan(fid, '%s', 'delimiter', '\t\n,', 'multipleDelimsAsOne', 1);
    val = dd{1};
catch
    val = [];
end
fclose(fid);
val = val(icatb_good_cells(val));
chk = cellfun('isempty', regexp(val, '\d+'));

inds = find(chk == 1);

compLabels = repmat(struct('name', '', 'value', []), 1, length(inds));
for nI = 1:length(inds)
    compLabels(nI).name = val{inds(nI)};
    if (nI == length(inds))
        endT = length(val);
    else
        endT = inds(nI + 1) - 1;
    end
    dd = str2num(char(val{inds(nI) + 1:endT}));
    compLabels(nI).value = dd(:)';
end

function files = listFiles(zipFile, outDir)

icatb_unzip(zipFile, outDir);

files = icatb_listFiles_inDir(outDir, '*.img;*.nii');

if (isempty(files))
    error(['No image files found in ', zipFile]);
end

files = icatb_fullFile('directory', outDir, 'files', files);

function removeFiles(outDir)

delete(fullfile(outDir, '*.*'));


