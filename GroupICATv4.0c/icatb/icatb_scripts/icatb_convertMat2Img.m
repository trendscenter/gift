function icatb_convertMat2Img(matFile, templateFile, pixdim)
% Converts .mat to image format and places the images in the directory
% where the .mat file resides.
%
% Inputs:
% 1. matFile - Full file path of the MAT file
% 2. templateFile - Full file path of the template image file that the
% resulting images will use.
% 3. pixdim - Pixel dimensions (Optional variable)


if ~exist('templateFile', 'var') | isempty(templateFile)
    % get the volume information from template file
    templateFile = which('nsingle_subj_T1_2_2_5.nii');
end


% get the file naming
[pathstr, file_name, extn] = fileparts(deblank(matFile));

if ~strcmpi(extn, '.mat')
    error('File should be in .mat format');
end

% load mat file
S = load(matFile);

% get the field names
fieldNames = fieldnames(S);

% check the variable field names
if length(fieldNames) > 1 | isempty(fieldNames)
    error('Only one variable should be there in the .mat file');
end

% retrieve the cdata
cdata = getfield(S, fieldNames{1});

% check the data is an array
if ~isnumeric(cdata)
    error('data should be in numeric format');
end

dataType = 'real';

size_cdata = size(cdata);

% check the size of the data
if length(size(cdata)) == 3
    size_cdata = [size(cdata, 1), size(cdata, 2), size(cdata, 3), 1];
elseif length(size(cdata)) == 2
    size_cdata = [size(cdata, 1), size(cdata, 2), 1, 1];
elseif length(size(cdata)) == 1
    size_cdata = [size(cdata, 1), 1, 1, 1];
elseif length(size(cdata)) > 4
    error('data dimensions cannot be greater than 4');
end

% volume information
V = icatb_get_vol_nifti(templateFile);

V = V(1);

% dimension
V.dim(1:3) = [size_cdata(1), size_cdata(2), size_cdata(3)];

dim = V.dim(1:3);

% check the existence of variable
if exist('pixdim', 'var')

    % origin
    origin = (1 + double(dim)) / 2;

    % calculate offset
    off = -pixdim.*origin;

    % create a transformation matrix
    V.mat = [pixdim(1) 0 0 off(1); 0 pixdim(2) 0 off(2); 0 0 pixdim(3) off(3); ...
        0 0 0 1];

end
% end for calculating transformation matrix

imageType = 'real';

msgString = ['Writing out a new set of images ... '];
disp(msgString);
numFiles = size_cdata(4);
% loop over time points
for ii  = 1:numFiles
    % show the status of writing images
    if ii == 1
        statusH = icatb_statusBar('init', numFiles , 'Writing Images', '', '');
        statusH = icatb_statusBar('increment', 1);
    elseif ii > 1 & ii < numFiles
        statusH = icatb_statusBar('increment', 1);
    end
    % print string
    %disp(['Writing image ', num2str(ii)]);
    % fileNaming
    if ii < 10
        fileNaming = ['_00', num2str(ii)];
    elseif ii >= 10 & ii < 100
        fileNaming = ['_0', num2str(ii)];
    else
        fileNaming = ['_', num2str(ii)];
    end
    file = file_name;
    % file name
    V.fname = fullfile(pathstr, [file, fileNaming, '.img']);
    % real data
    icatb_write_vol(V, squeeze(cdata(:, :, :, ii)));
end

statusH = icatb_statusBar('increment', 1);
% finish status bar
icatb_statusBar('finished');
msgString = ['The new set of images are stored in the directory ', pathstr];
disp(msgString);
