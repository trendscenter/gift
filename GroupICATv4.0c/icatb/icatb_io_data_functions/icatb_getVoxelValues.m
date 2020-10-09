function [voxelValues] = icatb_getVoxelValues(voxel, files, structFile, dataType, complexInfo, writeType, ...
    zipFileName, files_in_zip)
% get voxel values
%
% Input:
% 1. voxel - voxel location
% 2. files - analyze images
% 3. structFile - structural file
%
% Output:
% voxelValues - voxel values

% get spm files path

if ~exist('dataType', 'var')
    dataType = 'real';
end

if ~exist('complexInfo', 'var')
    complexInfo = [];
end

if ~exist('writeType', 'var')
    writeType = 'read';
end

if ~exist('zipFileName', 'var')
    zipFileName = {};
end

pathstr_comp = fileparts(deblank(files(1, :)));
% Unzip files
if ~isempty(zipFileName)
    icatb_unzip(zipFileName, pathstr_comp);
end

% get the file naming
allfiles = icatb_get_complex_files_naming(files, dataType, complexInfo, writeType);

clear files;


if isstruct(allfiles)
    % first set of files (real or magnitude)
    voxelValuesFirst = getVoxelValues(voxel, str2mat(allfiles.first), structFile);

    % first set of files (imaginary or phase)
    voxelValuesSecond = getVoxelValues(voxel, str2mat(allfiles.second), structFile);

    % get only the magnitude part
    if strcmpi(complexInfo.complexType, 'real&imaginary') & strcmpi(writeType, 'read')
        voxelValues = abs(complex(voxelValuesFirst, voxelValuesSecond));
    else
        voxelValues = voxelValuesFirst;
    end

else
    % get the voxel values
    voxelValues = getVoxelValues(voxel, allfiles, structFile);

end

% zip the images back
if ~isempty(zipFileName)
    fileN = str2mat(files_to_zip);
    %fileN = icatb_fullFile('directory', pathstr_comp, 'files', fileN);
    icatb_delete_file_pattern(fileN, pathstr_comp);
end


function [voxelValues] = getVoxelValues(voxel, files, structFile)

%numComp = icatb_get_countTimePoints(files);

% different time points
numComp = zeros(1, size(files, 1));

% loop over files
for ii = 1:size(files, 1)
    % get the count for the time points
    numComp(ii) = icatb_get_countTimePoints(deblank(files(ii, :)));
end
% end loop over files

% get the file numbers
%file_numbers = [1:1:numComp];

% get the structural volume
structVol = icatb_get_vol_nifti(structFile);

% slice definition
[parameters] = icatb_get_slice_def(structVol, 'axial');

% get the required parameters
transform = parameters.transform; slicedef = parameters.slicedef;

slices = parameters.slices;

% permuted order
permuteOrder = parameters.permuteOrder;

clear parameters;

% Initialise img structure
img = repmat(struct('V', []), sum(numComp), 1);

volCount = 0;
% loop over files
for ii = 1:size(files, 1)
    tempV = icatb_spm_vol(deblank(files(ii, :)));
    for nn = 1:length(tempV)
        volCount = volCount + 1;
        img(volCount).V = tempV(nn);
    end
    clear tempV;
end
% end loop over files

% form grid
% get coordinates for plane
X = 1;Y = 2;Z = 3;
dims = slicedef;
% get the xmm and ymm
xmm = dims(X,1):dims(X,2):dims(X,3);
ymm = dims(Y,1):dims(Y,2):dims(Y,3);

% get the slices for the corresponding voxel
xmm = xmm(voxel(1));
ymm = ymm(voxel(2));
% get the slice only for
zmm = slices(voxel(3));

[y, x] = meshgrid(ymm, xmm');
vdims = [length(xmm), length(ymm), length(zmm)];

% no of slices, and panels (an extra for colorbars)
nslices = vdims(Z);
minnpanels = nslices + 2;

% cycle through slices displaying images
nvox = prod(vdims(1:2));

images = zeros(sum(numComp), vdims(1), vdims(2));

ixyzmm = [x(:)'; y(:)'; ones(1,nvox)*zmm; ones(1, nvox)];
% loop over images
for jj = 1 : sum(numComp)
    % to voxel space of image
    vixyz = (transform*img(jj).V(1).mat) \ ixyzmm;
    % raw data
    i1 = icatb_spm_sample_vol(img(jj).V(1), vixyz(X, :), vixyz(Y, :), vixyz(Z, :), [1, nan]);
    % assign the images
    images(jj) = i1;
end

checkNan = isnan(images);
% set the nan indices to zero
images(checkNan) = 0;

% Initialise voxel values
voxelValues = zeros(1, sum(numComp));

% get the voxel values
for ii = 1:sum(numComp)
    voxelValues(ii) = images(ii);
end