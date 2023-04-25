function [images, coords, HInfo, slices, text_left_right] = icatb_MontageImMemResize(struImBgTemplate, struImComponent, anatomicalView, slices, spmPath)
% function that takes structural file and component files as input and
% returns the interpolated component images
% This function belongs to icatb_MontageImMem.m
%
% Input:
% 1. struImBgTemplate - structural volume
% 2. struImComponent - fullfile path for the component files
% 3. anatomicalView - slice plane
% 4. slices - slices in mm
%
% Output:
% 1. images - interpolated image
% 2. coords - Real world coordinates
% 3. HInfo - header info
% 4. slices - slices in mm
% 5. Text for left-right - size equal to the number of components by 2


% Note: image dimension will be the same as structural

%slicedef - slice def other than selected

%slices - selected slice plane

%loop over slices

%loop over images

% %t the current image
%
% get the transform (anatomical view) and multiply with the transformation matrix
% for each slice
%
% evaluate spm sample vol

% Initiations
icatb_defaults;
global FLIP_ANALYZE_IMAGES;
global INTERP_VAL;

interp_val = INTERP_VAL;
if (isempty(interp_val))
    interp_val = 0;
end

% check the existence of the vars
if ~exist('struImBgTemplate', 'var')
    error('Structural volume is missing');
end

if ~exist('struImComponent', 'var')
    error('component file is missing');
end

% By default anatomical view is axial
if ~exist('anatomicalView', 'var')
    anatomicalView = 'axial';
else
    anatomicalView = lower(anatomicalView);
end

if ~exist('slices', 'var')
    slices = [];
end

% only 1 component needed for this plot 042523
numComp = 1;
% Initialise img structure
img = repmat(struct('V', []), numComp + 1, 1);
img(1).V=struImBgTemplate;
img(2).V=struImComponent;

% voxel size
VOX = double(struImBgTemplate(1).private.hdr.pixdim(2:4));

% make the voxel sizes positive
VOX = abs(VOX);

% slice definition
[parameters] = icatb_get_slice_def(struImBgTemplate, anatomicalView);

% get the required parameters
transform = parameters.transform; slicedef = parameters.slicedef;

if isempty(slices)
    slices = parameters.slices;
end

% permuted order
permuteOrder = parameters.permuteOrder;

% form grid
% get coordinates for plane
X = 1;Y = 2;Z = 3;
dims = slicedef;
xmm = dims(X,1):dims(X,2):dims(X,3);
ymm = dims(Y,1):dims(Y,2):dims(Y,3);
zmm = slices;

[y, x] = meshgrid(ymm, xmm');
vdims = [length(xmm), length(ymm), length(zmm)];

% no of slices, and panels (an extra for colorbars)
nslices = vdims(Z);
minnpanels = nslices + 2;

% cycle through slices displaying images
nvox = prod(vdims(1:2));

images = zeros(numComp + 1, vdims(1), vdims(2), vdims(3));
coords = zeros(vdims(1), vdims(2), vdims(3), 3);

% loop over slices
for ii = 1:nslices
    ixyzmm = [x(:)'; y(:)'; ones(1,nvox)*zmm(ii); ones(1, nvox)];
    % loop over images
    for jj = 1 : numComp + 1
        % to voxel space of image
        vixyz = (transform*img(jj).V(1).mat) \ ixyzmm;
        % raw data
        i1 = icatb_spm_sample_vol(img(jj).V(1), vixyz(X, :), vixyz(Y, :), vixyz(Z, :), [interp_val, nan]);
        i1 = reshape(i1, vdims(1:2));
        
        % assign the images
        images(jj, :, :, ii) = i1;
        
        % Convert to real world coordinates
        if jj == 1
            RCP = [vixyz(X, :); vixyz(Y, :); vixyz(Z, :)];
            RCP(4, :) = 1;
            tempCoords = (img(jj).V(1).mat(1:3, :)*RCP)';
            tempCoords = reshape(tempCoords, [vdims(1:2), 3]);
            coords(:, :, ii, :) = tempCoords;
            clear tempCoords;
            clear RCP;
        end
        % End for converting real world coordinates
    end
end

% permute the order
images = permute(images, permuteOrder);

% Real world coordinates
coords = permute(coords, [ceil(permuteOrder(2:end)-1), 4]);

% set the nan indices to zero
images((isfinite(images) == 0)) = 0;

%save Header info in structure
HInfo = struct(...
    'DIM', [size(images, 2) size(images, 3) size(images, 4) ], ...
    'V',    img(1).V,...
    'VOX', VOX);

% left right from read image
text_left_right = returnLeftRightText(img(2).V.mat(1, 1));



function text_left_right = returnLeftRightText(m)
% Return left right text

if  sign(m(1, 1)) == -1
    text_left_right = 'RL';
else
    text_left_right = 'LR';
end
