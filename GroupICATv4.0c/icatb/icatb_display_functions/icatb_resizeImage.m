function [images, coords, HInfo, slices, text_left_right] = icatb_resizeImage(structVol, compFile, anatomicalView, slices, file_numbers, spmPath)
% function that takes structural file and component files as input and
% returns the interpolated component images
%
% Input:
% 1. structVol - structural volume
% 2. compFile - fullfile path for the component files
% 3. anatomicalView - slice plane
% 4. slices - slices in mm
% 5. file_numbers - components needed
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



icatb_defaults;
global FLIP_ANALYZE_IMAGES;


% check the existence of the vars
if ~exist('structVol', 'var')
    error('Structural volume is missing');
end

if ~exist('compFile', 'var')
    error('component file is missing');
end

if ~exist('file_numbers', 'var')
    file_numbers = [];
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

compFile = icatb_rename_4d_file(compFile);

if (isempty(file_numbers))
    file_numbers = (1:size(compFile, 1));
end

compFile = compFile(file_numbers, :);

% get the number of components by reading the dimensions (handles nifti
% also)
numComp = length(file_numbers);

% voxel size
VOX = double(structVol(1).private.hdr.pixdim(2:4));

% make the voxel sizes positive
VOX = abs(VOX);

% slice definition
[parameters] = icatb_get_slice_def(structVol, anatomicalView);

% get the required parameters
transform = parameters.transform; slicedef = parameters.slicedef;

if isempty(slices)
    slices = parameters.slices;
end

% permuted order
permuteOrder = parameters.permuteOrder;

% get the extension
[path_comp, file_comp, extn_comp] = fileparts(deblank(compFile(1, :)));
extn_comp = icatb_parseExtn(extn_comp);

% Initialise img structure
img = repmat(struct('V', []), numComp + 1, 1);

% check the image extension
if (strcmpi(extn_comp, '.img') || strcmpi(extn_comp, '.nii'))
    % analyze convention
    % loop over image vols
    for ii = 1:length(img)
        if ii == 1
            img(ii).V = structVol;
        else
            % get the volume from spm vols
            img(ii).V = icatb_spm_vol(deblank(compFile(ii - 1, :)));
            
            if ii == 2
                text_left_right = returnLeftRightText(img(ii).V.mat(1, 1));
            end
            
        end
    end
    % end loop over image vols
    
    % elseif strcmpi(extn_comp, '.nii')
    %     % nifti extension
    %     % loop over image vols
    %     for ii = 1:length(img)
    %         if ii == 1
    %             img(ii).V = structVol;
    %         else
    %             % get the volume from spm vol nifti
    %             img(ii).V = icatb_get_vol_nifti([compFile, ',', num2str(file_numbers(ii-1))]);
    %
    %             if ii == 2
    %                 text_left_right = returnLeftRightText(img(ii).V.mat(1, 1));
    %             end
    %
    %         end
    %     end
    %     % end loop over image vols
else
    % print error
    error('Unknown image extensions');
end
% end for checking image extensions


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
        i1 = icatb_spm_sample_vol(img(jj).V(1), vixyz(X, :), vixyz(Y, :), vixyz(Z, :), [1, nan]);
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



function text_left_right = returnLeftRightText(m)
% Return left right text

if  sign(m(1, 1)) == -1
    text_left_right = 'RL';
else
    text_left_right = 'LR';
end
