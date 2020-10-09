function [V, HInfo] = icatb_returnHInfo(P)
% Return Volume and Header information of a vector of 3D analyze images,
% one 4D nifti or one 4D analyze image
%
% Input:
%
% 1. P - file names as a character array (Pass either a nifti, 4D analyze
% file or a vector of 3D analyze files)
%
% Output:
%
% 1. V - volume information
% 2. HInfo - header information


if size(P, 1) > 1
    P = deblank(P(1, :));
end

V = icatb_get_vol_nifti(P);

% voxel size
VOX = double(V(1).private.hdr.pixdim(2:4));
% make the voxel sizes positive
VOX = abs(VOX);

%save Header info in structure
HInfo = struct('DIM', V(1).dim(1:3), 'V', V, 'VOX', VOX);