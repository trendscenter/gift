function realWorldCoord = icatb_voxel_to_real(V, voxelCoord)
% Voxel to real world coordinates

tmp = V(1).mat;
realWorldCoord = tmp(1:3,:)*[voxelCoord(:); 1];

realWorldCoord = realWorldCoord(:)';