function voxelCoord = icatb_real_to_voxel(V, realWorldCoord)
% Real to voxel coordinates

tmp = V(1).mat;
is = inv(tmp);
voxelCoord = is(1:3,1:3)*realWorldCoord(:) + is(1:3,4);

voxelCoord = voxelCoord(:)';