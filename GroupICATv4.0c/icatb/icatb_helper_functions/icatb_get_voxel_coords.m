function xyz = icatb_get_voxel_coords(dims)
%% Get voxel coords
%

xdim = dims(1);
ydim = dims(2);
zdim = dims(3);

[xords, yords] = ndgrid(1:xdim, 1:ydim);
xords = xords(:)'; yords = yords(:)';
CrPl    = 1:zdim;
zords   = CrPl(:)*ones(1,xdim*ydim);
xyz   = [repmat(xords,1,numel(CrPl)); ...
    repmat(yords,1,numel(CrPl)); ...
    reshape(zords',1,[])];