function [b,nn,mm,DIM] = icatb_get_montage(a,voxelSize,gridsize)
% [b,nn,mm,DIM] = icatb_get_montage(a,gridsize)
%--------------------------------------------------------------------------
% CREATED BY: Eric Egolf 
% LAST MODIFIED: 12-29-03
% ABOUT: Makes montage of 3D image
%
% #########################################################################
% 
% USAGE:
%
% [b] = icatb_get_montage(a)
%   INFO: makes montage of 3D image
%   PARAMETERS:
%       a = 4D (x,y,1,z) image
%
%   OUTPUT:
%   b = 2D image that is a montage of 3D image
%       ex: if if had a 3D image with 5 slices in the z direction
%       this function would return a 2d image with a montage of each of the
%       slices
% #############################################################
% 
% LICENSING:
% 
% 
%------------------------------------------------------
if(~exist('voxelSize','var'))
    voxelSize = [1 1 1];    
end
if (~exist('gridsize','var')),
   gridsize = -999;
end;

siz = [size(a,1) size(a,2) size(a,4)];
nn = sqrt(prod(siz))/(siz(2));
mm = (siz(3))/nn;
if(~all(voxelSize==[1 1 1]))
    siz2=[siz(1)*voxelSize(2) siz(2)*voxelSize(1) siz(3)];
    nn=sqrt(prod(siz2))/(siz2(2));
    mm=(siz2(3))/nn;
end
if (ceil(nn)-nn) < (ceil(mm)-mm),
   nn = ceil(nn); mm = ceil(siz(3)/nn);
else
   mm = ceil(mm); nn = ceil(siz(3)/mm);
end

if (gridsize ~= -999),
   nn = gridsize(1);
   mm = gridsize(2);
end;

b = a(1,1); % to inherit type 
% b(1,1) = 0; % from a

b = repmat(b, [mm*siz(1), nn*siz(2), size(a,3), 1]);

rows = 1:siz(1); cols = 1:siz(2);
for i=0:mm-1,
   for j=0:nn-1,
      k = j+i*nn+1;
      if k<=siz(3),
         b(rows+i*siz(1),cols+j*siz(2),:) = a(:,:,:,k);
      end
   end
end
DIM(1)=size(a,1);
DIM(2)=size(a,2);