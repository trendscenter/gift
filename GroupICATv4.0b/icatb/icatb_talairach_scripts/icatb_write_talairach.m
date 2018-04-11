function  icatb_write_talairach(data, dist, thresh, xdim, ydim, zdim, xlen, ylen, zlen, origin, out3coord, out3coordZ)
%function [] = write_talairach(data,dist,thresh,xdim,ydim,zdim, ...
%   xlen,ylen,zlen,origin,out3coord,out3coordZ)
%data = can be binary or valued...results will be sorted by this value
%dist = distance between contiguous voxels
%thresh = threshold
%xdim,ydim,zdim....'nuf said
%xlen,ylen,zlen...(voxel size in mm)
%origin...e.g. [27 38 11]
%out3coord = 'Driving_Red.txt'
%out3coordZ = 'Driving_Red_Z.txt'
%Author: Vince Calhoun, 1 June 2001

tmp = reshape(data,xdim,ydim,zdim);
%thresh = 2;
%dist = 3;%mm
%b = (derivative3d(tmp.*(tmp>thresh)));
%ind = find((tmp>thresh)&(abs(b)<.005));
ind = find((tmp>thresh));
[Y I] = sort(tmp(ind));
ind = ind(flipud(I));
[x y z] = ind2sub([xdim ydim zdim],ind);

%xlen = 3;ylen = 3;zlen = 5;
xlen2 = xlen^2;ylen2 = ylen^2;zlen2 = zlen^2;
clear x2 y2 z2;
x2(1) = x(1);y2(1) = y(1);z2(1) = z(1);cnt = 1;
for j = 1:length(x),
   flag = 1;
   for k = 1:length(x2),
      if (sqrt(xlen2*(x(j)-x2(k))^2+ylen2*(y(j)-y2(k))^2+zlen2*(z(j)-z2(k))^2) < dist),
         flag = 0;%don't count...too close
      end;
   end;
   if (flag),
      cnt = cnt+1;
      x2(cnt) = x(j);y2(cnt) = y(j);z2(cnt) = z(j);
   end;
end;
%[x2' y2' z2']

clear x2t y2t z2t;
for j = 1:length(x2);
   x2t(j) = (x2(j)-origin(1))*xlen;
   y2t(j) = (y2(j)-origin(2))*ylen;
   z2t(j) = (z2(j)-origin(3))*zlen;
   Zval(j) = tmp(x2(j),y2(j),z2(j));
end;
a = round(icatb_mni2tal([x2t' y2t' z2t']));

fid1 = fopen(out3coord,'wt');
fid2 = fopen(out3coordZ,'wt');
for j = 1:size(a,1),
   fprintf(fid1,'%i %i %i\n', [a(j,1);a(j,2);a(j,3)]);
   fprintf(fid2,'%i %i %i %d MNI: %i %i %i\n', ...
      [a(j,1);a(j,2);a(j,3);Zval(j);x2t(j);y2t(j);z2t(j)]);
end;
fclose(fid1);
fclose(fid2);

