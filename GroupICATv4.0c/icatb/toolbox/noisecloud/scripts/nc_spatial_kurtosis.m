function feature = nc_spatial_kurtosis(spatialMap);
  % FEATURE NAME: Kurtosis;
      % We take a component image, sum the voxel values quadrupled, and divide by the ;
 % number of voxels, subtract 3.  We then normalize the value by linear scaling transform of abs(ln(kurtosis));
 % This does the same thing as doing kurtosis(curr_network(:));
      ICkurt = sum(power(spatialMap(:),4)) / length(spatialMap(:));
 feature = abs(log(ICkurt));
      end