function feature = nc_spatial_skewness(spatialMap);
      % FEATURE 2/118: Skewness: remove mean of IC, normalize variance, cube,;
     % divide by number of voxels.  This returns same as doing;
     % skewness(curr_network(:));
     ICnorm = (spatialMap(:) - mean(spatialMap(:))) / std(spatialMap(:));
     ICskew = sum(power(ICnorm(:),3)) / length(ICnorm(:));
     feature = abs(log(ICskew));
      end