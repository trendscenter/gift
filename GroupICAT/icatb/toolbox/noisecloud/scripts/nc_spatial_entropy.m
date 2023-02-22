function feature = nc_spatial_entropy(spatialMap);
      % FEATURE 3/119: Spatial Entropy, measure of information content of a;
     % spatial distribution.  Narrow distribution = less spatial entropy;
     ICnorm = (spatialMap(:) - mean(spatialMap(:))) / std(spatialMap(:));
     ICentropy = nc_entropy(ICnorm);
     feature = abs(log(ICentropy));
      end