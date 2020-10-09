function feature = nc_spatial_nodedist(spatialMap);
          % Features include:;
     % avg_distance_btw_10_highest_connected_nodes,average_Z_score_10_highest_connected_nodes;
      % FEATURE 4/120: Average Distance between 10 nodes with highest;
     % connectivity (10 local max).  I would think that "small world networks";
     % would have a smaller average distance.  First find local max;
     [Maxima,MaxPos,~,~] = MinimaMaxima3D(spatialMap,1,0,10,0);
     % Now calculate pairwise euclidian distances, sum, and divide by number;
     % of unique distances;
     pairwise_dist = nc_dist(MaxPos');
     avg_dist = (sum(sum(pairwise_dist))/2) / ((length(pairwise_dist(:)) - size(MaxPos,1)) / 2);
     feature(1) = max(avg_dist,0);
              % FEATURE 5/121: Average Z Score 10 nodes highest connectivity;
     feature(2) = max(mean(Maxima),0);
         end