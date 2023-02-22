function feature = nc_temporal_entropy(T)
    % Feature 23: Temporal Entropy;
   feature = abs(log(nc_entropy(T)));
  end