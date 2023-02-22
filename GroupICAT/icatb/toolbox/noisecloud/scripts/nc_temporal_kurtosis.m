function feature = nc_temporal_kurtosis(T);
      % Feature 17: Kurtosis of timeseries;
     % Kurtosis Kurtosis is a measure of how outlier-prone a distribution is;
      % The kurtosis of the normal distribution is 3. Distributions that are more;
      % outlier-prone than the normal distribution have kurtosis greater than 3;
      % distributions that are less outlier-prone have kurtosis less than 3;
     feature = nc_kurtosis(T);
  end