function feature = nc_temporal_skewness(T);
      % Feature 16: skewness of timeseries;
     % First calculate skewness. Skewness is a measure of the asymmetry of ;
     % the data around the sample mean. If skewness is negative, the data are ;
     % spread out more to the left of the mean than to the right. If skewness is;
     % positive, the data are spread out more to the right. The skewness of the;
      % normal distribution (or any perfectly symmetric distribution) is zero.;
     feature = nc_skewness(T);
  end