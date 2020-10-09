function feature = nc_temporal_ac(T);
      % Features include;
     % one_lag_temporal_ac,two_lag_temporal_ac,three_lag_temporal_ac,four_lag_temporal_ac,five_lag_temporal_ac;
      % Feature 18-22: One,two,three,four,five lag auto-correlation;
     [ACF,dd1,dd2] = nc_autocorr(T,5);
     feature(1) = ACF(2);
     feature(2) = ACF(3);
     feature(3) = ACF(4);
     feature(4) = ACF(5);
     feature(5) = ACF(6);
  end