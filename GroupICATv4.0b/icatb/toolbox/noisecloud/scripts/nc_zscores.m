function z = nc_zscores(a)
%% Z-scores
%

a = detrend(a, 0);
z = a*diag(1./std(a));