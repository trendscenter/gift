function k = nc_kurtosis(x)
% Kurtosis

x = detrend(x, 0);
s2 = mean(x.^2);
m4 = mean(x.^4);
k = (m4 ./ (s2.^2));