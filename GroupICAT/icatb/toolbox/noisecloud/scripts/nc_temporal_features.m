function f = nc_temporal_features(T, convertToZ)
%% Compute temporal features
%

if (~exist('convertToZ', 'var'))
    convertToZ = 'no';
end

T = T(:);
if (strcmpi(convertToZ, 'yes'))
    T = T ./ std(T);
end


tmp = nc_temporal_ac(T);
f{1} = tmp(:)';
tmp = nc_temporal_bins(T);
f{2} = tmp(:)';
tmp = nc_temporal_dynrange(T);
f{3} = tmp(:)';
tmp = nc_temporal_entropy(T);
f{4} = tmp(:)';
tmp = nc_temporal_highfreqnoise(T);
f{5} = tmp(:)';
tmp = nc_temporal_kurtosis(T);
f{6} = tmp(:)';
f{7} = nc_temporal_mean(T);
tmp = nc_temporal_peaks(T);
f{8} = tmp(:)';
tmp = nc_temporal_psd(T);
f{9} = tmp(:)';
tmp = nc_temporal_skewness(T);
f{10} = tmp(:)';
tmp = nc_temporal_energyratio(T);
f{11} = tmp(:)';

f = [f{:}];

