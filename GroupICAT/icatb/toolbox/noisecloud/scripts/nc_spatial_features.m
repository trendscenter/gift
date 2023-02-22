function f = nc_spatial_features(sm, convertToZ)
%% Compute spatial map features
%

if (~exist('convertToZ', 'var'))
    convertToZ = 'no';
end

sm(isfinite(sm) == 0) = 0;
if (strcmpi(convertToZ, 'yes'))
    mask = abs(sm) > eps;
    tmp = sm(mask);
    sm(mask) = tmp ./ std(tmp);
end

tmp = nc_spatial_kurtosis(sm);
f{1} = tmp(:)';
tmp = nc_spatial_aal(sm);
f{2} = tmp(:)';
tmp = nc_spatial_degreecluster(sm);
f{3} = tmp(:)';
tmp = nc_spatial_nodedist(sm);
f{4} = tmp(:)';
tmp = nc_spatial_entropy(sm);
f{5} = tmp(:)';
tmp = nc_spatial_mirror(sm);
f{6} = tmp(:)';
tmp = nc_spatial_skewness(sm);
f{7} = tmp(:)';
tmp = nc_spatial_tissues(sm);
f{8} = tmp(:)';

f = [f{:}];



