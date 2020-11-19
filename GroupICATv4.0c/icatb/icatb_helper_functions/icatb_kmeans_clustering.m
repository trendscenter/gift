function [IDXall, Call, SUMDall, Dall] = icatb_kmeans_clustering(data, num_clusters, opts)
%% K-means clustering
%
% Inputs:
% data - numeric data or cell array of filenames
% num_clusters - Number of clusters
% opts - K-means options
%
% Outputs:
%
% IDXall - Cluster indices of each point
% Call - K cluster centroid locations
% SUMDall - within-cluster sums of point-to-centroid distances
%

max_iter = 150;
try
    max_iter = opts.max_iter;
catch
end

kmeans_distance_method = 'City';
try
    kmeans_distance_method = opts.kmeans_distance_method;
catch
end

kmeans_num_replicates = 1;
try
    kmeans_num_replicates = opts.kmeans_num_replicates;
catch
end

try
    Cp = opts.Cp;
catch
end

try
    stat_opts = opts.stat_opts;
catch
end

variableToLoad = 'FNCdyn';
try
    variableToLoad = opts.variableToLoad;
catch
end


if (isnumeric(data))
    params = {'distance', kmeans_distance_method, 'Replicates', kmeans_num_replicates, 'MaxIter', max_iter, 'Display', 'iter', 'empty', 'drop'};
else
    params = {'Replicates', kmeans_num_replicates, 'MaxIter', max_iter, 'Display', 'iter'};
end

if (exist('Cp', 'var'))
    params(end + 1) = {'Start'};
    params(end + 1) = {Cp};
end

if (exist('stat_opts', 'var'))
    params(end + 1) = {'Options'};
    params(end + 1) = {stat_opts};
end


try
    if (isnumeric(data))
        [IDXall, Call, SUMDall, Dall] = kmeans(data, num_clusters, params{:});
    else
        [IDXall, Call, SUMDall, Dall] = kmeans_tall(data, num_clusters, variableToLoad, params);
    end
catch
    [IDXall, Call, SUMDall, Dall] = icatb_kmeans(data, num_clusters, params{:});
end


function  [IDXp, Cp, Sump, Dp]  = kmeans_tall(files, N, variableToLoad, params)
% Use k-means tall algorithm

Dp = [];

variableList = whos('-file', files{1});
chk = strmatch(variableToLoad, cellstr(char(variableList.name)), 'exact');

dims = [variableList(chk).size(1), length(files)];

ds = fileDatastore(files, 'ReadFcn', @(fname) getfield(load(fname), variableToLoad), 'UniformRead', true);
tdata = tall(ds);
[IDXp, Cp, Sump]  = kmeans(tdata, N, params{:});

IDXp = gather(IDXp);

IDXp = reshape(IDXp, dims);

IDXp = permute(IDXp, [2, 1]);

IDXp = IDXp(:);