function [gIdx, c, sumD] = icatb_compute_complex_kmeans(Corrected_Out_Distribution, coin, k_cluster, numRuns)
%% Compute complex k-means
%
% Inputs:
%
% 1. Corrected_Out_Distribution -
% 2. coin -
% 3. k_cluster - Number of clusters
% 4. numRuns - Number of runs
%
% Outputs:
% 1. gIdx - Indices
% 2. c - Centroids
% 3. sumD - Distances
%

if (~exist('numRuns', 'var'))
    numRuns = 1;
end

%%
incoi_mask_v = find(coin);
incoi_mask = coin;



%%
Corrected_Out_Distribution = permute(Corrected_Out_Distribution, [1 3 2]);


Corrected_Out_Distribution_v = reshape(Corrected_Out_Distribution, [], size(Corrected_Out_Distribution,3));
clearvars Corrected_Out_Distribution;

gIdx = cell(1, numRuns);
c = cell(1, numRuns);
sumD = cell(1, numRuns);

for run_num=1:numRuns
    disp(['Kmeans Run #', num2str(run_num) ' / ', num2str(numRuns)]);
    [gIdx{run_num},c{run_num}, sumD{run_num}] = icatb_complex_k_means(Corrected_Out_Distribution_v(:,:),k_cluster);
end

