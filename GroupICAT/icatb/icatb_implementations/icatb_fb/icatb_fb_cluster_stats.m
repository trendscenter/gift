function [cluster_sub_band_win] = icatb_fb_cluster_stats(dfncInfo, statsDir)
%% Compute dFNC correlations and save matrix of dimensions subjects x
% sessions x component pairs
%

%% Get params
outputDir = dfncInfo.outputDir;
post_process_file = fullfile(outputDir,  [dfncInfo.prefix, '_post_process.mat']);
load(post_process_file);
numOfSub = dfncInfo.userInput.numOfSub;
numOfSess = dfncInfo.userInput.numOfSess;

cluster_sub_band_win = zeros(length(dfncInfo.outputFiles),clusterInfo.numBands,clusterInfo.numWindows);
%% Loop over outputfiles
for n = 1:length(dfncInfo.outputFiles)
    % Sort each subject and window to its 
    % closest centroid for each frequency
    fname = fullfile(outputDir, dfncInfo.outputFiles{n});
    load(fname, 'FBdyn');
                
    % Dimensions
    [nBands, nWindows, nCorr] = size(FBdyn);
    nClusters = size(clusterInfo.Call, 1);
    
    % Reshape into observations x connectivity elements
    X = reshape(FBdyn, nBands*nWindows, nCorr);
    
    % Calculate distances to the 8 centroids
    D = sum(X.^2, 2) + sum(clusterInfo.Call.^2, 2)' - 2*(X*clusterInfo.Call');
    D = max(D, 0); % Avoid tiny negative values from numerical errors
    
    % Assign each observation to its nearest centroid
    [~, idx] = min(D, [], 2);
    
    % Restore original band and window dimensions
    cluster_sub_band_win(n,:,:) = reshape(idx, nBands, nWindows);            
    
    if (isfield(dfncInfo.userInput, 'selectedRegressors'))
        load(fname, 'task_connectivity');
        if (n == 1)
            dfncTaskConnectivity = zeros(n, size(task_connectivity, 2), size(task_connectivity, 1));
        end
        dfncTaskConnectivity(n, :, :) = task_connectivity';
    end
    
end
%% End of loop over outputfiles

nSubjects = numOfSub;
nClusters = size(clusterInfo.Dall_norm, 2);

% Verify dimensions
assert(size(clusterInfo.Dall_norm, 1) == ...
    nSubjects * nBands * nWindows);

% Convert existing fuzzy memberships into 4D matrix
members_sub_band_win_cluster = reshape( ...
    clusterInfo.Dall_norm, ...
    nWindows, nBands, nSubjects, nClusters);

members_sub_band_win_cluster = permute( ...
    members_sub_band_win_cluster, [3 2 1 4]);

cluster_stats_file = fullfile(statsDir,  [dfncInfo.prefix, '_cluster_stats_fb.mat']);
save(cluster_stats_file, 'cluster_sub_band_win', 'members_sub_band_win_cluster');

if (exist('dfncTaskConnectivity', 'var'))
    dfncTaskConnectivity = reshape(dfncTaskConnectivity, numOfSess, numOfSub, size(task_connectivity, 2), size(task_connectivity, 1));
    save(cluster_stats_file, 'dfncTaskConnectivity', '-append');
end
