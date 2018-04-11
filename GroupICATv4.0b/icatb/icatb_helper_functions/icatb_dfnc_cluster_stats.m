function [dfnc_corrs, state_vector_stats] = icatb_dfnc_cluster_stats(dfncInfo, statsDir, thresholdWindows)
%% Compute dFNC correlations and save matrix of dimensions subjects x
% sessions x component pairs
%

%% Get params
outputDir = dfncInfo.outputDir;
post_process_file = fullfile(outputDir,  [dfncInfo.prefix, '_post_process.mat']);
load (post_process_file);
numOfSub = dfncInfo.userInput.numOfSub;
numOfSess = dfncInfo.userInput.numOfSess;
M = numOfSub*numOfSess;
Nwin = length(clusterInfo.IDXall) / M;
aIND = reshape(clusterInfo.IDXall, M, Nwin);
num_comps = length(dfncInfo.comps);
num_comp_pairs = (num_comps*(num_comps - 1))/2;

if (~exist('thresholdWindows', 'var'))
    thresholdWindows = 1;
end


%% Initialise correlations to NaNs
dfnc_corrs = NaN*ones(M, num_comp_pairs, dfncInfo.postprocess.num_clusters);

%% Loop over outputfiles
for n = 1:length(dfncInfo.outputFiles)
    fname = fullfile(outputDir, dfncInfo.outputFiles{n});
    load(fname, 'FNCdyn');
    states = aIND(n, :);
    % Loop over clusters
    for nC = 1:dfncInfo.postprocess.num_clusters
        chk = find(states == nC);
        %if (~isempty(chk))
        if (length(chk) >= thresholdWindows)
            m = median(FNCdyn(chk, :));
            dfnc_corrs(n, :, nC) = m;
        end
    end
    % End of loop over clusters
end
%% End of loop over outputfiles

dfnc_corrs = reshape(dfnc_corrs, numOfSess, numOfSub, num_comp_pairs, dfncInfo.postprocess.num_clusters);


%% Compute state vector stats
aIND = aIND';
aFR = zeros(M, dfncInfo.postprocess.num_clusters);
aTM = zeros(M, dfncInfo.postprocess.num_clusters, dfncInfo.postprocess.num_clusters);
aMDT = zeros(M, dfncInfo.postprocess.num_clusters);
aNT = zeros(M, 1);
for ii = 1:M
    [FRii, TMii, MDTii, NTii] = icatb_dfnc_statevector_stats(aIND(:,ii), dfncInfo.postprocess.num_clusters);
    aFR(ii,:) = FRii;
    aTM(ii,:,:) = TMii;
    aMDT(ii, :) = MDTii;
    aNT(ii) = NTii;
end

%% Save info
state_vector_stats.frac_time_state = aFR;
state_vector_stats.mean_dwell_time = aMDT;
state_vector_stats.transition_matrix = aTM;
state_vector_stats.num_transitions = aNT;

cluster_stats_file = fullfile(statsDir,  [dfncInfo.prefix, '_cluster_stats.mat']);
save(cluster_stats_file, 'state_vector_stats', 'dfnc_corrs', 'thresholdWindows');
