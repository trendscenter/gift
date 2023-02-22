function icatb_display_dynamic_coherence(param_file)
%% Display dynamic coherence results
%
%

%% Select dynamic coherence param file
if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('title', 'Select dynamic coherence parameter file', 'typeEntity', 'file', 'typeSelection', 'single', 'filter', '*dyn_coh.mat');
    drawnow;
    if (isempty(param_file))
        error('Dynamic Coherence parameter file is not selected');
    end
end

drawnow;

load(param_file);

if (~exist('cohInfo', 'var'))
    error('Selected file is not a valid dynamic coherence parameter file');
end

outDir = fileparts(param_file);

if (isempty(outDir))
    outDir = pwd;
end


cd(outDir);

load(fullfile(outDir, [cohInfo.userInput.prefix, '.mat']));

numComps = length(cohInfo.comps);

outFile = fullfile(outDir, [cohInfo.prefix, '_results.mat']);

load(outFile);

cluster_metric = cellfun(@sum, clusterInfo.sumD);

[tval, tbest] = min(cluster_metric);

best_cluster_results = clusterInfo.IDX{tbest};
tcentroids = clusterInfo.centroids{tbest};
best_cluster_centroids = tcentroids(:,1:end);% +1i*tcentroids(:,end/2+1:end);
best_cluster_metric = tval;


icatb_plot_dyn_coh_clusters(best_cluster_centroids, best_cluster_results, coin, cohInfo.comps, cohInfo.TRN);

