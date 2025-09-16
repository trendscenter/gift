function icatb_postprocess_sdh(param_file)
%% Post-process spatial dynamics hierarchy
%


if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('title', 'Select spatial dynamics hierarchy file', 'typeEntity', 'file', 'typeSelection', ...
        'single', 'filter', '*sdh*info*mat');
    drawnow;
    
end

%% K-means params
if (~isstruct(param_file))
    load(param_file);
    
    outputDir = fileparts(param_file);
    if (isempty(outputDir))
        outputDir = pwd;
    end
    sdhInfo.outputDir = outputDir;
    
    dlg_title = 'Cluster Options';
    distance_opts = {'City', 'sqEuclidean', 'Hamming', 'Correlation', 'Cosine'};
    
    numParameters = 1;
    inputText(numParameters).promptString = 'Do You Want To Estimate Clusters?';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = {'No', 'Yes'};
    inputText(numParameters).dataType = 'string';
    inputText(numParameters).tag = 'est_clusters                                                                                                                ';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
    numParameters = numParameters + 1;
    inputText(numParameters).promptString = 'Enter number of clusters';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = '4';
    inputText(numParameters).dataType = 'numeric';
    inputText(numParameters).tag = 'num_clusters';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
    numParameters = numParameters + 1;
    inputText(numParameters).promptString = 'Enter maximum number of iterations';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = '150';
    inputText(numParameters).dataType = 'numeric';
    inputText(numParameters).tag = 'max_iter';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString = 'Select distance method';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = distance_opts;
    inputText(numParameters).dataType = 'string';
    inputText(numParameters).tag = 'kmeans_distance_method';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString = 'Number of times to repeat the clustering';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = '10';
    inputText(numParameters).dataType = 'numeric';
    inputText(numParameters).tag = 'kmeans_num_replicates';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
    
    answers = icatb_inputDialog('inputtext', inputText, 'Title', dlg_title, 'handle_visibility', 'on');
    
    est_clusters = answers{1};
    num_clusters = answers{2};
    max_iter = answers{3};
    kmeans_distance_method = answers{4};
    kmeans_num_replicates = answers{5};
    
else
    sdhInfo = param_file;
    num_clusters = 4;
    try
        num_clusters = sdhInfo.postprocess.num_clusters;
    catch
    end
    
    max_iter = 150;
    try
        max_iter = sdhInfo.postprocess.max_iter;
    catch
    end
    
    kmeans_distance_method = 'City';
    try
        kmeans_distance_method = sdhInfo.postprocess.kmeans_distance_method;
    catch
    end
    
    kmeans_num_replicates = 10;
    try
        kmeans_num_replicates = sdhInfo.postprocess.kmeans_num_replicates;
    catch
    end
    
    est_clusters = 'no';
    try
        est_clusters = sdhInfo.postprocess.est_clusters;
    catch
    end
    
end


%% get vars
outputDir = sdhInfo.outputDir;
cd(outputDir);
mask_ind = sdhInfo.mask_ind;
tempV = sdhInfo.V;
subjectFiles = sdhInfo.subjectFiles;
prefix = sdhInfo.prefix;
comp_network_names = sdhInfo.comp_network_names;

if (length(num_clusters) == 1)
    num_clusters = repmat(num_clusters, 1, size(comp_network_names, 1));
end

sdhInfo.postprocess.num_clusters = num_clusters;
sdhInfo.postprocess.max_iter = max_iter;
sdhInfo.postprocess.kmeans_distance_method = kmeans_distance_method;
sdhInfo.postprocess.kmeans_num_replicates = kmeans_num_replicates;
sdhInfo.postprocess.est_clusters = est_clusters;


numOfSub = sdhInfo.numOfSub;
numOfSess = sdhInfo.numOfSess;
time_points = sdhInfo.time_points;

outputFiles = cell(1, size(comp_network_names, 1));

%% Load data-sets for each functional domain
for nFD = 1:size(comp_network_names, 1)
    
    fname = comp_network_names{nFD, 1};
    comps = comp_network_names{nFD, 2};
    disp(['Loading domain ', fname, ' ...']);
    endTp = 0;
    data = zeros(sum(time_points), length(mask_ind));
    SP = cell(numOfSub*numOfSess, 1);
    for countD = 1:numOfSub*numOfSess
        startTp = endTp + 1;
        endTp = endTp + time_points(countD);
        tmp = load(fullfile(outputDir, subjectFiles{countD}));
        tmpDat = tmp.tc(:, comps)*tmp.ic(comps, :);
        data(startTp:endTp, :) = tmpDat;
        DEV = std(tmpDat, [], 2);
        [xmax, imax, xmin, imin] = icatb_extrema(DEV); % find the extrema
        pIND = sort(imax);
        SP{countD} = tmpDat(pIND, :);
    end
    
    SPflat = cell2mat(SP);
    
    if (~strcmpi(est_clusters, 'yes'))
        current_no_clusters = num_clusters(nFD);
    else
        cluster_estimate_results = icatb_optimal_clusters(SPflat, min([max(size(SPflat)), 10]), 'method', 'all', 'cluster_opts', ...
            {'Replicates', kmeans_num_replicates, 'Distance', kmeans_distance_method, ...
            'MaxIter', max_iter}, 'num_tests', kmeans_num_replicates, 'display', 1);
        current_no_clusters = 0;
        for Ntests = 1:length(cluster_estimate_results)
            current_no_clusters = current_no_clusters + (cluster_estimate_results{Ntests}.K(1));
        end
        current_no_clusters = ceil(current_no_clusters/length(cluster_estimate_results));
        clusterInfo.cluster_estimate_results = cluster_estimate_results;
        disp(['Number of estimated clusters used is mean of all tests: ', num2str(current_no_clusters)]);
        fprintf('\n');
    end
    
    
    try
        [IDXp, Cp, SUMDp, Dp] = kmeans(SPflat, current_no_clusters, 'distance', kmeans_distance_method, 'Replicates', ...
            kmeans_num_replicates, 'MaxIter', max_iter, 'Display', 'iter', 'empty', 'drop');
    catch
        [IDXp, Cp, SUMDp, Dp] = icatb_kmeans(SPflat, current_no_clusters, 'distance', kmeans_distance_method, 'Replicates', ...
            kmeans_num_replicates, 'MaxIter', max_iter, 'Display', 'iter', 'empty', 'drop');
    end
    
    % Save subsampled data
    %clusterInfo.SP = SPflat;
    clusterInfo.Cp = Cp;
    
    clear IDXp SUMDp Dp
    
    %  if (strcmpi(methodInfo.name, 'k-means'))
    try
        [IDXp, Call, SUMDp, Dp] = kmeans(data, current_no_clusters, 'distance', kmeans_distance_method, 'Replicates',  1, 'MaxIter', max_iter, ...
            'Display', 'iter', 'empty', 'drop', 'Start', Cp);
    catch
        [IDXp, Call, SUMDp, Dp] = icatb_kmeans(data, current_no_clusters, 'distance', kmeans_distance_method, 'Replicates',1, 'MaxIter', max_iter, ...
            'Display', 'iter', 'empty', 'drop', 'Start', Cp);
    end
    
    clusterInfo.IDXall = IDXp;
    clusterInfo.Call = Call;
    clusterInfo.SUMDp = SUMDp;
    clusterInfo.Dp = Dp;
    
    domainInfo.clusterInfo = clusterInfo;
    domainInfo.num_clusters = current_no_clusters;
    clusterFiles = cell(1, current_no_clusters);
    for nC = 1:current_no_clusters
        tmpV = tempV;
        tmpV.n(1) = 1;
        tmp_fname = [prefix, '_', fname, '_cluster_', icatb_returnFileIndex(nC), '.nii'];
        clusterFiles{nC} = tmp_fname;
        tmp3D = zeros(tmpV(1).dim(1:3));
        tmp3D(mask_ind) = Call(nC, :);
        tmpV.fname = fullfile(outputDir, tmp_fname);        
        icatb_write_vol(tmpV, tmp3D);
    end
    
    domainInfo.clusterFiles = clusterFiles;
    
    outputFiles{nFD} = [sdhInfo.prefix, '_results_domain_', fname, '.mat'];
    
    
    icatb_parSave(fullfile(outputDir, outputFiles{nFD}), {domainInfo}, {'domainInfo'});
    
end


sdhInfo.postprocess.outputFiles = outputFiles;

fname = fullfile(outputDir, sdhInfo.param_file);
disp(['Saving parameters info in ', fname]);
save(fname, 'sdhInfo');
fprintf('\n\n');

