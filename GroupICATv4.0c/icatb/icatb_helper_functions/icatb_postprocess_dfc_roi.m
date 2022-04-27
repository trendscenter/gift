function icatb_postprocess_dfc_roi(param_file)
%% Post-process dfc ROI correlations
%


if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('title', 'Select dFC ROI param file', 'typeEntity', 'file', 'typeSelection', ...
        'single', 'filter', '*dfc*roi.mat');
    drawnow;
    
end

%% K-means params
if (~isstruct(param_file))
    load(param_file);
    
    outputDir = fileparts(param_file);
    if (isempty(outputDir))
        outputDir = pwd;
    end
    dfcRoiInfo.outputDir = outputDir;
    
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
    dfcRoiInfo = param_file;
    num_clusters = 4;
    try
        num_clusters = dfcRoiInfo.postprocess.num_clusters;
    catch
    end
    
    max_iter = 150;
    try
        max_iter = dfcRoiInfo.postprocess.max_iter;
    catch
    end
    
    kmeans_distance_method = 'City';
    try
        kmeans_distance_method = dfcRoiInfo.postprocess.kmeans_distance_method;
    catch
    end
    
    kmeans_num_replicates = 10;
    try
        kmeans_num_replicates = dfcRoiInfo.postprocess.kmeans_num_replicates;
    catch
    end
    
    est_clusters = 'no';
    try
        est_clusters = dfcRoiInfo.postprocess.est_clusters;
    catch
    end
    
end


%% get vars
tempVol = dfcRoiInfo.userInput.V;
analysisType = dfcRoiInfo.userInput.analysisType;
outputDir = dfcRoiInfo.outputDir;
cd(outputDir);
subjectFiles = dfcRoiInfo.outputFiles;
dfc_prefix = [dfcRoiInfo.userInput.prefix, '_dfc_roi'];
roi_info = dfcRoiInfo.userInput.roi_info;
roi_labels = roi_info{1};
roi_indices = roi_info{2};
roi_labels_sel = roi_labels(roi_indices);

if (strcmpi(analysisType, 'roi-voxel'))
    if (length(num_clusters) == 1)
        num_clusters = repmat(num_clusters, 1, length(roi_indices));
    end
    outputFiles = cell(1, length(roi_indices));
else
    num_clusters = num_clusters(1);
    outputFiles = cell(1, 1);
end

dfcRoiInfo.postprocess.num_clusters = num_clusters;
dfcRoiInfo.postprocess.max_iter = max_iter;
dfcRoiInfo.postprocess.kmeans_distance_method = kmeans_distance_method;
dfcRoiInfo.postprocess.kmeans_num_replicates = kmeans_num_replicates;
dfcRoiInfo.postprocess.est_clusters = est_clusters;
masks_y = dfcRoiInfo.masks_y;

numOfDataSets = dfcRoiInfo.userInput.numOfDataSets;


%% Load data-sets for each functional domain
for nFD = 1:length(outputFiles)
    
    if (strcmpi(analysisType, 'roi-roi'))
        disp('Loading dfC correlations ...');
    else
        current_label = roi_labels_sel{nFD};
        disp(['Loading dFC correlations ', current_label, ' ...']);
    end
    data = cell(numOfDataSets, 1);
    SP = cell(numOfDataSets, 1);
    for countD = 1:numOfDataSets
        tmp = load(fullfile(outputDir, subjectFiles{countD}));
        if (iscell(tmp.FNCdyn))
            tmpDat = tmp.FNCdyn{nFD};
            tmpDat = tmpDat(:, masks_y{nFD});
        else
            tmpDat = tmp.FNCdyn;
        end
        DEV = std(tmpDat, [], 2);
        [xmax, imax, xmin, imin] = icatb_extrema(DEV); % find the extrema
        pIND = sort(imax);
        SP{countD} = tmpDat(pIND, :);
        data{countD} = tmpDat;
    end
    
    SPflat = cell2mat(SP);
    
    if (~strcmpi(est_clusters, 'yes'))
        if (strcmpi(analysisType, 'roi-voxel'))
            current_no_clusters = num_clusters(nFD);
        else
            current_no_clusters = num_clusters;
        end
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
    
    %data = cell2mat(data);
    data = cat(3, data{:});
    data = permute(data, [3, 1, 2]);
    data = reshape(data, size(data, 1)*size(data, 2), size(data, 3));
    
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
    clusterInfo.num_clusters = current_no_clusters;
    
    % Save state-wise windowed corrs and subject-wise state vectors
    [clusterInfo.corrs_states, clusterInfo.states] = getStateCorr(dfcRoiInfo, clusterInfo); 
    
    if (strcmpi(analysisType, 'roi-voxel'))
        clusterFiles = cell(1, current_no_clusters);
        for nC = 1:current_no_clusters
            tmpV = tempVol;
            tmpV.n(1) = 1;
            tmp_fname = [dfc_prefix, '_', 'label_', icatb_returnFileIndex(nFD), '_cluster_', icatb_returnFileIndex(nC), '.nii'];
            clusterFiles{nC} = tmp_fname;
            tmp3D = zeros(tmpV(1).dim(1:3));
            tmp3D(masks_y{nFD}) = Call(nC, :);
            tmpV.fname = fullfile(outputDir, tmp_fname);
            icatb_write_vol(tmpV, tmp3D);
        end
        clusterInfo.clusterFiles = clusterFiles;
    end
    
    
    outputFiles{nFD} = [dfcRoiInfo.prefix, '_cluster_results_', icatb_returnFileIndex(nFD), '.mat'];
    
    
    icatb_parSave(fullfile(outputDir, outputFiles{nFD}), {clusterInfo}, {'clusterInfo'});
    
end


dfcRoiInfo.postprocess.outputFiles = outputFiles;
dfcRoiInfo.postprocess.roi_labels_sel = roi_labels_sel;

fname = fullfile(outputDir, dfcRoiInfo.param_file);
disp(['Saving parameters info in ', fname]);
save(fname, 'dfcRoiInfo');
fprintf('\n\n');


function [corrs, states] = getStateCorr(dfcRoiInfo, clusterInfo)
% corrs - subjects x sessions x clusters
% states - subjects x sessions x windows
outDir = dfcRoiInfo.outputDir;

Nwin = size(clusterInfo.IDXall, 1)/(dfcRoiInfo.userInput.filesInfo.numOfSess*dfcRoiInfo.userInput.numOfDataSets);

% subjects x sessions x windows
states = reshape(clusterInfo.IDXall, dfcRoiInfo.userInput.filesInfo.numOfSess, dfcRoiInfo.userInput.numOfDataSets, Nwin);
states = permute(states, [2, 1, 3]);

bInit = true;
nVoxOpt = 1;
nVoxOptMax = 2; %temporary variable for while below
while nVoxOpt <= nVoxOptMax %loop for voxel option
    count = 0;
    clear subCorrs; %voxel option needs several correlation matrices, generating sub level
    for nSub = 1:dfcRoiInfo.userInput.numOfDataSets
        for nSess = 1:dfcRoiInfo.userInput.filesInfo.numOfSess
            count = count + 1;
            fn = fullfile(outDir, dfcRoiInfo.outputFiles{count});
            load(fn);
            tmp = squeeze(states(nSub, nSess, :));
            if iscell(FNCdyn) 
                % It is voxel option
                bVox = true;
                if bInit
                    % Initiates corrs with number of matrices
                    nVoxOptMax = length(FNCdyn); %number of cell strings needed
                    corrs=cell(nVoxOptMax,1);
                    bInit = false;
                end
                for nState = 1:dfcRoiInfo.postprocess.num_clusters
                    inds = tmp == nState;
                    subCorrs{nSub, nSess, nState} = FNCdyn{nVoxOpt}(inds, :);
                end
            else
                % It is regular (non voxel option) and only one cell of matrices
                bVox = false;
                clear bInit ; %not needed for non voxel
                for nState = 1:dfcRoiInfo.postprocess.num_clusters
                    inds = tmp == nState;
                    corrs{nSub, nSess, nState} = FNCdyn(inds, :);
                end
                nVoxOptMax = 0; % looping not needed for non-voxel
            end
        end
    end
    if bVox
        corrs{nVoxOpt} = subCorrs; %save several matrices for voxel opt.
    end
    nVoxOpt = nVoxOpt + 1;
end