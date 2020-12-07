function icatb_postprocess_spatial_chronnectome(param_file)
%% Postprocessing on spatial chronnectome data
%
% Inputs:
% 1. param_file - ICA Parameter file or spatial chronnectome parameter file
%
%

if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('title', 'Select spatial chronnectome Parameter File', 'typeEntity', 'file', 'typeSelection', ...
        'single', 'filter', '*schronn.mat');
    drawnow;
    
end


kmeans_init = 'none';
use_tall_array = 'no';
est_clusters = 'no';

if (~isstruct(param_file))
    load(param_file);
    
    outputDir = fileparts(param_file);
    if (isempty(outputDir))
        outputDir = pwd;
    end
    schronnInfo.outputDir = outputDir;
    
    answers = post_process_spatial_chronnectome;
    
    num_bins = answers.num_bins;
    num_intervals = answers.num_intervals;
    tmap_transition_matrix_threshold = answers.tmap_transition_matrix_threshold;
    est_clusters = answers.est_clusters;
    num_clusters = answers.num_clusters;
    max_iter = answers.max_iter;
    kmeans_distance_method = answers.kmeans_distance_method;
    kmeans_num_replicates = answers.kmeans_num_replicates;
    
    try
        kmeans_init = answers.kmeans_init;
    catch
    end
    
    try
        use_tall_array = answers.use_tall_array;
    catch
    end
    
else
    schronnInfo = param_file;
    
    num_bins = 10;
    num_intervals = 1;
    tmap_transition_matrix_threshold = 1.5;
    
    try
        num_bins = schronnInfo.postprocess.num_bins;
    catch
    end
    
    try
        num_intervals = schronnInfo.postprocess.num_intervals;
    catch
    end
    
    
    try
        tmap_transition_matrix_threshold = schronnInfo.postprocess.tmap_transition_matrix_threshold;
    catch
    end
    
    
    num_clusters = 4;
    try
        num_clusters = schronnInfo.postprocess.num_clusters;
    catch
    end
    
    max_iter = 150;
    try
        max_iter = schronnInfo.postprocess.max_iter;
    catch
    end
    
    kmeans_distance_method = 'City';
    try
        kmeans_distance_method = schronnInfo.postprocess.kmeans_distance_method;
    catch
    end
    
    kmeans_num_replicates = 10;
    try
        kmeans_num_replicates = schronnInfo.postprocess.kmeans_num_replicates;
    catch
    end
    
    try
        kmeans_init = schronnInfo.postprocess.kmeans_init;
    catch
    end
    
    try
        use_tall_array = schronnInfo.postprocess.use_tall_array;
    catch
    end
    
    try
        est_clusters = schronnInfo.postprocess.est_clusters;
    catch
    end
    
end

outputDir = schronnInfo.outputDir;
numComps = size(schronnInfo.outputFiles, 2);

if (length(num_clusters) == 1)
    num_clusters = repmat(num_clusters, 1, numComps);
end

schronnInfo.postprocess.num_clusters = num_clusters;
schronnInfo.postprocess.max_iter = max_iter;
schronnInfo.postprocess.kmeans_distance_method = kmeans_distance_method;
schronnInfo.postprocess.kmeans_num_replicates = kmeans_num_replicates;
schronnInfo.postprocess.num_bins = num_bins;
schronnInfo.postprocess.num_intervals = num_intervals;
schronnInfo.postprocess.tmap_transition_matrix_threshold = tmap_transition_matrix_threshold;

comps = schronnInfo.comps;
results_files = cell(1, numComps);

stat_opts = [];
try
    stat_opts = statset('UseParallel', 1, 'display', 'final');
catch
end


variableToLoad = 'dynamic_coupling_maps';


for nComps = 1:numComps
    
    SP = cell(size(schronnInfo.outputFiles, 1), 1);
    %k1_peaks = zeros(1, size(schronnInfo.outputFiles, 1));
    
    disp(['Computing kmeans on component ', num2str(comps(nComps)), ' ...']);
    data = cell(size(schronnInfo.outputFiles, 1), 1);
    for nSub = 1:size(schronnInfo.outputFiles, 1)
        tmpFName = fullfile(outputDir, schronnInfo.outputFiles{nSub, nComps});
        tmp = load(tmpFName);
        if (nSub == 1)
            dims = size(tmp.dynamic_coupling_maps);
        end
        if (strcmpi(use_tall_array, 'no'))
            data{nSub} = tmp.(variableToLoad);
        else
            data{nSub} = tmpFName;
        end
        DEV = std(tmp.dynamic_coupling_maps, [], 2);
        [xmax, imax, xmin, imin] = icatb_extrema(DEV); % find the extrema
        pIND = sort(imax);
        %k1_peaks(nSub) = length(pIND);
        tmp = tmp.(variableToLoad);
        SP{nSub} = tmp(pIND, :);
    end
    
    SPflat = cell2mat(SP);
    
    dims = [size(schronnInfo.outputFiles, 1), dims];
    
    
    if (~strcmpi(est_clusters, 'yes'))
        current_no_clusters = num_clusters(nComps);
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
    
    if (~strcmpi(kmeans_init, 'none'))
        
        [IDXp, Cp, SUMDp, Dp] = icatb_kmeans_clustering(SPflat, current_no_clusters, schronnInfo.postprocess);
        
        % Save subsampled data
        %clusterInfo.SP = SPflat;
        clusterInfo.Cp = Cp;
        
    end
    
    clear IDXp SUMDp Dp
    
    if (strcmpi(use_tall_array, 'no'))
        data = cat(3, data{:});
        data = permute(data, [3, 1, 2]);
        %dims = [size(data, 1), size(data, 2), size(data, 3)];
        data = reshape(data, dims(1)*dims(2), dims(3));
    end
    
    
    post_process_opts = schronnInfo.postprocess;
    post_process_opts.variableToLoad = variableToLoad;
    
    if (~strcmpi(kmeans_init, 'none'))
        
        post_process_opts.Cp = Cp;
        post_process_opts.kmeans_num_replicates = 1;
        
        [IDXp, Call, SUMDp, Dp] = icatb_kmeans_clustering(data, current_no_clusters, post_process_opts);
        
    else
        
        [IDXp, Call, SUMDp, Dp] = icatb_kmeans_clustering(data, current_no_clusters, post_process_opts);
        
    end
    
    clusterInfo.IDX = reshape(IDXp, dims(1:2));
    clusterInfo.IDXall = IDXp;
    clusterInfo.Call = Call;
    clusterInfo.SUMD = SUMDp;
    clusterInfo.Dp = Dp;
    clusterInfo.dims = dims;
    
    
    disp('... Saving post-processing results ...');
    
    post_process_results_dir = [schronnInfo.prefix, '_postprocess_results'];
    if (exist(fullfile(outputDir, post_process_results_dir), 'dir') ~= 7)
        mkdir(outputDir, post_process_results_dir);
    end
    
    results_file = fullfile(post_process_results_dir, [schronnInfo.prefix, '_postprocess_comp_', icatb_returnFileIndex(comps(nComps)), '.mat']);
    
    clusterFiles = cell(1, schronnInfo.postprocess.num_clusters(nComps));
    for nC = 1:size(Call, 1)
        tmpV = schronnInfo.HInfo.V(1);
        tmpV.n(1) = 1;
        tmp_fname = fullfile(post_process_results_dir, [schronnInfo.prefix, '_comp_', icatb_returnFileIndex(comps(nComps)), ...
            '_cluster_', icatb_returnFileIndex(nC), '.nii']);
        clusterFiles{nC} = tmp_fname;
        tmp3D = zeros(tmpV(1).dim(1:3));
        tmp3D(schronnInfo.mask_ind) = Call(nC, :);
        tmpV.fname = fullfile(outputDir, tmp_fname);
        icatb_write_vol(tmpV, tmp3D);
    end
    
    clusterInfo.clusterFiles = clusterFiles;
    
    if (strcmpi(use_tall_array, 'no'))
        data = reshape(data, dims);
    end
    
    
    subject_cluster_files = cell(1, current_no_clusters);
    for nc = 1:current_no_clusters
        
        subject_states = zeros([prod(tmpV.dim(1:3)), dims(1)]);
        for n = 1:dims(1)
            if (strcmpi(use_tall_array, 'no'))
                tmp = squeeze(data(n, :, :));
            else
                load(data{n}, variableToLoad);
                tmp = dynamic_coupling_maps;
                clear dynamic_coupling_maps;
            end
            idx = clusterInfo.IDX(n, :);
            inds = find(idx == nc);
            
            if (~isempty(inds))
                subject_states(schronnInfo.mask_ind, n) = mean(tmp(inds, :));
            end
            
        end
        
        tmp_dat = subject_states(schronnInfo.mask_ind, :)';
        tmp_dat(tmp_dat == 0) = NaN;
        [~, tmp_dat] = icatb_ttest(tmp_dat, 0);
        tmp_dat(isfinite(tmp_dat) == 0) = 0;
        tmp_dat = tmp_dat >= tmap_transition_matrix_threshold;
        
        if (nc == 1)
            tmap_all = tmp_dat;
        else
            tmap_all = tmp_dat | tmap_all;
        end
        
        
        subject_states = reshape(subject_states, [tmpV.dim(1:3), dims(1)]);
        subject_cluster_files{nc} = fullfile(post_process_results_dir, [schronnInfo.prefix, '_subject_comp_', icatb_returnFileIndex(comps(nComps)), '_states_', ...
            icatb_returnFileIndex(nc), '.nii']);
        icatb_write_nifti_data(fullfile(outputDir, subject_cluster_files{nc}), repmat(tmpV, dims(1), 1), subject_states);
        clear subject_states;
        
    end
    
    clusterInfo.subject_cluster_files = subject_cluster_files;
    
    
    
    %% Compute state vector stats
    aIND =  clusterInfo.IDX';
    aFR = zeros(dims(1), current_no_clusters);
    aTM = zeros(dims(1), current_no_clusters, current_no_clusters);
    aMDT = zeros(dims(1), current_no_clusters);
    aNT = zeros(dims(1), 1);
    for ii = 1:dims(1)
        [FRii, TMii, MDTii, NTii] = icatb_dfnc_statevector_stats(aIND(:,ii), current_no_clusters);
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
    
    clusterInfo.state_vector_stats = state_vector_stats;
    
    
    % Spatial_transition matrix
    glcm = zeros(dims(1), num_bins, num_bins);
    glcms_Contrast = zeros(1, dims(1));
    glcms_Correlation = glcms_Contrast;
    glcms_Energy = glcms_Contrast;
    glcms_Homogeneity = glcms_Contrast;
    
    for n = 1:dims(1)
        if (strcmpi(use_tall_array, 'no'))
            tmp = squeeze(data(n, :, :));
        else
            tmpDc = load(data{n}, variableToLoad);
            tmp = tmpDc.(variableToLoad);
        end
        tmp = tmp(:, tmap_all)';
        [dumvar2, SI] = graycomatrix(tmp, 'NumLevels', num_bins, 'Offset', [0, num_intervals]);
        glcm(n, :, :) = 100*dumvar2/(sum(dumvar2(:)) + eps);
        
        dumvar3 = graycoprops(dumvar2);
        glcms_Contrast(n) = dumvar3.Contrast;
        glcms_Correlation(n) = dumvar3.Correlation;
        glcms_Energy(n) = dumvar3.Energy;
        glcms_Homogeneity(n) = dumvar3.Homogeneity;
    end
    
    spatial_trans_stats.spatial_transition_matrix = glcm;
    spatial_trans_stats.contrast = glcms_Contrast;
    spatial_trans_stats.energy = glcms_Energy;
    spatial_trans_stats.correlation = glcms_Correlation;
    spatial_trans_stats.homogenity = glcms_Homogeneity;
    
    clusterInfo.spatial_trans_stats = spatial_trans_stats;
    
    save(fullfile(outputDir, results_file), 'clusterInfo');
    
    results_files{nComps} = results_file;
    
end

schronnInfo.postprocess.results_files = results_files;



%% Compute spatial transition matrices


save(fullfile(outputDir, [schronnInfo.prefix, '.mat']), 'schronnInfo');

disp('Done');
fprintf('\n');