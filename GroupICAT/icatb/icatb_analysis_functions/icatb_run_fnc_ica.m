function icatb_run_fnc_ica(sesInfo)
%% Run FNC ICA

icatb_defaults;
global NUM_RUNS_GICA;

if ~exist('sesInfo', 'var')
    sesInfo = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Parameter File', 'filter', '*ica*param*.mat');
    if isempty(sesInfo)
        error('Parameter file is not selected for analysis');
    end
end


if ischar(sesInfo)
    fileN = sesInfo;
    clear sesInfo;
    load(fileN);
end

if ~strcmpi(sesInfo.userInput.modality, 'fnc')
    tmpToolboxName = 'GIFT';
    if strcmpi(sesInfo.userInput.modality, 'eeg')
        tmpToolboxName = 'EEGIFT';
    elseif strcmpi(sesInfo.userInput.modality, 'smri')
        tmpToolboxName = 'SBM';
    end
    error(['Analysis selected is not FNC ICA. Please use ', tmpToolboxName, ' toolbox.']);
end

sesInfo.outputDir = sesInfo.userInput.pwd;
sesInfo.prefix = sesInfo.userInput.prefix;
sesInfo.numComp = sesInfo.userInput.numComp;
sesInfo.algorithm = sesInfo.userInput.algorithm;
sesInfo.which_analysis = sesInfo.userInput.which_analysis;
sesInfo.scaleType = sesInfo.userInput.scaleType;
sesInfo.modality = sesInfo.userInput.modality;

try
    sesInfo.icasso_opts = sesInfo.userInput.icasso_opts;
catch
end

try
    sesInfo.mst_opts = sesInfo.userInput.mst_opts;
catch
end

try
    sesInfo.cross_isi_opts = sesInfo.userInput.cross_isi_opts;
catch
end

if (sesInfo.which_analysis == 2)
    icasso_opts = struct('sel_mode', 'randinit', 'num_ica_runs', max([2, NUM_RUNS_GICA]));
    if isfield(sesInfo, 'icasso_opts')
        icasso_opts = sesInfo.icasso_opts;
    end
    sesInfo.icasso_opts = icasso_opts;
end


% performing batch analysis
output_LogFile = fullfile(sesInfo.outputDir, [sesInfo.userInput.prefix, '_results.log']);

% Print output to a file
diary(output_LogFile);

disp('Starting Analysis ');
fprintf('\n');

numComp = sesInfo.numComp;

contrast_vector = [];
try
    contrast_vector = sesInfo.userInput.contrast_vector;
catch
end

fnc_variable_mat_file = 'fnc_corrs_all';
try
    fnc_variable_mat_file = sesInfo.userInput.fnc_variable_mat_file;
catch
end


tic;

%% Load FNC data
fnc_file_name = fullfile(sesInfo.outputDir, [sesInfo.prefix, '_fnc_data.mat']);

numSessions = length(sesInfo.userInput.dataInfo);
dataInfo = sesInfo.userInput.dataInfo;
for nSess = 1:numSessions
    sessName = dataInfo(nSess).name;
    files =  cellstr(dataInfo(nSess).files);
    if (nSess == 1)
        input_data_file_patterns = cell(size(files, 1), size(files, 2));
    end
    input_data_file_patterns(:, nSess) = cellstr(files);
end

inputData.contrast_vector = contrast_vector;
inputData.fnc_variable_mat_file = fnc_variable_mat_file;
inputData.input_data_file_patterns = input_data_file_patterns;

fnc_matrix = icatb_fnc_input(inputData);
save(fnc_file_name, 'fnc_matrix');

fnc_matrix = icatb_mat2vec(fnc_matrix);

sesInfo.numOfSub = size(fnc_matrix, 1);

disp('Running data reduction step ...');

fnc_matrix = fnc_matrix';

%% Run PCA
[whitesig, dewhiteM, Lambda, V, whiteM] = icatb_calculate_pca(fnc_matrix, numComp, 'whiten', true);
pca_mat_file = fullfile(sesInfo.outputDir, [sesInfo.prefix, '_pca_r1-1.mat']);
save(pca_mat_file, 'whitesig', 'dewhiteM', 'Lambda', 'V', 'whiteM');

disp('Done data reduction step');

%% Calculate ICA
disp('Calculate ICA ...');
[sesInfo, A, W, icasig] = runICA(sesInfo, whitesig');
ica_mat_file = fullfile(sesInfo.outputDir, [sesInfo.prefix, '_ica.mat']);
save(ica_mat_file, 'A', 'W', 'icasig');
disp('Done calculate ICA');

%% Save Back-reconstructed output
disp('Saving back-reconstructed output ...');

ica_loadings = dewhiteM*A;
br_mat_file = fullfile(sesInfo.outputDir, [sesInfo.prefix, '_ica_br1.mat']);
compSet.ic = icasig;
compSet.tc = ica_loadings;
save(br_mat_file, 'compSet');
disp('Done back-reconstruction');

%% Scale components
if ~strcmpi(sesInfo.scaleType, 'no')
    disp('Converting to z-scores ...');
    ica_loadings = icatb_zscore(ica_loadings);
    icasig = icatb_zscore(icasig')';
end

loadings_fname = fullfile(sesInfo.outputDir, [sesInfo.prefix, '_group_loading_coeff_.txt']);
save(loadings_fname, 'ica_loadings', '-ascii');

ic = icatb_vec2mat(icasig);
tc = ica_loadings;

%% Save calibration step
sc_mat_file = fullfile(sesInfo.outputDir, [sesInfo.prefix, '_ica_c1-1.mat']);
save(sc_mat_file, 'ic', 'tc');

sesInfo.isInitialized = 1;

disp('Saving parameter file ...');
param_file = fullfile(sesInfo.outputDir, [sesInfo.prefix, '_ica_parameter_info.mat']);
save(param_file, 'sesInfo');

% Use tic and toc instead of cputime
t_end = toc;

disp(['Time taken to run the analysis is ', num2str(t_end), ' seconds']);

fprintf('\n');

disp(['All the analysis information is stored in the file ', output_LogFile]);

disp('Finished with Analysis');
fprintf('\n');

diary('off');


function [sesInfo, A, W, icasig] = runICA(sesInfo, data)
%% Run ICA

icatb_defaults;
global NUM_RUNS_GICA;

if isempty(NUM_RUNS_GICA)
    NUM_RUNS_GICA = 1;
end

if NUM_RUNS_GICA < 1
    numRuns = 1;
else
    numRuns = ceil(NUM_RUNS_GICA);
end

outputDir = sesInfo.outputDir;

sesInfo.num_runs_gica = numRuns;

which_analysis = sesInfo.which_analysis;
numOfIC = sesInfo.numComp;

icaAlgo = icatb_icaAlgorithm; % available ICA algorithms

algoVal = sesInfo.algorithm; % algorithm index

% selected ICA algorithm
algorithmName = deblank(icaAlgo(algoVal, :));

ICA_Options = {};
try
    ICA_Options = sesInfo.userInput.ICA_Options;
catch
end

try
    icasso_opts = sesInfo.icasso_opts;
catch
end


if (which_analysis == 1)
    
    fprintf('\n');
    disp(['Number of times ICA will run is ', num2str(numRuns)]);
    
    % Loop over number of runs
    for nRun = 1:numRuns
        
        fprintf('\n');
        
        disp(['Run ', num2str(nRun), ' / ', num2str(numRuns)]);
        
        fprintf('\n');
        
        % Run ICA
        [icaAlgo, W, A, icasig] = icatb_icaAlgorithm(algorithmName, data, ICA_Options);
        
        %A = dewhiteM*pinv(W);
        
        if (nRun == 1)
            icasig2 = zeros(numRuns, size(icasig, 1), size(icasig, 2));
        end
        
        if (nRun > 1),
            
            rho = zeros(size(icasig, 1), size(icasig, 1));
            for k1 = 1:size(icasig, 1)
                for k2 = 1:size(icasig, 1)
                    rho(k1, k2) = icatb_corr2(flatrow(icasig(k1, :)), flatrow(icasig2(1, k2, :)));
                end
            end
            % rho = (rho1+rho2)/2;
            
            Y = zeros(1, size(icasig, 1));
            I = Y;
            Ys = Y;
            
            for k = 1:size(icasig,1)
                [Y(k) I(k)] = max(abs(rho(:,k)));
                Ys(k) = rho(I(k),k);%get signed correlation
                rho(I(k), k) = 0;
            end;
            
            %reorder and force to be positively correlated
            icasig = sign(repmat(Ys', 1, size(icasig,2))).*icasig(I,:);
            A = sign(repmat(Ys, size(A,1),1)).*A(:,I);
            
        end
        
        % store icasig and A
        icasig2(nRun, :, :) = icasig;
        A2(nRun, :, :) = A;
        
    end
    % end loop over number of runs
    
    if numRuns > 1
        icasig = squeeze(mean(icasig2));
        A = squeeze(mean(A2));
        clear W;
        W = pinv(A);
    end
    
    clear icasig2;
    clear A2;
    
    
elseif (which_analysis == 2)
    % ICASSO
    
    
    %%%%% Calculate PCA and Whitening matrix %%%%%
    % PCA
    [V, Lambda] = icatb_v_pca(data, 1, numOfIC, 0, 'transpose', 'yes');
    
    % Whiten matrix
    [w, White, deWhite] = icatb_v_whiten(data, V, Lambda, 'transpose');
    
    clear V Lambda;
    
    sR = icatb_icassoEst(icasso_opts.sel_mode, data, icasso_opts.num_ica_runs, 'numOfPC', numOfIC, 'algoIndex', sesInfo.algorithm, ...
        'dewhiteM', deWhite, 'whiteM', White, 'whitesig', w, 'icaOptions', ICA_Options);
    
    clear data w deWhite White;
    
    
    
    %%%% Visualization %%%%%%
    
    sR = icassoExp(sR);
    
    %%% Visualization & returning results
    %%% Allow to disable visualization
    disp(['Launch Icasso visualization supposing ', num2str(numOfIC), ' estimate-clusters.']);
    disp('Show demixing matrix rows.');
    icassoShow(sR, 'L', numOfIC, 'estimate', 'demixing');
    
    disp(['Launch Icasso visualization supposing ', num2str(numOfIC), ' estimate-clusters.']);
    disp('Show IC source estimates (default), reduce number of lines');
    disp('Collect results.');
    iq = icassoShow(sR, 'L', numOfIC, 'colorlimit', [.8 .9]);
    
    try
        minClusterSize = icasso_opts.min_cluster_size;
    catch
        minClusterSize = 2;
    end
    
    try
        maxClusterSize = icasso_opts.max_cluster_size;
    catch
        maxClusterSize = icasso_opts.num_ica_runs;
    end
    
    if (minClusterSize <= 1)
        minClusterSize = 2;
    end
    
    if (minClusterSize > icasso_opts.num_ica_runs)
        minClusterSize = icasso_opts.num_ica_runs;
    end
    
    [metric_Q, A, W, icasig] = getStableRunEstimates(sR, minClusterSize, maxClusterSize);
    
    icassoResultsFile = fullfile(sesInfo.outputDir, [sesInfo.userInput.prefix, '_icasso_results.mat']);
    try
        if (sesInfo.write_analysis_steps_in_dirs)
            icassoResultsFile = fullfile(outputDir, [sesInfo.userInput.prefix, '_ica_files', filesep, 'icasso_results.mat']);
        end
    catch
    end
    
    icatb_save(icassoResultsFile, 'iq', 'A', 'W', 'icasig', 'sR', 'algorithmName', 'metric_Q');
    
    clear sR;
    
    
elseif (which_analysis == 4)
    % Cross isi
    WR = zeros(size(data, 1), size(data, 1), sesInfo.cross_isi_opts.num_ica_runs);
    parfor nRI = 1:sesInfo.cross_isi_opts.num_ica_runs
        fprintf('\n\n%s\n\n',['Randomization using ', algorithmName, ': Round ' num2str(nRI) '/' ...
            num2str(sesInfo.cross_isi_opts.num_ica_runs)]);
        [icaAlgoxx, WR(:, :, nRI), Axx, icxx] = icatb_icaAlgorithm(algorithmName, data, ICA_Options);
    end
    
    [bestRun, W, A, icasig, cross_isi] = RunSelection_crossISIidx(WR, data);
    
    
    criResultsFile = fullfile(outputDir, [sesInfo.userInput.prefix, '_cross_isi_results.mat']);
    
    icatb_save(criResultsFile, 'A', 'W', 'icasig', 'WR', 'algorithmName', 'bestRun', 'cross_isi');
    
else
    % MST
    
    icasigR = icatb_parMST(sesInfo, algorithmName, sesInfo.mst_opts.num_ica_runs, '', data);
    
    [corrMetric, W, A, icasig, bestRun] = icatb_bestRunSelection(icasigR, data);
    
    clear data;
    
    mstResultsFile = fullfile(outputDir, [sesInfo.userInput.prefix, '_mst_results.mat']);
    
    
    icatb_save(mstResultsFile, 'A', 'W', 'icasig', 'icasigR', 'algorithmName', 'corrMetric', 'bestRun');
    
    clear icasigR;
    
end


function [metric_Q, A, W, icasig, stableRun] = getStableRunEstimates(sR, minClusterSize, maxClusterSize)
%% Get stable run based on code by Sai Ma. Stable run estimates will be used instead of centrotype
%

% number of runs and ICs
numOfRun = length(sR.W);
numOfIC = size(sR.W{1},1);

% Get the centrotype for each cluster and Iq
index2centrotypes = icassoIdx2Centrotype(sR,'partition', sR.cluster.partition(numOfIC,:));
Iq = icassoStability(sR, numOfIC, 'none');

% Find IC index  within each cluster
partition = sR.cluster.partition(numOfIC, :);
clusterindex = cell(1, numOfIC);
for i = 1:numOfIC
    temp = (partition == i);
    clusterindex{i} = sR.index(temp, :);
    clear temp;
end

% Compute stability metric for each run within each cluster
eachRun = zeros(numOfRun, numOfIC);
qc = 0; % num of qualified clusters
for i = 1:numOfIC
    thisCluster = (clusterindex{i}(:,1))';
    clusterSize = length(clusterindex{i});
    if ((clusterSize >= minClusterSize) && (clusterSize <= maxClusterSize) && (Iq(i)>=0.7))
        qc = qc + 1;
        for k = 1:numOfRun
            thisRun = find(thisCluster == k);
            ICindex = (clusterindex{i}(thisRun,1)-1)*numOfIC + clusterindex{i}(thisRun,2);
            if ~isempty(thisRun)
                eachRun(k,i) = max(sR.cluster.similarity(index2centrotypes(i),ICindex'));
            end
            clear thisRun ICindex;
        end
    end
    clear thisCluster clusterSize;
end

%% Find stable run
metric_Q = sum(eachRun,2)/qc;
[dd, stableRun] = max(metric_Q);

%% Get stable run estimates
W = sR.W{stableRun};
clusters_stablerun = partition((stableRun - 1)*numOfIC + 1 : stableRun*numOfIC);
[dd, inds] = sort(clusters_stablerun);
W = W(inds, :);
A = pinv(W);
icasig = W*sR.signal;

