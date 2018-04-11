function [sesInfo] = icatb_calculateICA(sesInfo, statusHandle)

% ICA is performed on the reduced data set
% Input: sesInfo - structure containing all parameters necessary for group
% ica analysis


if ~exist('sesInfo','var')
    %P=icatb_spm_get(1,'*.mat','Select Parameter File');
    [P] = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Parameter File', 'filter', '*param*.mat');
    if isempty(P)
        error('Parameter file is not selected for analysis');
    end
    [pathstr,fileName]=fileparts(P);
    outputDir = pathstr;
    sesInfo.outputDir = outputDir;
    % Make sure parameter file exists
    load(P);
    if(~exist('sesInfo','var'))
        %         infoCell{1} = P;
        %         icatb_error('The selected file does not contain sesInfo variable', infoCell);
        error(['The selected file ', P, ' does not contain the sesInfo variable']);
    end
else
    outputDir = sesInfo.outputDir;
end

if ~exist('statusHandle', 'var')
    statusHandle = [];
end

if sesInfo.isInitialized == 0
    error('Parameter file has not been initialized');
end


icatb_defaults;
global INDIVIDUAL_ICA_INDEX;
global GROUP_ICA_INDEX;
global WRITE_COMPLEX_IMAGES;
global NUM_RUNS_GICA;
global OPEN_DISPLAY_GUI;


which_analysis = 1;
if (isfield(sesInfo, 'which_analysis'))
    which_analysis = sesInfo.which_analysis;
end

if (which_analysis == 2)
    icasso_opts = struct('sel_mode', 'randinit', 'num_ica_runs', max([2, NUM_RUNS_GICA]));
    if isfield(sesInfo, 'icasso_opts')
        icasso_opts = sesInfo.icasso_opts;
    end
end


icaAlgo = icatb_icaAlgorithm; % available ICA algorithms

algoVal = sesInfo.algorithm; % algorithm index

% selected ICA algorithm
algorithmName = deblank(icaAlgo(algoVal, :));


disp(' ');
disp('---------------------------------------------------------------------');


if (strcmpi(algorithmName, 'moo-icar'))
    algorithmName = 'gig-ica';
end

if (~strcmpi(algorithmName, 'iva-gl') && ~strcmpi(algorithmName, 'iva-l') && ~strcmpi(algorithmName, 'gig-ica') && ~strcmpi(algorithmName, 'constrained ica (spatial)'))
    if (which_analysis == 1)
        disp('STARTING GROUP ICA STEP ');
    elseif (which_analysis == 2)
        disp('STARTING GROUP ICA STEP USING ICASSO');
    else
        disp('STARTING GROUP ICA STEP USING MST');
    end
elseif (strcmpi(algorithmName, 'iva-gl') || strcmpi(algorithmName, 'iva-l'))
    disp('STARTING GROUP IVA STEP');
else
    disp(['STARTING ', upper(algorithmName)]);
end
disp('---------------------------------------------------------------------');


% Modality type
[modalityType, dataTitle, compSetFields] = icatb_get_modality;

if isempty(NUM_RUNS_GICA)
    NUM_RUNS_GICA = 1;
end

if NUM_RUNS_GICA < 1
    numRuns = 1;
else
    numRuns = ceil(NUM_RUNS_GICA);
end

sesInfo.num_runs_gica = numRuns;


%% Open parallel mode
parallel_info.mode = 'serial';
parallel_info.num_workers = 4;

try
    parallel_info = sesInfo.userInput.parallel_info;
catch
end


toolboxNames = ver;
parallelCluster = ~isempty(find(strcmpi(cellstr(char(toolboxNames.Name)), 'parallel computing toolbox') == 1));

if (strcmpi(parallel_info.mode, 'parallel'))
    statusHandle = [];
end

if ~isfield(sesInfo.userInput, 'dataType')
    dataType = 'real';
else
    dataType = sesInfo.userInput.dataType;
end

sesInfo.dataType = dataType;

%number of components to extract
numOfIC = sesInfo.numComp;
% data = reshape(data,xdim*ydim*zdim,size(data,2));
mask_ind = sesInfo.mask_ind;

ICA_Options = {};
% get ica options
if (isfield(sesInfo, 'ICA_Options'))
    ICA_Options = sesInfo.ICA_Options;
else
    if (isfield(sesInfo.userInput, 'ICA_Options'))
        ICA_Options = sesInfo.userInput.ICA_Options;
    end
end

% convert to cell
if isempty(ICA_Options)
    if ~iscell(ICA_Options)
        ICA_Options = {};
    end
end

if (strcmpi(algorithmName, 'iva-l') && ~isempty(ICA_Options))
    % Run IVA-L only once if GPCA is used as initialization
    chk = strmatch('type', lower(ICA_Options(1:2:end)), 'exact');
    if (isempty(chk))
        error('Initialization type is not passed');
    end
    
    initType = ICA_Options{2*chk};
    
    if (strcmpi(initType, 'gpca'))
        which_analysis = 1;
    end
end

appDataName = 'gica_waitbar_app_data';

if ~isempty(statusHandle)
    % get the status handles
    statusData = getappdata(statusHandle, appDataName);
    statusData.perCompleted = statusData.perCompleted + statusData.unitPerCompleted;
    setappdata(statusHandle, appDataName, statusData);
    set(statusHandle, 'name', [num2str(round(statusData.perCompleted*100)), '% analysis done']); waitbar(statusData.perCompleted, statusHandle);
end


useTemporalICA = 0;
if (strcmpi(modalityType, 'fmri'))
    try
        useTemporalICA = strcmpi(sesInfo.group_ica_type, 'temporal');
    catch
    end
    
    if (useTemporalICA)
        if (strcmpi(algorithmName, 'iva-gl') || strcmpi(algorithmName, 'iva-l') || strcmpi(algorithmName, 'gig-ica') || strcmpi(algorithmName, 'constrained ica (spatial)') || ...
                strcmpi(algorithmName, 'semi-blind infomax'))
            error(['Temporal ica cannot be run using algorithm ', algorithmName]);
        end
        disp('Using temporal ica ...');
    else
        if (~strcmpi(algorithmName, 'iva-gl') && ~strcmpi(algorithmName, 'iva-l'))
            disp('Using spatial ica ...');
        end
    end
end

if (~strcmpi(algorithmName, 'gig-ica') && ~strcmpi(algorithmName, 'constrained ica (spatial)'))
    data = icatb_getDataForICA(sesInfo, algorithmName, statusHandle);
else
    parICAReference(sesInfo, algorithmName, parallel_info, statusHandle);
    
    msgStr = ['DONE CALCULATING ', upper(algorithmName)];
    
    disp('---------------------------------------------------------------------');
    disp(msgStr);
    disp('---------------------------------------------------------------------');
    disp('');
    
    if (~isempty(statusHandle))
        
        % get the status handles
        statusData = getappdata(statusHandle, appDataName);
        statusData.perCompleted = statusData.perCompleted + statusData.unitPerCompleted;
        setappdata(statusHandle, appDataName, statusData);
        set(statusHandle, 'name', [num2str(round(statusData.perCompleted*100)), '% analysis done']); waitbar(statusData.perCompleted, statusHandle);
        
    end
    
    return;
    
end



if strcmpi(algorithmName, 'semi-blind infomax')
    ICA_Options = icatb_sbica_options(ICA_Options, dewhiteM);
    sesInfo.userInput.ICA_Options = ICA_Options;
end

if (~strcmpi(algorithmName, 'iva-gl') && ~strcmpi(algorithmName, 'iva-l'))
    
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
        
        if (useTemporalICA)
            temporal_icasig = icasig;
            clear icasig;
        end
        
        
    elseif (which_analysis == 2)
        % ICASSO
        
        if (strcmpi(parallel_info.mode, 'serial') || parallelCluster)
            
            %%%%% Calculate PCA and Whitening matrix %%%%%
            % PCA
            [V, Lambda] = icatb_v_pca(data, 1, numOfIC, 0, 'transpose', 'yes');
            
            % Whiten matrix
            [w, White, deWhite] = icatb_v_whiten(data, V, Lambda, 'transpose');
            
            clear V Lambda;
            
            if strcmpi(parallel_info.mode, 'serial')
                sR = icatb_icassoEst(icasso_opts.sel_mode, data, icasso_opts.num_ica_runs, 'numOfPC', numOfIC, 'algoIndex', sesInfo.algorithm, ...
                    'dewhiteM', deWhite, 'whiteM', White, 'whitesig', w, 'icaOptions', ICA_Options);
            else
                % use PCT
                sR = icatb_parIcassoEst_cluster(icasso_opts.sel_mode, data, icasso_opts.num_ica_runs, 'numOfPC', numOfIC, 'algoIndex', sesInfo.algorithm, ...
                    'dewhiteM', deWhite, 'whiteM', White, 'whitesig', w, 'icaOptions', ICA_Options);
            end
            
            clear data w deWhite White;
            
        else
            % Run ICASSOn on matlab multiple sessions
            sR = parallel_icassoEst(sesInfo);
        end
        
        %%%% Visualization %%%%%%
        
        sR = icassoExp(sR);
        
        %%% Visualization & returning results
        %%% Allow to disable visualization
        if OPEN_DISPLAY_GUI
            disp(['Launch Icasso visualization supposing ', num2str(numOfIC), ' estimate-clusters.']);
            disp('Show demixing matrix rows.');
            icassoShow(sR, 'L', numOfIC, 'estimate', 'demixing');
            
            disp(['Launch Icasso visualization supposing ', num2str(numOfIC), ' estimate-clusters.']);
            disp('Show IC source estimates (default), reduce number of lines');
            disp('Collect results.');
            iq = icassoShow(sR, 'L', numOfIC, 'colorlimit', [.8 .9]);
        else
            iq = icassoResult(sR, numOfIC);
        end
        
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
        
        icassoResultsFile = fullfile(outputDir, [sesInfo.userInput.prefix, '_icasso_results.mat']);
        try
            if (sesInfo.write_analysis_steps_in_dirs)
                icassoResultsFile = fullfile(outputDir, [sesInfo.userInput.prefix, '_ica_files', filesep, 'icasso_results.mat']);
            end
        catch
        end
        
        if (~useTemporalICA)
            icatb_save(icassoResultsFile, 'iq', 'A', 'W', 'icasig', 'sR', 'algorithmName', 'metric_Q');
        else
            temporal_icasig = icasig;
            clear icasig;
            icatb_save(icassoResultsFile, 'iq', 'A', 'W', 'temporal_icasig', 'sR', 'algorithmName', 'metric_Q');
        end
        
        clear sR;
        
    else
        % MST
        if (strcmpi(parallel_info.mode, 'serial') || parallelCluster)
            if (strcmpi(parallel_info.mode, 'serial'))
                icasigR = icatb_parMST(sesInfo, algorithmName, sesInfo.mst_opts.num_ica_runs, '', data);
            else
                icasigR = icatb_parMST_cluster(sesInfo, algorithmName, sesInfo.mst_opts.num_ica_runs, data);
            end
        else
            icasigR = parallel_mstEst(sesInfo, algorithmName, sesInfo.mst_opts.num_ica_runs);
        end
        
        [corrMetric, W, A, icasig, bestRun] = icatb_bestRunSelection(icasigR, data);
        
        clear data;
        
        mstResultsFile = fullfile(outputDir, [sesInfo.userInput.prefix, '_mst_results.mat']);
        try
            if (sesInfo.write_analysis_steps_in_dirs)
                mstResultsFile = fullfile(outputDir, [sesInfo.userInput.prefix, '_ica_files', filesep, 'mst_results.mat']);
            end
        catch
        end
        
        if (~useTemporalICA)
            icatb_save(mstResultsFile, 'A', 'W', 'icasig', 'icasigR', 'algorithmName', 'corrMetric', 'bestRun');
        else
            temporal_icasig = icasig;
            clear icasig;
            icatb_save(mstResultsFile, 'A', 'W', 'temporal_icasig', 'icasigR', 'algorithmName', 'corrMetric', 'bestRun');
        end
        
        clear icasigR;
        
    end
    
else
    
    if (which_analysis == 1)
        % IVA
        [i, W] = icatb_icaAlgorithm(algorithmName, data, ICA_Options);
        
        icasig = zeros(sesInfo.numComp, length(mask_ind));
        for n = 1:size(W, 3)
            tmp = squeeze(W(:, :, n)*data(:, :, n));
            icasig = icasig + tmp;
        end
        
        icasig = icasig / size(W, 3);
        
    else
        
        if (which_analysis == 2)
            numRuns = icasso_opts.num_ica_runs;
            disp('ICASSO is not implemented when using IVA algorithm. Using MST instead ...');
        else
            % MST
            numRuns =  sesInfo.mst_opts.num_ica_runs;
        end
        
        % Run IVA several times using MST
        if (strcmpi(parallel_info.mode, 'serial') || parallelCluster)
            if (strcmpi(parallel_info.mode, 'serial'))
                icasigR = icatb_parMST(sesInfo, algorithmName, numRuns, '', data);
            else
                icasigR = icatb_parMST_cluster(sesInfo, algorithmName, numRuns, data);
            end
        else
            icasigR = parallel_mstEst(sesInfo, algorithmName, numRuns);
        end
        
        [corrMetric, W, A, icasig, bestRun] = icatb_bestRunSelection(icasigR, data);
        icasig = squeeze(mean(icasig, 3));
        
        mstResultsFile = fullfile(outputDir, [sesInfo.userInput.prefix, '_mst_results.mat']);
        try
            if (sesInfo.write_analysis_steps_in_dirs)
                mstResultsFile = fullfile(outputDir, [sesInfo.userInput.prefix, '_ica_files', filesep, 'mst_results.mat']);
            end
        catch
        end
        
        if (~useTemporalICA)
            icatb_save(mstResultsFile, 'A', 'W', 'icasig', 'icasigR', 'algorithmName', 'corrMetric', 'bestRun');
        else
            temporal_icasig = icasig;
            clear icasig;
            icatb_save(mstResultsFile, 'A', 'W', 'temporal_icasig', 'icasigR', 'algorithmName', 'corrMetric', 'bestRun');
        end
        
        clear icasigR;
        
    end
    
end


clear data;

fprintf('\n');

if ~isempty(statusHandle)
    
    % get the status handles
    statusData = getappdata(statusHandle, appDataName);
    statusData.perCompleted = statusData.perCompleted + statusData.unitPerCompleted;
    setappdata(statusHandle, appDataName, statusData);
    set(statusHandle, 'name', [num2str(round(statusData.perCompleted*100)), '% analysis done']); waitbar(statusData.perCompleted, statusHandle);
end

clear data;clear sphere;clear signs;clear bias;clear lrates;

sesInfo.numComp = size(W, 1);
numOfIC = sesInfo.numComp;

% save in matlab format
icaout = [sesInfo.ica_mat_file, '.mat'];
icaout = fullfile(outputDir, icaout);

if (strcmpi(modalityType, 'fmri'))
    if (~useTemporalICA)
        skew = zeros(1, numOfIC);
        if (~strcmpi(algorithmName, 'iva-gl') && ~strcmpi(algorithmName, 'iva-l'))
            [A, W, icasig, skew] = changeSignOfComponents(A, icasig);
        end
    else
        icatb_save(icaout, 'temporal_icasig');
        [A, W, icasig, temporal_icasig, skew] = backReconTemporalICAComps(sesInfo, A, temporal_icasig, compSetFields);
    end
elseif (strcmpi(modalityType, 'smri'))
    
    [A, W, icasig, skew] = changeSignOfComponents(A, icasig);
    
end


% enforce dimensions to save in analyze format
if(size(icasig,1)>size(icasig,2))
    icasig = icasig';
end

drawnow;


if (exist('skew', 'var'))
    icatb_save(icaout, 'W', 'icasig', 'mask_ind', 'skew');
else
    icatb_save(icaout, 'W', 'icasig', 'mask_ind');
end

if (exist('temporal_icasig', 'var'))
    icatb_save(icaout, 'temporal_icasig', '-append');
end

msgStr = 'DONE CALCULATING GROUP IVA';

if (~strcmpi(algorithmName, 'iva-gl') && ~strcmpi(algorithmName, 'iva-l'))
    msgStr = 'DONE CALCULATING GROUP ICA';
    icatb_save(icaout, 'A', '-append');
    
    % naming of complex images
    [sesInfo, complexInfo] = icatb_name_complex_images(sesInfo, 'write');
    
    % global variable necessary for defining the real&imag or magnitude&phase
    WRITE_COMPLEX_IMAGES = sesInfo.userInput.write_complex_images;
    
    
    if strcmpi(modalityType, 'fmri')
        
        if isfield(sesInfo, 'zipContents')
            sesInfo = rmfield(sesInfo, 'zipContents');
        end
        
        sesInfo.zipContents.zipFiles = {};
        sesInfo.zipContents.files_in_zip.name = {};
        
        if(sesInfo.numOfSub ==1 & sesInfo.numOfSess==1)
            % complex data
            if ~isreal(icasig)
                % convert to complex data structure
                icasig = complex_data(icasig);
                A = complex_data(A);
            end
            
            outfile = sesInfo.icaOutputFiles(1).ses(1).name(1,:);
            [fileNames, zipfilename, files_in_zip] = icatb_saveICAData(outfile, icasig, A, mask_ind, numOfIC, sesInfo.HInfo, ...
                sesInfo.dataType, complexInfo, outputDir);
        else
            outfile = sesInfo.aggregate_components_an3_file;
            % complex data
            if ~isreal(icasig)
                % convert to complex data
                icasig = complex_data(icasig);
                A = complex_data(A);
            end
            [aggFileNames, zipfilename, files_in_zip] = icatb_saveICAData(outfile, icasig, A, mask_ind, numOfIC, sesInfo.HInfo, ...
                sesInfo.dataType, complexInfo, outputDir);
            drawnow;
        end
        
        %[pp, fileName] = fileparts(sesInfo.userInput.param_file);
        %save(fullfile(outputDir, [fileName, '.mat']), 'sesInfo');
        
        % store the zip filenames to a structure
        sesInfo.zipContents.zipFiles{1} = zipfilename;
        sesInfo.zipContents.files_in_zip(1).name = files_in_zip;
        
    end
    
end

[pp, fileName] = fileparts(sesInfo.userInput.param_file);
drawnow;

icatb_save(fullfile(outputDir, [fileName, '.mat']), 'sesInfo');

disp('---------------------------------------------------------------------');
disp(msgStr);
disp('---------------------------------------------------------------------');
disp('');

if ~isempty(statusHandle)
    
    % get the status handles
    statusData = getappdata(statusHandle, appDataName);
    statusData.perCompleted = statusData.perCompleted + statusData.unitPerCompleted;
    setappdata(statusHandle, appDataName, statusData);
    set(statusHandle, 'name', [num2str(round(statusData.perCompleted*100)), '% analysis done']); waitbar(statusData.perCompleted, statusHandle);
end


function [data] = flatrow(data)

data = data(:);


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



function [A, W, icasig, temporal_icasig, skew] = backReconTemporalICAComps(sesInfo, A, temporal_icasig, compSetFields)
%% Project timecourses on to data to get spatial maps
%

outputDir = sesInfo.outputDir;
disp('Projecting ICA timecourses on to the data to get spatial maps ...');

%% Compute back-reconstructed components
icasig = icatb_parBackReconstruct(sesInfo, (1:length(sesInfo.inputFiles)), 1);

% Compute skewness of the components
[A, W, icasig, skew] = changeSignOfComponents(A, icasig);


%% Correct sign and save back-reconstructed MAT files
skewMatrix = diag(sign(skew));

temporal_icasig = skewMatrix*temporal_icasig;

for nF = 1:length(sesInfo.inputFiles)
    
    subFile = [sesInfo.back_reconstruction_mat_file, num2str(nF), '.mat'];
    load(fullfile(outputDir, subFile), 'compSet');
    
    tc = compSet.(compSetFields{2});
    
    if (sesInfo.numOfSub*sesInfo.numOfSess == 1)
        tc = skewMatrix*tc;
    else
        tc = tc*skewMatrix;
    end
    
    compSet.(compSetFields{2}) = tc;
    compSet.(compSetFields{1}) = skewMatrix*compSet.(compSetFields{1});
    
    msgString = ['-saving back reconstructed ica data for set ', num2str(nF),' -> ',subFile];
    disp(msgString);
    drawnow;
    icatb_save(fullfile(outputDir, subFile), 'compSet');
    
end


function [A, W, icasig, skew] = changeSignOfComponents(A, icasig)
% Change sign of components

numOfIC = size(icasig, 1);
skew = zeros(1, numOfIC);
disp('Using skewness of the distribution to determine the sign of the components ...');
%force group images to be positive
for compNum = 1:numOfIC
    v = icatb_recenter_image(icasig(compNum, :));
    skew(compNum) = icatb_skewness(v) + eps;
    clear v;
    if (sign(skew(compNum)) == -1)
        disp(['Changing sign of component ',num2str(compNum)]);
        icasig(compNum, :) = icasig(compNum, :)*-1;
        A(:, compNum) = A(:, compNum)*-1;
    end
    
end

W = pinv(A);



function  icasig = parICAReference(sesInfo, algorithmName, parallel_info, statusHandle)
%% ICA with reference
%


[dd, paramFile, extn] = fileparts(sesInfo.userInput.param_file);
paramFile = fullfile(sesInfo.outputDir, [paramFile, extn]);

parallelMode = parallel_info.mode;
num_workers = 4;
try
    num_workers = parallel_info.num_workers;
catch
end
toolboxNames = ver;
parallelCluster = ~isempty(find(strcmpi(cellstr(char(toolboxNames.Name)), 'parallel computing toolbox') == 1));
dataSetsToRun = (1:sesInfo.numOfSub*sesInfo.numOfSess);
backReconFiles = cellstr(strcat(fullfile(sesInfo.outputDir, sesInfo.back_reconstruction_mat_file), icatb_numberToString(dataSetsToRun(:)), '.mat'));

if (strcmpi(parallelMode, 'serial'))
    %% Serial
    icatb_parICAReference(sesInfo, dataSetsToRun, algorithmName, statusHandle);
else
    if (parallelCluster)
        %% Use PCT
        icatb_parICAReference_cluster(sesInfo, dataSetsToRun, algorithmName);
    else
        %% Run in background mode
        cValue = clock;
        % remove a minute
        cValue(end-1) = max([0, cValue(end-1) - 1]);
        currentTime = datenum(cValue);
        giftPath = fileparts(which('gift.m'));
        dummyScriptPath = fullfile(giftPath, 'icatb_parallel_files', 'icatb_dummyScript.m');
        totalSubjectsIn = length(dataSetsToRun);
        if (totalSubjectsIn < num_workers)
            num_workers = totalSubjectsIn;
        end
        increments = ceil(totalSubjectsIn/num_workers);
        eW = 0;
        for nF = 1:num_workers
            sW = eW + 1;
            eW = eW + increments;
            eW = min([totalSubjectsIn, eW]);
            tmpDataSetsToRun = dataSetsToRun(sW:eW);
            % Run separate matlab sessions in background mode
            % (Dummyscript i.e., no text required to run in linux OS)
            commandStr = ['!matlab -nodesktop -nosplash -r "addpath(genpath(''', giftPath, '''));icatb_parICAReference(''', paramFile,  ''',[', ...
                num2str(tmpDataSetsToRun), '],''', algorithmName, ''');exit" < "', dummyScriptPath, '" &'];
            eval(commandStr);
        end
        
        icatb_waitForTaskCompletion(backReconFiles, currentTime);
    end
end

% Modality type
[modalityType, dataTitle, compSetFields] = icatb_get_modality;

countD = 0;
for nF = dataSetsToRun
    countD = countD + 1;
    load(backReconFiles{nF});
    tmp = compSet.(compSetFields{1});
    if (countD == 1)
        icasig = zeros(size(tmp));
    end
    icasig = icasig + tmp;
end
icasig = icasig./length(dataSetsToRun);
icatb_save(fullfile(sesInfo.outputDir, [sesInfo.ica_mat_file, '.mat']), 'icasig');



function sR = parallel_icassoEst(sesInfo, sel_mode, num_ica_runs)
%% Do ICASSO Estimation in parallel mode (sessions)
%

icatb_defaults;
global NUM_RUNS_GICA;


cValue = clock;
% remove a minute
cValue(end-1) = max([0, cValue(end-1) - 1]);
currentTime = datenum(cValue);

%% Initialize varaibles
num_workers = 4;
try
    num_workers = sesInfo.parallel_info.num_workers;
catch
end

if (~exist('sel_mode', 'var'))
    sel_mode = 'randinit';
    try
        sel_mode = sesInfo.icasso_opts.sel_mode;
    catch
    end
end

if (~exist('num_ica_runs', 'var'))
    num_ica_runs = max([2, NUM_RUNS_GICA]);
    try
        num_ica_runs = sesInfo.icasso_opts.num_ica_runs;
    catch
    end
end

runs = (1:num_ica_runs);

if (num_workers > num_ica_runs)
    num_workers = num_ica_runs;
end

[ddd, paramFile, extn] = fileparts(sesInfo.userInput.param_file);
paramFile = fullfile(sesInfo.outputDir, [paramFile, extn]);
giftPath = fileparts(which('gift.m'));
dummyScriptPath = fullfile(giftPath, 'icatb_parallel_files', 'icatb_dummyScript.m');
increments = ceil(num_ica_runs/num_workers);
icassoFiles = cellstr(strcat(fullfile(sesInfo.outputDir, [sesInfo.userInput.prefix, '_icasso_inc_']), icatb_numberToString((1:num_workers)'), '.mat'));

%% Run ICASSO in separate sessions
eW = 0;
for nF = 1:num_workers
    sW = eW + 1;
    eW = eW + increments;
    eW = min([num_ica_runs, eW]);
    tmpRun = length(runs(sW:eW));
    % Run separate matlab sessions in background mode
    % (Dummyscript i.e., no text required to run in linux OS)
    commandStr = ['!matlab -nodesktop -nosplash -r "addpath(genpath(''', giftPath, '''));try;rng(''shuffle'');catch;end;icatb_parIcassoEst(''', paramFile, ''', ''', sel_mode, ''',', ...
        num2str(tmpRun), ',''', icassoFiles{nF}, ''');exit" < "', dummyScriptPath, '" &'];
    eval(commandStr);
end

icatb_waitForTaskCompletion(icassoFiles, currentTime);

%% Gather output
% sR.A = {};
% sR.W = sR.A;
for nOut = 1:length(icassoFiles)
    p=load(icassoFiles{nOut});
    if (nOut == 1)
        sR = p.sR;
    else
        sR.A = [sR.A, p.sR.A];
        sR.W = [sR.W, p.sR.W];
    end
    delete(icassoFiles{nOut});%cleanup
end

n = size(sR.A{1}, 2);
index = cell(1, num_ica_runs);
for i = 1:num_ica_runs
    index{i} = [repmat(i, n, 1), (1:n)'];
end
sR.index = cat(1, index{:});




function icasigR = parallel_mstEst(sesInfo, algorithmName, num_ica_runs)
%% Do MST Estimation in parallel mode (sessions)
%

icatb_defaults;

cValue = clock;
% remove a minute
cValue(end-1) = max([0, cValue(end-1) - 1]);
currentTime = datenum(cValue);

%% Initialize varaibles
num_workers = 4;
try
    num_workers = sesInfo.parallel_info.num_workers;
catch
end

runs = (1:num_ica_runs);

if (num_workers > num_ica_runs)
    num_workers = num_ica_runs;
end

[ddd, paramFile, extn] = fileparts(sesInfo.userInput.param_file);
paramFile = fullfile(sesInfo.outputDir, [paramFile, extn]);
giftPath = fileparts(which('gift.m'));
dummyScriptPath = fullfile(giftPath, 'icatb_parallel_files', 'icatb_dummyScript.m');
increments = ceil(num_ica_runs/num_workers);
mstFiles = cellstr(strcat(fullfile(sesInfo.outputDir, [sesInfo.userInput.prefix, '_mst_inc_']), icatb_numberToString((1:num_workers)'), '.mat'));

%% Run MST in separate sessions
eW = 0;
for nF = 1:num_workers
    sW = eW + 1;
    eW = eW + increments;
    eW = min([num_ica_runs, eW]);
    tmpRun = length(runs(sW:eW));
    % Run separate matlab sessions in background mode
    % (Dummyscript i.e., no text required to run in linux OS)
    commandStr = ['!matlab -nodesktop -nosplash -r "addpath(genpath(''', giftPath, '''));try;rng(''shuffle'');catch;end;icatb_parMST(''', paramFile, ''', ''', algorithmName, ''',', ...
        num2str(tmpRun), ',''', mstFiles{nF}, ''');exit" < "', dummyScriptPath, '" &'];
    eval(commandStr);
end

icatb_waitForTaskCompletion(mstFiles, currentTime);

%% Gather output
for nOut = 1:length(mstFiles)
    p = load(mstFiles{nOut});
    tmp = [];
    try
        tmp = p.icasigR;
    catch
    end
    if (isempty(tmp))
        break;
    end
    if (nOut == 1)
        icasigR = tmp;
    else
        icasigR = [icasigR, tmp];
    end
    delete(mstFiles{nOut});
end




