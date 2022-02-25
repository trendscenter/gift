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


if (strcmpi(algorithmName, 'gig-ica'))
    algorithmName = 'moo-icar';
end

if (isempty(icatb_findstr(lower(algorithmName),'iva')) && ~strcmpi(algorithmName, 'moo-icar') && ~icatb_string_compare(algorithmName, 'constrained'))
    if (which_analysis == 1)
        disp('STARTING GROUP ICA STEP ');
    elseif (which_analysis == 2)
        disp('STARTING GROUP ICA STEP USING ICASSO');
    elseif (which_analysis == 3)
        disp('STARTING GROUP ICA STEP USING MST');
    else
        disp('STARTING GROUP ICA STEP USING Cross ISI');
    end
elseif (~isempty(icatb_findstr(lower(algorithmName),'iva')))
    %(strcmpi(algorithmName, 'iva-gl') || strcmpi(algorithmName, 'iva-l') || strcmpi(algorithmName, 'iva-l-sos'))
    disp('STARTING GROUP IVA STEP');
    if (strcmpi(algorithmName, 'iva-l-sos-adaptive'))
        disp('Using IVA-L-SOS-Adaptive algorithm ...');
        
    end
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
        if (~isempty(icatb_findstr(lower(algorithmName),'iva')) || strcmpi(algorithmName, 'moo-icar') || icatb_string_compare(algorithmName, 'constrained') || ...
                strcmpi(algorithmName, 'semi-blind infomax'))
            error(['Temporal ica cannot be run using algorithm ', algorithmName]);
        end
        disp('Using temporal ica ...');
    else
        if (isempty(icatb_findstr(lower(algorithmName),'iva')))
            disp('Using spatial ica ...');
        end
    end
end

if (~strcmpi(algorithmName, 'moo-icar') && ~icatb_string_compare(algorithmName, 'constrained'))
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
    load(fullfile(outputDir, [sesInfo.data_reduction_mat_file, '1-1']),'dewhiteM');
    ICA_Options = icatb_sbica_options(ICA_Options, dewhiteM);
    sesInfo.userInput.ICA_Options = ICA_Options;
end

if (isempty(icatb_findstr(lower(algorithmName),'iva')))
    
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
        try
            if (sesInfo.write_analysis_steps_in_dirs)
                criResultsFile = fullfile(outputDir, [sesInfo.userInput.prefix, '_ica_files', filesep, 'cross_isi_results.mat']);
            end
        catch
        end
        
        icatb_save(criResultsFile, 'A', 'W', 'icasig', 'WR', 'algorithmName', 'bestRun', 'cross_isi');
        
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
        
        if (useTemporalICA)
            which_analysis = 3;
        end
        
        if (which_analysis == 2)
            numRuns = icasso_opts.num_ica_runs;
            disp('ICASSO is not implemented when using IVA algorithm. Using MST instead ...');
        elseif (which_analysis == 3)
            % MST
            numRuns =  sesInfo.mst_opts.num_ica_runs;
        else
            % disp('Cross ISI is not implemented when using IVA algorithm. Using MST instead ...');
            numRuns =  sesInfo.cross_isi_opts.num_ica_runs;
        end
        
        if (which_analysis ~= 4)
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
            
            
        else
            % Cross ISI
            WR = cell(1, sesInfo.cross_isi_opts.num_ica_runs);
            parfor nRI = 1:sesInfo.cross_isi_opts.num_ica_runs
                fprintf('\n\n%s\n\n',['Randomization using ', algorithmName, ': Round ' num2str(nRI) '/' ...
                    num2str(sesInfo.cross_isi_opts.num_ica_runs)]);
                [icaAlgoxx, WR{nRI}, Axx, icxx] = icatb_icaAlgorithm(algorithmName, data, ICA_Options);
            end
            
            [selected_run, cross_isi_joint] = run_selection_crossJointISI(WR);
            
            bestRun = selected_run(1);
            [W, A, icasig]  = getIVASig(WR{bestRun}, data);
            
            criResultsFile = fullfile(outputDir, [sesInfo.userInput.prefix, '_cross_isi_results.mat']);
            try
                if (sesInfo.write_analysis_steps_in_dirs)
                    criResultsFile = fullfile(outputDir, [sesInfo.userInput.prefix, '_ica_files', filesep, 'cross_isi_results.mat']);
                end
            catch
            end
            
            icatb_save(criResultsFile, 'A', 'W', 'icasig', 'algorithmName', 'cross_isi_joint', 'bestRun');
        end
        
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
        if (isempty(icatb_findstr(lower(algorithmName),'iva')))
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

if (isempty(icatb_findstr(lower(algorithmName),'iva')))
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




function [consistentRunidx, W, A, icasig, crossISIW] = RunSelection_crossISIidx(W, data)
%% Select the most consistent run based on Cross ISI.
% Input:
% W - demixing matrices of different runs, with dimension as the number of components x the number of
% components x the number of runs.
% Outputs:
% consistentRun_idx - index of the run with the lowest average pairwise
% cross ISI values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code by Chunying Jia.
% Reference:
% Q. Long, C. Jia, Z. Boukouvalas, B. Gabrielson, D. Emge, and T. Adali, "CONSISTENT RUN SELECTION FOR INDEPENDENT COMPONENT ANALYSIS:
% APPLICATION TO FMRI ANALYSIS", 2018 IEEE International Conference on Acoustics, Speech and Signal Processing.
% Please contact chunyin1@umbc.edu or qunfang1@umbc.edu.
N_run = size(W,3);
crossISIW = zeros(N_run,1);

for i = 1:N_run
    crossISI_ij = 0;
    for j = 1:N_run
        if j==i
            continue;
        end
        crossISI_ij = crossISI_ij + calculate_ISI(W(:,:,j),inv(W(:,:,i)));
    end
    crossISIW(i,1) = crossISI_ij / (N_run -1);
end
[~,consistentRunidx] = min(crossISIW);


W = squeeze(W(:, :, consistentRunidx));
A = pinv(W);
icasig = W*data;

function ISI = calculate_ISI(W,A)

p = W*A;
p = abs(p);

N = size(p,1);

b1 = 0;
for i = 1:N
    a1 = 0;
    max_pij = max(p(i,:));
    for j = 1:N
        
        pij = p(i,j);
        a1 = pij/max_pij + a1;
    end
    b1 = a1 - 1 + b1;
end

b2 = 0;
for j = 1:N
    a2 = 0;
    max_pij = max(p(:,j));
    for i = 1:N
        
        pij = p(i,j);
        a2 = pij/max_pij + a2;
    end
    b2 = a2 - 1 + b2;
end

ISI = (b1+b2)/((N-1)*N*2);




function [selected_run, crossisi_jnt, crossisi_avg] = run_selection_crossJointISI(W)

% Input:
%       W : cell array of length equal to number of runs,
%           where each cell is a N x N x K matrix, where N is the number of
%           components, K is the number of datasets
% Output:
%       selected_run(1) : run selected using cross joint ISI
%       selected_run(2) : run selected using cross average ISI
%
% Code written by Suchita Bhinge (suchita1 at umbc.edu)
%
% References
%
%[1] Long, Q., C. Jia, Z. Boukouvalas, B. Gabrielson, D. Emge, and T. Adali.
%    "Consistent run selection for independent component analysis: Application
%    to fMRI analysis." IEEE International Conference on Acoustics, Speech and
%    Signal Processing (ICASSP), 2018.

MaxRuns = length(W);

for runs = 1 : MaxRuns
    for k = 1 : size(W{1},3)
        A{runs}(:,:,k) = inv(W{runs}(:,:,k));
    end
end
for runs = 1 : MaxRuns
    idx = setdiff(1:MaxRuns,runs);
    for j = 1: length(idx)
        %[isi(j),isiGrp(j)]=bss_isi(W{runs},A{idx(j)});
        [isi(j),isiGrp(j)]=bss_isi(W{idx(j)},A{runs}); % modified by Chunying Jia (chunyin1@umbc.edu)
    end
    crossisi_avg(runs) = mean(isi);clear isi
    crossisi_jnt(runs) = mean(isiGrp);clear isiGrp idx
end

[~,idx] = min(crossisi_jnt);
selected_run(1) = idx;
[~,idx] = min(crossisi_avg);
selected_run(2) = idx;
clear idx

function [isi,isiGrp,success,G]=bss_isi(W,A,s,Nuse)
% Non-cell inputs:
% isi=bss_isi(W,A) - user provides W & A where x=A*s, y=W*x=W*A*s
% isi=bss_isi(W,A,s) - user provides W, A, & s
%
% Cell array of matrices:
% [isi,isiGrp]=bss_isi(W,A) - W & A are cell array of matrices
% [isi,isiGrp]=bss_isi(W,A,s) - W, A, & s are cell arrays
%
% 3-d Matrices:
% [isi,isiGrp]=bss_isi(W,A) - W is NxMxK and A is MxNxK
% [isi,isiGrp]=bss_isi(W,A,s) - S is NxTxK (N=#sources, M=#sensors, K=#datasets)
%
% Measure of quality of separation for blind source separation algorithms.
% W is the estimated demixing matrix and A is the true mixing matrix.  It should be noted
% that rows of the mixing matrix should be scaled by the necessary constants to have each
% source have unity variance and accordingly each row of the demixing matrix should be
% scaled such that each estimated source has unity variance.
%
% ISI is the performance index given in Complex-valued ICA using second order statisitcs
% Proceedings of the 2004 14th IEEE Signal Processing Society Workshop, 2004, 183-192
%
% Normalized performance index (Amari Index) is given in Choi, S.; Cichocki, A.; Zhang, L.
% & Amari, S. Approximate maximum likelihood source separation using the natural gradient
% Wireless Communications, 2001. (SPAWC '01). 2001 IEEE Third Workshop on Signal
% Processing Advances in, 2001, 235-238.
%
% Note that A is p x M, where p is the number of sensors and M is the number of signals
% and W is N x p, where N is the number of estimated signals.  Ideally M=N but this is not
% guaranteed.  So if N > M, the algorithm has estimated more sources than it "should", and
% if M < N the algorithm has not found all of the sources.  This meaning of this metric is
% not well defined when averaging over cases where N is changing from trial to trial or
% algorithm to algorithm.

% Some examples to consider
% isi=bss_isi(eye(n),eye(n))=0
%
% isi=bss_isi([1 0 0; 0 1 0],eye(3))=NaN
%


% Should ideally be a permutation matrix with only one non-zero entry in any row or
% column so that isi=0 is optimal.

% generalized permutation invariant flag (default=false), only used when nargin<3
gen_perm_inv_flag=false;
success=true;

Wcell=iscell(W);
if nargin<2
    Acell=false;
else
    Acell=iscell(A);
end
if ~Wcell && ~Acell
    if ndims(W)==2 && ndims(A)==2
        if nargin==2
            % isi=bss_isi(W,A) - user provides W & A
            
            % Traditional Metric, user provided W & A separately
            G=W*A;
            [N,M]=size(G);
            Gabs=abs(G);
            if gen_perm_inv_flag
                % normalization by row
                max_G=max(Gabs,[],2);
                Gabs=repmat(1./max_G,1,size(G,2)).*Gabs;
            end
        elseif nargin==3
            % Equalize energy associated with each estimated source and true
            % source.
            %
            % y=W*A*s;
            % snorm=D*s; where snorm has unit variance: D=diag(1./std(s,0,2))
            % Thus: y=W*A*inv(D)*snorm
            % ynorm=U*y; where ynorm has unit variance: U=diag(1./std(y,0,2))
            % Thus: ynorm=U*W*A*inv(D)*snorm=G*snorm and G=U*W*A*inv(D)
            
            y=W*A*s;
            D=diag(1./std(s,0,2));
            U=diag(1./std(y,0,2));
            G=U*W*A/D; % A*inv(D)
            [N,M]=size(G);
            Gabs=abs(G);
        else
            error('Not acceptable.')
        end
        
        isi=0;
        for n=1:N
            isi=isi+sum(Gabs(n,:))/max(Gabs(n,:))-1;
        end
        for m=1:M
            isi=isi+sum(Gabs(:,m))/max(Gabs(:,m))-1;
        end
        isi=isi/(2*N*(N-1));
        isiGrp=NaN;
        success=NaN;
    elseif ndims(W)==3 && ndims(A)==3
        % IVA/GroupICA/MCCA Metrics
        % For this we want to average over the K groups as well as provide the additional
        % measure of solution to local permutation ambiguity (achieved by averaging the K
        % demixing-mixing matrices and then computing the ISI of this matrix).
        [N,M,K]=size(W);
        if M~=N
            error('This more general case has not been considered here.')
        end
        L=M;
        
        isi=0;
        GabsTotal=zeros(N,M);
        G=zeros(N,M,K);
        for k=1:K
            if nargin<=2
                % Traditional Metric, user provided W & A separately
                Gk=W(:,:,k)*A(:,:,k);
                Gabs=abs(Gk);
                if gen_perm_inv_flag
                    % normalization by row
                    max_G=max(Gabs,[],2);
                    Gabs=repmat(1./max_G,1,size(Gabs,2)).*Gabs;
                end
            else %if nargin==3
                % Equalize energy associated with each estimated source and true
                % source.
                %
                % y=W*A*s;
                % snorm=D*s; where snorm has unit variance: D=diag(1./std(s,0,2))
                % Thus: y=W*A*inv(D)*snorm
                % ynorm=U*y; where ynorm has unit variance: U=diag(1./std(y,0,2))
                % Thus: ynorm=U*W*A*inv(D)*snorm=G*snorm and G=U*W*A*inv(D)
                yk=W(:,:,k)*A(:,:,k)*s(:,:,k);
                Dk=diag(1./std(s(:,:,k),0,2));
                Uk=diag(1./std(yk,0,2));
                Gk=Uk*W(:,:,k)*A(:,:,k)/Dk;
                
                Gabs=abs(Gk);
            end
            G(:,:,k)=Gk;
            
            if nargin>=4
                Np=Nuse;
                Mp=Nuse;
                Lp=Nuse;
            else
                Np=N;
                Mp=M;
                Lp=L;
            end
            
            % determine if G is success by making sure that the location of maximum magnitude in
            % each row is unique.
            if k==1
                [~,colMaxG]=max(Gabs,[],2);
                if length(unique(colMaxG))~=Np
                    % solution is failure in strictest sense
                    success=false;
                end
            else
                [~,colMaxG_k]=max(Gabs,[],2);
                if ~all(colMaxG_k==colMaxG)
                    % solution is failure in strictest sense
                    success=false;
                end
            end
            
            GabsTotal=GabsTotal+Gabs;
            
            for n=1:Np
                isi=isi+sum(Gabs(n,:))/max(Gabs(n,:))-1;
            end
            for m=1:Mp
                isi=isi+sum(Gabs(:,m))/max(Gabs(:,m))-1;
            end
        end
        isi=isi/(2*Np*(Np-1)*K);
        
        Gabs=GabsTotal;
        if gen_perm_inv_flag
            % normalization by row
            max_G=max(Gabs,[],2);
            Gabs=repmat(1./max_G,1,size(Gabs,2)).*Gabs;
        end
        %       figure; imagesc(Gabs); colormap('bone'); colorbar
        isiGrp=0;
        for n=1:Np
            isiGrp=isiGrp+sum(Gabs(n,:))/max(Gabs(n,:))-1;
        end
        for m=1:Mp
            isiGrp=isiGrp+sum(Gabs(:,m))/max(Gabs(:,m))-1;
        end
        isiGrp=isiGrp/(2*Lp*(Lp-1));
    else
        error('Need inputs to all be of either dimension 2 or 3')
    end
elseif Wcell && Acell
    % IVA/GroupICA/MCCA Metrics
    % For this we want to average over the K groups as well as provide the additional
    % measure of solution to local permutation ambiguity (achieved by averaging the K
    % demixing-mixing matrices and then computing the ISI of this matrix).
    
    K=length(W);
    N=0; M=0;
    Nlist=zeros(K,1);
    for k=1:K
        Nlist(k)=size(W{k},1);
        N=max(size(W{k},1),N);
        M=max(size(A{k},2),M);
    end
    commonSources=false; % limits the ISI to first min(Nlist) sources
    if M~=N
        error('This more general case has not been considered here.')
    end
    L=M;
    
    % To make life easier below lets sort the datasets to have largest
    % dataset be in k=1 and smallest at k=K;
    [Nlist,isort]=sort(Nlist,'descend');
    W=W(isort);
    A=A(isort);
    if nargin > 2
        s=s(isort);
    end
    G=cell(K,1);
    isi=0;
    if commonSources
        minN=min(Nlist);
        GabsTotal=zeros(minN);
        Gcount=zeros(minN);
    else
        GabsTotal=zeros(N,M);
        Gcount=zeros(N,M);
    end
    for k=1:K
        if nargin==2
            % Traditional Metric, user provided W & A separately
            G{k}=W{k}*A{k};
            Gabs=abs(G{k});
            if gen_perm_inv_flag
                % normalization by row
                max_G=max(Gabs,[],2);
                Gabs=repmat(1./max_G,1,size(Gabs,2)).*Gabs;
            end
        elseif nargin>=3
            % Equalize energy associated with each estimated source and true
            % source.
            %
            % y=W*A*s;
            % snorm=D*s; where snorm has unit variance: D=diag(1./std(s,0,2))
            % Thus: y=W*A*inv(D)*snorm
            % ynorm=U*y; where ynorm has unit variance: U=diag(1./std(y,0,2))
            % Thus: ynorm=U*W*A*inv(D)*snorm=G*snorm and G=U*W*A*inv(D)
            yk=W{k}*A{k}*s{k};
            Dk=diag(1./std(s{k},0,2));
            Uk=diag(1./std(yk,0,2));
            G{k}=Uk*W{k}*A{k}/Dk;
            
            Gabs=abs(G{k});
        else
            error('Not acceptable.')
        end
        
        if commonSources
            Nk=minN;
            Gabs=Gabs(1:Nk,1:Nk);
        elseif nargin>=4
            commonSources=true;
            Nk=Nuse;
            minN=Nk;
        else
            Nk=Nlist(k);
        end
        
        if k==1
            [~,colMaxG]=max(Gabs(1:Nk,1:Nk),[],2);
            if length(unique(colMaxG))~=Nk
                % solution is a failure in a strict sense
                success=false;
            end
        elseif success
            if nargin>=4
                [~,colMaxG_k]=max(Gabs(1:Nk,1:Nk),[],2);
            else
                [~,colMaxG_k]=max(Gabs,[],2);
            end
            if ~all(colMaxG_k==colMaxG(1:Nk))
                % solution is a failure in a strict sense
                success=false;
            end
        end
        
        if nargin>=4
            GabsTotal(1:Nk,1:Nk)=GabsTotal(1:Nk,1:Nk)+Gabs(1:Nk,1:Nk);
        else
            GabsTotal(1:Nk,1:Nk)=GabsTotal(1:Nk,1:Nk)+Gabs;
        end
        Gcount(1:Nk,1:Nk)=Gcount(1:Nk,1:Nk)+1;
        for n=1:Nk
            isi=isi+sum(Gabs(n,:))/max(Gabs(n,:))-1;
        end
        for m=1:Nk
            isi=isi+sum(Gabs(:,m))/max(Gabs(:,m))-1;
        end
        isi=isi/(2*Nk*(Nk-1));
    end
    
    if commonSources
        Gabs=GabsTotal;
    else
        Gabs=GabsTotal./Gcount;
    end
    % normalize entries into Gabs by the number of datasets
    % contribute to each entry
    
    if gen_perm_inv_flag
        % normalization by row
        max_G=max(Gabs,[],2);
        Gabs=repmat(1./max_G,1,size(Gabs,2)).*Gabs;
    end
    isiGrp=0;
    
    if commonSources
        for n=1:minN
            isiGrp=isiGrp+sum(Gabs(n,1:minN))/max(Gabs(n,1:minN))-1;
        end
        for m=1:minN
            isiGrp=isiGrp+sum(Gabs(1:minN,m))/max(Gabs(1:minN,m))-1;
        end
        isiGrp=isiGrp/(2*minN*(minN-1));
    else
        for n=1:Nk
            isiGrp=isiGrp+sum(Gabs(n,:))/max(Gabs(n,:))-1;
        end
        for m=1:Nk
            isiGrp=isiGrp+sum(Gabs(:,m))/max(Gabs(:,m))-1;
        end
        isiGrp=isiGrp/(2*L*(L-1));
    end
    
else
    % Have not handled when W is cell and A is single matrix or vice-versa.  Former makes
    % sense when you want performance of multiple algorithms for one mixing matrix, while
    % purpose of latter is unclear.
end

return


function [W, A, SR]  = getIVASig(W, X)

A = zeros(size(W));
SR = zeros(size(X));

for n = 1:size(W, 3)
    SR(:, :, n) = squeeze(W(:, :, n)*X(:, :, n));
    A(:, :, n) = pinv(squeeze(W(:, :, n)));
end