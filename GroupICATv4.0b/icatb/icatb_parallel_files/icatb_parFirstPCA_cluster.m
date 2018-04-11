function icatb_parFirstPCA_cluster(sesInfo, dataSetsToRun, conserve_disk_space)
%% Run PCA in parallel mode
%

if (ischar(sesInfo))
    load(sesInfo);
end

j = 1;
data_reduction_mat_file = sesInfo.data_reduction_mat_file;
outputDir = sesInfo.outputDir;
numOfPC = sesInfo.reduction(j).numOfPCAfterReduction;
numOfSess = sesInfo.numOfSess;
modalityType = icatb_get_modality;
files = sesInfo.inputFiles;
mask_ind = sesInfo.mask_ind;
preproc_type = sesInfo.preproc_type;
precision = 'double';
max_iter = 1000;
tol = 1e-6;

group_pca_type = 'subject specific';
if isfield(sesInfo, 'group_pca_type')
    group_pca_type = sesInfo.group_pca_type;
end

doGrandMeanPCA = (strcmpi(group_pca_type, 'grand mean') && (sesInfo.numOfSub*sesInfo.numOfSess > 1));

LambdaV = [];
VV = [];
VM = [];

if (doGrandMeanPCA)
    pp = load(fullfile(sesInfo.outputDir, [sesInfo.data_reduction_mat_file, 'mean.mat']));
    VM = pp.VM;
end

%% Group PCA Options
gpca_opts = getPCAOpts(sesInfo);
pcaType = gpca_opts{1}.pcaType;
pca_opts = gpca_opts{1}.pca_opts;

try
    precision = pca_opts.precision;
catch
end

try
    tol = pca_opts.tolerance;
catch
end

try
    max_iter = pca_opts.max_iter;
catch
end


if (strcmpi(precision, 'single') && (tol < 1e-4))
    
    try
        pca_opts.tolerance = 1e-4;
    catch
    end
end


%% Loop over groups after concatenation
parfor nDataSet = 1:length(dataSetsToRun)
    
    tmpDataSetsToRun = dataSetsToRun;
    i = tmpDataSetsToRun(nDataSet);
    tmpFiles = files;
    
    pcasig = []; dewhiteM = []; whiteM = [];
    
    pcaout = [data_reduction_mat_file, num2str(j), '-', num2str(i), '.mat'];
    pcaout = fullfile(outputDir, pcaout);
    
    
    fprintf('\n');
    
    subNum = ceil(i/numOfSess);
    ses = mod(i-1, numOfSess) + 1;
    if (~strcmpi(modalityType, 'smri'))
        if (~doGrandMeanPCA)
            msg_string = ['--Doing pca on Subject #', num2str(subNum), ' Session #', num2str(ses)];
        else
            msg_string = ['-- Projecting Subject #', num2str(subNum), ' Session #', num2str(ses), ' to the eigen space of the mean of data-sets'];
        end
    else
        msg_string = 'Doing pca';
    end
    
    disp(msg_string);
    
    %% Load data
    data = preprocData(tmpFiles(i).name, mask_ind, preproc_type, precision);
    
    %% Project data on to the eigen space of the mean data
    if (doGrandMeanPCA)
        data = data*VM;
    end
    
    if (conserve_disk_space ~= 1)
        varsToSave = {'pcasig', 'dewhiteM', 'whiteM', 'Lambda', 'V'};
    else
        varsToSave = {'V', 'Lambda'};
    end
    
    if (~conserve_disk_space)
        [pcasig, dewhiteM, Lambda, V, whiteM] = icatb_calculate_pca(data, numOfPC, 'type', pcaType, 'remove_mean', 0, 'whiten', 1, 'pca_options', pca_opts);
    else
        [V, Lambda] = icatb_calculate_pca(data, numOfPC, 'type', pcaType, 'remove_mean', 0, 'whiten', 0, 'pca_options', pca_opts);
    end
    
    %% Apply the transformations to Eigen vectors, dewhitening
    % and whitening matrices
    if (doGrandMeanPCA)
        V = VM*V;
        if (exist('dewhiteM', 'var'))
            dewhiteM = VM*dewhiteM;
            whiteM = pinv(dewhiteM);
        end
    end
    
    %% Save PCA vars
    drawnow;
    
    if (conserve_disk_space == 1)
        Lambda = diag(Lambda);
        Lambda = Lambda(:)';
        LambdaV = [LambdaV;Lambda];
        VV = [VV;V];
    end
    
    if (conserve_disk_space ~= 1)
        icatb_parSave(pcaout, {pcasig, dewhiteM, whiteM, Lambda, V}, varsToSave);
    end
    
    
end

if (conserve_disk_space == 1)
    pcaout = [data_reduction_mat_file, num2str(j), '-1.mat'];
    pcaout = fullfile(outputDir, pcaout);
    icatb_parSave(pcaout, {VV, LambdaV}, {'V', 'Lambda'});
end


function data = preprocData(fileN, mask_ind, preProcType, precisionType)

% Load data
data = icatb_read_data(fileN, [], mask_ind, precisionType);

% Call pre-processing function
data = icatb_preproc_data(data, preProcType);

if (~strcmpi(preProcType, 'remove mean per timepoint'))
    % Remove mean per timepoint
    data = icatb_remove_mean(data, 1);
end


function pca_opts = getPCAOpts(sesInfo)
%% Get PCA opts
%

isOptsCell = 0;
try
    isOptsCell = iscell(sesInfo.pca_opts);
catch
end

if (~isOptsCell)
    
    if (isfield(sesInfo, 'pca_opts'))
        tmp_pca_opts = sesInfo.pca_opts;
    end
    
    %% Group PCA Options
    pcaType = 'standard';
    if (isfield(sesInfo, 'pcaType'))
        pcaType = sesInfo.pcaType;
    end
    
    sesInfo.pcaType = pcaType;
    sesInfo = icatb_check_pca_opts(sesInfo);
    pcaType = sesInfo.pcaType;
    tmp_pca_opts = sesInfo.pca_opts;
    
    pca_opts{1} = struct('pcaType', pcaType, 'pca_opts', tmp_pca_opts);
    
else
    
    pca_opts = sesInfo.pca_opts;
    for nP = 1:length(pca_opts)
        pca_opts{nP} = icatb_check_pca_opts(pca_opts{nP});
    end
    
end

if (length(pca_opts) ~= sesInfo.numReductionSteps)
    pca_opts = repmat(pca_opts(1), 1, sesInfo.numReductionSteps);
end
