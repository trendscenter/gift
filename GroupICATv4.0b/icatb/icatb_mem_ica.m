function varargout = icatb_mem_ica(sesInfo)
%% Calculate memory required for computing ICA/IVA.
%

icaAlgo = 'infomax';

precisionType = 'double';

%% Input parameters
voxels = 70000;
% Time points
time_points = 300;
% Number of subjects
numOfSub = 602;
% Number of sessions
numOfSess = 1;
% Data reduction steps
numDataReductionSteps = 2;
% Data reduction step numbers
numOfPC1 = 100;
numOfPC2 = 30;

%% End for input parameters



dispInfo = 1;
if (nargout > 0)
    dispInfo = 0;
end

%% Get info from session info
if (exist('sesInfo', 'var'))
    
    
    try
        voxels = length(sesInfo.userInput.mask_ind);
    catch
        dat = icatb_loadData(deblank(sesInfo.userInput.files(1).name(1, :)));
        voxels = size(dat, 1)*size(dat, 2)*size(dat, 3);
        clear dat;
    end
    
    modalityType = icatb_get_modality;
    
    if (~isfield(sesInfo.userInput, 'diffTimePoints'))
        if strcmpi(modalityType, 'fmri')
            % get the count for time points
            sesInfo.userInput.diffTimePoints = icatb_get_countTimePoints(sesInfo.userInput.files);
        else
            sesInfo.userInput.diffTimePoints = icatb_get_num_electrodes(sesInfo.userInput.files);
        end
    end
    
    time_points = max(sesInfo.userInput.diffTimePoints);
    numOfSub = sesInfo.userInput.numOfSub;
    numOfSess = sesInfo.userInput.numOfSess;
    numDataReductionSteps = sesInfo.userInput.numReductionSteps;
    if (numDataReductionSteps > 2)
        numDataReductionSteps = 2;
    end
    numOfPC1 = sesInfo.userInput.numOfPC1;
    numOfPC2 = sesInfo.userInput.numOfPC2;
    numOfPC3 = 0;
    %    numOfPC3 = sesInfo.userInput.numOfPC3;
    icaAlgo = sesInfo.userInput.algorithm;
    if (~ischar(icaAlgo))
        ica_types = cellstr(icatb_icaAlgorithm);
        icaAlgo = ica_types{icaAlgo};
    end
    
    try
        precisionType = sesInfo.userInput.pca_opts.precision;
    catch
        try
            precisionType = sesInfo.userInput.covariance_opts.precision;
        catch
        end
    end
    
end

if (numOfSub*numOfSess == 1)
    numDataReductionSteps = 1;
end

if (numDataReductionSteps >= 3)
    numDataReductionSteps = 2;
end
% End for checking

if (strcmpi(icaAlgo, 'iva-gl'))
    numDataReductionSteps = 1;
end


if (dispInfo)
    
    disp(['Number of voxels: ', num2str(voxels)]);
    disp(['Number of timepoints: ', num2str(time_points)]);
    disp(['Number of subjects: ', num2str(numOfSub)]);
    disp(['Number of sessions: ', num2str(numOfSess)]);
    disp(['Number of data reduction steps: ', num2str(numDataReductionSteps)]);
    
end

% Matlab version
matlab_version = icatb_get_matlab_version;
OSBIT = mexext;

precisionType = icatb_checkPrecision(precisionType, dispInfo);

minNumBytesVar = 8;
if strcmpi(precisionType, 'single')
    minNumBytesVar = 4;
end

if (dispInfo)
    
    disp(['Number of PC1: ', num2str(numOfPC1)]);
    if (numDataReductionSteps == 2)
        disp(['Number of PC2: ', num2str(numOfPC2)]);
    end
    
end

if (strcmpi(icaAlgo, 'moo-icar'))
    icaAlgo = 'gig-ica';
end

if strcmpi(icaAlgo, 'iva-gl')
    
    mems = compute_mem_requirements(1, voxels, time_points, time_points, numOfPC1, precisionType);
    chk = strcmpi(mems(:, 3), 'yes');
    mems = mems(chk, :);
    memory_GB = (voxels*numOfPC1*numOfSub*numOfSess*minNumBytesVar)/1024/1024/1024;
    
    pca_opts = repmat(struct('pca_type', 'standard', 'stack_data', 'yes', 'storage', 'full', 'max_mem', memory_GB), 1, size(mems, 1));
    for nOpts = 1:length(pca_opts)
        pca_opts(nOpts).pca_type = mems{nOpts, 1};
        pca_opts(nOpts).max_mem = max([memory_GB, mems{nOpts, 2}]);
        pca_opts(nOpts).precision = precisionType;
    end
    
elseif (strcmpi(icaAlgo, 'constrained ica (spatial)') || strcmpi(icaAlgo, 'gig-ica'))
    
    mems = compute_mem_requirements(1, voxels, time_points, time_points, numOfPC1, precisionType);
    chk = strcmpi(mems(:, 3), 'yes');
    mems = mems(chk, :);
    
    pca_opts = repmat(struct('pca_type', 'standard', 'stack_data', 'yes', 'storage', 'full', 'max_mem', []), 1, size(mems, 1));
    for nOpts = 1:length(pca_opts)
        pca_opts(nOpts).pca_type = mems{nOpts, 1};
        pca_opts(nOpts).max_mem = mems{nOpts, 2};
        pca_opts(nOpts).precision = precisionType;
    end
    
else
    
    numGroups = 10;
    try
        numGroups = sesInfo.userInput.pca_opts.numGroups;
    catch
    end
    numGroups = min([numGroups, numOfSub*numOfSess]);
    
    if (numDataReductionSteps == 1)
        p = time_points;
    else
        p = numOfPC1;
    end
    
    block_multiplier = 5;
    try
        block_multiplier = sesInfo.userInput.pca_opts.block_multiplier;
    catch
    end
    
    mems = compute_mem_requirements(numOfSub*numOfSess, voxels, time_points, numOfPC1, p, precisionType, block_multiplier, numGroups);
    
    pca_opts = repmat(struct('pca_type', 'standard', 'stack_data', 'yes', 'storage', 'full', 'max_mem', []), 1, size(mems, 1));
    for nOpts = 1:length(pca_opts)
        pca_opts(nOpts).pca_type = mems{nOpts, 1};
        pca_opts(nOpts).max_mem = mems{nOpts, 2};
        pca_opts(nOpts).stack_data = mems{nOpts, 3};
        pca_opts(nOpts).storage = mems{nOpts, 4};
        pca_opts(nOpts).precision = precisionType;
    end
    
end



if (dispInfo)
    fprintf('\n');
    for nOpts = 1:length(pca_opts)
        if (strcmpi(pca_opts(nOpts).pca_type, 'evd') || strcmpi(pca_opts(nOpts).pca_type, 'standard'))
            disp(['PCA Type: ', pca_opts(nOpts).pca_type, ', Stack datasets: ', pca_opts(nOpts).stack_data, ', Precision: ', precisionType, ' Covariance Storage: ', pca_opts(nOpts).storage]);
        else
            disp(['PCA Type: ', pca_opts(nOpts).pca_type, ', Stack datasets: ', pca_opts(nOpts).stack_data, ', Precision: ', precisionType]);
        end
        disp(['Approximate memory required for the analysis is ', num2str(pca_opts(nOpts).max_mem), ' GB']);
        fprintf('\n');
    end
end


if (nargout > 0)
    varargout{1} = pca_opts;
end



function mems = compute_mem_requirements(subjects, v, t, k, p, precision, l, g)
% subjects - no of subjects
% v - voxels
% t - Time points
% k - Number of components
% p - number of components retained in first pca step
% precision - double or single


if (~exist('p', 'var'))
    p = t;
end

if (~exist('precision', 'var'))
    precision = 'double';
end

if (strcmpi(precision, 'double'))
    numBytes = 8;
else
    numBytes = 4;
end

if (~exist('l', 'var'))
    l = 5;
end

if (~exist('g', 'var'))
    g = 10;
end

M = subjects;

%% Temporal concatenation
mem = ((v*M*p) + min([v^2, (M*p)^2])+ max([v,M*p])*k)*numBytes;

mem = mem/1024/1024/1024;


%% SVD
mems = {'SVD', mem, 'yes'};


%% EVD
mems(end+1, 1:4) = {'Standard', mem, 'yes', 'full'};


%% GIFT (Full covariance storage)
mem = (min([v^2, (M*p)^2])+ max([v,M*p])*k)*numBytes;

mem = mem/1024/1024/1024;
mems(end+1, 1:4) = {'Standard', mem, 'no', 'full'};


%% GIFT (Packed covariance storage)
mem = (min([v^2, (M*p)^2])+ max([v,M*p])*k)*numBytes;
mem = mem/1024/1024/1024;
mems(end+1, 1:4) = {'Standard', mem/2, 'no', 'packed'};


%% MPOWIT (stacked)
Ma = 2*(min([v, M*p])*l*k) + (l*k)^2;
Mb = min([v, M*p]*l*k)  + 3*(l*k)^2;
Mc = [];
mem = (v*M*p) + max([Ma, Mb, Mc]);
mem = mem*numBytes;

mem = mem/1024/1024/1024;
mems(end+1, 1:3) = {'MPOWIT', mem, 'yes'};

%% MPOWIT (un-stacked)

Ma = 2*(min([v, M*p])*l*k) + (l*k)^2;
Mb = min([v, M*p]*l*k)  + 3*(l*k)^2;
Mc = [];
mem = (v*p) + max([Ma, Mb, Mc]);
mem = mem*numBytes;

mem = mem/1024/1024/1024;
mems(end+1, 1:3) = {'MPOWIT', mem, 'no'};


%% STP
g1 = min([g, subjects]);
if (g1 ~= subjects)
    mem = (v*g1*p + (g1*p)^2 + 2*v*l*k + (g1*p*l*k) + 2*(l*k)^2);
    mem = mem*numBytes;
    mem = mem/1024/1024/1024;
else
    mem = mems{1, 2};
end

mems(end+1, 1:3) = {'STP', mem, 'no'};


%% EM PCA (stacked)
l=1;
Ma = 2*(min([v, M*p])*l*k) + (l*k)^2;
Mb = min([v, M*p]*l*k)  + 3*(l*k)^2;
Mc = [];
mem = (v*M*p) + max([Ma, Mb, Mc]);
mem = mem*numBytes;

mem = mem/1024/1024/1024;
mems(end+1, 1:3) = {'Expectation Maximization', mem, 'yes'};