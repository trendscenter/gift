function varargout = icatb_parCalculatePCA(files, numComp, varargin)
%% Principal component analysis.
% This code is based on:
%    I)    S. Rachakonda, R. F. Silva, J. Liu and V. D. Calhoun, "Memory Efficient PCA Methods for Large Group ICA", Frontiers in Neuroscience, 2016.
%    II)   S. Rachakonda, J. Liu and V. D. Calhoun. "Efficient Data Reduction in Group ICA Of fMRI Data", Proc. HBM, Seattle, WA. 2013.
%    III)  S. Roweis. "EM Algorithms for PCA and SPCA", in Neural Information Processing Systems (NIPS 1997)
%
% Available pca methods are eigen value decomposition (EVD), singular value decomposition (SVD), Expectation Maximization (empca), Multipower iteration (MPOWIT)
% and subsampled time pca (STP).
%
% Inputs:
% 1. files - Data could be specified as a 2D numeric array or file names
% must be passed in a cell array. Available formats are .nii, .img and
% .mat. When using .MAT extension it is advised to specify the variable to
% load using vartoload parameter.
% 2. numComp - Number of components.
% Optional parameters must be passed in pairs
%   a. type - Options available are standard or evd, svd, empca, mpowit, stp and
%   best. If best option is selected, best pca method is selected based on
%   the data size and maximum RAM available.
%   b. pca_options - PCA options specific to the method. Available options
%   must be specified in a data structure. Please see icatb_pca_options(pcaType) for more details.
%   c. mask - Specify binary mask for brain data. If you want to include all
%   the values, leave it as empty.
%   d. whiten - Whitening
%   e. preproc_type - Pre-process data. Available options are Remove Mean
%   Per Timepoint, Remove Mean Per Voxel, Intensity Normalization, Variance
%   Normalization and None.
%   f. verbose - Print info. Options are 0 and 1.
%   g. vartoload - Variable name in MAT file.
%   h. maxRAM - Available RAM in GB.
%
% Outputs:
% a. No whitening:
%   V - Eigen vectors of dimensions columns x numComp
%   Lambda - Diagonal matrix of eigen values sorted in ascending order.
% b. Whitening:
%   whitesig - PCA components of dimensions rows x numComp
%   dewhiteM - Dewhitening matrix of dimensions columns x numComp
%   Lambda - Diagonal matrix of eigen values sorted in ascending order.
%   V - Eigen vectors of dimensions columns x numComp
%   whiteM - Whitening matrix of dimensions numComp x rows
%
% Examples:
% a. Numeric 2D array:
%   data = randn(3000, 200);
%   numComp = 5;
%   [V, Lambda] = icatb_parCalculatePCA(data, numComp, 'whiten', false);
%   [whitesig, dewhiteM, Lambda, V, whiteM] = icatb_parCalculatePCA(data, numComp, 'whiten', true);
%
% b. Files in NII format
% files = {'C:\ns1.nii', 'C:\ns2.nii', 'C:\ns3.nii'};
% mask = 'myMask.nii';
% [V, Lambda] = icatb_parCalculatePCA(files, numComp, 'mask', mask, 'whiten', false);
%
% c. Files in MAT format
% files = {'C:\ns1.mat', 'C:\ns2.mat', 'C:\ns3.mat'};
% [V, Lambda] = icatb_parCalculatePCA(files, numComp, 'whiten', false, 'vartoload', 'pcasig');
%
% d. Best PCA:
% files = {'C:\ns1.nii', 'C:\ns2.nii', 'C:\ns3.nii'};
% mask = 'myMask.nii';
% [V, Lambda] = icatb_parCalculatePCA(files, numComp, 'mask', mask, 'whiten', false, 'type', 'best', 'maxram', 4);

%% Initialize variables
removeMean = 1;
pcaType = 'standard';
whitenMatrix = 1;
verbose = true;
mask = [];
maxRAM = 4;
preproc_type = 'none';
varToLoad = 'pcasig';
minDimSize = 10000;
chunks = 5000;

for i = 1:2:length(varargin)
    if (strcmpi(varargin{i}, 'type'))
        pcaType = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'pca_options')
        pca_options = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'mask'))
        mask = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'whiten')
        whitenMatrix = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'verbose')
        verbose = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'remove_mean')
        removeMean = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'preproc_type')
        preproc_type = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'vartoload')
        varToLoad = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'maxram')
        maxRAM = varargin{i + 1};
    end
end

if (~isnumeric(files))
    if (~iscell(files))
        files = {files};
    end
end

%% Load mask
if (~isempty(mask))
    if (ischar(mask))
        mask = icatb_rename_4d_file(mask);
        mask = icatb_loadData(deblank(mask(1, :)));
    end
    mask = find(mask ~= 0);
end

%% Get dimensions of the data
if (~exist('dims', 'var'))
    dims = getDims(files, mask, varToLoad);
end

if ~isnumeric(files)
    [pp, fN, extn] = fileparts(deblank(files{1}(1,:)));
end

df = min(dims);

rows = dims(1);
cols = dims(2);

if (numComp > min(dims))
    error(['Number of components (', num2str(numComp), ') selected is greater than the least dimension of the data (', num2str(df), ')']);
end


%% Number of subjects
numSubjects = 1;
if (~isnumeric(files))
    numSubjects = length(files);
end

precisionType = 'double';
try
    precisionType = pca_options.precision;
catch
end

numGroups = 10;
try
    numGroups = pca_options.numGroups;
catch
    
end
numGroups = min([numGroups, numSubjects]);

block_multiplier = 5;
try
    block_multiplier = pca_options.block_multiplier;
catch
end

interComp = min([df, 500]);
try
    interComp = pca_options.num_comp;
catch
end

tolerance = 1e-5;
try
    tolerance = pca_options.tolerance;
catch
end

if (~exist('pcaType', 'var') || strcmpi(pcaType, 'best'))
    %% Determine best pca method based on maximum RAM available
    
    mems = compute_mem_requirements(numSubjects, dims(1), dims(2), numComp, dims(2), precisionType, block_multiplier, numGroups);
    chkStackData = strcmpi('yes', mems(:, 3));
    chkMEM = [mems{:, 2}] < maxRAM;
    
    includeEVD = ones(length(mems(:, 1)), 1);
    chk = (strcmpi('standard', mems(:, 1)) | strcmpi('svd', mems(:, 1)));
    if (min(dims) <= minDimSize)
        includeEVD(chk) = 1;
    else
        includeEVD(chk) = 0;
    end
    
    includeEVD = (includeEVD == 1);
    
    selectedPCAOptions = mems(chkStackData(:) & chkMEM(:) & includeEVD(:), :);
    
    if (isempty(selectedPCAOptions))
        %% Use MPOWIT instead
        includeMPOWIT = strcmpi('mpowit', mems(:, 1));
        if (~isnumeric(files))
            chkStackData = strcmpi('no', mems(:, 3));
        else
            chkStackData = strcmpi('yes', mems(:, 3));
        end
        selectedPCAOptions = mems(chkStackData(:) & includeMPOWIT(:), :);
    end
    
    selectedPCAOptions = selectedPCAOptions(1, :);
    pcaType = selectedPCAOptions{1, 1};
    pca_options.precision = precisionType;
    pca_options.numGroups = numGroups;
    pca_options.block_multiplier = block_multiplier;
    pca_options.storage = selectedPCAOptions{1, 4};
    pca_options.precision = precisionType;
    pca_options.stack_data = selectedPCAOptions{1, 3};
    pca_options.num_comp = interComp;
    chkPCA.pcaType = pcaType;
    chkPCA.pca_opts = pca_options;
    chkPCA = icatb_check_pca_opts(chkPCA);
    pca_options = chkPCA.pca_opts;
end


% pca_options.numGroups = numGroups;
% pca_options.block_multiplier = block_multiplier;

if (isnumeric(files) && strcmpi(pcaType, 'stp'))
    %disp('Numeric data is passed. Using MPOWIT instead of STP ...');
    pcaType = 'evd';
    pca_options.stack_data = 'yes';
end


if (~exist('pca_options', 'var'))
    pca_options = icatb_pca_options(pcaType);
end

pca_options = icatb_pca_options(pcaType, pca_options, 'off');

stack_data = 'no';
try
    stack_data = pca_options.stack_data;
catch
end

if (strcmpi(pcaType, 'svd'))
    stack_data = 'yes';
end


%% Use mpowit instead of empca when pca is run in unstacked mode
if (strcmpi(stack_data, 'no'))
    if (strcmpi(pcaType, 'expectation maximization') || strcmpi(pcaType, 'empca'))
        pcaType = 'mpowit';
        pca_options = icatb_pca_options(pcaType, pca_options, 'off');
    end
end

if (isnumeric(files) || strcmpi(stack_data, 'yes'))
    %% Load all data
    
    if (isnumeric(files))
        data = files;
        clear files;
    else
        data = loadData(files, mask, preproc_type, precisionType, varToLoad);
    end
    
    if (removeMean)
        data = icatb_remove_mean(data, 0);
    end
    
    
    switch (lower(pcaType))
        
        case {'evd', 'standard'}
            %% Eigen value decomposition
            if (verbose)
                disp('Using Eigen Value Decomposition ...');
            end
            %eig_solver = pca_options.eig_solver;
            [V, Lambda] = compute_evd(data, numComp, chunks, verbose);
            
        case 'svd'
            %% Singular value decomposition
            [V, Lambda] = icatb_svd(data, numComp, 'solver', pca_options.solver, 'verbose', verbose);
            
        case 'mpowit'
            %% MPOWIT
            if (verbose)
                disp('Using Multi-power iteration ...');
            end
            Vi = randn(df, min([pca_options.block_multiplier*numComp, df]));
            tol = pca_options.tolerance;
            if strcmpi(class(data), 'single')
                if (tol < 1e-4)
                    tol = 1e-4;
                end
            end
            [V, Lambda, whitesig, final_iter] = multipowit(data, numComp, 'Vi', Vi, 'max_iter', pca_options.max_iter, 'tolerance', tol, 'verbose', verbose, ...
                'svd', 1);
        case {'empca', 'expectation maximization'}
            %% Expectation maximization
            Vi = randn(size(data, 2), numComp)';
            tol = pca_options.tolerance;
            if strcmpi(class(data), 'single')
                if (tol < 1e-4)
                    tol = 1e-4;
                end
            end
            [V, Lambda] = icatb_calculate_em_pca(data, Vi, 'tolerance', tol, 'max_iter', pca_options.max_iter, 'verbose', verbose);
    end
    
    checkEig(Lambda);
    
    if (whitenMatrix)
        [whiteM, dewhiteM] = get_pca_info(V, Lambda);
        if (~strcmpi(pcaType, 'mpowit'))
            whitesig = data*whiteM';
        end
        varargout{1} = whitesig;
        varargout{2} = dewhiteM;
        varargout{3} = Lambda;
        varargout{4} = V;
        varargout{5} = whiteM;
    else
        varargout{1} = V;
        varargout{2} = Lambda;
    end
    
    if (strcmpi(pcaType, 'mpowit'))
        varargout{end+1} = final_iter;
    end
    
else
    %% Load one data-set at a time
    
    switch (lower(pcaType))
        
        case {'evd', 'standard'}
            %% Eigen value decomposition
            if (verbose)
                disp('Using Eigen Value Decomposition ...');
            end
            eig_solver = pca_options.eig_solver;
            if (rows > cols)
                cov_m = computeCovInTimeUnStacked(files, mask, preproc_type, precisionType, varToLoad, dims, chunks, pca_options.storage);
                [V, Lambda] = icatb_eig_symm(cov_m, dims(2), 'num_eigs', numComp, 'eig_solver', eig_solver, 'verbose', verbose);
                checkEig(Lambda);
                if (whitenMatrix)
                    [whiteM, dewhiteM] = get_pca_info(V, Lambda);
                    whitesig = applyTransform(files, mask, preproc_type, precisionType, varToLoad, dims, whiteM');
                    varargout{1} = whitesig;
                    varargout{2} = dewhiteM;
                    varargout{3} = Lambda;
                    varargout{4} = V;
                    varargout{5} = whiteM;
                else
                    varargout{1} = V;
                    varargout{2} = Lambda;
                end
            else
                [whitesig, Lambda, V]  = computeEvdInVoxelsUnStacked(files, numComp, mask, preproc_type, precisionType, varToLoad, dims, chunks, pca_options.storage, numGroups);
                checkEig(Lambda);
                if (whitenMatrix)
                    [whiteM, dewhiteM] = get_pca_info(V, Lambda);
                    varargout{1} = whitesig;
                    varargout{2} = dewhiteM;
                    varargout{3} = Lambda;
                    varargout{4} = V;
                    varargout{5} = whiteM;
                else
                    varargout{1} = V;
                    varargout{2} = Lambda;
                end
            end
            
        case 'stp'
            %% STP
            if (verbose)
                disp('Using STP ...');
            end
            [whitesig, S0, Lambda, V] = parstp(files, numComp, dims, 'mask', mask, 'preproc_type', preproc_type, 'precision', precisionType, 'varToLoad', varToLoad, ...
                'numGroups', numGroups, 'interComp', interComp);
            checkEig(Lambda);
            if (whitenMatrix)
                [whiteM, dewhiteM] = get_pca_info(V, Lambda);
                varargout{1} = whitesig;
                varargout{2} = dewhiteM;
                varargout{3} = Lambda;
                varargout{4} = V;
                varargout{5} = whiteM;
            else
                varargout{1} = V;
                varargout{2} = Lambda;
            end
            
        case {'mpowit', 'empca', 'expectation maximization'}
            %% MPOWIT
            if (verbose)
                disp('Using Multi-power iteration ...');
            end
            [whitesig, Lambda, final_iter, V] = mpowit_pca2(files, numComp, dims, 'mask', mask, 'preproc_type', preproc_type, 'precision', precisionType, 'varToLoad', varToLoad, ...
                'numGroups', numGroups, 'block_multiplier', pca_options.block_multiplier, 'tolerance', tolerance);
            checkEig(Lambda);
            if (whitenMatrix)
                [whiteM, dewhiteM] = get_pca_info(V, Lambda);
                varargout{1} = whitesig;
                varargout{2} = dewhiteM;
                varargout{3} = Lambda;
                varargout{4} = V;
                varargout{5} = whiteM;
                varargout{6} = final_iter;
            else
                varargout{1} = V;
                varargout{2} = Lambda;
                varargout{3} = final_iter;
            end
    end
    
    
end



function [V, Lambda] = compute_evd(data, numComp, chunks, verbose)
%% EVD
%

[rows, cols] = size(data);

df = min(size(data));

if (~exist('verbose', 'var'))
    verbose = 1;
end

eig_solver = 'all';
if (df > 500)
    eig_solver = 'selective';
end

if (~exist('chunks', 'var'))
    chunks = df;
end

if (rows >= cols)
    cov_m = computeCovInTime(data, chunks);
    [V, Lambda] = icatb_eig_symm(cov_m, size(cov_m, 1), 'num_eigs', numComp, 'eig_solver', eig_solver, 'verbose', verbose);
else
    cov_m = computeCovInVoxels(data, chunks);
    [V, Lambda] = icatb_eig_symm(cov_m, size(cov_m, 1), 'num_eigs', numComp, 'eig_solver', eig_solver, 'verbose', verbose);
    V = (pinv(V)*data)';
    V = V* diag(1./sqrt(sum(V.^2)));
end


function varargout = parstp(files, numComp, dims, varargin)
%% Subsampled Time PCA based on three stage pca method
%
% Inputs:
% 1. files - File names in cell array. Length of cell array is equal to
% number of subjects
% 2. numComp - Number of components extracted from data
% 3. interComp - Number of intermediate components retained in each group.
% 4. numGroups - Number of subjects in each group
%
% Outputs:
% S - Principal Components
% S0 - Intermediate principal components
% V - Eigen vectors
% Lamda - Eigen values


preproc_type = 'none';
varToLoad = '';
precisionType = 'double';
interComp = numComp;
numGroups = 10;
mask = [];

for n = 1:2:length(varargin)
    if (strcmpi(varargin{n}, 'intercomp'))
        interComp = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'numgroups'))
        numGroups = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'mask'))
        mask = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'preproc_type'))
        preproc_type = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'precision'))
        precisionType = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'vartoload'))
        varToLoad = varargin{n + 1};
    end
end

if (~iscell(files))
    error('File names must be passed in cell array');
end

if (isempty(mask))
    mask = (1:dims(1));
end

NumWorkers = 4;

try
    clusterInfo = parcluster;
    NumWorkers = clusterInfo.NumWorkers;
catch
end

NumWorkers = min([NumWorkers, length(files)]);
subsPerGroup = ceil(length(files)/NumWorkers);

S = cell(1, NumWorkers);

parfor n = 1:NumWorkers
    tmpFiles = files;
    startIndex = (n - 1)*subsPerGroup + 1;
    endIndex = n*subsPerGroup;
    endIndex = min([endIndex, length(files)]);
    
    [dd, S{n}] = stp(tmpFiles(startIndex:endIndex), interComp, dims, 'numgroups', numGroups, 'mask', mask, 'preproc_type', preproc_type, 'precision', precisionType, 'vartoload', varToLoad);
    
end

S = [S{:}];
V = compute_evd(S, min([interComp, size(S, 2)]));
S = S*V;

S0 = S;
S = S(:, end-numComp+1:end);

%% Get eigen values and whiten pca components
Lambda = diag(diag(S'*S)/(size(S, 1) - 1));
S =  S*diag(1./std(S));

varargout{1} = S;
varargout{2} = S0;
varargout{3} = Lambda;


%% Compute eigen vectors
if (nargout == 4)
    V = cell(1, length(files));
    parfor nF = 1:length(files)
        dat = loadData(files(nF), mask, preproc_type, precisionType, varToLoad);
        V{nF} = (pinv(S)*dat)';
    end
    V = cat(1, V{:});
    varargout{4} = V;
end


function varargout = stp(files, numComp, dims, varargin)
%% Subsampled Time PCA based on three stage pca method
%
% Inputs:
% 1. files - File names in cell array. Length of cell array is equal to
% number of subjects
% 2. numComp - Number of components extracted from data
% 3. interComp - Number of intermediate components retained in each group.
% 4. numGroups - Number of subjects in each group
%
% Outputs:
% S - Principal Components
% S0 - Intermediate principal components
% V - Eigen vectors
% Lamda - Eigen values


preproc_type = 'none';
varToLoad = '';
precisionType = 'double';
interComp = numComp;
numGroups = 10;
mask = [];

for n = 1:2:length(varargin)
    if (strcmpi(varargin{n}, 'intercomp'))
        interComp = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'numgroups'))
        numGroups = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'mask'))
        mask = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'preproc_type'))
        preproc_type = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'precision'))
        precisionType = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'vartoload'))
        varToLoad = varargin{n + 1};
    end
end

if (~iscell(files))
    error('File names must be passed in cell array');
end

if (isempty(mask))
    mask = (1:dims(1));
end

numSubjects = length(files);
numGroups = min([numSubjects, numGroups]);


nLoops = ceil(length(files)/numGroups);

e = 0;

for n = 1:nLoops
    
    s = e + 1;
    e = e + numGroups;
    e = min([e, length(files)]);
    
    inds = (s:e);
    
    dat = cell(1, length(inds));
    for nS = 1:length(inds)
        pcasig = loadData(files(inds(nS)), mask, preproc_type, precisionType, varToLoad);
        dat{nS} = pcasig;
    end
    
    clear pcasig;
    dat = [dat{:}];
    
    if (n == 1)
        
        V = compute_evd(dat, min([interComp, size(dat, 2)]));
        S = dat*V;
        
    else
        
        V = compute_evd(dat, min([interComp, size(dat, 2)]));
        S2 = dat*V;
        
        dat = [S, S2];
        V = compute_evd(dat, min([interComp, size(dat, 2)]));
        S = dat*V;
        
    end
    
    
end

clear dat V S2;

if (numComp > size(S, 2))
    numComp = size(S, 2);
end

S0 = S;
S = S(:, end-numComp+1:end);

%% Get eigen values and whiten pca components
Lambda = diag(diag(S'*S)/(size(S, 1) - 1));
S =  S*diag(1./std(S));

varargout{1} = S;
varargout{2} = S0;
varargout{3} = Lambda;

%% Compute eigen vectors
if (nargout == 4)
    V = zeros(dims(2), size(S, 2));
    endT = 0;
    for nF = 1:length(files)
        dat = loadData(files(nF), mask, preproc_type, precisionType, varToLoad);
        startT = endT + 1;
        endT = endT + size(dat, 2);
        V(startT:endT, :) = (pinv(S)*dat)';
    end
    varargout{4} = V;
end


function varargout = multipowit(data, numComp, varargin)
%% Principal Component Analysis Using Multipower Iteration
%
% Inputs:
%
% 1. data - Data  of dimensions M x N
% 2. numComp - No. of principal components
% 3. varargin - Pass name value in pairs. Names are as follows
%   a. max_iter - Maximum no. of iterations
%   b. tolerance - Stopping criteria tolerance.
%   c. Vi - Initial eigen vectors of dimension K x numComp where K is the
%   smallest dimension of the data.
%
% Outputs:
%
% 1. V - Eigen vectors of dimension K x numComp where K is the smallest
% dimension of the data.
% 2. D - Eigen values in a diagonal matrix of dimensions numComp x numComp
% 3. U - Eigen maps of dimension J x numComp where J is the largest
% dimension of the data.
% 4. n - No. of iterations required to converge.
%


%% Initialise variables
max_iter = 1000;

verbose = 1;
isCov = 0;

%% Parse arguments
for ii = 1:2:length(varargin)
    if (strcmpi(varargin{ii}, 'max_iter'))
        max_iter = varargin{ii + 1};
    elseif (strcmpi(varargin{ii}, 'tolerance'))
        tol = varargin{ii + 1};
    elseif (strcmpi(varargin{ii}, 'Vi'))
        Vn = varargin{ii + 1};
    elseif (strcmpi(varargin{ii}, 'verbose'))
        verbose = varargin{ii + 1};
    elseif (strcmpi(varargin{ii}, 'svd'))
        useSvd = varargin{ii + 1};
    elseif (strcmpi(varargin{ii}, 'iscov'))
        isCov = varargin{ii + 1};
    end
end

%% Dimensions check
[num_rows, num_cols] = size(data);

useTranspose = 0;
if (num_rows > num_cols)
    useTranspose = 1;
end


minDims = min([num_rows, num_cols]);

if (numComp > minDims)
    error(['No. of components cannot be greater than degrees of freedom in the data (', num2str(minDims), ')']);
end


if (numComp >= ceil(0.6*minDims))
    %% Use EVD if number of components to be extracted is atleast 60% of minimum dimensions or degrees of freedom
    
    [V, D] = icatb_calculate_pca(data, numComp, 'type', 'evd', 'verbose', 0, 'whiten', 0);
    
    varargout{1} = V;
    varargout{2} = D;
    
    if (nargout > 2)
        U = data*V;
        U = U*diag(1./std(U));
        varargout{3} = U;
        varargout{4} = 1;
    end
    
    return;
    
end



if (~exist('Vn', 'var'))
    Vn = randn(minDims, min([5*numComp, minDims]));
    if (~isreal(data))
        Vn = complex (randn(minDims, min([5*numComp, minDims])), randn(minDims, min([5*numComp, minDims])));
    end
end

if ((size(Vn, 1) ~= num_rows) && (size(Vn, 1) ~= num_cols))
    error('Please check the dimensions of Vi as it should match no. of rows or columns');
end

if (size(Vn, 1) == num_rows)
    useTranspose = 0;
end

if (size(Vn, 2) > minDims)
    Vn = Vn(:, 1:minDims);
end


%% Use svd for orthogonalizing eigen vectors. SVD is preferred for complex valued data.
if (~exist('useSvd', 'var'))
    if (isreal(data))
        useSvd = 0;
    else
        useSvd = 1;
    end
end

if (~exist('max_iter', 'var'))
    max_iter = 1000;
end

Vn = ortho(Vn, useSvd);

D0 = zeros(numComp, 1);

isDataSingle = strcmpi(class(data), 'single');


%% Loop over iterations
for n = 1:max_iter
    
    if (~isCov)
        
        % Project eigen vectors on to the data
        if (useTranspose)
            EI = (data*Vn);
        else
            EI = (Vn'*data)';
        end
        
        if (isDataSingle)
            EI = double(EI);
        end
        
        D = eig(covMatrix(EI, num_rows), 'nobalance');
        D = sort(D);
        D = D(end-numComp+1:end);
        
        if (~exist('tol', 'var'))
            if (n == 1)
                tol = max(size(data))*eps(max(D));
            end
        end
        
        err = max(abs(D - D0));
        
        if (err < 10*tol)
            break;
        end
        
        
        if (useTranspose)
            Vn = (EI'*data)';
        else
            Vn = data*EI;
        end
        
        if (isDataSingle)
            Vn = double(Vn);
        end
        
    else
        
        Vn = data*Vn;
        if (isDataSingle)
            Vn = double(Vn);
        end
        D = eig(covMatrix(Vn, num_rows), 'nobalance');
        D = sort(D);
        D = D(end-numComp+1:end);
        
        if (~exist('tol', 'var'))
            if (n == 1)
                tol = max(size(data))*eps(max(D));
            end
        end
        
        err = max(abs(D - D0));
        
        if (err < 10*tol)
            Vn = ortho(Vn, useSvd);
            break;
        end
        
        
    end
    
    Vn = ortho(Vn, useSvd);
    
    if (verbose)
        disp(['Iteration: ', num2str(n), ' Error: ', num2str(err)]);
    end
    
    D0 = D;
    
end


if (~isCov)
    
    if (useTranspose)
        cov_m = data*Vn;
    else
        cov_m = (Vn'*data)';
    end
    
    
    if (~useSvd)
        [V, D] = eig(covMatrix(cov_m, num_rows), 'nobalance');
    else
        [dd, D, V] = svd(cov_m, 0);
        D = diag(D);
        D = D ./ sqrt(size(cov_m, 1) - 1);
        [D, inds] = sort(D.^2);
        V = V(:, inds);
        D = diag(D);
        clear dd;
    end
    
    Vn = Vn*V;
    
    V = Vn(:, end-numComp+1:end);
    
    V = norm2(V);
    
    D = D(end-numComp+1:end, end-numComp+1:end);
    
    
else
    
    [V, D] = eig(covMatrix(data*Vn, num_rows), 'nobalance');
    Vn = Vn*V;
    V = Vn(:, end-numComp+1:end);
    V = norm2(V);
    D = V'*data*V;
    
end

%if (nargout >= 3)
if (useTranspose)
    U = data*V;
else
    U = (V'*data)';
end
%end

U = norm2(U);

if (useTranspose)
    
    varargout{1} = V;
    varargout{2} = D;
    U = U*diag(1./std(U));
    varargout{3} = U;
    varargout{4} = n;
    
else
    
    varargout{1} = U;
    varargout{2} = D;
    V = V*diag(1./std(V));
    varargout{3} = V;
    varargout{4} = n;
    
end

if (verbose)
    fprintf('Done\n\n');
end



function varargout = mpowit_pca2(files, numComp, dims, varargin)
%% MPOWIT (un-stacked mode)
%

preproc_type = 'none';
varToLoad = '';
precisionType = 'double';
numGroups = 10;
mask = [];
max_iter = 20;
tol = 1e-6;
block_multiplier = 10;

for n = 1:2:length(varargin)
    if (strcmpi(varargin{n}, 'numgroups'))
        numGroups = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'mask'))
        mask = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'preproc_type'))
        preproc_type = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'precision'))
        precisionType = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'vartoload'))
        varToLoad = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'block_multiplier'))
        block_multiplier = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'tolerance'))
        tol = varargin{n + 1};
    end
end


if ((tol < 1e-4) && strcmpi(precisionType, 'single'))
    tol = 1e-4;
end

numGroups = min([numGroups, length(files)]);

interComp = min([dims, block_multiplier*numComp]);

% disp('Computing initial PCA subspace using STP to minimize number of data-loads...');
% [dd, S, D] = parstp(files, numComp, dims, 'mask', mask, 'preproc_type', preproc_type, 'precision', precisionType, 'varToLoad', varToLoad, ...
%     'numGroups', numGroups, 'interComp', interComp);
%
% Dold = diag(D);
%

try
    disp('Computing initial PCA subspace using STP to minimize number of data-loads...');
    [dd, S, D] = parstp(files, numComp, dims, 'mask', mask, 'preproc_type', preproc_type, 'precision', precisionType, 'varToLoad', varToLoad, ...
        'numGroups', numGroups, 'interComp', interComp);
    Dold = diag(D);
catch
    disp(['STP initialization failed due to ', lasterr]);
    disp('Using random initialization instead of STP initialization ...');
    S = ortho(randn(dims(1), interComp), 1);
    Dold = zeros(numComp, 1);
end

S = norm2(S);

for i = 1:max_iter
    
    if (i > 1)
        S = ortho(S, 1);
    end
    
    [S, cov_m] = doSum(files, numGroups, S, mask, preproc_type, precisionType, varToLoad);
    
    if (strcmpi(precisionType, 'single'))
        S = double(S);
        cov_m = double(cov_m);
    end
    
    [V, D] = eig(cov_m, 'nobalance');
    S = S*V;
    D = diag(D);
    D = D(end-numComp+1:end);
    err = norm(D - Dold);
    
    disp(['Iter: ', num2str(i), ' err: ', num2str(err)]);
    
    if (err < 10*tol)
        break;
    end
    
    Dold = D;
    
end

S = S(:, end-numComp+1:end);
D = diag(D);

varargout{1} = S;
varargout{2} = D;
varargout{3} = i;

%% Compute eigen vectors
if (nargout == 4)
    V = zeros(dims(2), size(S, 2));
    endT = 0;
    for nF = 1:length(files)
        dat = loadData(files(nF), mask, preproc_type, precisionType, varToLoad);
        startT = endT + 1;
        endT = endT + size(dat, 2);
        V(startT:endT, :) = (pinv(S)*dat)';
    end
    varargout{4} = V;
end


function dims = getDims(files, mask, varToLoad)
% Get dimensions of the matrix


if (~isnumeric(files))
    
    [dd, fn, extn] = fileparts(deblank(files{1}(1, :)));
    
    if (strcmpi(extn, '.mat'))
        tp = 0;
        for nF = 1:length(files)
            cn = files{nF};
            tmp = whos('-file', cn);
            if (~isempty(varToLoad))
                chk = find(strcmpi(varToLoad, cellstr(char(tmp.name))) == 1);
            else
                chk = 1;
            end
            if (isempty(chk))
                error(['Variable ', varToLoad, ' does not exist in file ', cn]);
            end
            chkdims = tmp(chk).size;
            if (nF == 1)
                voxels = chkdims(1);
            else
                if (voxels ~= chkdims(1))
                    error(['No of rows in file ', cn, ' do not match with the first file']);
                end
            end
            tp = tp + chkdims(2);
        end
        if (~isempty(mask))
            mask_ind = find(mask);
            dims = [length(mask_ind), tp];
        else
            dims = [voxels, tp];
        end
    else
        if (~exist('mask', 'var') || isempty(mask))
            error('Mask does not exist');
        end
        mask_ind = find(mask);
        tp = 0;
        for nF = 1:length(files)
            tmp = icatb_rename_4d_file(files{nF});
            tp = tp + size(tmp, 1);
        end
        dims = [length(mask_ind), tp];
    end
    
else
    dims = [size(files, 1), size(files, 2)];
end



function data = loadData(files, mask_ind, preProcType, precisionType, varToLoad)
%% Load data and preprocess if specified
%

data = cell(1, length(files));

verbose = 0;
if (length(files) > 1)
    verbose = 1;
end

for nF = 1:length(files)
    if (verbose)
        disp(['Loading data-set ', num2str(nF), ' ...']);
    end
    tmp = icatb_read_data(files{nF}, [], mask_ind, precisionType, varToLoad);
    tmp = squeeze(tmp);
    size_d = size(tmp);
    tmp = reshape(tmp, prod(size_d(1:end-1)), size_d(end));
    if (~strcmpi(preProcType, 'none'))
        %Call pre-processing function
        tmp = icatb_preproc_data(tmp, preProcType);
        if (~strcmpi(preProcType, 'remove mean per timepoint'))
            %Remove mean per timepoint
            tmp = icatb_remove_mean(tmp, 1);
        end
    end
    data{nF} = tmp;
    
end

data = [data{:}];

function cov_m = computeCovInVoxels(data, chunks)
%% Compute covariance in voxels dimension
%

[voxels, tp] = size(data);

df = max([1, voxels - 1]);

if (icatb_get_matlab_version < 2013)
    
    cov_m = zeros(voxels, voxels);
    
    nLoops = ceil(voxels/chunks);
    
    endTp = 0;
    for nL = 1:nLoops
        startTp = endTp + 1;
        endTp = endTp + chunks;
        endTp = min([endTp, voxels]);
        cov_m(:, startTp:endTp) = data*data(startTp:endTp, :)';
    end
    
else
    cov_m = data*data';
end

cov_m = cov_m./df;


function cov_m = computeCovInTime(data, chunks)
%% Compute covariance in time dimension
%

[voxels, tp] = size(data);

df = max([1, voxels - 1]);

if (icatb_get_matlab_version < 2013)
    
    cov_m = zeros(tp, tp);
    
    nLoops = ceil(tp/chunks);
    
    endTp = 0;
    for nL = 1:nLoops
        startTp = endTp + 1;
        endTp = endTp + chunks;
        endTp = min([endTp, tp]);
        cov_m(startTp:endTp, :) = data(:, startTp:endTp)'*data;
    end
    
else
    cov_m = data'*data;
end

cov_m = cov_m./df;



function cov_m = computeCovPacked(cov_m, A, num_blocks, df)

[M, N] = size(A);

if (~exist('df', 'var'))
    df = N - 1;
end

%num_blocks = 1000;
%[M, N] = size(A);

nLoops = ceil(N/num_blocks);
endB = 0;
inds2 = 0;
B = A';
for nL = 1:nLoops
    startB = endB + 1;
    endB = endB + num_blocks;
    endB = min([endB, N]);
    tmp = B*A(:, startB:endB) / df;
    tmp = tmp(startB:end, :);
    inds = find(tril(ones(size(tmp))) ~= 0);
    inds2 = inds2 + (1:length(inds));
    cov_m(inds2) = cov_m(inds2) + tmp(inds);
    inds2 = max(inds2);
    clear tmp;
end


function cov_m = computeCovFull(cov_m, A, num_blocks, df)

[M, N] = size(A);

if (~exist('df', 'var'))
    df = M - 1;
end

nLoops = ceil(N/num_blocks);
endB = 0;
B = A';

for nL = 1:nLoops
    startB = endB + 1;
    endB = endB + num_blocks;
    endB = min([endB, N]);
    tmp = B*A(:, startB:endB) / df;
    cov_m(:, startB:endB) = cov_m(:, startB:endB) + tmp;
    clear tmp;
end


function cov_m = computeCovInTimeUnStacked(files, mask_ind, preProcType, precisionType, varToLoad, dims, chunks, storage)
%% Compute covariance in time dimension

nLoops = ceil(dims(1)/chunks);
numPairs = length(files)*(length(files) - 1)/2;

if ((nLoops + 1) < numPairs)
    
    if (isempty(mask_ind))
        mask_ind = (1:dims(1));
    end
    
    meanData = zeros(1, dims(2));
    if (~strcmpi(preProcType, 'none'))
        endT = 0;
        for nF = 1:length(files)
            tmp = icatb_read_data(files{nF}, [], mask_ind, precisionType, varToLoad);
            tmp = squeeze(tmp);
            size_d = size(tmp);
            tmp = reshape(tmp, prod(size_d(1:end-1)), size_d(end));
            startT = endT + 1;
            endT = endT + size(tmp, 2);
            if (~strcmpi(preProcType, 'none'))
                %Call pre-processing function
                tmp = icatb_preproc_data(tmp, preProcType);
                meanData(:, startT:endT) = mean(tmp);
            end
        end
    end
    
    files = char(files);
    endLoop = 0;
    if (strcmpi(storage, 'full'))
        cov_m = zeros(dims(2), dims(2), precisionType);
    else
        cov_m = zeros(dims(2)*(dims(2) + 1)/2, 1, 'single');
    end
    for nF = 1:nLoops
        
        startLoop = endLoop + 1;
        endLoop = endLoop + chunks;
        endLoop = min([dims(1), endLoop]);
        
        tmp = icatb_read_data(files, [], mask_ind(startLoop:endLoop), precisionType, varToLoad);
        tmp = reshape(tmp, size(tmp, 1), size(tmp, 2)*size(tmp, 3));
        if (~strcmpi(preProcType, 'none'))
            %Call pre-processing function
            tmp = icatb_preproc_data(tmp, preProcType);
        end
        
        tmp = rmBaseLine(tmp, meanData);
        
        if (strcmpi(storage, 'full'))
            cov_m = computeCovFull(cov_m, tmp, chunks, (dims(1) - 1));
        else
            cov_m = computeCovPacked(cov_m, tmp, chunks, (dims(1) - 1));
        end
        
    end
    
else
    
    if (strcmpi(storage, 'full'))
        
        cov_m = zeros(dims(2), dims(2), precisionType);
        
        currentRowEndIndex = 0;
        for i = 1:length(files)
            data1 = loadData(files(i), mask_ind, preProcType, precisionType, varToLoad);
            currentRowStartIndex = currentRowEndIndex + 1;
            currentRowEndIndex = currentRowEndIndex + size(data1, 2);
            rows = (currentRowStartIndex:currentRowEndIndex);
            currentColEndIndex = 0;
            for j = 1:length(files)
                data2 = loadData(files(j), mask_ind, preProcType, precisionType, varToLoad);
                currentColStartIndex = currentColEndIndex + 1;
                currentColEndIndex = currentColEndIndex + size(data2, 2);
                cols = (currentColStartIndex:currentColEndIndex);
                cov_m(rows, cols) = data1'*data2/(dims(1) - 1);
                cov_m(cols, rows) = data2'*data1/(dims(1) - 1);
            end
        end
        
    else
        
        cov_m = zeros(dims(2)*(dims(2) + 1)/2, 1, 'single');
        currentLength = dims(2);
        numElem = currentLength:-1:1;
        
        %% Diagonal elements
        diags = [1, numElem(1:end-1)];
        diags = cumsum(diags);
        countEndDiag = 0;
        for i = 1:length(files)
            data1 = loadData(files(i), mask_ind, preProcType, precisionType, varToLoad);
            PCs = (size(data1, 2):-1:1);
            countStartDiag = countEndDiag + 1;
            countEndDiag = countEndDiag + size(data1, 2);
            tempElem = diags(countStartDiag:countEndDiag);
            
            countSize = 0;
            for j = i:length(files)
                data2 = loadData(files(j), mask_ind, preProcType, precisionType, varToLoad);
                tmpCov = data2'*data1/(dims(1) - 1);
                if (i ~= j)
                    temp_numpcs = repmat(PCs, size(data2, 2), 1);
                    temp_inds = repmat((1:size(data2, 2))', 1, size(data1, 2));
                    temp_diags = repmat(tempElem, size(data2, 2), 1);
                    inds = temp_inds + temp_diags + temp_numpcs - 1 + countSize;
                    countSize = countSize + size(data2, 2);
                else
                    tril_inds = find(tril(ones(size(tmpCov))) ~= 0);
                    tmpCov = tmpCov(tril_inds);
                    endTempInd = 0;
                    inds = zeros(1, length(tmpCov));
                    for nE = 1:length(PCs)
                        chkRows = tempElem(nE):tempElem(nE)+PCs(nE)-1;
                        startTempInd = endTempInd + 1;
                        endTempInd = endTempInd + length(chkRows);
                        inds(startTempInd:endTempInd) = chkRows;
                    end
                end
                cov_m(inds) = tmpCov;
                
            end
        end
        
    end
    
end


function [S, Lambda, V] = computeEvdInVoxelsUnStacked(files, numComp, mask, preproc_type, precisionType, varToLoad, dims, chunks, pca_storage, numSubjectsInGroup)
%% Compute covariance in voxel dimension
%

voxels = dims(1);

if (strcmpi(pca_storage, 'full'))
    cov_m = zeros(voxels, voxels);
else
    cov_m = zeros(voxels*(voxels + 1)/2, 1, 'single');
end

%% Reduce number of loops by grouping subjects
numSubjectsInGroup = min([numSubjectsInGroup, length(files)]);
numGroups = ceil(length(files)/numSubjectsInGroup);

e = 0;
df = voxels - 1;
for nG = 1:numGroups
    s = e + 1;
    e = e + numSubjectsInGroup;
    e = min([e, length(files)]);
    
    data = loadData(files(s:e), mask, preproc_type, precisionType, varToLoad);
    
    if (strcmpi(pca_storage, 'full'))
        cov_m = cov_m + ((data*data')/df);
    else
        cov_m = computeCovPacked(cov_m, data', chunks, df);
    end
    
end

drawnow;

%% Compute evd in voxel dimension
[S, Lambda] = icatb_eig_symm(cov_m, voxels, 'num_eigs', numComp, 'eig_solver', 'selective');

S = S*diag(1./std(S));


e = 0;
V = zeros(dims(2), numComp);
endTp = 0;
for nG = 1:numGroups
    
    s = e + 1;
    e = e + numSubjectsInGroup;
    e = min([e, length(files)]);
    data = loadData(files(s:e), mask, preproc_type, precisionType, varToLoad);
    startTp = endTp + 1;
    endTp = endTp + size(data, 2);
    V(startTp:endTp, :) = (pinv(S)*data)';
    
end

V = norm2(V);


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

function [whiteningMatrix, dewhiteningMatrix] = get_pca_info(V, Lambda)
%% Get Whitening and de-whitening matrix
%
% Inputs:
% 1. V - Eigen vectors
% 2. Lambda - Eigen values diagonal matrix
%
% Outputs:
% 1. whiteningMatrix - Whitening matrix
% 2. dewhiteningMatrix - Dewhitening matrix
%


whiteningMatrix = sqrtm(Lambda) \ V';
dewhiteningMatrix = V * sqrtm(Lambda);



function X = rmBaseLine(X, orig_mean)
%% Remove baseline

try
    if (~isempty(which('bsxfun.m')))
        %% Use bsxfun if possible
        X = bsxfun(@minus, X, orig_mean);
    else
        tempVar = repmat(orig_mean, size(X, 1), 1);
        X = X - tempVar;
        clear tempVar;
    end
catch
    clear tempVar;
    %% Remove mean of sample one at a time
    for nM = 1:size(X, 2)
        X(:, nM) = X(:, nM) - orig_mean(nM);
    end
end



function whitesig = applyTransform(files, mask, preProcType, precisionType, varToLoad, dims, dewhiteM)
%% Apply transformation matrix
%

whitesig = zeros(dims(1), size(dewhiteM, 2), precisionType);

endTp = 0;
for nF = 1:length(files)
    tmp = icatb_read_data(files{nF}, [], mask, precisionType, varToLoad);
    tmp = squeeze(tmp);
    tmp = reshape(tmp, size(tmp, 1), size(tmp, 2)*size(tmp, 3));
    startTp = endTp + 1;
    endTp = endTp + size(tmp, 2);
    if (~strcmpi(preProcType, 'none'))
        %Call pre-processing function
        tmp = icatb_preproc_data(tmp, preProcType);
        if (~strcmpi(preProcType, 'remove mean per timepoint'))
            %Remove mean per timepoint
            tmp = icatb_remove_mean(tmp, 1);
        end
    end
    whitesig = whitesig + tmp*dewhiteM(startTp:endTp, :);
end


function data = ortho(data, useSvd)
%% Orthogonalise data
%

if (~useSvd)
    cov_m = covMatrix(data);
    [V1, D1] = eig(cov_m);
    V1 = V1(:, end:-1:1);
    data = data*V1;
    data = norm2(data);
else
    [data, dd] = svd(data, 0);
end


function cov_m = covMatrix(data, l)
%% Covariance matrix
%

if (~exist('l', 'var'))
    l = size(data, 1);
end

l = max([1, l - 1]);
cov_m = (data'*data)/l;


function data = norm2(data)
%% Norm2
%

if (isreal(data))
    data =  data*diag(1./sqrt(sum(data.^2)));
else
    data =  data*diag(1./sqrt(sum(abs(data).^2)));
end


function [S, cov_m] = doSum(files, numSubjectsInGroup, Sn, mask, preproc_type, precisionType, varToLoad)

numGroups = ceil(length(files)/numSubjectsInGroup);

S = zeros(size(Sn));
cov_m = zeros(size(Sn, 2), size(Sn, 2));

parfor n = 1:numGroups
    
    startT = (n - 1)*numSubjectsInGroup + 1;
    endT = n*numSubjectsInGroup;
    endT = min([length(files), endT]);
    inds = (startT:endT);
    
    count = 0;
    
    data = cell(1, length(inds));
    for nD = inds
        count = count + 1;
        pcasig = loadData(files(nD), mask, preproc_type, precisionType, varToLoad);
        data{count} = pcasig;
    end
    
    data = [data{:}];
    b = (Sn'*data)';
    S = S + data*b;
    
    cov_m = cov_m + b'*b/(size(data, 1) - 1);
    
end


function checkEig(Lambda)
%% Check eigen values
%

if (numel(Lambda) ~= length(Lambda))
    Lambda = diag(Lambda);
end

L = length(find(Lambda <= eps(class(Lambda))));

if (L > 1)
    error([num2str(L), ' eigen values are less than or equal to machine precision']);
elseif (L == 1)
    error('One of the eigen values is less than or equal to machine precision');
end