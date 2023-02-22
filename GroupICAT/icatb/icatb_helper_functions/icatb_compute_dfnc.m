function corrInfo = icatb_compute_dfnc(X, TR, varargin)
%% Compute dFNC using sliding window approach. Currently, gaussian window is implemented.
%
% Inputs:
% 1. X - Timecourses of dimensions time by components
% 2. TR - TR in seconds
% 3. varargin - Arguments passed in pairs
%   a. Y - If Y is empty, sliding window correlation is computed on X. Y is of
%       dimensions Time by voxels.
%   b. time_points - Specify all the subjects timepoints and also specify minTR
%   for accuarte interpolation.
%   c. minTR - Min TR across subjects
%   d. preprocess_params - Pre-process params.
%   preprocess_params.detrend_no = 3; % options are 0, 1, 2 and 3.
%   preprocess_params.doDespike = 'yes'; % Options are no and yes.
%   preprocess_params.tc_filter = 0.15; % Specify in HZ. you can also use
%   bandpass.
%   e. modelTC - Specify model timecourse for task-related dfnc.
%   f. windowing_params - Windowing params like wsize, window_alpha.
%   % windowing_params.wsize = 30; % specify window size
%   % windowing_params.window_alpha = 3; % specify alpha.
%
% Output:
% corrInfo.X - Pre-processed X
% corrInfo.Y - Pre-processed Y if any
% corrInfo.FNCdyn - Dynamic windowed correlations.
%

icatb_defaults;
global DETRENDNUMBER;
global DFNC_DEFAULTS;

if (numel(X) == length(X))
    X = X(:);
end


current_tp = size(X, 1);

covariateInfo = [];
modelTC = [];

Y = [];

%% Loop over varargin
for n = 1:2:length(varargin)
    if (strcmpi(varargin{n}, 'mintr'))
        minTR = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'minTP'))
        minTP = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'preprocess_params'))
        preprocess_params = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'windowing_params'))
        windowing_params = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'modelTC'))
        modelTC = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'covariateInfo'))
        covariateInfo = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'y'))
        Y = varargin{n + 1};
    end
end

% if (~exist('time_points', 'var') || isempty(time_points))
%     time_points = current_tp;
% end

if (~exist('minTR', 'var') || isempty(minTR))
    minTR = min(TR);
end

if (~isempty(Y))
    if (numel(Y) == length(Y))
        Y = Y(:);
    end
end

%% Initialize variables required
interpFactor = TR/minTR;
[num, denom] = rat(interpFactor);
%minTP = min(ceil(time_points*num /denom)); % TRuncate timepoints to least length of the subjects in
detrend_no = DETRENDNUMBER;
doDespike = 'no';
tc_filter = 0;
method = 'none';
wsize = 30;
window_alpha = 3;
window_type = 'gaussian';
num_L1_repetitions = 10;
initial_lambdas = (0.1:.03:.40);

% run task dfnc
task_dfnc = exist('modelTC', 'var') && ~isempty(modelTC);

% Pre-process params
if (exist('preprocess_params', 'var'))
    try
        detrend_no = preprocess_params.detrend_no;
    catch
    end
    
    try
        doDespike = preprocess_params.doDespike;
    catch
    end
    
    try
        tc_filter = preprocess_params.tc_filter;
    catch
    end
end

% Windowing params
if (exist('windowing_params', 'var'))
    try
        method = windowing_params.method;
    catch
    end
    try
        wsize = windowing_params.wsize;
    catch
    end
    try
        window_alpha = windowing_params.window_alpha;
    catch
    end
    
    try
        num_L1_repetitions = windowing_params.num_L1_repetitions;
    catch
    end
    
    try
        window_type = windowing_params.window_type;
    catch
    end
end

if (strcmpi(method, 'shared trajectory'))
    task_dfnc = 0;
end


%% PRE-PROCESS
if (~isempty(detrend_no))
    disp(['Detrending time series with value ', num2str(detrend_no), ' ...']);
    X = icatb_detrend(X, 1, [], detrend_no);
    if (~isempty(Y))
        Y = icatb_detrend(Y, 1, [], detrend_no);
    end
end

if (~isempty(covariateInfo))
    scansToInclude = [];
    try
        covariate_file = covariateInfo.file;
        disp('Variance associated with the covariates will be removed from the timecourses ... ');
        scansToInclude = covariateInfo.file_numbers;
    catch
    end
    
    X = regress_cov(X, covariate_file, scansToInclude);
    if (~isempty(Y))
        Y = regress_cov(Y, covariate_file, scansToInclude);
    end
end

% resample timeseries if needed
if (interpFactor > 1)
    X = resample(X, num, denom);
    if (~isempty(Y))
        Y = resample(Y, num, denom);
    end
    if (task_dfnc)
        modelTC = resample(modelTC, num, denom);
    end
    
end

if (strcmpi(doDespike, 'yes'))
    disp('Despiking timeseries ... ');
    X = icatb_despike_tc(X, minTR);
    if (~isempty(Y))
        Y = icatb_despike_tc(Y, minTR);
    end
end

% Filter timecourses
if (tc_filter > 0)
    disp('Filtering timeseries ... ');
    X = icatb_filt_data(X, minTR, tc_filter);
    if (~isempty(Y))
        Y = icatb_filt_data(Y, minTR, tc_filter);
    end
end


if (~exist('minTP', 'var'))
    minTP = size(X, 1);
end

X = X(1:minTP , :);
if (~isempty(Y))
    Y = Y(1:minTP, :);
end

if (task_dfnc)
    modelTC = modelTC(1:minTP, :);
end


%% Windowing related
window_step_size = 1;

try
    window_step_size = DFNC_DEFAULTS.step_size;
catch
end

if (~window_step_size)
    window_step_size = 1;
end

A = compute_sliding_window(minTP, window_alpha, wsize);

Nwin = minTP - wsize;

window_steps = (1:window_step_size:Nwin);

model_tcwin = zeros(Nwin, size(modelTC, 2));
XWin = cell(1, Nwin);
YWin = XWin;
for ii = 1:Nwin
    
    Ashift = circshift(A, round(-minTP/2) + round(wsize/2) + window_steps(ii));
    msk = icatb_window_mask(Ashift, window_steps(ii), window_alpha, wsize);
    
    try
        tcwin1 = bsxfun(@times, X(msk, :), Ashift(msk));
    catch
        tcwin1 = zeros(length(msk), size(X, 2));
        for nwin1 = 1:size(X, 2)
            tcwin1(:, nwin1) = X(msk, nwin1).* Ashift(msk);
        end
    end
    
    
    if (~isempty(Y))
        try
            tcwin2 = bsxfun(@times, Y(msk, :), Ashift(msk));
        catch
            tcwin2 = zeros(length(msk), size(Y, 2));
            for nwin2 = 1:size(Y, 2)
                tcwin2(:, nwin2) = Y(msk, nwin2).* Ashift(msk);
            end
        end
    end
    
    if (task_dfnc)
        % compute mean of the model TC within the window
        model_tcwin(ii, :) = mean(modelTC(msk, :).*Ashift(msk));
    end
    
    XWin{ii} = tcwin1;
    if (~isempty(Y))
        YWin{ii} = tcwin2;
    end
    
end

numComp1 = size(X, 2);
numComp2 = size(Y, 2);

if (isempty(Y))
    % Sliding window correlations on X
    
    FNCdyn = zeros(Nwin, numComp1*(numComp1 - 1)/2);
    
    if (strcmpi(method, 'none') || strcmpi(method, 'correlation'))
        % No L1
        for ii = 1:Nwin
            a = icatb_corr(squeeze(XWin{ii}));
            FNCdyn(ii, :) = icatb_mat2vec(a);
        end
        FNCdyn = icatb_r_to_z(FNCdyn);
        
    elseif (strcmpi(method, 'mutual information'))
        for ii = 1:Nwin
            a = icatb_compute_mi(squeeze(XWin{ii})');
            FNCdyn(ii, :) = icatb_mat2vec(a);
        end
    elseif (strcmpi(method, 'shared trajectory'))
        FNCdyn = calculate_shared_trajectory(X, wsize);
    else
        % glasso or L1 approach
        
        useMEX = 0;
        
        try
            GraphicalLassoPath([1, 0; 0, 1], 0.1);
            useMEX = 1;
        catch
        end
        
        disp('Using L1 regularisation ...');
        
        % L1 regularisation
        Pdyn = zeros(Nwin, numComp1*(numComp1 - 1)/2);
        
        Lambdas = zeros(num_L1_repetitions, length(initial_lambdas));
        
        fprintf('\t rep ')
        % Loop over no of repetitions
        for r = 1:num_L1_repetitions
            fprintf('%d, ', r)
            [trainTC, testTC] = split_timewindows(XWin, 1);
            trainTC = icatb_zscore(trainTC);
            testTC = icatb_zscore(testTC);
            [wList, thetaList] = computeGlasso(trainTC, initial_lambdas, useMEX);
            obs_cov = icatb_cov(testTC);
            L = cov_likelihood(obs_cov, thetaList);
            Lambdas(r, :) = L;
        end
        
        fprintf('\n')
        [mv, minIND] = min(Lambdas, [], 2);
        blambda = mean(initial_lambdas(minIND));
        fprintf('\tBest Lambda: %0.3f\n', blambda)
        best_lambda = blambda;
        
        % now actually compute the covariance matrix
        fprintf('\tWorking on estimating covariance matrix for each time window...\n')
        for ii = 1:Nwin
            [wList, thetaList] = computeGlasso(icatb_zscore(squeeze(XWin{ii})), blambda, useMEX);
            a = icatb_corrcov(wList);
            a = a - eye(numComp1);
            FNCdyn(ii, :) = icatb_mat2vec(a);
            InvC = -thetaList;
            r = (InvC ./ repmat(sqrt(abs(diag(InvC))), 1, numComp1)) ./ repmat(sqrt(abs(diag(InvC)))', numComp1, 1);
            r = r + eye(numComp1);
            Pdyn(ii, :) = icatb_mat2vec(r);
        end
        
        corrInfo.best_lambda = best_lambda;
        corrInfo.Pdyn = Pdyn;
        corrInfo.Lambdas = Lambdas;
        
        FNCdyn = icatb_r_to_z(FNCdyn);
        
    end
    
    
    if (task_dfnc)
        task_connectivity = zeros(size(model_tcwin, 2), size(FNCdyn, 2));
        for nRegress = 1:size(model_tcwin, 2)
            task_connectivity(nRegress, :) = icatb_r_to_z(icatb_corr(model_tcwin(:, nRegress), FNCdyn));
        end
        corrInfo.task_connectivity = task_connectivity;
    end
    
    
    
else
    % Correlation with X and Y
    
    FNCdyn = zeros(Nwin, size(XWin{ii}, 2), size(YWin{ii}, 2), class(Y));
    
    if (strcmpi(method, 'none') || strcmpi(method, 'correlation'))
        % No L1 approach
        for ii = 1:Nwin
            a = icatb_corr(squeeze(XWin{ii}), squeeze(YWin{ii}));
            FNCdyn(ii, :, :) = a;
        end
        FNCdyn = icatb_r_to_z(FNCdyn);
        
    elseif (strcmpi(method, 'mutual information'))
        % Mutual information
        for ii = 1:Nwin
            a = icatb_compute_mi(squeeze(XWin{ii})', squeeze(YWin{ii})');
            FNCdyn(ii, :, :) = a;
        end
    else
        % partial correlation approach
        for ii = 1:Nwin
            FNCdyn(ii, :, :) = icatb_partial_corr(squeeze(YWin{ii}), squeeze(XWin{ii}));
        end
        FNCdyn = icatb_r_to_z(FNCdyn);
    end
    
    
end

%% Output
corrInfo.X = X;
corrInfo.Y = Y;
corrInfo.FNCdyn = FNCdyn; % Windowed correlations


function c = compute_sliding_window(nT, win_alpha, wsize)
%% Compute sliding window
%

nT1 = nT;
if mod(nT, 2) ~= 0
    nT = nT + 1;
end

m = nT/2;
w = round(wsize/2);
gw = gaussianwindow(nT, m, win_alpha);
b = zeros(nT, 1);  b((m -w + 1):(m+w)) = 1;
c = conv(gw, b); c = c/max(c); c = c(m+1:end-m+1);
c = c(1:nT1);


function tc = regress_cov(tc, file_name, scansToInclude)
%% Regress covariates from timecourses
%

file_name = deblank(file_name);

X = icatb_load_ascii_or_mat(file_name);

if (~exist('scansToInclude', 'var') || isempty(scansToInclude))
    scansToInclude = (1:size(X, 1));
end

scansToInclude(scansToInclude > size(X, 1)) = [];

if (isempty(scansToInclude))
    error(['Please check file numbers specified for file ', file_name]);
end

X = icatb_zscore(X);

X = X(scansToInclude, :);

if (size(X, 1) ~= size(tc, 1))
    error(['Please check the timepoints in file ', file_name]);
end

betas = pinv(X)*tc;

% Remove variance associated with the covariates
tc = tc - X*betas;


function [wList, thetaList] = computeGlasso(tc, initial_lambdas, useMEX)
%% Compute graphical lasso


if (useMEX == 1)
    [wList, thetaList] = GraphicalLassoPath(tc, initial_lambdas);
else
    tol = 1e-4;
    maxIter = 1e4;
    S = icatb_cov(tc);
    thetaList = zeros(size(S, 1), size(S, 2), length(initial_lambdas));
    wList = thetaList;
    
    for nL = 1:size(wList, 3)
        [thetaList(:, :, nL), wList(:, :, nL)] = icatb_graphicalLasso(S, initial_lambdas(nL), maxIter, tol);
    end
    
end


function [trainTC, testTC] = split_timewindows(XWin, ntrain)
%% Split time windows into train and test for glasso approach

Nwin = length(XWin);

r = randperm(Nwin);

train_inds = r(1:ntrain);
test_inds = r(ntrain+1:end);

trainTC = XWin(train_inds);
testTC = XWin(test_inds);

trainTC = squeeze(cat(1, trainTC{:}));
testTC = squeeze(cat(1, testTC{:}));

function w = gaussianwindow(N,x0,sigma)
%% Gaussian window

x = 0:N-1;
w = exp(- ((x-x0).^2)/ (2 * sigma * sigma))';

function L = cov_likelihood(obs_cov, theta)
%% L = cov_likelihood(obs_cov, sigma)
% INPUT:
% obs_cov is the observed covariance matrix
% theta is the model precision matrix (inverse covariance matrix)
% theta can be [N x N x p], where p lambdas were used
% OUTPUT:
% L is the negative log-likelihood of observing the data under the model
% which we would like to minimize

nmodels = size(theta,3);

L = zeros(1,nmodels);
for ii = 1:nmodels
    % log likelihood
    theta_ii = squeeze(theta(:,:,ii));
    L(ii) = -log(det(theta_ii)) + trace(theta_ii*obs_cov);
end


function [Shared_traj] = calculate_shared_trajectory(tmp_dat,winsize)
%% Calculate shared trajectory
%
% Inputs:
% 1. tmp_dat - Time series T by component number
% 2. winsize - window size to calculate wap
%
% Outputs:
% wap - calulated wap for specified window size
%

diff_dat = diff(tmp_dat,1);
Shared_traj = cell(size(diff_dat,2));
for compn1=1:size(diff_dat,2)
    for compn2=(compn1+1):size(diff_dat,2)
        M_ang1 = atand(diff_dat(:,compn1)./diff_dat(:,compn2));
        M_ang2 = atand(diff_dat(:,compn2)./diff_dat(:,compn1));
        
        M_ang3 = [M_ang1 M_ang2];
        [~,min_id] = min(abs(M_ang3)');
        
        for n=1:length(M_ang1)
            tmp_val = [M_ang1(n) M_ang2(n)];
            [~,min_idx] = min(abs(tmp_val));
            diff_phase(n,1) = tmp_val(min_idx);
        end
        diff_mag = sqrt(diff_dat(:,compn1).^2+diff_dat(:,compn2).^2);
        %%
        if winsize>0
            tmp_wap = zeros(size(diff_dat,1)-winsize+1,1);
            for widx_s =1:(size(diff_dat,1)-winsize+1)
                widx = widx_s:(widx_s+winsize-1);
                tmp_wap(widx_s,1) = (diff_phase(widx)'*diff_mag(widx))/(sum(diff_mag(widx)));
            end
            Shared_traj{compn1,compn2} = tmp_wap;
        else
            Shared_traj{compn1,compn2} = diff_phase;
        end
    end
end
Shared_traj = cell2mat(m_mat2vec(Shared_traj')');
%end


function [vec, IND] = m_mat2vec(mat)
% [vec, IND] = mat2vec(mat)
% Returns the lower triangle of mat
% mat should be square [m x m], or if 3-dims should be [n x m x m]

if ndims(mat) == 2
    
    [n,m] = size(mat);
    if n ~=m
        error('mat must be square!')
    end
elseif ndims(mat) == 3
    
    [n,m,m2] = size(mat);
    if m ~= m2
        error('2nd and 3rd dimensions must be equal!')
    end
end

temp = ones(m);
%% find the indices of the lower triangle of the matrix
IND = find((temp-triu(temp))>0);
if ndims(mat) == 2
    vec = mat(IND);
else
    mat = reshape(mat, n, m*m2);
    vec = mat(:,IND);
end

%end

