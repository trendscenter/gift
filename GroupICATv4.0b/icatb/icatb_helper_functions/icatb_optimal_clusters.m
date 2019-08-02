function out = icatb_optimal_clusters(data, num_clusters, varargin)
%% Get optimal number of clusters using elbow method, gap statistic, silhouette, BIC, AIC, dunns.
%
% Inputs:
% 1. data - Observations of dimensions N x M
% 2. num_clusters - Number of clusters to evaluate
% Arguments passed in pairs
%   a. method - Options are gap, silhoutte, bic, aic and dunns
%   b. num_tests - Number of times data is generated when using gap method
%
% Outputs:
% 1. results.K - Optimal number of clusters
% 2. results.values - Test values
% 3. results.sem - Standard error of mean (gap stat only)
%

useGUI = 0;
if (~exist('data', 'var'))
    useGUI = 1;
end


%% Parse inputs
Nt = 10;
method = 'gap';
def_opts = {'Replicates', 10, 'Distance', 'sqEuclidean', 'MaxIter', 200, 'EmptyAction', 'singleton'};
displayResults = 0;
standardize_values = 0;
%num_clusters = 10;


if (useGUI)
    estParams = getParamsGUI;
    data = estParams.data;
    method = estParams.method;
    num_clusters = estParams.num_clusters;
    cluster_opts = {'Replicates', estParams.Replicates, 'Distance', estParams.Distance, 'MaxIter', estParams.MaxIter};
    Nt = estParams.num_tests;
    displayResults = 1;
end

for nV = 1:length(varargin)
    if (strcmpi(varargin{nV}, 'num_tests'))
        Nt = varargin{nV + 1};
    elseif (strcmpi(varargin{nV}, 'method'))
        method = varargin{nV + 1};
    elseif (strcmpi(varargin{nV}, 'cluster_opts'))
        cluster_opts = varargin{nV + 1};
    elseif (strcmpi(varargin{nV}, 'display'))
        displayResults = varargin{nV + 1};
    elseif (strcmpi(varargin{nV}, 'standardize'))
        standardize_values = varargin{nV + 1};
    end
end

opts = def_opts;

if (exist('cluster_opts', 'var'))
    def_names = lower(def_opts(1:2:end));
    field_names = lower(cluster_opts(1:2:end));
    for n = 1:length(field_names)
        chk = strmatch(field_names{n}, def_names, 'exact');
        if (~isempty(chk))
            opts{2*chk} = cluster_opts{2*n};
        end
    end
end

% Run in parallel if parpool exists
run_parallel = 0;
try
    pool = gcp;
    if (~isempty(pool))
        run_parallel = 1;
    end
catch
end

if (standardize_values)
    data = detrend(data, 0);
    data = data*diag(1./std(data));
end


if (length(num_clusters) == 1)
    KList = (1:num_clusters);
else
    KList = num_clusters;
end

idx = cell(1, length(KList));
C = idx;
sumD = idx;
if (~run_parallel)
    for j = 1:length(KList)
        [idx{j}, C{j}, sumD{j}, Dist{j}] = doKmeans(data, KList(j), opts);
    end
else
    parfor j = 1:length(KList)
        [idx{j}, C{j}, sumD{j}, Dist{j}] = doKmeans(data, KList(j), opts);
    end
end

if (strcmpi(method, 'all'))
    methodStr = {'elbow', 'gap', 'bic', 'aic', 'dunns', 'silhouette'};
else
    methodStr = cellstr(method);
end

clear method;

out = cell(1, length(methodStr));

for numMethod = 1:length(methodStr)
    
    method = methodStr{numMethod};
    if (strcmpi(method, 'elbow'))
        % Elbow method
        klist = KList;
        SSE = zeros(1, length(klist));
        R = zeros(1,length(klist));
        for k2 = 1:length(klist)
            [dd, R(k2)] = cluster_goodness(Dist{k2}, idx{k2});
            tmpSumD = sumD{k2};
            SSE(k2) = sum(tmpSumD(:));
        end
        
        % exclude values with only one cluster
        R(klist == 1) = [];
        klist(klist == 1) = [];
        
        [bestx, yfit] = fit_L_tocurve_area(klist,R);
        results.title = 'Elbow';
        results.klist = klist;
        results.values = R;
        results.fit = yfit;
        results.K = bestx;
        
    elseif (strcmpi(method, 'gap'))
        % Gap stat
        [K, values, sem] = computeGap(data, KList, Nt, run_parallel, opts, sumD);
        results.K = K;
        results.values = values;
        results.sem = sem;
        results.klist = KList;
        results.title = 'Gap';
    elseif (strcmpi(method, 'bic') || strcmpi(method, 'aic'))
        % BIC or AIC criteria
        [BIC, AIC] = compute_aic_bic(data, KList, idx, sumD);
        if (strcmpi(method, 'bic'))
            values = BIC;
            results.title = 'BIC';
        else
            values = AIC;
            results.title = 'AIC';
        end
        [dd, indsa]  = min(values, [], 2);
        results.K = KList(indsa(:)');
        results.values = values;
        results.klist = KList;
    elseif (strcmpi(method, 'dunns'))
        % Dunns
        klist = KList;
        idx2 = idx;
        idx2(klist == 1) = [];
        klist(klist == 1) = [];
        values = dunns_index(data, idx2);
        [dd, K] = max(values);
        K = klist(K);
        results.K = K;
        results.values = values;
        results.klist = klist;
        results.title = 'Dunns Index';
        
    else
        % Silhouette
        klist = KList;
        idx2 = idx;
        idx2(klist == 1) = [];
        klist(klist == 1) = [];
        [K, values] = computeSilh(data, klist, run_parallel, idx2);
        results.K = K;
        results.values = values;
        results.klist = klist;
        results.title = 'Silhouette';
    end
    
    out{numMethod} = results;
    
end

if (displayResults)
    
    num_rows = ceil(sqrt(length(out)));
    num_cols = ceil(length(out)/num_rows);
    
    sz = get(0, 'screensize');
    pos = [50, 50, sz(3) - 100, sz(4) - 100];
    figH = figure('name', 'Cluster Estimation', 'color', [1, 1, 1], 'position', pos);
    
    countN = 0;
    for nR = 1:num_rows
        for nC = 1:num_cols
            countN = countN + 1;
            if (countN > length(out))
                break;
            end
            sh = subplot(num_rows, num_cols, countN);
            res = out{countN};
            if (strcmpi(res.title, 'elbow'))
                plot(res.klist, res.values, 'parent', sh);
                hold on;
                plot(res.klist, res.fit, 'k', 'parent', sh);
                legend('Ratio of dispersion (Wthn/Bwn Groups)', 'Fit');
            elseif (strcmpi(res.title, 'gap'))
                try
                    errorbar(res.klist, res.values, res.sem, '-bo', 'linewidth', 1.5);
                catch
                end
            elseif (strcmpi(res.title, 'bic') || strcmpi(res.title, 'aic'))
                plotyy(res.klist, res.values(1, :), res.klist, res.values(2, :), 'parent', sh);
                legend('Loglikelihood', 'Euclidean');
            else
                plot(res.klist, res.values, '-bo', 'linewidth', 1.5, 'parent', sh);
            end
            axis(sh, 'tight');
            xlabel('Clusters', 'parent', sh);
            ylabel(res.title, 'parent', sh);
            if (strcmpi(res.title, 'bic') || strcmpi(res.title, 'aic'))
                title(['Optimal clusters using loglikehood = ', num2str(res.K(1)), ' & euclidean = ', num2str(res.K(2))], 'parent', sh);
            else
                title(['Optimal clusters  = ', num2str(res.K)], 'parent', sh);
            end
        end
    end
    
    lineH = findobj(figH, 'type', 'line');
    set(lineH, 'linewidth', 1.5);
    
end

if (nargout == 0)
    assignin('base', 'est_cluster_results', out);
end


%if (length(out) == 1)
%    out = out{1};
%end


function [optimal_clusters, gapvalues, sem] = computeGap(data, KList, Nt, run_parallel, opts, sumD)
%% Compute gap stat
%

%% Get expected value of logW
[pcaX, V] = doPCA(data);

rlogW = zeros(Nt, length(KList));
for n = 1:Nt
    X = getdata(pcaX, V);
    if (~run_parallel)
        for k = 1:length(KList)
            [dd, pp, sumd] = doKmeans(X, KList(k), opts);
            rlogW(n, k) = log(sum(sumd(:)));
        end
    else
        parfor k = 1:length(KList)
            [dd, pp, sumd] = doKmeans(X, KList(k), opts);
            rlogW(n, k) = log(sum(sumd(:)));
        end
    end
end

ElogW = mean(rlogW);
sem = std(rlogW, 1, 1)*sqrt(1 + (1/Nt));


%% Gap statistic and determine no. of optimal clusters
gapvalues = zeros(1, length(KList));
for j = 1:length(KList)
    sumd = sumD{j};
    %[dd, pp, sumd] = doKmeans(data, j, opts);
    logW = log(sum(sumd(:)));
    gapvalues(j) = ElogW(j) - logW;
end

chk = diff(gapvalues)./sem(2:end);
inds = find(chk <= 1);
optimal_clusters = min([length(KList), inds(:)']);



function [pcaX, V] = doPCA(data)
%% PCA
%

data = bsxfun(@minus, data, mean(data));
[U, S, V] = svd(data, 0);
pcaX = data*V;


function X = getdata(pcaX, V)
%% Generate uniform distribution
%
mn = min(pcaX);
mx = max(pcaX);
r = (mx - mn);
[rows, cols] = size(pcaX);

X = bsxfun(@plus, mn, bsxfun(@times, r, rand(rows, cols)));
X = X*V';


function  [dd, pp, sumd, D] = doKmeans(X, k, opts)
%% Kmeans
%

try
    [dd, pp, sumd, D] = kmeans(X, k, opts{:});
catch
    
    if (k == 1)
        dd = ones(size(X, 1), 1);
        pp = mean(X);
        data = X;
        X = bsxfun(@minus, X, pp);
        sumd = sum(sum(X.^2));
        D = data - repmat(pp, size(X, 1), 1);
        D = sum(D.^2, 2);
    else
        [dd, pp, sumd, D] = icatb_kmeans(X, k, opts{:});
    end
    
end


function [optimal_clusters, values] = computeSilh(data, klist, run_parallel, IDX)
%% Compute Silhouette
%

values = zeros(1, length(klist));

if (~run_parallel)
    
    for n = 1:length(klist)
        %[idx, cx, sumd] = doKmeans(data, n, opts);
        idx = IDX{n};
        silh = getsilh(data, idx, run_parallel);
        silh(isfinite(silh) == 0) = [];
        values(n) = mean(silh);
    end
    
else
    
    parfor n = 1:length(klist)
        idx = IDX{n};
        %[idx, cx, sumd] = doKmeans(data, n, opts);
        silh = getsilh(data, idx, run_parallel);
        silh(isfinite(silh) == 0) = [];
        values(n) = mean(silh);
    end
    
end

[dd, optimal_clusters] = max(values);
optimal_clusters = klist(optimal_clusters);


function silh = getsilh(X, idx, run_parallel)
%% Silhouette values using squared euclidean

k = length(unique(idx));
n = size(X, 1);
mbrs = (repmat(1:k, n, 1) == repmat(idx, 1, k));
count = histc(idx(:)', 1:k);
avgDWithin = repmat(NaN, n, 1);
avgDBetween = repmat(NaN, n, k);

if (~run_parallel)
    for j = 1:n
        mbrs_tmp = mbrs;
        count_tmp = count;
        distj = sum(bsxfun(@minus, X, X(j,:)).^2, 2);
        for i = 1:k
            if i == idx(j)
                avgDWithin(j) = sum(distj(mbrs_tmp(:,i))) ./ max(count_tmp(i)-1, 1);
            else
                avgDBetween(j, i) = sum(distj(mbrs_tmp(:,i))) ./ count_tmp(i);
            end
        end
    end
else
    parfor j = 1:n
        mbrs_tmp = mbrs;
        count_tmp = count;
        distj = sum(bsxfun(@minus, X, X(j,:)).^2, 2);
        for i = 1:k
            if i == idx(j)
                avgDWithin(j) = sum(distj(mbrs_tmp(:,i))) ./ max(count_tmp(i)-1, 1);
            else
                avgDBetween(j, i) = sum(distj(mbrs_tmp(:,i))) ./ count_tmp(i);
            end
        end
    end
end

% Calculate the silhouette values
minavgDBetween = min(avgDBetween, [], 2);
silh = (minavgDBetween - avgDWithin) ./ max(avgDWithin, minavgDBetween);


function [BIC, AIC] = compute_aic_bic(X, KList, idx, sumD)
%% Compute BIC and AIC (Information from two cluster stats and stats.stackexchange.com)
%
%

if (numel(X) == length(X))
    X = X(:);
end

%
% if (~exist('K', 'var'))
%     K = min([min(size(X)), 10]);
% end

BIC = zeros(2, length(KList));
AIC = zeros(2, length(KList));

for k = 1:length(KList)
    [BIC(:, k), AIC(:, k)] = aic_bic(X, KList(k), idx{k}, sumD{k});
end

function [BIC, AIC] = aic_bic(X, K, idx, sumd)
%% AIC and BIC
%

N = size(X, 1);
P = size(X, 2);

%idx = doKmeans(X, K, opts);


Nc = zeros(1, K);
for n = 1:K
    chk = (idx == n);
    Nc(n) = length(find(chk==1));
    x = X(chk, :);
    if (length(chk) > 1)
        tmp = var(x, 0, 1);
    else
        tmp = var(x, 1, 1);
    end
    
    if (n == 1)
        Vc = zeros(length(tmp), K);
    end
    
    Vc(:, n) = tmp';
end

V = var(X, 0, 1);
V = repmat(V', 1, K);


LL = -Nc.* sum( log(Vc + V)/2);
D1 = -2*sum(LL);
D2 = sum(sumd(:));
mj = K*(2*P);
BIC(1) = D1 + mj * log(N);
AIC(1) = D1 + 2*mj;

BIC(2) = D2 + mj * log(N);
AIC(2) = D2 + 2*mj;


function di = dunns_index(data, IDX)
%% Compute dunns index
%

distances = icatb_pdist(data, 'all', 'euclidean');
distances = icatb_vec2mat(distances);
%K = length(IDX);
di = zeros(1, length(IDX));
for ii = 1:length(IDX)
    idx = IDX{ii};
    num_clusters = length(unique(idx));
    num = [];
    denom = [];
    for n = 1:num_clusters
        x = (idx == n);
        y = (idx ~= n);
        tmp = distances(x, y);
        num = min([num; tmp(:)]);
        tmp = distances(x, x);
        denom = max([denom; tmp(:)]);
    end
    
    di(ii) = num/denom;
    
end

%klist = 2:K;


function [bestx, yfit] = fit_L_tocurve_area(x,y)
%[bestx, F] = fit_L_tocurve_area(x,y, plotit)
%
% if nargin < 3
%     plotit = 1;
% end
x = x(:);
y = y(:);
%%
% P = [xbp m1 m2 b]
options = optimset('TolFun',1e-14, 'TolX', 1e-14, 'MaxFunEvals', 100000, 'MaxIter', 10000);
P0(1) = x(5);
P0(2) = (y(4)-y(2))/(x(4)-x(2));
P0(3) = (y(end)-y(end-1))/(x(end)-x(end-1));
P0(4) = y(5);

LB = [x(2)    -Inf -Inf  min(y)];
UB = [x(end-1) 0     0   max(y)];

[PF,RESNORM,RESIDUAL,EXITFLAG,OUTPUT] = icatb_lsqcurvefit(@Lcurve,P0,x,y, LB, UB, options);

bestx = ceil(PF(1));

yfit = Lcurve(PF, x);


function Z = Lcurve(P, XDATA)
% P = [xbp m1 m2 b]

XDATA = XDATA-P(1);

% curve 1
Y1 = P(2)*XDATA + P(4);
% curve 2
Y2 = P(3)*XDATA + P(4);

Z = zeros(size(Y1));
Z(XDATA < 0) = Y1(XDATA < 0);
Z(XDATA >= 0) = Y2(XDATA >= 0);


function [SSE, R] = cluster_goodness(D, IDX)
% Get SSE and dispersion ratio

[nobs, k] = size(D);
Din = zeros(1, k);
Dout = zeros(1, k);

for ii = 1:k
    Din(ii) = sum(D(IDX == ii,ii).^2); % dispersion in cluster
    Dout(ii) = sum(D(IDX ~= ii,ii).^2); % sum of squared distances
end
SSE = sum(Din);
R = mean(Din./Dout);



function handles = getParamsGUI
%% Get params
%

param_file = icatb_selectEntry('title', 'Select dFNC Parameter File/2D Matrix File', 'typeEntity', 'file', 'typeSelection', 'single', 'filter', '*dfnc.mat;*.mat');

drawnow;

if (isempty(param_file))
    error('Data is not selected');
end

methodsStr = {'Elbow', 'Gap', 'BIC', 'AIC', 'Dunns', 'Silhouette'};

index = icatb_listdlg('PromptString', 'Select method', 'SelectionMode','multiple', 'ListString', methodsStr, 'movegui', 'center', 'windowStyle', 'modal', ...
    'title_fig', 'Cluster estimate method');
if (isempty(index))
    error('Method is not selected');
end

methodsStr = methodsStr(index);
handles.method = methodsStr;

dlg_title = 'Select cluster options';

numParameters = 1;

inputText(numParameters).promptString = 'Enter number of clusters to be evaluated';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = '10';
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'num_clusters';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;

numParameters = numParameters + 1;

inputText(numParameters).promptString = 'Enter maximum number of iterations';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = '200';
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'MaxIter';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;

distance_opts =  {'City', 'sqEuclidean', 'Hamming', 'Correlation', 'Cosine'};
chk = strmatch('sqeuclidean', lower(distance_opts), 'exact');
if (isempty(chk))
    chk = 1;
end


numParameters = numParameters + 1;

inputText(numParameters).promptString = 'Select distance method';
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerString = distance_opts;
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'Distance';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = chk;

numParameters = numParameters + 1;

inputText(numParameters).promptString = 'Number of times to repeat the clustering';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = '10';
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'Replicates';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;

numParameters = numParameters + 1;

inputText(numParameters).promptString = 'Number of reference data-sets for computing gap';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = '20';
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'num_tests';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;


answers = icatb_inputDialog('inputtext', inputText, 'Title', dlg_title, 'handle_visibility', 'on');

if (~isempty(answers))
    
    handles.num_clusters = answers{1};
    handles.MaxIter = answers{2};
    handles.Distance = answers{3};
    handles.Replicates = answers{4};
    handles.num_tests = answers{5};
else
    error('input window was quit');
end


isDFNC = 0;
try
    d = whos('-file', param_file);
    chk = strmatch('dfncInfo', cellstr(char(d.name)),'exact');
    if (~isempty(chk))
        isDFNC = 1;
    end
catch
end

disp('Loading data ...');

if (isDFNC)
    outputDir = fileparts(param_file);
    if (isempty(outputDir))
        outputDir = pwd;
    end
    load(param_file);
    % Load data
    for nR = 1:length(dfncInfo.outputFiles)
        
        current_file = fullfile(outputDir,  dfncInfo.outputFiles{nR});
        load(current_file, 'FNCdyn');
        if (nR == 1)
            data = zeros(length(dfncInfo.outputFiles), size(FNCdyn, 1), size(FNCdyn, 2));
        end
        
        data(nR, :, :) = FNCdyn;
        
    end
    data = reshape(data, size(data, 1)*size(data, 2), size(data, 3));
else
    data = icatb_load_ascii_or_mat(param_file);
end

handles.data = data;


