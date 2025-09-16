function meta_states_info = icatb_dfnc_meta_state_analysis(FNCdynflat, num_comps, varargin)
%% Meta state analysis
%
% Inputs:
% 1. FNCdynflat - Matrix of dimensions subjects x windows x component pairs
% 2. num_comps - Number of components/clusters
% 3. varargin - Arguments passed in pairs
% a. 'method' - Meta state method used. Options are 'k-means', 'pca',
% 'tica', 'sica' and 'all'. If you select 'all', all the meta state methods
% are used.
% b. 'dmethod' - Distance method used in k-means
% c. 'kmeans_max_iter' - Maximum no. of iterations used in kmeans
% d. 'num_ica_runs' - Number of times ICA is run when using temporal or
% spatial ICA.
% e. 'ica_algorithm' - ICA algorithm used in temporal or spatial ICA.
%
% Output:
% meta_states_info - Data structure containing meta state info
%

%% Initialise vars
num_ica_runs = 1;
kmeans_max_iter = 150;
dmethod = 'City';
ica_algorithm = 'infomax';
method = 'all';
kmeans_num_replicates = 5;

for n = 1:2:length(varargin)
    if (strcmpi(varargin{n}, 'dmethod'))
        dmethod = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'kmeans_max_iter'))
        kmeans_max_iter = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'num_ica_runs'))
        num_ica_runs = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'ica_algorithm'))
        ica_algorithm = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'method'))
        method = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'replicates'))
        kmeans_num_replicates = varargin{n + 1};
    end
end

M = size(FNCdynflat, 1);
Nwin = size(FNCdynflat, 2);
num_pairs = size(FNCdynflat, 3);
FNCdynflat = reshape(FNCdynflat, M*Nwin,  num_pairs);


% subjects x windows x component pairs
FNCdynflat_nomean = FNCdynflat;
FNCdynflat = icatb_remove_mean(FNCdynflat);

%% PCA
if (~strcmpi(method, 'k-means'))
    [S_pca, DD] = icatb_calculate_pca(FNCdynflat, num_comps, 'type', 'evd', 'whiten', 0, 'verbose', 1, 'preproc_type', 'none', 'remove_mean', 0);
    TC_pca = FNCdynflat*S_pca;
    if (strcmpi(method, 'pca') || strcmpi(method, 'all'))
        meta_states_info.pca.S = S_pca;
        % reshape it to Windows x components x subjects for quantile discretization
        tmp = reshape(TC_pca, M, Nwin, num_comps); tmp = permute(tmp, [2, 3, 1]);
        meta_states_info.pca.A = TC_pca;
        meta_states_info.pca.Tcs = Qrtl(tmp);
        meta_states_info.pca.summary = getSummaryMetaStates(meta_states_info.pca.Tcs);
        clear tmp;
    end
end

%% Temporal ICA
if (strcmpi(method, 'tica') || strcmpi(method, 'all'))
    TC_pca = TC_pca*diag(1./std(TC_pca)); % whiten the matrix
    [Wt, At, TC_tica] = run_stability_analysis(ica_algorithm, TC_pca', num_ica_runs);
    TC_tica = TC_tica';
    S_tica = (pinv(TC_tica)*FNCdynflat)';
    meta_states_info.tica.A = TC_tica;
    meta_states_info.tica.S = S_tica;
    tmp = reshape(TC_tica, M, Nwin, num_comps); tmp = permute(tmp, [2, 3, 1]);
    % signed quantile discretization
    meta_states_info.tica.Tcs = Qrtl(tmp);
    meta_states_info.tica.summary = getSummaryMetaStates(meta_states_info.tica.Tcs);
    clear TC_tica S_tica tmp;
end

%% Spatial ICA
if (strcmpi(method, 'sica') || strcmpi(method, 'all'))
    S_pca = S_pca*diag(1./std(S_pca)); % Whiten the matrix
    [Ws, As, S_sica] = run_stability_analysis(ica_algorithm, S_pca', num_ica_runs);
    S_sica = S_sica';
    TC_sica = (pinv(S_sica)*FNCdynflat')';
    meta_states_info.sica.A = TC_sica;
    meta_states_info.sica.S = S_sica;
    tmp = reshape(TC_sica, M, Nwin, num_comps); tmp = permute(tmp, [2, 3, 1]);
    % signed quantile discretization
    meta_states_info.sica.Tcs = Qrtl(tmp);
    meta_states_info.sica.summary = getSummaryMetaStates(meta_states_info.sica.Tcs);
    clear tmp TC_sica S_sica;
end

%% Kmeans
if (strcmpi(method, 'k-means') || strcmpi(method, 'all'))
    try
        [IDXall, Call, SUMDall, Dall] = kmeans(FNCdynflat_nomean, num_comps, 'distance', dmethod, 'Replicates', kmeans_num_replicates, 'Display', 'iter', 'MaxIter', kmeans_max_iter, ...
            'empty', 'drop');
    catch
        [IDXall, Call, SUMDall, Dall] = icatb_kmeans(FNCdynflat_nomean, num_comps, 'distance', dmethod, 'Replicates', kmeans_num_replicates, 'Display', 'iter', 'MaxIter', kmeans_max_iter, ...
            'empty', 'drop');
    end
    wfncs = reshape(FNCdynflat_nomean, M, Nwin , size(FNCdynflat_nomean, 2));
    wfncs = permute(wfncs, [2, 3, 1]);
    Tcs = kMeansTCs(wfncs, Call', dmethod);
    meta_states_info.kmeans.S = Call';
    meta_states_info.kmeans.IDXall = IDXall;
    meta_states_info.kmeans.SUMDall = SUMDall;
    meta_states_info.kmeans.Dall = Dall;
    meta_states_info.kmeans.Tcs = Qrtl(Tcs);
    meta_states_info.kmeans.summary = getSummaryMetaStates(meta_states_info.kmeans.Tcs);
end

function meta_summary = getSummaryMetaStates(Tcs)
% Get summary from meta states
%
% Number of unique rows (NUM_STATES)
% Change in states (CHANGE_STATES)
% MAX L1 distance between rows (STATE_SPAN)
% SUM of L1 distances between successive states (TOTAL_DIST)

num_states = zeros(size(Tcs, 3), 1);
change_states = zeros(size(Tcs, 3), 1);
state_span = zeros(size(Tcs, 3), 1);
total_dist = zeros(size(Tcs, 3), 1);
for nS = 1:size(Tcs, 3)
    tmp = squeeze(Tcs(:, :, nS));
    d1 = icatb_pdist(tmp, 'all', 1);
    d2 = icatb_pdist(tmp, 'successive', 1);
    idx = unique(tmp, 'rows');
    num_states(nS) = size(idx, 1);
    change_states(nS) = sum(d2 ~= 0); %length(find(max(abs(tmp_diff')) ~= 0));
    state_span(nS) = max(d1);
    total_dist(nS) = sum(d2);
end

freq_levels_info = compute_freq_levels(Tcs);

meta_summary.num_states = num_states;
meta_summary.change_states = change_states;
meta_summary.state_span = state_span;
meta_summary.total_dist = total_dist;
meta_summary.freq_levels_info = freq_levels_info;



function [W, A, icasig] = run_stability_analysis(ica_algorithm, data, num_ica_runs)
% Stability analysis
icasigR = cell(1, num_ica_runs);
fprintf('\n');
disp(['Number of times ICA will run is ', num2str(num_ica_runs)]);
for nRun = 1:length(icasigR)
    fprintf('\n');
    disp(['Run ', num2str(nRun), ' / ', num2str(num_ica_runs)]);
    fprintf('\n');
    [dd1, W, A, icasigR{nRun}]  = icatb_icaAlgorithm(ica_algorithm, data);
end
clear dd1 dd2 dd3;

if (num_ica_runs > 1)
    clear W A;
    [corrMetric, W, A, icasig, bestRun] = icatb_bestRunSelection(icasigR, data);
else
    icasig = icasigR{1};
end


function y = vals2quantile(x, p)

p = 100.*p;
n = length(x);
x = sort(x,1);
q = [0 100*(0.5:(n-0.5))./n 100]';
xx = [x(1,:); x(1:n,:); x(n,:)];
y = interp1q(q,xx,p(:));
y = reshape(y, size(p));


function Y = Qrtl(X)

% if length(size(X))<3 | min(size(X))<2 | sum(X(:)<0)==0
%     error('this is specized function for nonnegative 3D arrays')
% end

Y=NaN.*ones(size(X));

posXinds=find(X>=0);
posX=X(posXinds);
posCutoffs=[0,vals2quantile(posX,[0.25,0.5,0.75]),max(posX(:))+1];
posCutoffs=repmat(posCutoffs,length(posX),1);
posX=repmat(posX,1,5);
pos=posX-posCutoffs;
pos=sign(pos);
pos=diff(pos');
[pos_Qrts,Fj_pos]=find(abs(pos)>=1);
pos_Qrts=pos_Qrts(unique(Fj_pos,'First'));
Y(posXinds)=pos_Qrts;

negXinds = find(X<0);
if (~isempty(negXinds))
    negX=X(negXinds); negX=abs(negX);
    negCutoffs=[0,vals2quantile(abs(negX),[0.25,0.5,0.75]),max(negX(:))+1];
    negCutoffs=repmat(negCutoffs,length(negX),1);
    negX=repmat(negX,1,5);
    neg=negX-negCutoffs;
    neg=sign(neg);
    neg=diff(neg');
    [neg_Qrts,Fj_neg]=find(abs(neg)>=1);
    neg_Qrts=neg_Qrts(unique(Fj_neg,'First'));
    Y(negXinds)=-neg_Qrts;
end

function TCs = kMeansTCs(wFNCs, KmeansCentroids, opts1)
%wFNCs should be numWindows by numCorrelations by numSubjects
%KmeansCentroids should be numCorrelations by numClusters


% if nargin<3
%     opts1=1;
% end
% if nargin<4
%     if length(opts1)==0
%         opts1=1;
%     end
%     opts2=1;
% end

if (~exist('opts1', 'var'))
    opts1 = 'sqeuclidean';
end

eps=0.0001;
nClusts=size(KmeansCentroids,2);

%distance_opts =  {'City', 'sqEuclidean', 'Hamming', 'Correlation', 'Cosine'};

Ds = zeros(size(wFNCs, 1), nClusts, size(wFNCs,3));
for c=1:nClusts
    CtrdRep = repmat(squeeze(KmeansCentroids(:, c))', [size(wFNCs, 1), 1, size(wFNCs, 3)]);
    %ctrdrep=repmat(squeeze(KmeansCentroids(:,c)),[1,size(wFNCs,3)]);
    %     CtrdRep = zeros(size(wFNCs, 1), size(ctrdrep, 1), size(ctrdrep, 2));
    %     for w=1:size(wFNCs,1)
    %         CtrdRep(w,:,:)=ctrdrep;
    %     end
    %     clear ctrdrep
    if (strcmpi(opts1, 'city') || strcmpi(opts1, 'cityblock'))
        % L1 metric
        Ds(1:size(wFNCs,1),c,1:size(wFNCs,3))=squeeze(sum(abs(wFNCs - CtrdRep),2));
    elseif (strcmpi(opts1, 'sqeuclidean'))
        % L2 metric
        Ds(1:size(wFNCs,1),c,1:size(wFNCs,3))=squeeze((sum((wFNCs - CtrdRep).^2,2)));
    elseif (strcmpi(opts1, 'correlation'))
        % Correlation metric
        C1 = num2cell(wFNCs, 2);
        C2 = num2cell(CtrdRep, 2);
        dd = cellfun(@icatb_corr2, C1, C2, 'UniformOutput', false);
        Ds(1:size(wFNCs, 1), c, 1:size(wFNCs, 3)) = cell2mat(dd);
        %Ds(1:size(wFNCs,1),c,1:size(wFNCs,3))= icatb_corr(squeeze(wFNCs(:, );
    elseif (strcmpi(opts1, 'cosine'))
        % Cosine metric
        func_h = @(xa, ya) sum(xa(:).*ya(:))/(sqrt(sum(xa(:).*xa(:))*sum(ya(:).*ya(:))));
        C1 = num2cell(wFNCs, 2);
        C2 = num2cell(CtrdRep, 2);
        dd = cellfun(func_h, C1, C2, 'UniformOutput', false);
        Ds(1:size(wFNCs, 1), c, 1:size(wFNCs, 3)) = cell2mat(dd);
        % Ds(1:size(wFNCs, 1), c, 1:size(wFNCs, 3)) = 1- ((wFNC.C_i)/sqrt((wFNC.wFNC)(C_i.C_i))]/2;
        %Ds = 1;
    else
        % hamming
        if (~isempty(which('pdist2.m')))
            func_h = @(xa, ya) pdist2(xa(:)', ya(:)', 'hamming');
        else
            func_h = @(xa, ya) sum(xor(xa(:), ya(:)))/length(xa);
        end
        
        C1 = num2cell(wFNCs, 2);
        C2 = num2cell(CtrdRep, 2);
        dd = cellfun(func_h, C1, C2, 'UniformOutput', false);
        Ds(1:size(wFNCs, 1), c, 1:size(wFNCs, 3)) = cell2mat(dd);
        
    end
end

%if opts2==1
if (strcmpi(opts1, 'city') || strcmpi(opts1, 'cityblock') || strcmpi(opts1, 'sqeuclidean'))
    denom = squeeze(sum(Ds,2));
    TCs = zeros(size(wFNCs, 1), nClusts, size(wFNCs,3));
    for c=1:nClusts
        TCs(:,c,:)=squeeze(Ds(:,c,:))./denom;
    end
    TCs = (1 - TCs);
elseif (strcmpi(opts1, 'correlation'))
    TCs = (1 + Ds)/2;
elseif (strcmpi(opts1, 'cosine'))
    TCs = (1 - Ds)/2;
else
    % hamming
    TCs = (1 - Ds);
end
%elseif opts2==2
%   denom=Ds+eps.*ones(size(Ds));
%   TCs=1./denom;
%end


function levels_info = compute_freq_levels(TCs)
%% Get frequency for each level, component and subject
%

Levels = unique(TCs(:));
Levels = Levels(:)';
num_windows = size(TCs, 1);
FreqCompLevels = zeros(size(TCs, 2), length(Levels), size(TCs, 3));

for sub = 1:size(TCs, 3)
    for comp = 1:size(TCs, 2)
        for lev = 1:length(Levels)
            FreqCompLevels(comp, lev, sub) = sum(TCs(:, comp, sub) == Levels(lev))/num_windows;
        end
    end
end


levels_info.levels = Levels;
levels_info.FreqCompLevels = FreqCompLevels;