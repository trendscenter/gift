function icatb_post_process_spatial_dfnc(param_file, varargin)
%% Post process spatial dfnc
%


%% Initialize vars
num_clusters = 3;
kmeans_max_iter = 150;
dmethod = 'city';
threshold_type = 'none';
threshold = 0.05;
transition_threshold = 0.35;
num_permutations = 20;

for n = 1:2:length(varargin)
    if (strcmpi(varargin{n}, 'comps'))
        comps = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'num_clusters'))
        num_clusters = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'dmethod'))
        dmethod = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'max_iter'))
        kmeans_max_iter = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'threshold_type'))
        threshold_type = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'threshold'))
        threshold = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'transition_threshold'))
        transition_threshold = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'num_permutations'))
        num_permutations = varargin{n + 1};
    end
end

%% Select sdFNC file
if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('title', 'Select spatial dFNC Parameter File', 'typeEntity', 'file', 'typeSelection', 'single', 'filter', '*_sdfnc.mat');
    drawnow;
    if (isempty(param_file))
        error('ICA/sdFNC parameter file is not selected');
    end
end

drawnow;

[outputDir, fN, extn] = fileparts(param_file);
if (isempty(outputDir))
    outputDir = pwd;
end
param_file = fullfile(outputDir, [fN, extn]);

load (param_file);

if (~exist('sdfncInfo', 'var'))
    error('Selected file is not a valid spatial dfnc file');
end

if (~isfield(sdfncInfo, 'resultFiles'))
    error('Please run spatial dfnc in order to do post-processing');
end


%% Open GUI to get postprocess params.
outputFiles = sdfncInfo.resultFiles;
outputDir = sdfncInfo.userInput.outputDir;
%numComp = sdfncInfo.numComp;
numOfSub = sdfncInfo.numOfSub;
numOfSess = sdfncInfo.numOfSess;
prefix = sdfncInfo.userInput.prefix;
windows = sdfncInfo.windows;

cd (outputDir);

compFiles = icatb_rename_4d_file(fullfile(outputDir, sdfncInfo.compFiles));

useGUI = 0;
if (~exist('comps', 'var'))
    useGUI = 1;
    comps = [];
end

defaultParams = struct('comps', comps, 'num_clusters', num_clusters, 'max_iter', kmeans_max_iter, 'dmethod', dmethod, 'threshold_type', threshold_type, ...
    'threshold', threshold, 'transition_threshold', transition_threshold, 'num_permutations', num_permutations);

params = defaultParams;
fNames = fieldnames(defaultParams);
try
    params = sdfncInfo.postprocess.params;
    for i = 1:length(fNames)
        if (~isfield(params, fNames{i}))
            params.(fNames{i}) = defaultParams.(fNames{i});
        end
    end
catch
end

if (useGUI)
    params = getParams(compFiles, params);
end

drawnow;

sdfncInfo.postprocess.params = params;
icatb_save(param_file, 'sdfncInfo');

talpha = params.threshold;

%% Post-process file
post_process_file = [prefix, '_sdfnc_post_process.mat'];
sdfncInfo.postprocess.file = post_process_file;
post_process_file = fullfile(outputDir, post_process_file);

comps = params.comps;
numComp = length(comps);

%% Compute FNC metrics
disp('Computing FNC metrics using mutual information between spatial components ...');
MI = zeros(numComp, numComp, numOfSub, numOfSess, windows);
for nSub = 1:size(outputFiles, 1)
    for nSess = 1:size(outputFiles, 2)
        fname = fullfile(outputDir, outputFiles{nSub, nSess});
        load(fname, 'ic');
        ic = ic(comps, :, :);
        MI(:, :, nSub, nSess, :) = compute_mi(ic);
        clear ic;
    end
end
icatb_save(post_process_file, 'MI');

drawnow;


%% Get FC metrics for each group
% Average sessions or runs
MI_subjects = squeeze(mean(MI, 4));
% Standard deviation across windows
std_MI_subjects = std(MI_subjects + eps, 0, 4);

groupMIvals = cell(1, sdfncInfo.numGroups);
MIvals = cell(1, sdfncInfo.numGroups);
MIStdVals = MIvals;

eS = 0;
for nG = 1:sdfncInfo.numGroups
    sS = eS + 1;
    eS = eS + length(sdfncInfo.userInput.group(nG).val);
    val = std_MI_subjects(:, :, sS:eS);
    MIStdVals{nG} = val;
    MIvals{nG} = MI_subjects(:, :, sS:eS, :);
end

icatb_save(post_process_file, 'MIStdVals', '-append');

%% KL divergence between windows
disp('Computing KL divergence metrics for each group between consecutive windows ...');
KL = computeKLD(MIvals);
icatb_save(post_process_file, 'KL', '-append');

%% Markov modeling
% Perform markov modeling to get finite states:
%  Use each component MI's across window and subjects to compute states.
disp('Computing k-means on each component ...');
clusterInfo = sdfnStateVecStats(MIvals, params);
icatb_save(post_process_file, 'clusterInfo', '-append');

%% Compute t-tests for each window and median test of std values across windows
if (sdfncInfo.numGroups > 1)
    
    groupNames = cellstr(char(sdfncInfo.userInput.group.name));
    group_combinations = nchoosek(1:sdfncInfo.numGroups, 2);
    groupCombNames = cell(1, size(group_combinations, 1));
    
    
    median_test_results = NaN(length(groupCombNames), numComp, numComp, 2);
    for nC = 1:size(group_combinations, 1)
        g1Ind = group_combinations(nC, 1);
        g2Ind = group_combinations(nC, 2);
        combName = [groupNames{g1Ind}, ' vs ', groupNames{g2Ind}];
        groupCombNames{nC} = combName;
        
        me_pvalues = NaN*ones(numComp, numComp);
        diff_values = me_pvalues;
        
        %% Median
        for m = 1:numComp
            for n = m + 1:numComp
                a = squeeze(MIStdVals{g1Ind}(m, n, :));
                b = squeeze(MIStdVals{g2Ind}(m, n, :));
                [p, h, stats] = icatb_ranksum(a, b, 'alpha', talpha);
                if (h == 1)
                    me_pvalues(m, n) = p;
                    me_pvalues(n, m) = p;
                    diff_values(m, n) = median(a) - median(b);
                    diff_values(n, m) =  diff_values(m, n);
                    %median_test_results(nC, m, n, 1) = p;
                    %median_test_results(nC, m, n, 2) = median(a) - median(b);
                    % median_test_results(nC, m, n, 3) = median(b);
                end
            end
        end
        
        median_test_results(nC, :, :, 1) = me_pvalues;
        median_test_results(nC, :, :, 2) = diff_values;
        
        
        ttest_results = NaN(length(groupCombNames), windows, numComp, numComp, 4);
        
        %% T-tests at each window
        for nwin = 1:windows
            
            t_values = NaN*ones(numComp, numComp);
            p_values = t_values;
            
            for m = 1:numComp
                for n = m + 1:numComp
                    a = squeeze(MIvals{g1Ind}(m, n, :, nwin));
                    b = squeeze(MIvals{g2Ind}(m, n, :, nwin));
                    [p_value, t_value] = icatb_ttest2(a, b, 0);
                    if (p_value < talpha)
                        t_values(m, n) = t_value;
                        t_values(n, m) = t_values(m, n);
                        p_values(m, n) = p_value;
                        p_values(n, m) = p_values(m, n);
                        %ttest_results(nC, nwin, m, n, 1) = t_value;
                        %ttest_results(nC, nwin, m, n, 2) = p_value;
                        ttest_results(nC, nwin, m, n, 3) = mean(a);
                        ttest_results(nC, nwin, n, m, 3) =  ttest_results(nC, nwin, m, n, 3);
                        ttest_results(nC, nwin, m, n, 4) = mean(b);
                        ttest_results(nC, nwin, n, m, 4) = ttest_results(nC, nwin, m, n, 4);
                    end
                end
            end
            
            good_p_vals = getSignificantVals(p_values, params.threshold, params.threshold_type);
            t_values(good_p_vals == 0) = NaN;
            p_values(good_p_vals == 0) = NaN;
            
            ttest_results(nC, nwin, :, :, 1) = t_values;
            ttest_results(nC, nwin, :, :, 2) = p_values;
            
            clear t_values p_values t_value p_value;
            
        end
        
    end
    
    
    icatb_save(post_process_file, 'groupCombNames', 'median_test_results', 'ttest_results', '-append');
    
    
    
end





%% Compute group differences at each window
% Do permutation testing and compute KLD for adjacent windows



function MI = compute_mi(icc)
%% Compute FNC metrics using mutual information
%

order = size(icc, 1);
windows = size(icc, 3);

MI = zeros(order, order, windows);
for k = 1:windows
    temp = squeeze(icc(:, :, k));
    temp_fnc = icatb_compute_mi(temp);
    %     temp_fnc = zeros(order, order);
    %     for m = 1:order
    %         for n = m+1:order
    %             %MI(m, n, k) = sqrt(1-exp(-2*mutualinfo(icatb_zscore(temp(m, :)), icatb_zscore(temp(n, :)))));
    %             temp_fnc(m, n) = sqrt(1-exp(-2*mutualinfo(icatb_zscore(temp(m, :)), icatb_zscore(temp(n, :)))));
    %             temp_fnc(n, m) = temp_fnc(m, n);
    %         end
    %     end
    MI(:, :, k) = temp_fnc;
    clear temp temp_fnc;
end

% function KL = computeKLD(data, numPermutations, fractionTrials)
% %% KL divergence
% % data - comps x comps x groups
%
%
% nosub = size(data, 2);
% windows = size(data, 3);
%
% frac1 = round(nosub*fractionTrials);
%
% %% KL metric
% KL = zeros(numPermutations, frac1, windows - 1);
%
% %% Loop over permutations
% for p = 1:numPermutations
%
%     I1 = randperm(nosub);
%     I1 = sort(I1(1:frac1));
%
%     % Compute KLD metric
%     for nwin = 1:windows - 1
%         for nSub = 1:frac1
%             temp1 = squeeze(data(:, I1(nSub), nwin));
%             temp2 = squeeze(data(:, I1(nSub), nwin + 1));
%             f_w1 = ksdensity(temp1);
%             f_w2 = ksdensity(temp2);
%             temp = KLdist(f_w1, f_w2);
%             KL(p, nSub, nwin) = 1 - exp(-temp);
%         end
%     end
%
% end
% %% End of loop over permutations


function KL = computeKLD(MIvals)
%% KL divergence

numGroups = length(MIvals);
KL = cell(1, numGroups);

[numComp, b, numOfSub, windows] = size(MIvals{1});

%% KL metric
for nG = 1:numGroups
    mi_vals = MIvals{nG};
    kl_vals = zeros(numOfSub, windows - 1);
    for nSub = 1:numOfSub
        for w = 1:windows - 1
            temp1 = squeeze(mi_vals(:, :, nSub, w));
            temp2 = squeeze(mi_vals(:, :, nSub, w + 1));
            temp1 = icatb_mat2vec(temp1);
            temp2 = icatb_mat2vec(temp2);
            %f_w1 = ksdensity(temp1);
            %f_w2 = ksdensity(temp2);
            f_w1 = icatb_kde(temp1);
            f_w2 = icatb_kde(temp2);
            val = KLdist(f_w1, f_w2);
            kl_vals(nSub, w) = 1 - exp(-val);
        end
    end
    KL{nG} = kl_vals;
end


function clusterInfo = sdfnStateVecStats(MIvals, params)
% average MI state for each group
% counts for each cluster state for each group
% Transition matrix for each group
%

comps = params.comps;
clusterno = params.num_clusters;
kmeans_max_iter = params.max_iter;
kmeans_dmethod = params.dmethod;
numComp = length(params.comps);
num_permutations = params.num_permutations;
transition_threshold = params.transition_threshold;

numGroups = length(MIvals);
clusterInfo = cell(1, numComp);

totalSubjects = 0;
for i = 1:length(MIvals)
    totalSubjects = totalSubjects + size(MIvals{i}, 3);
end

for nComp = 1:numComp
    
    compsToInclude = (1:numComp);
    compsToInclude(nComp) = [];
    
    disp(['Performing k-means clustering on component ', num2str(comps(nComp)), ' ...']);
    
    e = 0;
    subjects = zeros(1, numGroups);
    for nG = 1:length(MIvals)
        mi = MIvals{nG};
        mi = squeeze(mi(nComp, compsToInclude, :, :));
        mi = permute(mi, [3, 2, 1]);
        windowno = size(mi, 1);
        nosub = size(mi, 2);
        subjects(nG) = nosub;
        s = e + 1;
        e = e + nosub*windowno;
        mi = reshape(mi, windowno*nosub, size(mi, 3));
        if (nG == 1)
            X = zeros(totalSubjects*windowno, length(compsToInclude));
        end
        X(s:e, :) = mi;
        clear mi;
    end
    
    try
        [IDXall, Call, SUMDall, Dall] = kmeans(X, clusterno, 'distance', kmeans_dmethod, 'Replicates', 1, 'Display', 'iter', 'MaxIter', kmeans_max_iter, ...
            'empty', 'drop');
    catch
        [IDXall, Call, SUMDall, Dall] = icatb_kmeans(X, clusterno, 'distance', kmeans_dmethod, 'Replicates', 1, 'Display', 'iter', 'MaxIter', kmeans_max_iter, ...
            'empty', 'drop');
    end
    
    e = 0;
    frequency_states = zeros(numGroups, clusterno);
    sums_groups = cell(1, numGroups);
    transition_matrix = cell(1, numGroups);
    state_mean = cell(numGroups, clusterno);
    state_sem = cell(numGroups, clusterno);
    
    for nG = 1:numGroups
        s = e + 1;
        e = e + subjects(nG)*windowno;
        idx = IDXall(s:e);
        [transition_matrix{nG}, sums_groups{nG}] = getTM(idx, subjects(nG), windowno, clusterno, num_permutations, transition_threshold);
        for cn = 1:clusterno
            clusterInds = find(idx == cn);
            frequency_states(nG, cn)   = length(clusterInds);
            if (length(clusterInds) > 0)
                state_mean{nG, cn} = mean(X(clusterInds, :));
                state_sem{nG, cn} = std(X(clusterInds, :))./sqrt(subjects(nG));
            end
        end
    end
    
    out.frequency_states = frequency_states;
    out.transition_matrix = transition_matrix;
    out.sum = sums_groups;
    out.state_mean = state_mean;
    out.state_sem = state_sem;
    out.IDXall = IDXall;
    clusterInfo{nComp} = out;
    clear out;
    
end


function [TM, sum_hc] = getTM(idx, nosub, windowno, clusterno, totalpermuteno, transThreshold)

sum_hc = zeros(clusterno);
TM = sum_hc;

for n = 1:totalpermuteno
    selectsub = round(nosub*0.7);
    Ih  = rand_int(1, nosub, selectsub, 1, 1);
    temp1 = sort(Ih);
    hc_c= zeros(length(temp1), windowno);
    for j = 1:length(temp1)
        hc_c(j,:) = idx(windowno*(temp1(j)-1)+1:windowno*temp1(j));
    end
    hc_tm = zeros(clusterno);
    for k = 1:selectsub
        for i = 1:windowno-1
            hc_tm(hc_c(k,i),hc_c(k,i+1)) = hc_tm(hc_c(k,i),hc_c(k,i+1)) + 1;
        end
    end
    for l = 1:clusterno
        if sum(hc_tm(l,:))>0
            hc_tm(l,:) = hc_tm(l,:)/sum(hc_tm(l,:));
        end
    end
    
    [row,col] = find(hc_tm >= transThreshold);
    for m = 1:size(row,1)
        sum_hc(row(m),col(m)) = sum_hc(row(m),col(m)) + 1;
    end
    TM = TM + hc_tm;
end

TM = TM/totalpermuteno;

function params = getParams(compFiles, defaults)
%% Get params
%

% num_clusters = 6;
% kmeans_max_iter = 150;
% dmethod = 'city';
% comps = (1:size(compFiles, 1));
% list_string = num2str((1:length(comps))');
% threshold_type = 'none';
% threshold = 0.05;
%
% for n = 1:2:length(varargin)
%     if (strcmpi(varargin{n}, 'comps'))
%         comps = varargin{n + 1};
%     elseif (strcmpi(varargin{n}, 'num_clusters'))
%         num_clusters = varargin{n + 1};
%     elseif (strcmpi(varargin{n}, 'method'))
%         dmethod = varargin{n + 1};
%     elseif (strcmpi(varargin{n}, 'max_iter'))
%         kmeans_max_iter = varargin{n + 1};
%     elseif (strcmpi(varargin{n}, 'threshold_type'))
%         threshold_type = varargin{n + 1};
%     elseif (strcmpi(varargin{n}, 'threshold'))
%         threshold = varargin{n + 1};
%     end
% end

handles = sdfnc_post_process({compFiles, defaults});

waitfor(handles);

appName = 'getParamsAppData';
if (isappdata(0, appName))
    params = getappdata(0, appName);
    rmappdata(0, appName);
else
    error('Parameters are not selected');
end

%
% handles = guidata(handles);
%
% params.comps = get(handles.comp, 'userdata');
% params.num_clusters = str2num(get(handles.num_clusters, 'string'));
% params.max_iter = str2num(get(handles.max_iter, 'string'));
% opts = cellstr(get(handles.dmethod, 'string'));
% val = get(handles.dmethod, 'value');
% params.d_method = lower(opts{val});
%
% opts = cellstr(get(handles.threshold_type, 'string'));
% val = get(handles.threshold_type, 'value');
% params.threshold_type = lower(opts{val});
% params.threshold = str2num(get(handles.threshold, 'string'));


function good_inds = getSignificantVals(ps, thresh, threshdesc)

ps = ps(:);

if (strcmpi(threshdesc, 'fdr'))
    p_masked = icatb_fdr(ps, thresh);
else
    p_masked = thresh;
end

good_inds = ps < p_masked;

%pF(good_inds) = ps(good_inds);




function dist = KLdist(p,q)
idx = find(q==0);
qq = q;
qq(idx) = eps;
pp=p/sum(p);
qq=qq/sum(qq);

tem = 0;
for j = 1:length(pp)
    m=pp(j);n=qq(j);
    if m~=0
        temp =  m*log(m/n);
        tem = tem + temp;
    end
end

dist = tem;
