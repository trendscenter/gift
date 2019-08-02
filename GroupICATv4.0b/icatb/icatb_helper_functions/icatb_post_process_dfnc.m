function dfncInfo = icatb_post_process_dfnc(param_file)
%% Post process dFNC
%


icatb_defaults;

%% Select dFNC file
showGUI = 0;
if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('title', 'dFNC Parameter File', 'typeEntity', 'file', 'typeSelection', 'single', 'filter', '*dfnc.mat');
    drawnow;
    if (isempty(param_file))
        error('dFNC parameter file is not selected');
    end
    showGUI = 1;
end

if (ischar(param_file))
    load(param_file);
    outputDir = fileparts(param_file);
else
    dfncInfo = param_file;
    outputDir = dfncInfo.outputDir;
end

if (isempty(outputDir))
    outputDir = pwd;
end

cd (outputDir);

outputDir = pwd;


%% Initialise vars
TR = dfncInfo.TR;

TR = min(TR);

num_clusters = 6;
kmeans_max_iter = 150;
dmethod = 'city';
ica_algorithm = 'infomax';
ica_algorithms_list = {'Infomax', 'Fast ICA', 'Erica', 'Simbec', 'Evd', 'Jade Opac', ...
    'Amuse', 'SDD ICA', 'Radical ICA', 'Combi', 'ICA-EBM', 'ERBM'};
distance_opts =  {'City', 'sqEuclidean', 'Hamming', 'Correlation', 'Cosine'};
kmeans_num_replicates = 5;
num_tests_est_clusters = 10;

try
    num_clusters = dfncInfo.postprocess.num_clusters;
    kmeans_max_iter = dfncInfo.postprocess.kmeans_max_iter;
    dmethod = dfncInfo.postprocess.dmethod;
    kmeans_num_replicates = dfncInfo.postprocess.kmeans_num_replicates;
catch
end


try
    num_tests_est_clusters = dfncInfo.postprocess.num_tests_est_clusters;
catch
end


try
    ica_algorithm = dfncInfo.postprocess.ica.algorithm;
catch
end

num_ica_runs = 5;
try
    num_ica_runs = dfncInfo.postprocess.ica.num_ica_runs;
catch
end

if (isnumeric(ica_algorithm))
    ica_algorithm = lower(ica_algorithms_list{ica_algorithm});
end

%icaIndex = strmatch(ica_algorithm, lower(ica_algorithms_list), 'exact');
regressCovFile = '';
try
    regressCovFile = dfncInfo.postprocess.regressCovFile;
catch
end

num_comps_ica = num_clusters;

try
    num_comps_ica = dfncInfo.postprocess.ica.num_comps;
catch
end

estimate_clusters = 'no';

try
    estimate_clusters = dfncInfo.postprocess.estimate_clusters;
catch
end

meta_method = 'k-means';

try
    meta_method = dfncInfo.postprocess.meta_method;
catch
end

if (showGUI)
    
    guiInputs = struct('estimate_clusters', estimate_clusters, 'num_clusters', num_clusters, 'kmeans_max_iter', kmeans_max_iter, 'dmethod', dmethod, ...
        'ica_algorithm', ica_algorithm, 'num_ica_runs', num_ica_runs, 'num_comps_ica', num_comps_ica, 'regressCovFile', regressCovFile, ...
        'meta_method', meta_method, 'kmeans_num_replicates', kmeans_num_replicates, 'num_tests_est_clusters', num_tests_est_clusters);
    
    results = post_process_dfnc(guiInputs);
    
    if (isempty(results))
        error('Post-processing parameters are not selected');
    end
    
    field_Names = fieldnames(results);
    for n = 1:length(field_Names)
        tmp = results.(field_Names{n});
        if (isempty(tmp))
            tmp = '';
        end
        if (isnumeric(tmp))
            tmp = num2str(tmp);
        else
            tmp = ['''', tmp, ''''];
        end
        str = [field_Names{n}, '=', tmp, ';'];
        eval(str);
    end
    
    %     numParameters = 1;
    %
    %     inputText(numParameters).promptString = 'Select number of clusters/components to compute k-means/ICA on dynamic FNC correlations';
    %     inputText(numParameters).uiType = 'edit';
    %     inputText(numParameters).answerString = num2str(num_clusters);
    %     inputText(numParameters).dataType = 'numeric';
    %     inputText(numParameters).tag = 'num_clusters';
    %     inputText(numParameters).enable = 'on';
    %
    %     numParameters = numParameters + 1;
    %
    %     inputText(numParameters).promptString = 'Select max iterations for computing k-means';
    %     inputText(numParameters).uiType = 'edit';
    %     inputText(numParameters).answerString = num2str(kmeans_max_iter);
    %     inputText(numParameters).dataType = 'numeric';
    %     inputText(numParameters).tag = 'kmeans_max_iter';
    %     inputText(numParameters).enable = 'on';
    %
    %     numParameters = numParameters + 1;
    %
    %     opts = char('City', 'sqEuclidean', 'Hamming', 'Correlation', 'Cosine');
    %     vals = find(strcmpi(opts, dmethod) == 1);
    %     if (isempty(vals))
    %         vals = 1;
    %     end
    %     vals = vals(1);
    %     inputText(numParameters).promptString = 'Select distance method for computing k-means';
    %     inputText(numParameters).uiType = 'popup';
    %     inputText(numParameters).answerString = opts;
    %     inputText(numParameters).dataType = 'string';
    %     inputText(numParameters).tag = 'dmethod';
    %     inputText(numParameters).enable = 'on';
    %     inputText(numParameters).value = vals;
    %
    %
    %     numParameters = numParameters + 1;
    %
    %     inputText(numParameters).promptString = 'Select ICA algorithm in the specified list to do meta state analysis';
    %     inputText(numParameters).uiType = 'popup';
    %     inputText(numParameters).answerString = ica_algorithms_list;
    %     inputText(numParameters).dataType = 'string';
    %     inputText(numParameters).tag = 'ica_algorithm';
    %     inputText(numParameters).enable = 'on';
    %     inputText(numParameters).value = icaIndex;
    %
    %
    %     numParameters = numParameters + 1;
    %
    %     inputText(numParameters).promptString = 'Enter number of ICA runs (MST algorithm)';
    %     inputText(numParameters).uiType = 'edit';
    %     inputText(numParameters).answerString = num2str(num_ica_runs);
    %     inputText(numParameters).dataType = 'numeric';
    %     inputText(numParameters).tag = 'num_ica_runs';
    %     inputText(numParameters).enable = 'on';
    %
    %
    %     % Input dialog box
    %     answer = icatb_inputDialog('inputtext', inputText, 'Title', 'Select Post-processing options', 'handle_visibility',  'on');
    %
    %     if (isempty(answer))
    %         error('Post-processing options are not selected');
    %     end
    %
    %     drawnow;
    %
    %     for n = 1:length(answer)
    %         tmp = answer{n};
    %         if (isnumeric(tmp))
    %             tmp = num2str(tmp);
    %         else
    %             tmp = ['''', tmp, ''''];
    %         end
    %         str = [inputText(n).tag, '=', tmp, ';'];
    %         eval(str);
    %     end
    
end

dfncInfo.postprocess.estimate_clusters = estimate_clusters;
dfncInfo.postprocess.num_clusters = num_clusters;
dfncInfo.postprocess.kmeans_max_iter = kmeans_max_iter ;
dfncInfo.postprocess.dmethod = dmethod;
dfncInfo.postprocess.kmeans_num_replicates = kmeans_num_replicates;
dfncInfo.postprocess.num_tests_est_clusters = num_tests_est_clusters;
dfncInfo.postprocess.meta_method = meta_method;
dfncInfo.postprocess.ica.algorithm = lower(ica_algorithm);
dfncInfo.postprocess.ica.num_ica_runs = num_ica_runs;
dfncInfo.postprocess.ica.num_comps = num_comps_ica;
dfncInfo.postprocess.regressCovFile = regressCovFile;

%% Regress covariates
%if (~showGUI)
%chkCov = icatb_questionDialog('title', 'Regress covariates?',  'textbody', 'Do you want to remove the variance associated with the covariates from dFNC correlations?');
%else
chkCov = 0;
%if (isfield(dfncInfo.postprocess, 'regressCovFile'))
% regressCovFile = dfncInfo.postprocess.regressCovFile;
if (~isempty(regressCovFile))
    chkCov = 1;
end
%end
%end

dist_val = strmatch(lower(dmethod), lower(distance_opts), 'exact');
dmethod = distance_opts{dist_val};

if (chkCov)
    % if (showGUI)
    %     regressCovFile = icatb_selectEntry('title', 'Select covariates file (continuous covariate or reduced model from mancova fnc results) ...', 'typeEntity', 'file', ...
    %         'typeSelection', 'single', 'filter', '*mancovan*results*fnc.mat;*txt');
    % end
    drawnow;
    checkUNI = 0;
    try
        chkVarNames = whos('-file', regressCovFile);
        checkUNI = ~isempty(strmatch('UNI', char (chkVarNames.name), 'exact'));
    catch
    end
    
    if (checkUNI)
        load(regressCovFile);
        regressCovX = UNI.stats{1}.X;
        allTerms = UNI.stats{1}.Terms;
        chkMeanIndex = [];
        for nT = 1:length(allTerms)
            if all(allTerms{nT} == 0)
                chkMeanIndex = nT;
                break;
            end
        end
        
        if (~isempty(chkMeanIndex))
            regressCovX(:, chkMeanIndex) = [];
        end
        
        if (isempty(regressCovX))
            error('No covariates found from reduced model.');
        end
    else
        regressCovX = icatb_load_ascii_or_mat(regressCovFile);
        % remove mean of the covariate (assuming its continuous)
        regressCovX = detrend(regressCovX, 0);
    end
    
    if (size(regressCovX, 1) ~= dfncInfo.userInput.numOfSub)
        error(['Please check the dimensions of covariates. Number of rows must match the number of subjects ', num2str(dfncInfo.userInput.numOfSub), ')']);
    end
    
    doRegressCovariates(dfncInfo, regressCovX);
    
end

M = length(dfncInfo.outputFiles);

SP = cell(M,1);
k1_peaks = zeros(1,M);

disp('Computing spectra on FNC correlations ...');

%% Frequency analysis
for nR = 1:length(dfncInfo.outputFiles)
    
    current_file = fullfile(outputDir,  dfncInfo.outputFiles{nR});
    load(current_file, 'FNCdyn');
    
    if (nR == 1)
        Nwin = size(FNCdyn, 1);
        Fs = 1/TR;
        nfft = 2^(nextpow2(Nwin)+1); % pad to the next power of 2
        df=Fs/nfft;
        f=0:df:Fs;    % all possible frequencies
        f=f(1:nfft);
        fpass = [0.0, 1/(2*TR)];
        findx=find(f>=fpass(1) & f<=fpass(end)); % just one-half of the spectrum
        f = f(findx);
        dfncInfo.freq = f;
        FNCdynflat = zeros(M, Nwin, size(FNCdyn, 2));
    end
    
    
    FNCdynflat(nR, :, :) = FNCdyn;
    
    DEV = std(FNCdyn, [], 2);
    [xmax, imax, xmin, imin] = icatb_extrema(DEV); % find the extrema
    pIND = sort(imax);
    k1_peaks(nR) = length(pIND);
    SP{nR} = FNCdyn(pIND, :);
    
    %% Compute spectra
    FNCdyn = FNCdyn - repmat(mean(FNCdyn, 1), Nwin, 1);
    FNCdyn = FNCdyn.*repmat(icatb_hamming(Nwin), 1, size(FNCdyn, 2));
    S = fft(FNCdyn, nfft, 1)/Fs;
    S = S(findx,:);
    spectra_fnc = sqrt(S.*conj(S));
    
    icatb_save(current_file, 'spectra_fnc', '-append');
    
    %tmp_amp = squeeze(trapz(f, spectra_fnc));
    tmp_amp = squeeze(std(spectra_fnc));
    tmp_cm = squeeze(mean(spectra_fnc.*repmat(f(:), [1, size(spectra_fnc, 2)])));
    tmp_cm = tmp_cm ./ mean(spectra_fnc);
    
    if (nR == 1)
        FNCamp = zeros(size(tmp_amp));
        FNCcm = zeros(size(tmp_cm));
    end
    
    FNCamp = FNCamp + tmp_amp;
    FNCcm = FNCcm + tmp_cm;
    
    clear FNCdyn;
end

FNCamp = FNCamp / M;
FNCcm = FNCcm / M;


post_process_file = fullfile(outputDir,  [dfncInfo.prefix, '_post_process.mat']);
icatb_save(post_process_file, 'FNCamp', 'FNCcm');

clear FNCamp FNCcm;

fprintf('\n');

disp('Computing k-means on FNC correlations ...');

%% Cluster
SPflat = cell2mat(SP);

clear SP;

FNCdynflat = reshape(FNCdynflat, M*Nwin , size(FNCdynflat, 3));

% Determine optmial number of clusters
if (strcmpi(estimate_clusters, 'yes'))
    
    cluster_estimate_results = icatb_optimal_clusters(SPflat, min([max(size(SPflat)), 10]), 'method', 'all', 'cluster_opts', {'Replicates', kmeans_num_replicates, 'Distance', dmethod, ...
        'MaxIter', kmeans_max_iter}, 'num_tests', num_tests_est_clusters, 'display', 1);
    
    %     disp('Using gap method ...');
    %     gap_stats = icatb_optimal_clusters(SPflat, min([max(size(SPflat)), 10]), 'method', 'gap', 'num_tests', num_tests_est_clusters);
    %     if (gap_stats.K == 1)
    %         warning('Number of clusters estimated is 1 using gap statistic. Setting number of estimated clusters to 2.');
    %         gap_stats.K = 2;
    %         % error('Number of clusters estimated is 1 using gap statistic. Use specified number of clusters.');
    %     end
    %     disp('Using Silhouette method ...');
    %     silh_stats = icatb_optimal_clusters(SPflat, min([max(size(SPflat)), 10]), 'method', 'silhouette');
    %
    %     %num_clusters = gap_stats.K;
    %
    %     if (~isempty(which('errorbar.m')))
    %         figH = figure('name', 'Cluster Estimation', 'color', [1, 1, 1]);
    %
    %         sh = subplot(2, 1, 1);
    %         errorbar((1:length(gap_stats.values)), gap_stats.values, gap_stats.sem, '-bo');
    %         xlabel('Clusters', 'parent', sh);
    %         ylabel('Gap values', 'parent', sh);
    %         axis(sh, 'tight');
    %         title(['Optimal clusters  = ', num2str(gap_stats.K)], 'parent', sh);
    %
    %         sh = subplot(2, 1, 2);
    %         plot(silh_stats.klist, silh_stats.values, '-go');
    %         xlabel('Clusters', 'parent', sh);
    %         ylabel('Silhouette', 'parent', sh);
    %         axis(sh, 'tight');
    %         title(['Optimal clusters  = ', num2str(silh_stats.K)], 'parent', sh);
    %         drawnow;
    %
    %     end
    
    %num_clusters = ceil(mean([gap_stats.K, silh_stats.K]));
    num_clusters = 0;
    for Ntests = 1:length(cluster_estimate_results)
        num_clusters = num_clusters + (cluster_estimate_results{Ntests}.K(1));
    end
    num_clusters = ceil(num_clusters/length(cluster_estimate_results));
    %     cluster_estimate_results.num_clusters = num_clusters;
    %     cluster_estimate_results.gap_stats = gap_stats;
    %     cluster_estimate_results.silh_stats = silh_stats;
    disp(['Number of estimated clusters used in dFNC standard analysis is mean of all tests: ', num2str(num_clusters)]);
    fprintf('\n');
    
    dfncInfo.postprocess.num_clusters = num_clusters;
    
end



try
    [IDXp, Cp, SUMDp, Dp] = kmeans(SPflat, num_clusters, 'distance', dmethod, 'Replicates', kmeans_num_replicates, 'MaxIter', kmeans_max_iter, 'Display', 'iter', 'empty', 'drop');
catch
    [IDXp, Cp, SUMDp, Dp] = icatb_kmeans(SPflat, num_clusters, 'distance', dmethod, 'Replicates', kmeans_num_replicates, 'MaxIter', kmeans_max_iter, 'Display', 'iter', 'empty', 'drop');
end

% Save subsampled data
clusterInfo.SP = SPflat;
clear SPflat IDXp SUMDp;

%% Use as a starting point to cluster all the data
try
    [IDXall, Call, SUMDall, Dall] = kmeans(FNCdynflat, num_clusters, 'distance', dmethod, 'Replicates', 1, 'Display', 'iter', 'MaxIter', kmeans_max_iter, ...
        'empty', 'drop', 'Start', Cp);
catch
    [IDXall, Call, SUMDall, Dall] = icatb_kmeans(FNCdynflat, num_clusters, 'distance', dmethod, 'Replicates', 1, 'Display', 'iter', 'MaxIter', kmeans_max_iter, ...
        'empty', 'drop', 'Start', Cp);
end

clusterInfo.Call = Call;
clusterInfo.Cp = Cp;
clusterInfo.IDXall = IDXall;
clusterInfo.SUMDall = SUMDall;
clusterInfo.Dall = Dall;

[clusterInfo.corrs_states, clusterInfo.states] = getStateCorr(dfncInfo, clusterInfo);

icatb_save(post_process_file, 'clusterInfo', '-append');

if (exist('cluster_estimate_results', 'var'))
    icatb_save(post_process_file, 'cluster_estimate_results', '-append');
end


%% Meta state analysis
FNCdynflat = reshape(FNCdynflat, M, Nwin, size(FNCdynflat, 2));
meta_states_info = icatb_dfnc_meta_state_analysis(FNCdynflat, num_comps_ica, 'dmethod', dmethod, 'kmeans_max_iter', kmeans_max_iter, 'num_ica_runs', num_ica_runs, ...
    'ica_algorithm', lower(ica_algorithm), 'method', meta_method, 'replicates', kmeans_num_replicates);

icatb_save(post_process_file, 'meta_states_info', '-append');


param_file = fullfile(outputDir, [dfncInfo.prefix, '.mat']);
icatb_save(param_file, 'dfncInfo');

disp('Done');

fprintf('\n');


function doRegressCovariates(dfncInfo, regressCovX)
%% Regress variance from the dFNC correlations
%

disp('Regressing variance associated with the covariates from dfnc correlations ...');

%% Loop over sessions
for nSess = 1:dfncInfo.userInput.numOfSess
    
    %% Load subject dfnc correlations and form matrix subjects x (number of windows x component pairs)
    for nSub = 1:dfncInfo.userInput.numOfSub
        results_file = [dfncInfo.userInput.prefix, '_sub_', icatb_returnFileIndex(nSub), '_sess_', icatb_returnFileIndex(nSess), '_results.mat'];
        load(fullfile(dfncInfo.outputDir, results_file));
        FNCdyn = FNCdyn(:)';
        if (nSub == 1)
            FNC_corrs = zeros(dfncInfo.userInput.numOfSub, length(FNCdyn));
        end
        FNC_corrs(nSub, :) = FNCdyn;
    end
    
    %% Regress variance of covariates from the dFNC correlations
    betas = pinv(regressCovX)*FNC_corrs;
    FNC_corrs = FNC_corrs - regressCovX*betas;
    
    
    %% Write variance removed dFNC correlations
    for nSub = 1:dfncInfo.userInput.numOfSub
        results_file = [dfncInfo.userInput.prefix, '_sub_', icatb_returnFileIndex(nSub), '_sess_', icatb_returnFileIndex(nSess), '_results.mat'];
        load(fullfile(dfncInfo.outputDir, results_file));
        [Nwin, comp_pairs] = size(FNCdyn);
        FNCdyn_old = FNCdyn;
        clear FNCdyn;
        FNCdyn = FNC_corrs(nSub, :);
        FNCdyn = reshape(FNCdyn, Nwin, comp_pairs);
        icatb_save(fullfile(dfncInfo.outputDir, results_file), 'FNCdyn', 'FNCdyn_old', '-append');
    end
    
    clear FNC_corrs;
    
end
%% End of loop over sessions

fprintf('Done\n');



function [corrs, states] = getStateCorr(dfncInfo, clusterInfo)
% corrs - subjects x sessions x clusters
% states - subjects x sessions x windows
outDir = dfncInfo.outputDir;

Nwin = size(clusterInfo.IDXall, 1)/(dfncInfo.userInput.numOfSess*dfncInfo.userInput.numOfSub);

% subjects x sessions x windows
states = reshape(clusterInfo.IDXall, dfncInfo.userInput.numOfSess, dfncInfo.userInput.numOfSub, Nwin);
states = permute(states, [2, 1, 3]);

count = 0;
for nSub = 1:dfncInfo.userInput.numOfSub
    for nSess = 1:dfncInfo.userInput.numOfSess
        count = count + 1;
        fn = fullfile(outDir, dfncInfo.outputFiles{count});
        load(fn);
        tmp = squeeze(states(nSub, nSess, :));
        for nState = 1:dfncInfo.postprocess.num_clusters
            inds = find(tmp == nState);
            corrs{nSub, nSess, nState} = FNCdyn(inds, :);
        end
    end
end