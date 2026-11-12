function dfncInfo = icatb_fb_post_process(param_file, showGUI)

%% Post process dFNC
%

icatb_defaults;
global DFNC_DEFAULTS;
use_meta_state = 0;
try
    use_meta_state = DFNC_DEFAULTS.META_STATE;
catch
end

use_tall_array = 'no';

try
    use_tall_array = DFNC_DEFAULTS.use_tall_array;
catch
end

%% Select dFNC file
if (~exist('showGUI', 'var'))
    showGUI = 0;
end
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
    try
        outputDir = dfncInfo.outputDir;
    catch
    end
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
s_dist_method = 'sqEuclidean';
ica_algorithm = 'infomax';
ica_algorithms_list = {'Infomax', 'Fast ICA', 'Erica', 'Simbec', 'Evd', 'Jade Opac', ...
    'Amuse', 'SDD ICA', 'Radical ICA', 'Combi', 'ICA-EBM', 'ERBM'};
distance_opts =  {'City', 'sqEuclidean', 'Hamming', 'Correlation', 'Cosine'};
kmeans_num_replicates = 30;
num_tests_est_clusters = 10;
ref_chk_ena_stateguided = 1; % May be used in case batch is run without var
tag_edt_stateguided_numcomps = 6; % May be used in case batch is run without var
b_chk_ena_stateguided = 0; % May be used in case batch is run without var
b_chk_ena_statebased = 1;% May be used in case batch is run without var
b_chk_ena_metastate = 0;% May be used in case batch is run without var

try
    num_clusters = dfncInfo.postprocess.num_clusters;
    kmeans_max_iter = dfncInfo.postprocess.kmeans_max_iter;
    s_dist_method = dfncInfo.postprocess.dmethod;
    kmeans_num_replicates = dfncInfo.postprocess.kmeans_num_replicates;
    ref_chk_ena_stateguided = dfncInfo.postprocess.ref_chk_ena_stateguided;
catch
end

try
    num_tests_est_clusters = dfncInfo.postprocess.num_tests_est_clusters;
catch
end

try
    use_tall_array = dfncInfo.postprocess.use_tall_array;
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

kmeans_start = 'default';

try

    kmeans_start = dfncInfo.postprocess.kmeans_start;

catch

end

backProject = 0;

try

    if (isfield(dfncInfo.postprocess, 'Cp'))

        backProject = 1;

        Cp = dfncInfo.postprocess.Cp;

        if (size(Cp, 1) ~= num_clusters)

            error('No of rows in the initial centroids passed doesn''t match the number of clusters entered');

        end

    end

catch

end

if (showGUI)
    s_dist_method = 'sqEuclidean'; %Default
    fuzzival = 1.5; %Default
    guiInputs = struct('estimate_clusters', estimate_clusters, 'num_clusters', num_clusters, ...
        'kmeans_max_iter', kmeans_max_iter, 's_dist_method', s_dist_method, 'ica_algorithm', ica_algorithm, ...
        'num_ica_runs', num_ica_runs, 'num_comps_ica', num_comps_ica, ...
        'fuzzival', fuzzival, 'regressCovFile', regressCovFile, 'meta_method', meta_method, ...
        'kmeans_num_replicates', kmeans_num_replicates, ...
        'num_tests_est_clusters', num_tests_est_clusters, 'kmeans_start', kmeans_start, ...
        'use_tall_array', use_tall_array, ...
        'tag_edt_stateguided_numcomps', tag_edt_stateguided_numcomps, ...
        'tag_chk_ena_stateguided', b_chk_ena_stateguided, ...
        'tag_chk_ena_statebased', b_chk_ena_statebased, ...
        'tag_chk_ena_metastate', b_chk_ena_metastate);

    results = icatb_fb_gui_postproc(guiInputs);

    if (isempty(results))

        error('Post-processing parameters are not selected');

    end

    field_Names = fieldnames(results);

    for n = 1:length(field_Names)

        tmp = results.(field_Names{n});

        if (isempty(tmp))

            tmp = '';

        end

        %         if (isnumeric(tmp))

        %             tmp = num2str(tmp);

        %         else

        %             tmp = ['''', tmp, ''''];

        %         end

        %str = [field_Names{n}, '=', tmp, ';'];

        str = [field_Names{n}, ' = tmp;'];

        eval(str);

    end

end

% State Guided

dfncInfo.postprocess.b_chk_ena_stateguided = b_chk_ena_stateguided;

if b_chk_ena_stateguided

    if showGUI

        dfncInfo.postprocess.ref_chk_ena_stateguided = ref_chk_ena_stateguided;

        dfncInfo.postprocess.tag_edt_stateguided_numcomps = str2num(tag_edt_stateguided_numcomps);

    %else if batch script was run these variables are already set

    end

end

% State Based

dfncInfo.postprocess.b_chk_ena_statebased = b_chk_ena_statebased;

if b_chk_ena_statebased

    dfncInfo.postprocess.estimate_clusters = estimate_clusters;

    dfncInfo.postprocess.num_clusters = num_clusters;

end

% MEta state

dfncInfo.postprocess.b_chk_ena_metastate = b_chk_ena_metastate;

dfncInfo.postprocess.dmethod = s_dist_method; %







   

   



% Regress

dfncInfo.postprocess.regressCovFile = regressCovFile;

%

dfncInfo.postprocess.num_tests_est_clusters = num_tests_est_clusters;

dfncInfo.postprocess.kmeans_max_iter = kmeans_max_iter ;

dfncInfo.postprocess.kmeans_num_replicates = kmeans_num_replicates;

dfncInfo.postprocess.kmeans_start = kmeans_start;

if (strcmpi(kmeans_start, 'user Input'))

    estimate_clusters = 'no';

    if (~isempty(Cp))

        disp('Using initial centroids as input in the standard dfnc ...');

        fprintf('\n');

    end

end

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

dist_val = strmatch(lower(s_dist_method), lower(distance_opts), 'exact');

s_dist_method = distance_opts{dist_val};

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

    fbcRegressCovariates(dfncInfo, regressCovX);

end

% -------------------------------------------------------------------------
% FBC-specific processing
% FBdyn is [bands x time x connections] in each output file.
% The clustering feature dimension is connections ONLY.
% -------------------------------------------------------------------------
if ~strcmpi(use_tall_array, 'no')
    error('FBC:UnsupportedTall', 'FBC version currently requires use_tall_array = no.');
end
if b_chk_ena_stateguided
    error('FBC:UnsupportedStateGuided', ...
        'State-guided ICA requires a band-aware implementation; disable it for FBC.');
end
if b_chk_ena_metastate
    error('FBC:UnsupportedMetaState', ...
        'Meta-state and TV-dFNC processing are not implemented for FBC; disable meta-state.');
end

M = numel(dfncInfo.outputFiles);
if M ~= dfncInfo.userInput.numOfSub * dfncInfo.userInput.numOfSess
    error('FBC:FileCount', 'outputFiles count must equal subjects * sessions.');
end
if M == 0
    error('FBC:NoFiles', 'dfncInfo.outputFiles is empty.');
end

needPeaks = strcmpi(estimate_clusters, 'yes') || ...
    strcmpi(kmeans_start, 'subject exemplars');
SP = cell(M, 1);
FBflat = cell(M, 1);
numBands = [];
numWindows = [];
numConn = [];

fprintf('Loading FBC and computing within-band spectra ...\n');
for nR = 1:M
    current_file = fullfile(outputDir, dfncInfo.outputFiles{nR});
    loaded = load(current_file, 'FBdyn');
    if ~isfield(loaded, 'FBdyn')
        error('FBC:MissingData', 'Missing FBdyn in %s', current_file);
    end
    FBdyn = loaded.FBdyn;
    if ndims(FBdyn) ~= 3
        error('FBC:Dimensions', 'FBdyn must be [bands x windows x connections] in %s', current_file);
    end
    [B, T, C] = size(FBdyn);
    if isempty(numBands)
        numBands = B;
        numWindows = T;
        numConn = C;
        Fs = 1/TR;
        nfft = 2^(nextpow2(T)+1);
        ffull = (0:nfft-1)' * (Fs/nfft);
        freqInd = find(ffull <= Fs/2);
        dfncInfo.freq = ffull(freqInd)';
    elseif B ~= numBands || T ~= numWindows || C ~= numConn
        error('FBC:InconsistentDimensions', 'Inconsistent FBdyn dimensions in %s', current_file);
    end
    if ~all(isfinite(FBdyn(:)))
        error('FBC:NonFinite', 'FBdyn has NaN/Inf values in %s', current_file);
    end

    % Each band contributes T independent observations of length C.
    % Row order: time changes fastest within a band, then band.
    Xr = reshape(permute(FBdyn, [2 1 3]), T*B, C);
    FBflat{nR} = Xr;

    if needPeaks
        peakBlocks = cell(B,1);
        for b = 1:B
            Xb = reshape(FBdyn(b,:,:), T, C);
            DEV = std(Xb, [], 2);
            [~, imax] = icatb_extrema(DEV);
            pIND = sort(imax(:));
            if isempty(pIND)
                [~, pIND] = max(DEV);
            end
            peakBlocks{b} = Xb(pIND, :);
        end
        SP{nR} = vertcat(peakBlocks{:});
    end

    % Spectra are calculated along TIME for each band independently.
    % No post-filterbank detrending or filtering is applied to clustering data.
    spectra_fbc = zeros(numel(freqInd), C, B);
    h = icatb_hamming(T);
    h = h(:);
    for b = 1:B
        Xb = reshape(FBdyn(b,:,:), T, C);
        Xb = bsxfun(@minus, Xb, mean(Xb,1));
        Xb = bsxfun(@times, Xb, h);
        Y = fft(Xb, nfft, 1) / Fs;
        spectra_fbc(:,:,b) = abs(Y(freqInd,:));
    end
    icatb_save(current_file, 'spectra_fbc', '-append');
    amp = reshape(std(spectra_fbc, 0, 1), C, B);
    cm = zeros(C,B);
    for b = 1:B
        Sb = spectra_fbc(:,:,b);
        denom = mean(Sb,1);
        numer = mean(bsxfun(@times, Sb, dfncInfo.freq(:)),1);
        cm(:,b) = (numer ./ max(denom,eps))';
    end
    if nR == 1
        FNCamp = zeros(C,B);
        FNCcm = zeros(C,B);
    end
    FNCamp = FNCamp + amp;
    FNCcm = FNCcm + cm;
end
FNCamp = FNCamp / M;
FNCcm = FNCcm / M;
post_process_file = fullfile(outputDir, [dfncInfo.prefix, '_post_process.mat']);
if exist(post_process_file,'file')
    icatb_save(post_process_file, '-nocompression', '-append', 'FNCamp', 'FNCcm');
else
    icatb_save(post_process_file, '-nocompression', 'FNCamp', 'FNCcm');
end

dfncInfo.postprocess.num_bands = numBands;
dfncInfo.postprocess.num_windows = numWindows;
dfncInfo.postprocess.num_connections = numConn;

if b_chk_ena_statebased
    fprintf('Clustering FBC across all subjects, sessions, times and bands ...\n');
    if needPeaks
        SPflat = vertcat(SP{:});
        if isempty(SPflat)
            error('FBC:NoExemplars', 'No subject exemplars were found.');
        end
    end
    if strcmpi(estimate_clusters, 'yes')
        cluster_estimate_results = icatb_optimal_clusters(SPflat, ...
            min([size(SPflat,1), 10]), 'method', 'all', ...
            'cluster_opts', {'Replicates', kmeans_num_replicates, ...
            'Distance', lower(s_dist_method), 'MaxIter', kmeans_max_iter}, ...
            'num_tests', num_tests_est_clusters, 'display', 1);
        K = zeros(numel(cluster_estimate_results),1);
        for j = 1:numel(K)
            K(j) = cluster_estimate_results{j}.K(1);
        end
        num_clusters = ceil(mean(K));
        dfncInfo.postprocess.num_clusters = num_clusters;
    end

    % Faghiri et al. use common centroids across all frequency bands.
    % The peak exemplars are an optional GIFT initialization, not part of
    % the published FBC clustering procedure.
    if strcmpi(kmeans_start, 'subject exemplars')
        try
            [~, Cp] = kmeans(SPflat, num_clusters, 'distance', lower(s_dist_method), ...
                'Replicates', kmeans_num_replicates, ...
                'MaxIter', kmeans_max_iter, 'EmptyAction', 'drop');
        catch
            [~, Cp] = icatb_kmeans(SPflat, num_clusters, 'distance', lower(s_dist_method), ...
                'Replicates', kmeans_num_replicates, ...
                'MaxIter', kmeans_max_iter, 'empty', 'drop');
        end
    elseif strcmpi(kmeans_start, 'default')
        Cp = [];
    elseif strcmpi(kmeans_start, 'user input')
        if ~exist('Cp','var') || isempty(Cp)
            error('FBC:MissingCentroids', 'User-input centroids Cp are missing.');
        end
    else
        error('FBC:KmeansStart', 'Unknown kmeans_start: %s', kmeans_start);
    end

    X = vertcat(FBflat{:});
    post_process_opts = struct('max_iter', kmeans_max_iter, ...
        'kmeans_num_replicates', kmeans_num_replicates, ...
        'kmeans_distance_method', lower(s_dist_method), 'Cp', Cp);
    if ~isempty(Cp)
        post_process_opts.kmeans_num_replicates = 1;
    end
    if ~backProject
        [IDXall, Call, SUMDall, Dall] = ...
            icatb_kmeans_clustering(X, num_clusters, post_process_opts);
    else
        [Call, IDXall, Dall, SUMDall] = fbcGetCenters(X, Cp);
    end
    IDXall = IDXall(:);
    if numel(IDXall) ~= M*T*B
        error('FBC:AssignmentCount', 'Unexpected number of kmeans labels.');
    end

    % params needed are fuzzyness param & 'sqeuclidean' & 'euclidean'
    U = zeros(size(Dall));    
    d_fuzzy = results.fuzzival;
    pow_transform = 1 / (d_fuzzy - 1); % For 'sqeuclidean' (MATLAB kmeans default):
    % Identify windows with zero distance to any centroid
    zeroRows = any(Dall == 0, 2);
    % Calculate memberships for nonzero-distance windows
    Dn = Dall(~zeroRows, :);
    % Numerically stable inverse-distance weighting
    logWeights = -pow_transform * log(Dn);
    logWeights = logWeights - max(logWeights, [], 2);
    weights = exp(logWeights);
    U(~zeroRows, :) = weights ./ sum(weights, 2);
    % Handle exact matches to one or more centroids
    Z = double(Dall(zeroRows, :) == 0);
    U(zeroRows, :) = Z ./ sum(Z, 2);

    clusterInfo = struct();
    clusterInfo.Dall_norm = U;
    clusterInfo.Call = Call;
    clusterInfo.Cp = Cp;
    clusterInfo.IDXall = IDXall;
    clusterInfo.SUMDall = SUMDall;
    clusterInfo.Dall = Dall;
    if needPeaks
        clusterInfo.SP = SPflat;
    end
    clusterInfo.numBands = B;
    clusterInfo.numWindows = T;
    clusterInfo.numConnections = C;
    clusterInfo.labelOrder = 'time x band x file (subject/session)';
    % statesFile: [time x band x file]; file order = dfncInfo.outputFiles.
    clusterInfo.statesFile = reshape(IDXall, T, B, M);
    % membershipFile: [file x band x cluster]
    membershipFile = zeros(M,B,num_clusters);
    fractionFile = zeros(M,B,num_clusters);
    for r = 1:M
        for b = 1:B
            labels = clusterInfo.statesFile(:,b,r);
            for k = 1:num_clusters
                membershipFile(r,b,k) = sum(labels == k);
            end
            fractionFile(r,b,:) = membershipFile(r,b,:) / T;
        end
    end
    clusterInfo.membershipFile = membershipFile;
    clusterInfo.fractionRateFile = fractionFile;
    clusterInfo.bandMembership = reshape(sum(membershipFile,1), B, num_clusters);
    clusterInfo.bandFractionRate = clusterInfo.bandMembership / (M*T);

    % Subject/session indexing is assigned from file names, not assumed
    % from a particular ordering of outputFiles.
    % When filenames follow GIFT's sub_XX_sess_YY convention, also expose
    % arrays indexed explicitly by subject and session.
    [fileSub, fileSess] = fbcFileIndices(dfncInfo.outputFiles, ...
        dfncInfo.userInput.numOfSub, dfncInfo.userInput.numOfSess);
    clusterInfo.fileSubject = fileSub;
    clusterInfo.fileSession = fileSess;
    clusterInfo.states = [];
    clusterInfo.membership = [];
    clusterInfo.fractionRate = [];
    if all(fileSub > 0) && all(fileSess > 0)
        clusterInfo.states = zeros(dfncInfo.userInput.numOfSub, ...
            dfncInfo.userInput.numOfSess, T, B);
        clusterInfo.membership = zeros(dfncInfo.userInput.numOfSub, ...
            dfncInfo.userInput.numOfSess, B, num_clusters);
        clusterInfo.fractionRate = zeros(size(clusterInfo.membership));
        for r = 1:M
            clusterInfo.states(fileSub(r),fileSess(r),:,:) = ...
                reshape(clusterInfo.statesFile(:,:,r),1,1,T,B);
            clusterInfo.membership(fileSub(r),fileSess(r),:,:) = ...
                reshape(membershipFile(r,:,:),1,1,B,num_clusters);
            clusterInfo.fractionRate(fileSub(r),fileSess(r),:,:) = ...
                reshape(fractionFile(r,:,:),1,1,B,num_clusters);
        end
    end
    clusterInfo.corrs_states = cell(M,B,num_clusters);
    for r = 1:M
        for b = 1:B
            Xb = FBflat{r}((b-1)*T+(1:T),:);
            for k = 1:num_clusters
                clusterInfo.corrs_states{r,b,k} = ...
                    Xb(clusterInfo.statesFile(:,b,r) == k,:);
            end
        end
    end
    clusterInfo.fileNames = dfncInfo.outputFiles;
    icatb_save(post_process_file, 'clusterInfo', '-append');
    if exist('cluster_estimate_results','var')
        icatb_save(post_process_file, 'cluster_estimate_results', '-append');
    end
end

dfncInfo.postprocess.dmethod = s_dist_method;
dfncInfo.postprocess.num_clusters = num_clusters;
dfncInfo.postprocess.fuzzival = results.fuzzival;
dfncInfo.postprocess.b_continuous_kmeans = results.b_continuous_kmeans;
param_file = fullfile(outputDir, [dfncInfo.prefix, '.mat']);
icatb_save(param_file, 'dfncInfo');
fprintf('FBC post-processing done.\n');
end

function fbcRegressCovariates(dfncInfo, regressCovX)
% Regress subject-level covariates from each band/time/connection element.
% Subject files must have identical FBdyn dimensions across subjects.
for nSess = 1:dfncInfo.userInput.numOfSess
    files = cell(dfncInfo.userInput.numOfSub,1);
    for nSub = 1:dfncInfo.userInput.numOfSub
        files{nSub} = fullfile(dfncInfo.outputDir, sprintf('%s_sub_%s_sess_%s_results.mat', ...
            dfncInfo.userInput.prefix, icatb_returnFileIndex(nSub), ...
            icatb_returnFileIndex(nSess)));
        q = load(files{nSub}, 'FBdyn');
        if ~isfield(q,'FBdyn')
            error('FBC:MissingData', 'Missing FBdyn in %s', files{nSub});
        end
        if nSub == 1
            sz = size(q.FBdyn);
            Y = zeros(dfncInfo.userInput.numOfSub, numel(q.FBdyn));
        elseif ~isequal(size(q.FBdyn),sz)
            error('FBC:Dimensions', 'Inconsistent dimensions in %s', files{nSub});
        end
        Y(nSub,:) = q.FBdyn(:)';
    end
    Y = Y - regressCovX*(pinv(regressCovX)*Y);
    for nSub = 1:dfncInfo.userInput.numOfSub
        q = load(files{nSub}, 'FBdyn');
        FBdyn_old = q.FBdyn;
        FBdyn = reshape(Y(nSub,:), sz);
        icatb_save(files{nSub}, 'FBdyn', 'FBdyn_old', '-append');
    end
end
end

function [Cn, IDXall, Dall, SumD] = fbcGetCenters(X, C)
% Retain the original GIFT centroid back-projection formulation.
xx = pinv(C') * X';
Cn = pinv(xx') * X;
Dall = zeros(size(X,1),size(Cn,1));
for k = 1:size(Cn,1)
    Dall(:,k) = sqrt(sum(bsxfun(@minus,X,Cn(k,:)).^2,2));
end
[~,IDXall] = min(Dall,[],2);
SumD = [];
end

function [sub, sess] = fbcFileIndices(names, nSub, nSess)
% Extract explicit indices from GIFT-style result filenames.
sub = zeros(numel(names),1);
sess = zeros(numel(names),1);
for r = 1:numel(names)
    tok = regexp(names{r}, '_sub_(\d+)_sess_(\d+)_results\.mat$', 'tokens', 'once');
    if isempty(tok)
        sub(:) = 0;
        sess(:) = 0;
        return;
    end
    sub(r) = str2double(tok{1});
    sess(r) = str2double(tok{2});
end
if any(sub < 1 | sub > nSub | sess < 1 | sess > nSess) || ...
        numel(unique(sub + nSub*(sess-1))) ~= numel(names)
    sub(:) = 0;
    sess(:) = 0;
end
end
