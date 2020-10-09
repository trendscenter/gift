function icatb_spatial_chronn_stats(schronnInfo)
%% Compute statistics on the transition matrices and subject states
%

outputDir = schronnInfo.outputDir;
results_files = schronnInfo.postprocess.results_files;
statsInfo = schronnInfo.postprocess.statsInfo;

statsDir = [schronnInfo.userInput.prefix, '_stats'];
if (exist(fullfile(schronnInfo.outputDir, statsDir), 'dir') ~= 7)
    mkdir(schronnInfo.outputDir, statsDir);
end

schronnInfo.postprocess.statsInfo = statsInfo;


try
    statsInfo = getDesignMat(statsInfo);
catch
end

stats_files = cell(1, length(results_files));

% loop over results files
for nR = 1:length(results_files)
    
    in_file = fullfile(outputDir, results_files{nR});
    disp(['Loading file file ', in_file]);
    load(in_file);
    
    
    % spatial transition matrix
    sm = clusterInfo.spatial_trans_stats.spatial_transition_matrix;
    sm_dims = size(sm);
    sm = reshape(sm, schronnInfo.userInput.numOfSess, schronnInfo.userInput.numOfSub, size(sm, 2)*size(sm, 3));
    sm(abs(sm) < eps) = NaN;
    sm = squeeze(icatb_nanmean(sm, 1));
    
    spatial_trans_matrix_stats = run_model(sm, statsInfo);
    spatial_trans_matrix_stats.dims = sm_dims;
    
    % dwell time
    mdwt = clusterInfo.state_vector_stats.mean_dwell_time;
    mdwt(abs(mdwt) < eps) = NaN;
    mdwt = reshape(mdwt, schronnInfo.userInput.numOfSess, schronnInfo.userInput.numOfSub, size(mdwt, 2));
    mdwt = squeeze(icatb_nanmean(mdwt, 1));
    mean_dwell_time_stats = run_model(mdwt, statsInfo);
    
    % subject clusters
    subject_cluster_stats = cell(1, length(clusterInfo.subject_cluster_files));
    for nState = 1:length(clusterInfo.subject_cluster_files)
        
        tmpf = fullfile(schronnInfo.outputDir, clusterInfo.subject_cluster_files{nState});
        tmp_dat = icatb_read_data(tmpf, [], schronnInfo.mask_ind)';
        tmp_dat(abs(tmp_dat) < eps) = NaN;
        tmp_dat = reshape(tmp_dat, schronnInfo.userInput.numOfSess, schronnInfo.userInput.numOfSub, size(tmp_dat, 2));
        tmp_dat = squeeze(icatb_nanmean(tmp_dat, 1));
        subject_cluster_stats{nState} = run_model(tmp_dat, statsInfo);
        
    end
    
    
    [~, statsfN, extn] = fileparts(results_files{nR});
    out_file = fullfile(statsDir, [statsfN, '_stats', extn]);
    disp(['... saving file ', fullfile(schronnInfo.outputDir, out_file)]);
    save(fullfile(schronnInfo.outputDir, out_file), 'spatial_trans_matrix_stats', 'mean_dwell_time_stats', 'subject_cluster_stats');
    stats_files{nR} = out_file;
    
end


%% Save params
param_file = fullfile(schronnInfo.outputDir, [schronnInfo.prefix, '.mat']);
schronnInfo.postprocess.stats_files = stats_files;
save(param_file, 'schronnInfo');

disp('Done');
fprintf('\n');


%% write out results
resultsFile = 'icatb_display_stats_spatial_chronnectome';
outDir = fullfile(schronnInfo.outputDir, [schronnInfo.prefix, '_stats_display']);
if (exist(outDir, 'dir') ~= 7)
    mkdir (outDir);
end
assignin('base', 'schronnInfo', schronnInfo);
opts.codeToEvaluate = 'icatb_display_stats_spatial_chronnectome(schronnInfo);';
opts.outputDir = outDir;
opts.showCode = false;
opts.useNewFigure = false;
opts.format = 'html';
opts.createThumbnail = true;
if (strcmpi(opts.format, 'pdf'))
    opts.useNewFigure = false;
end

publish(resultsFile, opts);

close all;


function statsInfo = getDesignMat(statsInfo)
%% Get design matrix
%

covInfo = statsInfo.cov;
numOfSub = statsInfo.numOfSub;

design_matrix = zeros(numOfSub, length(covInfo));
%chkNan = double(isnan(design_matrix));

for i = 1:length(covInfo)
    cval = covInfo(i).value;
    ctype = covInfo(i).type;
    if (strcmpi(ctype, 'categorical'))
        %% Convert Categorical to numeric value
        if (isnumeric(cval))
            cval = cval(:);
            cval = num2str(cval);
        end
        cval = cellstr(cval);
        cval = cval(:)';
        
        unique_cval = lower(unique(cval));
        
        cvalue = zeros(length(cval), 1);
        
        for nU = 1:length(unique_cval)
            if (~strcmpi(unique_cval{nU}, 'nan'))
                tmpVal = nU-1;
            else
                tmpVal = NaN;
            end
            cvalue(strcmpi(cval, unique_cval{nU})) = tmpVal;
        end
        
        covInfo(i).labels = {unique_cval};
        
    else
        %% Continuous variables
        if (~isnumeric(cval))
            cvalue = str2num(str2mat(cval));
        else
            cvalue = cval;
        end
        
        cvalue = cvalue + eps;
        
        if (~isempty(covInfo(i).transformation))
            eval(['cvalue =',  covInfo(i).transformation, '(cvalue);']);
        end
        %        cvalue = detrend(cvalue(:), 0);
        
        covInfo(i).labels = covInfo(i).name;
        covInfo(i).labels = {[covInfo(i).transformation, '(', covInfo(i).labels, ')']};
    end
    
    if (i == 1)
        chkNan = double(isnan(cvalue));
    else
        chkNan = chkNan | double(isnan(cvalue));
    end
    
    design_matrix(:, i) = cvalue;
    
end

chkNan = (chkNan == 1);
good_sub_inds = (1:numOfSub);


%% Include good subjects only
%if (~isempty(chkNan))

good_sub_inds(chkNan) = [];
design_matrix(chkNan, :) = [];

chkContinuous = find(strcmpi('categorical', cellstr(char(covInfo.type))) == 0);

if (~isempty(chkContinuous))
    design_matrix(:, chkContinuous) = detrend(design_matrix(:, chkContinuous), 0);
end

for i = 1:length(covInfo)
    covInfo(i).value(chkNan) = [];
end

% Initial design matrix
statsInfo.X0 = design_matrix;

Stepwise_options = {};

%% Design including covariates
group_ind = strmatch('categorical', lower(cellstr(char(covInfo.type))), 'exact');
cov_ind = strmatch('continuous', lower(cellstr(char(covInfo.type))), 'exact');

group_ind = group_ind(:)';
cov_ind = cov_ind(:)';

groups = design_matrix(:, group_ind);
covariates = design_matrix(:, cov_ind);

terms = { 0 };
X     = ones(size(design_matrix, 1), 1);

%% Categorical covariates
if (~isempty(groups))
    [ x, t ] = mG2X(groups, 0, Stepwise_options);
    terms    = cat(2, terms, t);
    X        = cat(2, X, x);
end
%% Continuous covariates
if (~isempty(covariates))
    [ x, t ] = mC2X(covariates, size(groups, 2), Stepwise_options);
    terms    = cat(2, terms, t);
    X        = cat(2, X, x);
end

cov_names = cellstr(char(covInfo.name));
onames = cov_names([group_ind, cov_ind]);

for ii = 2:length(terms) %skip the constant term
    if length(terms{ii}) == 1
        %main effect
        term_names{ii-1} =  onames{terms{ii}};
    else
        %interaction term
        term_names{ii-1} =  [onames{terms{ii}(1)} '_X_' onames{terms{ii}(2)}];
    end
end

statsInfo.terms = terms;
statsInfo.X = X;
statsInfo.cov = covInfo;
statsInfo.good_sub_inds = good_sub_inds;
statsInfo.regressors = term_names;


function s = get_contrast_label(s, tname, allnames, alllabels, addTag)
%% Get contrast labels
%

if (~exist('addTag', 'var'))
    addTag = 0;
end

prefixV = '';
if (addTag)
    prefixV = [tname, '_'];
end

if length(s.Term) == 1 %main effect
    varIND = find(strcmp(tname, allnames));
    labels = alllabels{varIND};
    if (~iscell(labels))
        labels = {labels};
    end
    for jj = 1:size(s.Levels,1)
        if s.Levels(jj,2) == 0 && length(labels) == 1
            s.Contrast{jj} = [prefixV, labels{1}];
        else
            s.Contrast{jj} = [prefixV, '(' labels{s.Levels(jj,1)+1} ') - (' labels{s.Levels(jj,2)+1} ')'];
            
        end
    end
    
else %interaction
    term1_end = strfind(tname, '_X_')-1;
    term2_start = term1_end + 4;
    clear varIND
    varIND(1) = find(strcmp(tname(1:term1_end), allnames));
    varIND(2) = find(strcmp(tname(term2_start:end), allnames));
    labels_1 = alllabels{varIND(1)};
    labels_2 = alllabels{varIND(2)};
    
    if (~iscell(labels_1))
        labels_1 = {labels_1};
    end
    
    if (~iscell(labels_2))
        labels_2 = {labels_2};
    end
    
    
    for jj = 1:size(s.Levels,1)
        if length(labels_1)*length(labels_2) == 1 % both are continuous variables
            s.Contrast{jj} = [prefixV, '(' labels_1{1} ') X (' labels_2{1} ')'];
        elseif length(labels_1) == 1 || length(labels_2) == 1 %one continuous, one categorical
            
            if length(labels_1) > length(labels_2)
                temp = labels_1;
                labels_1 = labels_2;
                labels_2 = temp;
            end
            
            s.Contrast{jj} = [prefixV, '(' labels_1{1} ') X [(' labels_2{s.Levels(jj,1)+1} ') - (' labels_2{s.Levels(jj,2)+1} ')]'];
        elseif length(labels_1)*length(labels_2) == 4 %categorical, two by two
            s.Contrast{jj} = [prefixV, '[(' labels_1{s.Levels(jj,1)+1} ') - (' labels_1{s.Levels(jj,2)+1} ')] X [(' ...
                labels_2{s.Levels(jj,1)+1} ') - (' labels_2{s.Levels(jj,2)+1} ')]'];
        elseif length(labels_1)*length(labels_2) == 6 %categorical, two by three
            if length(labels_1) > length(labels_2)
                temp = labels_1;
                labels_1 = labels_2;
                labels_2 = temp;
            end
            s.Contrast{jj} = [prefixV, '[(' labels_1{2} ') - (' labels_1{1} ')] X [(' labels_2{s.Levels(jj,1)+1} ') - (' labels_2{s.Levels(jj,2)+1} ')]'];
            
        else
            s.Contrast{jj} = tname;
            %             combinations = [];
            %             for j = 1 : length(labels_1)
            %                 for k = 1 :length(labels_2)
            %                     combinations(end + 1, :) = [ j k ];
            %                 end
            %             end
        end
    end
    
end



function UNI = run_model(data, statsInfo)


univariate_tests = statsInfo.univariate_tests;
test_names = {};
try
    test_names = univariate_tests(:,1);
    cov_names = cellstr(char(statsInfo.cov.name));
    cov_labels = [statsInfo.cov.labels];
    terms = statsInfo.terms;
catch
end

t_u = cell(1, length(test_names) + 1);
p_u = t_u;
stats_u = t_u;


[ t_u{1}, p_u{1}, stats_u{1}] = icatb_nan_mT(data, ones(size(data, 1), 1), [], 0, {'verbose'});
stats_u{1} = get_contrast_label(stats_u{1}, {'one_sample_ttest'}, {'one_sample_ttest'}, {'ttest'});

for nUnivInfo = 1:length(test_names)
    
    %[dd, termNo] = intersect(mancovanInfo.regressors, univInfo(nUnivInfo).name);
    %termNo = terms{termNo};
    
    termNo = ismember(statsInfo.regressors, test_names{nUnivInfo});
    termNo = find(termNo == 1);
    termNo = termNo(:)';
    
    ia = termNo;
    try
        other_regressors = cellstr(char(univariate_tests{nUnivInfo, 2}));
        %[dd, ia] = intersect(mancovanInfo.regressors, other_regressors);
        ia = ismember(statsInfo.regressors, other_regressors);
        ia = find(ia == 1);
        ia = ia(:)';
        ia = sort([termNo, ia]);
    catch
    end
    
    X_reduced = statsInfo.X(:, [1, ia + 1]);
    allTerms = terms([1, ia + 1]);
    %termNo = terms{termNo + 1};
    termNo = unique(cell2mat(terms(termNo + 1)));
    
    fprintf('Working on term %d of %d\n', nUnivInfo, length(test_names))
    [ t_u{nUnivInfo+1}, p_u{nUnivInfo+1}, stats_u{nUnivInfo+1}] = icatb_nan_mT(data, X_reduced, allTerms, termNo, {'verbose'});
    addTag = 1;
    chkNames = strmatch(lower(test_names{nUnivInfo}), lower(cov_names), 'exact');
    if (~isempty(chkNames))
        if (strcmpi(statsInfo.cov(chkNames(1)).type, 'continuous'))
            addTag = 0;
        end
    end
    stats_u{nUnivInfo+1} = get_contrast_label(stats_u{nUnivInfo+1}, test_names{nUnivInfo}, cov_names, cov_labels, addTag);
    
end

tests = cell(length(test_names) + 1, 1);
tests{1} = 'ttest';
if (~isempty(test_names))
    tests(2:end) = test_names;
end

UNI.t_u = t_u;
UNI.p_u = p_u;
UNI.stats_u = stats_u;
UNI.tests = tests;