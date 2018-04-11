function mancovanInfo = icatb_mancovan_full_design(mancovanInfo, interactions)
%%  Get full design including interaction terms
%

desCriteria = 'mancova';

try
    desCriteria = mancovanInfo.userInput.designCriteria;
catch
end

if (~strcmpi(desCriteria, 'mancova'))
    
    mancovanInfo.prefix = mancovanInfo.userInput.prefix;
    mancovanInfo.outputDir = mancovanInfo.userInput.outputDir;
    mancovanInfo.numOfSub = mancovanInfo.userInput.numOfSub;
    mancovanInfo.numOfSess = mancovanInfo.userInput.numOfSess;
    mancovanInfo.designCriteria = mancovanInfo.userInput.designCriteria;
    mancovanInfo.ttestOpts =  mancovanInfo.userInput.ttestOpts;
    
    if (strcmpi(desCriteria, 'one sample t-test'))
        X = ones(length(mancovanInfo.userInput.ttestOpts.t.val{1}), 1);
    elseif (strcmpi(desCriteria, 'two sample t-test'))
        N1 = length(mancovanInfo.userInput.ttestOpts.t.val{1});
        N2 = length(mancovanInfo.userInput.ttestOpts.t.val{2});
        X = ones(N1 + N2, 1);
        X(N1 + 1:end) = 0;
    else
        N1 = length(mancovanInfo.userInput.ttestOpts.t.val{1});
        N2 = length(mancovanInfo.userInput.ttestOpts.t.val{2});
        X = ones(N1, 1);
    end
    mancovanInfo.X = X;
    return;
end


% Check covariates
for nC = 1:length(mancovanInfo.userInput.cov)
    name = mancovanInfo.userInput.cov(nC).name;
    val = mancovanInfo.userInput.cov(nC).value;
    if (isempty(name) && isempty(val))
        error('Please specify name and value for each covariate in setup analysis');
    end
end

%% Model interactions
if (~exist('interactions', 'var'))
    mancovanInfo = icatb_mancovan_interactions(mancovanInfo);
else
    mancovanInfo = icatb_mancovan_interactions(mancovanInfo, interactions);
end


mancovanInfo.prefix = mancovanInfo.userInput.prefix;
mancovanInfo.outputDir = mancovanInfo.userInput.outputDir;
mancovanInfo.numOfSub = mancovanInfo.userInput.numOfSub;
mancovanInfo.cov = mancovanInfo.userInput.cov;
mancovanInfo.modelInteractions = mancovanInfo.userInput.modelInteractions;

design_matrix = zeros(mancovanInfo.numOfSub, length(mancovanInfo.cov));

for i = 1:length(mancovanInfo.cov)
    cval = mancovanInfo.cov(i).value;
    ctype = mancovanInfo.cov(i).type;
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
        
        mancovanInfo.cov(i).labels = {unique_cval};
        
    else
        %% Continuous variables
        if (~isnumeric(cval))
            cvalue = str2num(str2mat(cval));
        else
            cvalue = cval;
        end
        
        cvalue = cvalue + eps;
        
        if (~isempty(mancovanInfo.cov(i).transformation))
            eval(['cvalue =',  mancovanInfo.cov(i).transformation, '(cvalue);']);
        end
        %        cvalue = detrend(cvalue(:), 0);
        
        mancovanInfo.cov(i).labels = mancovanInfo.cov(i).name;
        mancovanInfo.cov(i).labels = {[mancovanInfo.cov(i).transformation, '(', mancovanInfo.cov(i).labels, ')']};
    end
    
    if (i == 1)
        chkNan = double(isnan(cvalue));
    else
        chkNan = chkNan | double(isnan(cvalue));
    end
    
    design_matrix(:, i) = cvalue;
    
end

chkNan = (chkNan == 1);
good_sub_inds = (1:mancovanInfo.numOfSub);


%% Include good subjects only
%if (~isempty(chkNan))

good_sub_inds(chkNan) = [];
design_matrix(chkNan, :) = [];

chkContinuous = find(strcmpi('categorical', cellstr(char(mancovanInfo.cov.type))) == 0);

if (~isempty(chkContinuous))
    design_matrix(:, chkContinuous) = detrend(design_matrix(:, chkContinuous), 0);
end

numOfSub = size(design_matrix, 1);
mancovanInfo.numOfSub = numOfSub;

for i = 1:length(mancovanInfo.cov)
    mancovanInfo.cov(i).value(chkNan) = [];
end

%end

mancovanInfo.good_sub_inds = good_sub_inds;

% Initial design matrix
mancovanInfo.X0 = design_matrix;

Stepwise_options = mancovanInfo.modelInteractions.types;

%% Design including covariates
group_ind = strmatch('categorical', lower(cellstr(char(mancovanInfo.cov.type))), 'exact');
cov_ind = strmatch('continuous', lower(cellstr(char(mancovanInfo.cov.type))), 'exact');

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

gglist = mancovanInfo.modelInteractions.gglist;
cclist = mancovanInfo.modelInteractions.cclist;
gclist = mancovanInfo.modelInteractions.gclist;

%% Categorical-Categorical interactions
if (size(groups, 2) > 1 && ~isempty(strmatch('group-group', Stepwise_options, 'exact')))
    [ x, t ] = mGG2X(groups, 0, Stepwise_options, gglist);
    terms    = cat(2, terms, t);
    X        = cat(2, X, x);
end

%% Continuous-Continuous interactions
if (size(covariates, 2) > 1 && ~isempty(strmatch('covariate-covariate', Stepwise_options, 'exact')))
    [ x, t ] = mCC2X(covariates, size(groups, 2), Stepwise_options, cclist);
    terms    = cat(2, terms, t);
    X        = cat(2, X, x);
end

%% Categorical-Continuous interactions
if (~isempty(groups) && ~isempty(covariates) && ~isempty(strmatch('group-covariate', Stepwise_options, 'exact')))
    [ x, t ] = mGC2X(groups, covariates, 0, size(groups, 2), Stepwise_options, gclist);
    terms    = cat(2, terms, t);
    X        = cat(2, X, x);
end

if (rcond(X'*X) < eps)
    error('Design matrix is highly ill-conditioned. Please remove redundant covariates');
end

%% Create labels for the columns
cov_names = cellstr(char(mancovanInfo.cov.name));
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


%% Variance inflation factor and leverage
X2 = detrend(X(:, 2:end), 0);
vif = diag(inv(icatb_corr(X2)));
leverage_design = diag(X2*inv(X2'*X2)*X2');

fprintf('\n');
disp('Variance inflation factor of covariates is/are: ');
disp([char(term_names), repmat(': ', length(term_names), 1), num2str(vif(:), '%0.6f')]);
fprintf('\n');

fprintf('\n');
if (any(vif > 5))
    warning('Variance inflation factor is greater than 5 in one or more covariates. Multi collinearity is high in the design matrix.');
end
fprintf('\n');
mancovanInfo.vif = vif;
mancovanInfo.leverage = leverage_design;
mancovanInfo.regressors = term_names;
mancovanInfo.X = X;
mancovanInfo.terms = terms;

