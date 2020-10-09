function [UNI_within, UNI_Tp, UNI_between, UNI_contrasts] = icatb_ranova(data, covariates, varargin)
%% Do repeated measures anova and report within subject effects followed by anova for between subject effects and do multiple comparisons on the catgeorical covariates
%
% Inputs
% 1. data - subjects x measures
% 2. covariates - data stucture containing fields like name, val and type
% for example covN(1).name = 'Gender'; covN(1).val = {'Male', 'Female', 'Female', 'Male'}; covN(1).type
% = 'categorical';
%
% vars - variable names followed by categorical values
% within - Withinsubject measures
%
% Outputs:
%
%
% Example:
%load repeatedmeas
% Age = (table2array(between(:,1)))
% Group = cellstr(table2array(between(:,3)))
% IQ = (table2array(between(:,2)))
% Gender = cellstr(table2array(between(:,4)))
% data = table2array(between(:,5:end));
% covN(1).name = 'Age';
% covN(1).val = Age;
% covN(1).type='continuous';
% covN(2).type='categorical';
% covN(2).name='IQ';
% covN(2).val=IQ
% covN(2).type='continuous';
% covN(3).type='categorical'
% covN(3).name='Group'
% covN(3).val=Group
% covN(4).val=Gender
% covN(4).name='Gender'
% covN(4).type='categorical'
% [UNI_within, UNI_Tp, UNI_between, UNI_contrasts] = icatb_ranova(data,covN, 'within',within);
%

UNI_Tp = [];
UNI_within = [];
UNI_between = [];
UNI_contrasts = [];

interactions = 1;

for n = 1:2:length(varargin)
    if (strcmpi(varargin{n}, 'within'))
        within = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'interactions'))
        interactions = varargin{n + 1};
    end
end

if (~exist('covariates', 'var'))
    covariates = [];
end

if (~exist('within', 'var'))
    within = table(repmat('A', size(data, 2), 1), num2str((1:size(data, 2))'), 'VariableNames', {'w1', 'w2'});
end


tableStruct = cell(1, size(data, 2) + length(covariates) + 2);

for n = 1:size(data, 2)
    tableStruct{n} = data(:, n);
end


for n = 1:length(covariates)
    tableStruct{n + size(data, 2)} = covariates(n).val;
end

yVar = cellstr(strcat('y', num2str((1:size(data, 2))')))';

tableStruct{size(data, 2) + length(covariates) + 1} = 'VariableNames';
if (~isempty(covariates))
    tableStruct{size(data, 2) + length(covariates) + 2} = [yVar, cellstr(char(covariates.name))'] ;
else
    tableStruct{size(data, 2) + length(covariates) + 2} = yVar ;
end

% Interaction term
t = table(tableStruct{:});
%timeMeasurements = table(within,'VariableNames',{'Time'});

timeMeasurements = within;

formulaY = ['y1-y', num2str(size(data, 2))];

if (~isempty(covariates))
    
    categ_covType = cellstr(char(covariates.type));
    inds = (strcmpi(categ_covType, 'categorical'));
    categ_cov = covariates(inds);
    categ_cov = cellstr(char(categ_cov.name));
    inds = find(inds==1);
    
    for nCov = 1:length(inds)
        tmpCategVal = categorical(covariates(inds(nCov)).val);
        if (isnumeric(tmpCategVal))
            tmpCategVal = cellstr(num2str(tmpCategVal(:)));
        end
        covariates(inds(nCov)).val = tmpCategVal;
    end
    
    inds = (strcmpi(categ_covType, 'continuous'));
    cont_cov = covariates(inds);
    cont_cov = cellstr(char(cont_cov.name));
    
    if (interactions)
        
        if (length(categ_cov) > 1)
            combinationsToChoose = nchoosek(1:length(categ_cov), 2);
            inter_terms = cell(size(combinationsToChoose, 1), 1);
            for n = 1:length(inter_terms)
                inter_terms{n} = [categ_cov{combinationsToChoose(n, 1)}, ':', categ_cov{combinationsToChoose(n, 2)}];
            end
            
            allTerms = [inter_terms;categ_cov;cont_cov];
        else
            allTerms = [categ_cov;cont_cov];
        end
        
    else
        
        allTerms = [categ_cov;cont_cov];
        
    end
    
    goodInds = icatb_good_cells(allTerms);
    allTerms = allTerms(goodInds);
    
    formulaStr = allTerms{1};
    for n = 2:length(allTerms)
        formulaStr = [formulaStr, '+', allTerms{n}];
    end
    
else
    
    formulaStr = '1';
    
end


%% Within subjects anova

rm = fitrm(t, [formulaY, '~', formulaStr],'WithinDesign',timeMeasurements);

tbl_time = ranova(rm);

if (~isempty(covariates))
    errorMsq = table2array(tbl_time(end, 3));
    errorDF = table2array(tbl_time(end, 2));
    tbl_time = tbl_time(2:end-1,:);
    rm2 = fitrm(t, [formulaY, '~1'],'WithinDesign',timeMeasurements);
    tbl_time2 = ranova(rm2);
    tbl_time2 = tbl_time2(1,:);
    timeMsq =  table2array(tbl_time2(1,3));
    fTime = timeMsq/errorMsq;
    pTime = 1 - fcdf(fTime, table2array(tbl_time2(1, 2)),  errorDF);
    tbl_time2.F = fTime;
    tbl_time2.pValue = pTime;
    tbl_time2.Properties.RowNames={'Time'};
    tbl_time = [tbl_time2; tbl_time];
else
    tbl_time = tbl_time(1,:);
    tbl_time.Properties.RowNames={'Time'};
end


for n = 1:size(tbl_time, 1)
    
    UNI_within(end+1).title = tbl_time.Properties.RowNames{n};
    UNI_within(end).p = tbl_time.pValue(n);
    UNI_within(end).F = tbl_time.F(n);
    
end


timeMeasurements2 = table2array(timeMeasurements);


tps = ((unique((timeMeasurements2(:, 2)))));
numTimePoints = length(tps);


timeMeasurements2 = cellstr(timeMeasurements2(:, 1));
[~, inds] = unique(timeMeasurements2);
inds = sort(inds);
conds = timeMeasurements2(inds);

if (length(conds) >= 2)
    
    comb_conds = nchoosek(1:length(conds), 2);
    
    UNI_Tp = repmat(struct('p', [], 't', []), 1, size(comb_conds, 1)*(1 + numTimePoints));
    countTp = 0;
    for n = 1:size(comb_conds, 1)
        
        condA = conds{comb_conds(n, 1)};
        condB = conds{comb_conds(n, 2)};
        
        condAInds = strmatch(condA, timeMeasurements2, 'exact');
        condBInds = strmatch(condB, timeMeasurements2, 'exact');
        
        
        for nTps = 1:length(condAInds)
            
            countTp = countTp + 1;
            x = data(:, condAInds(nTps));
            y = data(:, condBInds(nTps));
            
            titleStr = [condA, ' - ', condB, ' (Timepoint:', char(tps(nTps)), ')'];
            [h, pTp, c, statsTp] = ttest(x - y);
            
            UNI_Tp(countTp).p = pTp;
            UNI_Tp(countTp).t = statsTp.tstat;
            UNI_Tp(countTp).title = titleStr;
            UNI_Tp(countTp).tp = tps(nTps);
            
        end
        
        countTp = countTp + 1;
        x = sum(data(:, condAInds), 2);
        y = sum(data(:, condBInds), 2);
        [h, pTp, c, statsTp] = ttest(x - y);
        UNI_Tp(countTp).p = pTp;
        UNI_Tp(countTp).t = statsTp.tstat;
        titleStr = [condA, ' - ', condB, ' (Average Time)'];
        UNI_Tp(countTp).title = titleStr;
        
    end
    
    
end



%% Between subjects anova
if (~isempty(covariates))
    tbl_between = anova(rm);
    tbl_between = tbl_between(2:end-1, :);
    
    for n = 1:size(tbl_between, 1)
        
        UNI_between(end+1).title = char(tbl_between.Between(n));
        UNI_between(end).p = (tbl_between.pValue(n));
        UNI_between(end).F = tbl_between.F(n);
        
    end
    
    if (exist('categ_cov', 'var'))
        
        if (length(categ_cov) > 1)
            for n = 1:length(categ_cov)
                T = multcompare(rm, categ_cov{n});
                dd = getUniqueComparisons(T, categ_cov{n});
                UNI_contrasts = [UNI_contrasts, dd];
            end
        end
        
    end
    
    
end




function UNI = getUniqueComparisons(T, comparisonName)

A = T.([comparisonName, '_1']);
B = T.([comparisonName, '_2']);

inds = strcmpi(A, A{1});
B2 = B(inds);
A2 = [A(1);B2];

chkA = {};
for n = 1:length(A2)-1
    for m = n+1:length(A2)
        chkA{end + 1} = [A2{n}, '-', A2{m}];
    end
end

contrastNames = strcat(A, '-', B);

[~, inds] = intersect(contrastNames, chkA);

T = T(inds, :);

contrastNames = contrastNames(inds);
Tvalues = T.Difference./T.StdErr;
pValues = T.pValue;

UNI = [];

for n = 1:length(contrastNames)
    UNI(end+1).title = contrastNames{n};
    UNI(end).p = pValues(n);
    UNI(end).t = Tvalues(n);
end






