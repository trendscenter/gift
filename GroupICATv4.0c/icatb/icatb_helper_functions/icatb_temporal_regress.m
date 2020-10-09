function varargout = icatb_temporal_regress(param_file, selectedRegressors, displayFormat)
%% Compute multiple linear regression
%
% Inputs:
% 1. param_file - ICA parameter file (*ica*param*mat)
% 2. selectedRegressors - Selected regressor names

icatb_defaults;
global DETRENDNUMBER;
global PARAMETER_INFO_MAT_FILE;

filterP = ['*', PARAMETER_INFO_MAT_FILE, '*.mat'];

if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Parameter File', 'filter', filterP);
end

drawnow;

load(param_file);

outputDir = fileparts(param_file);

if isempty(outputDir)
    outputDir = pwd;
end

sesInfo.outputDir = outputDir;

if (isempty(sesInfo.userInput.designMatrix.name))
    tmpSPMFiles = icatb_selectEntry('typeEntity', 'file', 'title', 'Select SPM design matrix/matrices', 'filter', 'SPM.mat', 'typeSelection', 'multiple');
    if (isempty(tmpSPMFiles))
        error('Design matrix/matrices is/are not selected');
    end
    sesInfo.userInput.designMatrix.name = tmpSPMFiles;
    drawnow;
    save(param_file, 'sesInfo');
end

spm_files = cellstr(char(sesInfo.userInput.designMatrix.name));

if (length(spm_files) > sesInfo.numOfSub)
    error('Number of design matrix/matrices should not exceed the number of subjects');
end


drawnow;


regressor_names = cell(1, length(spm_files));
for nF = 1:length(spm_files)
    
    load(spm_files{nF});
    names = strtrim(regexprep(cellstr(char(SPM.xX.name)),'Sn\(\d+\)',''));
    regressor_names{nF} = char(lower(names));
    
end


regressor_names = cellstr(char(regressor_names));
[dd, inds] = unique(regressor_names);
regressor_names = regressor_names(sort(inds));

if (~exist('selectedRegressors', 'var'))
    dd = icatb_listdlg('promptstring', 'Select regressors of interest', 'liststring', regressor_names, 'selectionmode', 'multiple');
    selectedRegressors = regressor_names(dd);
end

drawnow;

if (isempty(selectedRegressors))
    error('Regressors are not selected');
end

if (isnumeric(selectedRegressors))
    selectedRegressors = regressor_names(selectedRegressors);
end

selectedRegressors = lower(selectedRegressors);

if (~exist('displayFormat', 'var'))
    displayFormat = questdlg('Which format do you want to use to display results?', 'Output Format', 'HTML', 'PDF', 'None', 'HTML');
end



betas = repmat(NaN, [sesInfo.numOfSess, sesInfo.numOfSub, sesInfo.numComp, length(selectedRegressors)]);
partial_correlations = zeros(sesInfo.numOfSess, sesInfo.numOfSub, sesInfo.numComp, length(selectedRegressors));

if (length(spm_files) == 1)
    spmInfo = load (deblank(spm_files{1}));
end

%% Fit model
ssr = zeros(1, sesInfo.numComp);
sse = zeros(1, sesInfo.numComp);

disp('Computing beta weights ...');

for nSub = 1:sesInfo.numOfSub
    
    if (length(spm_files) > 1)
        spmInfo = load (deblank(spm_files{nSub}));
    end
    
    for nSess = 1:sesInfo.numOfSess
        
        disp(['Subject ', num2str(nSub), ' Session ', num2str(nSess), ' ...']);
        
        TC = icatb_loadComp(sesInfo, (1:sesInfo.numComp), 'subjects', nSub, 'sessions', nSess, 'detrend_no', 0);
        TC = squeeze(TC);
        [modelX, model_inds] = getModelTimeCourse(TC, spmInfo, selectedRegressors, nSess);
        if (~isempty(model_inds))
            modelX = detrend(modelX(:, model_inds), 0);
            modelX = icatb_modelX(modelX, size(modelX, 1), DETRENDNUMBER);
            for nComp = 1:size(TC, 2)
                [tmpB, r2, residual] = icatb_regress(TC(:, nComp), modelX);
                betas(nSess, nSub, nComp, model_inds) = tmpB(1:length(model_inds));
                TC(:, nComp) = TC(:, nComp) - modelX(:, length(model_inds)+1:end)*tmpB(length(model_inds)+1:end);
                
                ssr(1, nComp) = ssr(1, nComp) + sum(residual.^2);
                sse(1, nComp) = sse(1, nComp) + sum(TC(:, nComp).^2);
                
                for nRegress = 1:length(model_inds)
                    currentTc = TC(:, nComp);
                    allInds = (1:length(model_inds));
                    allInds = find(allInds ~= nRegress);
                    currentTc = currentTc - modelX(:, allInds)*tmpB(allInds);
                    partial_correlations(nSess, nSub, nComp, model_inds(nRegress)) = icatb_corr2(currentTc, modelX(:, nRegress));
                end
                
                
            end
        end
    end
end

%% Store info in stats
R2 = 1 - ssr./sse;

[dd, inds] = sort(R2);
component_numbers = inds(end:-1:1);

betas = betas(:, :, component_numbers, :);
stats.component_numbers = component_numbers;
stats.betas = betas;
stats.R2 = R2(component_numbers);
stats.partial_correlations = partial_correlations(:, :, component_numbers, :);
stats.selectedRegressors = selectedRegressors;
stats.ttest.tstat = NaN(sesInfo.numOfSess, sesInfo.numComp, length(selectedRegressors));
stats.ttest.pval = stats.ttest.tstat;
stats.ttest.df = stats.ttest.pval;

% Compute one sample t-test values
for nSess = 1:sesInfo.numOfSess
    for nComp = 1:sesInfo.numComp
        for nRegress = 1:length(selectedRegressors)
            tmpDat = squeeze(betas(nSess, :, nComp, nRegress));
            tmpDat(isnan(tmpDat) == 1) = [];
            if (length(tmpDat) > 1)
                [stats.ttest.pval(nSess, nComp, nRegress), stats.ttest.tstat(nSess, nComp, nRegress), stats.ttest.df(nSess, nComp, nRegress)] = icatb_ttest(tmpDat);
            end
        end
    end
end

% Print R2 summary metrics
printR2(sesInfo, stats);

% Print partial correlations
printPartialCorr(sesInfo, stats);

disp(['Saving information in ', fullfile(sesInfo.outputDir, [sesInfo.userInput.prefix, '_temporal_regression.txt']), ' ...']);

sesInfo.temporal_stats_betas = stats;
save(param_file, 'sesInfo');

disp('Done');

if (~isempty(displayFormat) && ~strcmpi(displayFormat, 'none'))
    outDir = fullfile(outputDir, [sesInfo.userInput.prefix, '_temporal_sort_results']);
    opts.outputDir = outDir;
    opts.showCode = false;
    opts.useNewFigure = false;
    opts.format = lower(displayFormat);
    opts.createThumbnail = true;
    if (strcmpi(opts.format, 'pdf'))
        opt.useNewFigure = false;
    end
    temporal_stats = stats;
    assignin('base', 'temporal_stats', temporal_stats);
    opts.codeToEvaluate = 'icatb_disp_temp_regress_results(temporal_stats, 0);';
    disp('Generating temporal sorting reults summary. Please wait ....');
    drawnow;
    publish('icatb_disp_temp_regress_results', opts);
    close all;
    disp('Done');
    fprintf('\n');
    
    if (strcmpi(opts.format, 'html'))
        icatb_openHTMLHelpFile(fullfile(outDir, 'icatb_disp_temp_regress_results.html'));
    elseif (strcmpi(opts.format, 'pdf'))
        open(fullfile(outDir, 'icatb_disp_temp_regress_results.pdf'));
    end
    
end

if (nargout)
    varargout{1} = sesInfo;
end


function [modelX, model_inds] = getModelTimeCourse(TC, spmInfo, selectedRegressors, nSess)

modelX = repmat(NaN, size(TC, 1), length(selectedRegressors));


XX = spmInfo.SPM.xX.X;
names = cellstr(char(spmInfo.SPM.xX.name));

try
    spmInfo.SPM.nscan(nSess);
catch
    nSess = 1;
end

if (nSess == 1)
    startTp = 1;
else
    startTp = sum(spmInfo.SPM.nscan(1:nSess-1)) + 1;
end

endTp = sum(spmInfo.SPM.nscan(1:nSess));
model_inds = [];
for nR = 1:length(selectedRegressors)
    regressName = ['Sn(', num2str(nSess), ') ', selectedRegressors{nR}];
    inds = strmatch(lower(regressName), lower(spmInfo.SPM.xX.name), 'exact');
    if (~isempty(inds))
        model_inds = [model_inds, nR];
        modelX(:, nR) = XX(startTp:endTp, inds(1));
    end
end


function printR2(sesInfo, stats)

titlePrint = 'TEMPORAL SORTING OF COMPONENTS OF ALL DATA-SETS USING MULTIPLE REGRESSION CRITERIA:';

numPara = 1;
varStruct(numPara).tag = 'Component Number';
varStruct(numPara).value = stats.component_numbers(:);

numPara = numPara + 1;
varStruct(numPara).tag = 'Multiple Regression Value';
varStruct(numPara).value = stats.R2(:);

for nSub = 1:sesInfo.numOfSub
    for nSess = 1:sesInfo.numOfSess
        %numPara = numPara + 1;
        if (sesInfo.numOfSess == 1)
            currentTag = ['Subject ', num2str(nSub)];
        else
            currentTag = ['Subject ', num2str(nSub), ' Session ', num2str(nSess)];
        end
        for nR = 1:length(stats.selectedRegressors)
            numPara = numPara + 1;
            currentTagR = [currentTag, ' ', stats.selectedRegressors{nR}];
            varStruct(numPara).tag = currentTagR;
            varStruct(numPara).value = squeeze(stats.betas(nSess, nSub, :, nR));
        end
    end
end

icatb_printToFile(fullfile(sesInfo.outputDir, [sesInfo.userInput.prefix, '_temporal_regression.txt']), varStruct, titlePrint, 'row_wise');


function printPartialCorr(sesInfo, stats)

titlePrint = 'Partial correlation values of viewing set: All data-sets';

numPara = 1;
varStruct(numPara).tag = 'Component Number';
varStruct(numPara).value = stats.component_numbers(:);

for nSub = 1:sesInfo.numOfSub
    for nSess = 1:sesInfo.numOfSess
        %numPara = numPara + 1;
        if (sesInfo.numOfSess == 1)
            currentTag = ['Subject ', num2str(nSub)];
        else
            currentTag = ['Subject ', num2str(nSub), ' Session ', num2str(nSess)];
        end
        for nR = 1:length(stats.selectedRegressors)
            numPara = numPara + 1;
            currentTagR = [currentTag, ' ', stats.selectedRegressors{nR}];
            varStruct(numPara).tag = currentTagR;
            varStruct(numPara).value = squeeze(stats.partial_correlations(nSess, nSub, :, nR));
        end
    end
end

icatb_printToFile(fullfile(sesInfo.outputDir, [sesInfo.userInput.prefix, '_temporal_partial_corr.txt']), varStruct, titlePrint, 'row_wise');