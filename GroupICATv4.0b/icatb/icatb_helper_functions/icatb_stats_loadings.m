function icatb_stats_loadings(fileName)
%% Do stats on ICA loadings
% Input file name must be in image format (*.img or *.nii) or ascii file
%

%% Initialize params
icatb_defaults;
global UI_FS;

if (~exist('fileName', 'var'))
    fileName = icatb_selectEntry('title', 'Select ICA loading coefficents ...', 'typeEntity', 'file', 'typeSelection', 'single', 'filter', ...
        '*group*load*coeff*img;*group*load*coeff*nii');
end

if (isempty(fileName))
    error('ICA Loading coefficients file is not selected');
end

[outDir, fN, extn] = fileparts(fileName);

if (isempty(outDir))
    outDir = pwd;
end

fileName = fullfile(outDir, [fN, extn]);


loading_coeff = loadCoeff(fileName, extn);

drawnow;

desCriteriaOptions = {'One sample t-test', 'Anova', 'Multiple Regression'};
figureData.outputDir = outDir;
figureData.prefix = fN;
figureData.data = loading_coeff;
figureData.numOfSub = size(loading_coeff, 1);
figureData.cov = [];
figureData.contrasts = [];

clear loading_coeff;

%% Cleanup existing figures
figureTag = 'setup_sbm_stats_gui';
figHandle = findobj('tag', figureTag);
if (~isempty(figHandle))
    delete(figHandle);
end

%% Setup figure for GUI
InputHandle = icatb_getGraphics('SBM Stats GUI', 'normal', figureTag, 'on');
set(InputHandle, 'menubar', 'none');
set(InputHandle, 'userdata', figureData);

%% Plot controls

xOffset = 0.04; yOffset = 0.04;
popupTextHeight = 0.05; popupTextWidth = 0.45; yPos = 0.95;
popupTextPos = [xOffset, yPos - yOffset - popupTextHeight, popupTextWidth, popupTextHeight];

% Plot Text
designTextH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', ...
    'position', popupTextPos, 'String', 'Select design criteria', 'fontsize', UI_FS - 1, ...
    'horizontalalignment', 'center');

designTextH = icatb_wrapStaticText(designTextH);

popupTextPos = get(designTextH, 'position');

popupWidth = 0.32;
popupPos = popupTextPos;
popupPos(1) = popupPos(1) + popupPos(3) + xOffset;
popupPos(3) = popupWidth;

% Plot popup
designPopupH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'popup', ...
    'position', popupPos, 'String', desCriteriaOptions, 'fontsize', UI_FS - 1, 'tag', 'design_criteria', 'callback', {@designCallback, InputHandle});

controlWidth = 0.4;
listboxHeight = controlWidth; listboxWidth = controlWidth;
promptHeight = 0.05;
promptWidth = controlWidth;


%% Covariates listbox
promptPos = [0.5 - 0.5*controlWidth, popupTextPos(2) - popupTextPos(4) - 2*yOffset, promptWidth, promptHeight];

textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Add Covariates', 'tag', ...
    'prompt_covariates', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);
listboxXOrigin = promptPos(1);
listboxYOrigin = promptPos(2) - 2*yOffset - listboxHeight;
listboxPos = [listboxXOrigin, listboxYOrigin, listboxWidth, listboxHeight];
covNames = '';
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', listboxPos, 'string', covNames, 'tag', ...
    'cov', 'fontsize', UI_FS - 1, 'min', 0, 'max', 2, 'value', 1,  'callback', {@addCov, InputHandle});
addButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) + 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
removeButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) - 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', addButtonPos, 'string', '+', 'tag', 'add_cov_button', 'fontsize',...
    UI_FS - 1, 'callback', {@addCov, InputHandle});
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', removeButtonPos, 'string', '-', 'tag', 'remove_cov_button', 'fontsize',...
    UI_FS - 1, 'callback', {@removeCov, InputHandle});

calculateWidth = 0.16;
calculateHeight = 0.05;
calculatePos = [0.5 - 0.5*calculateWidth, listboxPos(2) - calculateHeight - 1.5*yOffset, calculateWidth, calculateHeight];
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', calculatePos, 'string', 'Calculate', ...
    'tag', 'calculate_button', 'fontsize', UI_FS - 1, 'callback', {@calculateCallback, InputHandle});


designCallback(designPopupH, [], InputHandle);


function addCov(hObject, event_data, figH)
%% Add Covariates
%

icatb_defaults;
global UI_FS;

figureTag = 'add_cov_mancovan';
covFigHandle = findobj('tag', figureTag);
if (~isempty(covFigHandle))
    delete(covFigHandle);
end

statsInfo = get(figH, 'userdata');

listH = findobj(figH, 'tag', 'cov');

covName = '';
transformationName = '';
covVals = '';
cov_type = 'continuous';
covTypes = {'Continuous', 'Categorical'};
covTypeVal = 1;

if (listH == hObject)
    if (~strcmpi(get(figH, 'selectionType'), 'open'))
        return;
    end
    val = get(listH, 'value');
    try
        covName = statsInfo.cov(val).name;
        covVals = statsInfo.cov(val).value;
        if (isnumeric(covVals))
            covVals = covVals(:);
            covVals = num2str(covVals);
        end
        covVals = cellstr(covVals);
        covVals = covVals(:)';
        transformationName = statsInfo.cov(val).transformation;
        cov_type = statsInfo.cov(val).type;
        covTypeVal = strmatch(lower(cov_type), lower(covTypes), 'exact');
    catch
    end
end

covFigHandle = icatb_getGraphics('Select Covariates', 'normal', figureTag);
set(covFigHandle, 'menubar', 'none');

promptWidth = 0.52; controlWidth = 0.32; promptHeight = 0.05;
xOffset = 0.02; yOffset = promptHeight; yPos = 0.9;
okWidth = 0.12; okHeight = promptHeight;

%% Covariate name and value
promptPos = [xOffset, yPos - 0.5*yOffset, promptWidth, promptHeight];
textH = icatb_uicontrol('parent', covFigHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter Covariate Name', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);

promptPos = get(textH, 'position');

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
editH = icatb_uicontrol('parent', covFigHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', covName, 'tag', 'cov_name', 'fontsize', UI_FS - 1);

%% Type of covariate (Continuous or categorical)
promptPos = [xOffset, editPos(2) - 0.5*promptHeight - yOffset, promptWidth, promptHeight];
textH = icatb_uicontrol('parent', covFigHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Select Type Of Covariate', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);

promptPos = get(textH, 'position');

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
covH = icatb_uicontrol('parent', covFigHandle, 'units', 'normalized', 'style', 'popup', 'position', editPos, 'string', {'Continuous', 'Categorical'}, 'tag', 'cov_type', 'fontsize', UI_FS - 1, ...
    'value', covTypeVal);

promptPos(2) = promptPos(2) - 0.5*promptHeight - yOffset;
promptPos(1) = 0.5 - 0.5*promptWidth;
promptPos(3) = promptWidth;

textH = icatb_uicontrol('parent', covFigHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Add Covariate', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);
promptPos = get(textH, 'position');

editWidth = promptWidth;
editHeight = 0.3;
editPos = [0.5 - 0.5*editWidth, promptPos(2) - yOffset - editHeight, editWidth, editHeight];
editH = icatb_uicontrol('parent', covFigHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', covVals, 'fontsize', UI_FS - 1, 'tag', 'cov_value', 'min', 0, 'max', 2, ...
    'callback', {@covValueCallback, figH});

cmenu = uicontextmenu;

set(editH, 'uicontextmenu', cmenu);

uimenu(cmenu, 'Label', 'Load File', 'callback', {@editContextMenuCallback, covFigHandle, figH});


%% Transformation name and value
promptPos(2) = editPos(2) - promptHeight - yOffset;
promptPos(1) = xOffset;
editPos = promptPos;
textH = icatb_uicontrol('parent', covFigHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter transformation function (like log, atanh). Leave it as empty if you don''t want to apply transformation.', 'fontsize', UI_FS - 1, ...
    'tag', 'prompt_cov_transformation');
icatb_wrapStaticText(textH);
promptPos = get(textH, 'position');
promptPos(2) = promptPos(2) + (editPos(4) - promptPos(4));
set(textH, 'position', promptPos);
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(2) = promptPos(2) + 0.5*promptPos(4) - 0.5*editPos(4);
editPos(3) = controlWidth;
transformH = icatb_uicontrol('parent', covFigHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', transformationName, 'tag', 'cov_transformation', 'fontsize', UI_FS - 1);

okPos = [0.5 - 0.5*okWidth, yOffset + 0.5*okHeight, okWidth, okHeight];
icatb_uicontrol('parent', covFigHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', okPos, 'string', 'Done', 'tag', 'cov_done', 'fontsize', UI_FS - 1, 'callback', ...
    {@setCovCallback, covFigHandle, figH});


set(findobj(covFigHandle, 'tag', 'cov_type'), 'callback', {@covTypeCallback, covFigHandle, transformH, textH});

covTypeCallback(findobj(covFigHandle, 'tag', 'cov_type'), [], covFigHandle, transformH, textH);


function editContextMenuCallback(hObject, event_data, handles, figH)
%% Context menu callback

statsInfo = get(figH, 'userdata');

txtFile = icatb_selectEntry('title', 'Select covariate file' , 'filter', '*.txt;*.asc', 'typeEntity', 'file', 'typeSelection', 'single');
drawnow;
covTypeH = findobj(handles, 'tag', 'cov_type');
opts = cellstr(get(covTypeH, 'string'));
covVal = get(covTypeH, 'value');

try
    val = icatb_mancovan_load_covariates(txtFile, opts{covVal}, statsInfo.numOfSub);
    covValueH = findobj(handles, 'tag', 'cov_value');
    set(covValueH, 'string', val);
catch
    icatb_errorDialog(lasterr, 'Covariate Selection');
end


function setCovCallback(hObject, event_data, covFigH, handles)
%% Set covariate name, value and type

statsInfo = get(handles, 'userdata');

covNameH = findobj(covFigH, 'tag', 'cov_name');
covValueH = findobj(covFigH, 'tag', 'cov_value');
covTypeH = findobj(covFigH, 'tag', 'cov_type');
covTransformationH = findobj(covFigH, 'tag', 'cov_transformation');

% Covariate name, value and type
cov_name = get(covNameH, 'string');
cov_value = get(covValueH, 'string');
opts = cellstr(get(covTypeH, 'string'));
val = get(covTypeH, 'value');
covType = lower(opts{val});
cov_transformation = get(covTransformationH, 'string');

try
    if (isempty(cov_name))
        error('Covariate name is not entered');
    end
    
    if (isempty(cov_value))
        error('Covariate vector is not entered');
    end
    
    if (length(statsInfo.cov) > 0)
        chk = strmatch(lower(cov_name), lower(cellstr(char(statsInfo.cov.name))), 'exact');
        if (~isempty(chk))
            ind = chk;
        end
    end
    
    if (~exist('ind', 'var'))
        ind = length(statsInfo.cov) + 1;
    end
    
    statsInfo.cov(ind).name = cov_name;
    statsInfo.cov(ind).value = deblank(cov_value(:)');
    statsInfo.cov(ind).transformation = lower(cov_transformation);
    statsInfo.cov(ind).type = covType;
    
    set(handles, 'userdata', statsInfo);
    
    covListH = findobj(handles, 'tag', 'cov');
    set(covListH, 'string', cellstr(char(statsInfo.cov.name)));
    delete(covFigH);
    
catch
    icatb_errorDialog(lasterr, 'Covariate Selection');
end


function removeCov(hObject, event_data, figH)
%% Remove Covariate
%

statsInfo = get(figH, 'userdata');
listH = findobj(figH, 'tag', 'cov');
val = get(listH, 'value');

if (~isempty(val))
    check = icatb_questionDialog('title', 'Remove Covariate', 'textbody', 'Do you want to remove the covariate from the list?');
    if (~check)
        return;
    end
end

try
    strs = cellstr(char(statsInfo.cov.name));
    statsInfo.cov(val) = [];
    strs(val) = [];
    set(listH, 'value', 1);
    set(listH, 'string', strs);
    set(figH, 'userdata', statsInfo);
catch
end

function covValueCallback(hObject, event_data, figH)
%% Covariate value callback
%

statsInfo = get(figH, 'userdata');
val = cellstr(get(hObject, 'string'));
inds = icatb_good_cells(val);
val = val(inds);
if (length(val) == 1)
    val = textscan(val{1}, '%s', 'delimiter', '\t,', 'multipleDelimsAsOne', 1);
    val = val{1};
end
eval('val2 = val;');
%eval(['val2 = [', [val{:}], '];']);

if (length(val2) ~= numel(val2))
    error('Covariate must be entered in a  vector');
end

if (length(val2) ~= statsInfo.numOfSub)
    error(['Covariate vector length must equal the no. of subjects (', num2str(statsInfo.numOfSub), ')']);
end

% if (isnumeric(val2))
%     val2 = val2(:);
%     val2 = num2str(val2);
% end

val2 = strtrim(cellstr(val2));

set(hObject, 'string', val2);


function loading_coeff = loadCoeff(fileName, extn)
% Load coefficients

if (strcmpi(extn, '.nii') || strcmpi(extn, '.img'))
    loading_coeff = icatb_loadData(fileName);
else
    loading_coeff = icatb_load_ascii_or_mat(fileName);
end

loading_coeff = squeeze(loading_coeff);

if (length(loading_coeff) == numel(loading_coeff))
    loading_coeff = loading_coeff(:);
end


function covTypeCallback(hObject, ed, handles, transformH, promptH)
%% Covariates Type callback
%

str = cellstr(get(hObject, 'string'));
str = str{get(hObject, 'value')};
set(promptH, 'visible', 'off');
set(promptH, 'enable', 'off');
set(transformH, 'visible', 'off');
set(transformH, 'enable', 'off');


function calculateCallback(hObject, ed, handles)
%% Calculate callback
%


try
    
    figData = get(handles, 'userdata');
    
    %% Design criteria
    designCriteriaH = findobj(handles, 'tag', 'design_criteria');
    designOptions = get(designCriteriaH, 'string');
    designCriteriaVal = get(designCriteriaH, 'value');
    designCriteria = designOptions{designCriteriaVal};
    
    if (isempty(figData.cov) && ~strcmpi(designCriteria, 'one sample t-test'))
        error('Please enter covariates in order to do stats');
    end
    
    data = figData.data;
    
    %% Covariates
    if (~strcmpi(designCriteria, 'one sample t-test'))
        covH = findobj(handles, 'tag', 'cov');
        covVal = get(covH, 'value');
        if (isempty(covVal))
            error('Please select covariates in order to do stats');
        end
        covInfo = figData.cov(covVal);
        covNames = cellstr(char(covInfo.name));
        covTypes = cellstr(char(covInfo.type));
        
        if (strcmpi(designCriteria, 'anova'))
            chkRegress = find(strcmpi(covTypes, 'continuous') == 1);
            if (~isempty(chkRegress))
                error('Please select categorical covariates only if you want to do anova');
            end
        else
            chkRegress = find(strcmpi(covTypes, 'categorical') == 1);
            if (~isempty(chkRegress))
                error('Please select continuous covariates only if you want to do multiple regression');
            end
        end
    end
    
    %% Model Type
    modelType = 'linear';
    if (strcmpi(designCriteria, 'anova'))
        if (length(covVal) > 1)
            chkModel = icatb_questionDialog('title', 'Model Interactions', 'textbody', 'Do You Want To Model Interactions?');
            if (chkModel)
                modelType = 'interaction';
            end
        end
        
        anovaParameters = cell(1, length(covVal));
        for nR = 1:length(anovaParameters)
            anovaParameters{nR} = covInfo(nR).value;
        end
    end
    
    
    set(handles, 'pointer', 'watch');
    
    helpStr = 'Computing stats on ICA loading coefficients ...';
    
    helpH = helpdlg(helpStr, 'Stats');
    
    disp(helpStr);
    
    outFileName = fullfile(figData.outputDir, [figData.prefix, '_stats_summary.txt']);
    
    
    if (~exist(outFileName, 'file'))
        
        %% Print title
        fid = fopen(outFileName, 'a+');
        
        if fid == -1
            error('Error:FileOpen', 'File %s can''t be opened', outFileName);
        end
        
        fprintf(fid, '%s\n', '.............................................');
        fprintf(fid, '%s\n', '       SUMMARY OF STATS ON ICA LOADINGS      ');
        fprintf(fid, '%s\n', '.............................................');
        fprintf(fid, '\n');
        fprintf(fid, '\n');
        
        fclose(fid);
        
    end
    
    componentNum = (1:size(data, 2));
    
    titlePrint = ['Design critera: ', designCriteria];
    
    if (strcmpi(designCriteria, 'one sample t-test'))
        %% T-test
        
        pValue = zeros(1, size(data, 2));
        tValue = pValue;
        
        fid = fopen(outFileName, 'a+');
        fprintf(fid, '%s\n\n', titlePrint);
        tbl = cell(size(data, 2) + 1, 3);
        tbl(1, :) = {'Component Number', 'p-value', 'T-value'};
        for nComp = 1:size(data, 2)
            tmpDat = data(:, nComp);
            tmpDat(isnan(tmpDat)) = [];
            [pValue(nComp), tValue(nComp)] = icatb_ttest(tmpDat);
            tbl(nComp + 1, :) = {num2str(nComp), num2str(pValue(nComp)), num2str(tValue(nComp))};
        end
        
        icatb_print_table(tbl, fid, 'a+', 0);
        fprintf(fid, '\n\n');
        fclose(fid);
        
    elseif (strcmpi(designCriteria, 'anova'))
        %% Anova
        
        
        sepStrs = repmat('.', 1, 30);
        
        fid = fopen(outFileName, 'a+');
        
        fprintf(fid, '%s\n\n', titlePrint);
        
        %% Print summary
        for nComp = 1:length(componentNum)
            
            fprintf(fid, '%s\n', sepStrs);
            fprintf(fid, '%s\n', ['Component ', num2str(nComp), ': ']);
            fprintf(fid, '%s\n', sepStrs);
            
            tmpDat = data(:, nComp);
            chkNans = find(isnan(tmpDat) == 1);
            tmpDat(chkNans) = [];
            
            anovaP = cell(1, length(anovaParameters));
            for nR = 1:length(anovaParameters)
                tmpAnovaR = anovaParameters{nR};
                tmpAnovaR(chkNans) = [];
                anovaP{nR} = tmpAnovaR;
            end
            
            %% Anova
            anovaR = icatb_anova(tmpDat, anovaP, 'model', modelType, 'var_names', covNames);
            tbl = icatb_anova_table(anovaR);
            
            clear anovaP tmpDat;
            
            icatb_print_table(tbl, fid, 'a+', 0);
            fprintf(fid, '\n');
            
        end
        %% End for printing summary
        
        fprintf(fid, '\n\n');
        
        fclose(fid);
        
    else
        %% Multiple regression
        
        partial_corr_names = strcat('Partial Corr (', covNames, ')');
        partial_corr_names = partial_corr_names';
        
        beta_names = strcat('Betas (', covNames, ')');
        beta_names = beta_names';
        headings = {'R-square', beta_names{:}, partial_corr_names{:}};
        
        modelX = zeros(figData.numOfSub, length(covInfo));
        
        for nM = 1:size(modelX, 2)
            modelX(:, nM) = str2num(char(covInfo(nM).value));
        end
        
        tbl = cell(2, length(headings));
        tbl(1, :) = headings;
        
        sepStrs = repmat('.', 1, 30);
        
        fid = fopen(outFileName, 'a+');
        
        fprintf(fid, '%s\n\n', titlePrint);
        
        %% Print summary
        for nComp = 1:length(componentNum)
            
            tmpDat = data(:, nComp);
            chkNans = find(isnan(tmpDat) == 1);
            tmpDat(chkNans) = [];
            tmpModelX = modelX;
            tmpModelX(chkNans, :) = [];
            tmpModelX = detrend(tmpModelX, 0);
            
            fprintf(fid, '%s\n', sepStrs);
            fprintf(fid, '%s\n', ['Component ', num2str(nComp), ': ']);
            fprintf(fid, '%s\n', sepStrs);
            [rSquare_stat, beta_weights, ModelIndices, otherIndices, linearRegress, removeTrend, tmpData, ...
                subject_partial_corr] = icatb_multipleRegression(tmpModelX, tmpDat, size(tmpModelX, 2), 1, ...
                length(tmpDat), 0);
            
            beta_weights = beta_weights(ModelIndices);
            betaStr = cellstr(num2str(beta_weights(:), '%0.6f'));
            partialCorrStr = cellstr(num2str(subject_partial_corr(:), '%0.6f'));
            tbl(2, :) = {num2str(rSquare_stat, '%0.6f'), betaStr{:}, partialCorrStr{:}};
            icatb_print_table(tbl, fid, 'a+', 0);
            fprintf(fid, '\n');
            
            clear tmpModelX tmpDat chkNans;
            
        end
        
        fprintf(fid, '\n\n');
        fclose(fid);
        
    end
    
    try
        delete(helpH);
    catch
    end
    
    disp(['Stats completed. Please see summary file ', outFileName]);
    disp('');
    disp('');
    
    set(handles, 'pointer', 'arrow');
    
catch
    
    icatb_errorDialog(lasterr, 'SBM Stats');
    
    msg = lasterror;
    
    if (exist('fid', 'var'))
        try
            fclose(fid);
        catch
        end
    end
    
    set(handles, 'pointer', 'arrow');
    rethrow(msg);
    
end

function designCallback(hObject, ed, handles)
%% Design callback
%


%% Design criteria
designCriteriaH = findobj(handles, 'tag', 'design_criteria');
designOptions = get(designCriteriaH, 'string');
designCriteriaVal = get(designCriteriaH, 'value');
designCriteria = designOptions{designCriteriaVal};

covH = findobj(handles, 'tag', 'cov');
addCovH = findobj(handles, 'tag', 'add_cov_button');
removeCovH = findobj(handles, 'tag', 'remove_cov_button');

% contrastListH = findobj(handles, 'tag', 'selContrasts');
% addContrastsH = findobj(handles, 'tag',  'addContrasts');
% removeContrastsH = findobj(handles, 'tag',  'removeContrasts');

if (strcmpi(designCriteria, 'one sample t-test'))
    
    set(covH, 'enable', 'off');
    set(addCovH, 'enable', 'off');
    set(removeCovH, 'enable', 'off');
    %     set(contrastListH, 'enable', 'off');
    %     set(addContrastsH, 'enable', 'off');
    %     set(removeContrastsH, 'enable', 'off');
    
else
    
    set(covH, 'enable', 'on');
    set(addCovH, 'enable', 'on');
    set(removeCovH, 'enable', 'on');
    %     set(contrastListH, 'enable', 'on');
    %     set(addContrastsH, 'enable', 'on');
    %     set(removeContrastsH, 'enable', 'on');
    
end