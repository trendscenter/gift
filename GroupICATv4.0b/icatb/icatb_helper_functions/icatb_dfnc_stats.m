function icatb_dfnc_stats(param_file, outputDir)
%% Do stats on dfnc correlations.
%

icatb_defaults;
global UI_FS;

if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('title', 'Select dFNC Parameter File', 'typeEntity', 'file', 'typeSelection', ...
        'single', 'filter', '*dfnc.mat');
    drawnow;
end

if (isempty(param_file))
    error('Please select dFNC parameter file');
end

if (ischar(param_file))
    load(param_file);
    dfncOutputDir = fileparts(param_file);
    if (isempty(dfncOutputDir))
        dfncOutputDir = pwd;
    end
    if (~exist('dfncInfo', 'var'))
        error(['Selected file ', param_file, ' is not a valid dFNC parameter file']);
    end
else
    dfncInfo = param_file;
    dfncOutputDir = dfncInfo.userInput.outputDir;
    clear param_file;
end

drawnow;

%% Check post process file
post_process_file = fullfile(dfncOutputDir,  [dfncInfo.prefix, '_post_process.mat']);
if (~exist(post_process_file, 'file'))
    error('Please run postprocess step in order to use dfnc cluster stats');
end


load (post_process_file, 'clusterInfo');
Nwin = length(clusterInfo.IDXall) / length(dfncInfo.outputFiles);

if (~exist('outputDir', 'var'))
    outputDir = icatb_selectEntry('title', 'Select stats output directory ...', 'typeEntity', 'directory');
end

drawnow;

if (isempty(outputDir))
    outputDir = pwd;
end

defaultWindowSize = min([Nwin, 10]);

% %% Compute cluster stats
% cluster_stats_file = fullfile(outputDir,  [dfncInfo.prefix, '_cluster_stats.mat']);
% % open input dialog box
% prompt = {['Enter threshold in windows to compute median of correlations across windows (Max is ', num2str(Nwin), ')']};
% dlg_title = 'Threshold in windows';
% num_lines = 1;
% def = {num2str(defaultWindowSize)};
% % save the file with the file name specified
% answer = icatb_inputdlg2(prompt, dlg_title, num_lines, def);
% if (isempty(answer))
%     error('Threshold in windows is not selected');
% end
% thresholdWindows = str2num(answer{1});
% disp(['Atleast ', num2str(thresholdWindows), ' windows must be present in each state to compute median of dFNC correlations']);
% disp('Using only median of dFNC correlations for each cluster state...');
% icatb_dfnc_cluster_stats(dfncInfo, outputDir, thresholdWindows);
% drawnow;
%
% load (cluster_stats_file, 'dfnc_corrs');
%
% drawnow;


desCriteriaOptions = {'One sample t-test', 'Two sample t-test', 'Paired T-test'};
figureData.outputDir = outputDir;
figureData.prefix = dfncInfo.prefix;
figureData.dfncOutputDir = dfncOutputDir;
figureData.Nwin = Nwin;
%figureData.cluster_stats_file = cluster_stats_file;
figureData.numOfSub = dfncInfo.userInput.numOfSub;
figureData.numOfSess = dfncInfo.userInput.numOfSess;
figureData.cov = [];
figureData.contrasts = [];
figureData.threshdesc = 'none';
figureData.p_thresh = 0.05;
figureData.dfncInfo = dfncInfo;

clear loading_coeff;

%% Cleanup existing figures
figureTag = 'setup_dfnc_stats_gui';
figHandle = findobj('tag', figureTag);
if (~isempty(figHandle))
    delete(figHandle);
end

%% Setup figure for GUI
InputHandle = icatb_getGraphics('dFNC Stats GUI', 'normal', figureTag, 'on');
set(InputHandle, 'menubar', 'none');
set(InputHandle, 'userdata', figureData);


defaultsmenu = uimenu('parent', InputHandle, 'label', 'Stats-Defaults', 'tag', 'dfnc_stats_defaults');
set(defaultsmenu, 'callback', {@defaultsCallback, InputHandle});

displaymenu = uimenu('parent', InputHandle, 'label', 'HTML Report', 'tag', 'dfnc_stats_defaults_html');
set(displaymenu, 'callback', {@displayCallback, InputHandle});


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

promptHeight = 0.05;
promptWidth = popupTextWidth;

promptPos = [xOffset, popupPos(2) - popupPos(4) - 2*yOffset, promptWidth, promptHeight];
promptH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', ...
    'position', promptPos, 'String', ['Enter threshold in windows (Max = ', num2str(Nwin), ')'], 'fontsize', UI_FS - 1, ...
    'horizontalalignment', 'center');
promptH = icatb_wrapStaticText(promptH);

popupPos = get(promptH, 'position');
popupPos(1) = popupPos(1) + popupPos(3) + xOffset;
popupPos(3) = 0.12;
popupPos(4) = promptHeight;
editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', 'position', popupPos, 'String', num2str(defaultWindowSize), 'fontsize', ...
    UI_FS - 1, 'tag', 'threshold_windows');


% controlWidth = 0.4;
% listboxHeight = controlWidth; listboxWidth = controlWidth;
%
% %% Covariates listbox
% promptPos = [0.5 - 0.5*controlWidth, popupTextPos(2) - popupTextPos(4) - 2*yOffset, promptWidth, promptHeight];
%
% textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Add Covariates', 'tag', ...
%     'prompt_covariates', 'fontsize', UI_FS - 1);
% icatb_wrapStaticText(textH);
% listboxXOrigin = promptPos(1);
% listboxYOrigin = promptPos(2) - 2*yOffset - listboxHeight;
% listboxPos = [listboxXOrigin, listboxYOrigin, listboxWidth, listboxHeight];
% covNames = '';
% icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', listboxPos, 'string', covNames, 'tag', ...
%     'cov', 'enable', 'off', 'fontsize', UI_FS - 1, 'min', 0, 'max', 2, 'value', 1,  'callback', {@addCov, InputHandle});
% addButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) + 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
% removeButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) - 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
% icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', addButtonPos, 'string', '+', 'tag', 'add_cov_button', 'fontsize',...
%     UI_FS - 1, 'enable', 'off', 'callback', {@addCov, InputHandle});
% icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', removeButtonPos, 'string', '-', 'tag', 'remove_cov_button', 'fontsize',...
%     UI_FS - 1, 'enable', 'off', 'callback', {@removeCov, InputHandle});
%

calculateWidth = 0.16;
calculateHeight = 0.05;
yPos = 0.16;
calculatePos = [0.5 - 0.5*calculateWidth, yPos, calculateWidth, calculateHeight];
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', calculatePos, 'string', 'Calculate', ...
    'tag', 'calculate_button', 'fontsize', UI_FS - 1, 'callback', {@calculateCallback, InputHandle});


% calculateWidth = 0.16;
% calculateHeight = 0.05;
% calculatePos = [0.75 - 0.5*calculateWidth, listboxPos(2) - calculateHeight - 1.5*yOffset, calculateWidth, calculateHeight];
% icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', calculatePos, 'string', 'Display', ...
%     'tag', 'display_button', 'fontsize', UI_FS - 1, 'callback', {@displayCallback, InputHandle}, 'tooltipstring', 'Open HTML in web browser ...');

%designCallback(designPopupH, [], InputHandle);


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

figData = get(handles, 'userdata');

%% Compute median of correlations and save it
cluster_stats_file = fullfile(figData.outputDir,  [figData.prefix, '_cluster_stats.mat']);
figData.cluster_stats_file = cluster_stats_file;
threshH = findobj(handles, 'tag', 'threshold_windows');
thresholdWindows = str2num(get(threshH, 'string'));
if (isempty(thresholdWindows))
    error('Enter a valid number for threshold');
end
thresholdWindows = ceil(thresholdWindows);
disp(['Atleast ', num2str(thresholdWindows), ' windows must be present in each state to compute median of dFNC correlations']);
disp('Using only median of dFNC correlations for each cluster state...');
icatb_dfnc_cluster_stats(figData.dfncInfo, figData.outputDir, thresholdWindows);


stats_log = fullfile(figData.outputDir, [figData.prefix, '_stats.log']);
diary (stats_log);
try
    
    %% Design criteria
    designCriteriaH = findobj(handles, 'tag', 'design_criteria');
    designOptions = get(designCriteriaH, 'string');
    designCriteriaVal = get(designCriteriaH, 'value');
    designCriteria = designOptions{designCriteriaVal};
    
    load (figData.cluster_stats_file);
    
    if (isempty(figData.cov) && strcmpi(designCriteria, 'anova'))
        error('Please enter covariates in order to do stats');
    end
    
    
    %data = figData.data;
    
    %% Covariates
    if (strcmpi(designCriteria, 'anova'))
        covH = findobj(handles, 'tag', 'cov');
        covVal = get(covH, 'value');
        if (isempty(covVal))
            error('Please select covariates in order to do stats');
        end
        covInfo = figData.cov(covVal);
        covNames = cellstr(char(covInfo.name));
        covTypes = cellstr(char(covInfo.type));
        chkRegress = find(strcmpi(covTypes, 'continuous') == 1);
        if (~isempty(chkRegress))
            error('Please select categorical covariates only if you want to do anova');
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
    
    cluster_stats_directory = figData.outputDir;
    %     if (~exist(cluster_stats_directory, 'dir'))
    %         mkdir(cluster_stats_directory);
    %     end
    
    try
        groupNames = figData.ttestOpts.t.name;
    catch
        error('Select datasets for computing t-tests');
    end
    
    try
        groupVals = figData.ttestOpts.t.val;
    catch
        error('Select datasets for computing t-tests');
    end
    
    
    disp('....................... Computing stats on dFNC correlations ................................');
    disp('');
    disp('');
    
    if (strcmpi(designCriteria, 'one sample t-test'))
        
        disp('Selected design criteria is one sample t-test');
        disp('');
        
        %% T-test
        if (size(dfnc_corrs, 1) > 1)
            dfnc_corrs = squeeze(mean(dfnc_corrs, 1));
        else
            dfnc_corrs = squeeze(dfnc_corrs);
        end
        
        subInds = [figData.ttestOpts.t.val{1}];
        dfnc_corrs = dfnc_corrs(subInds, :, :);
        
        numClusters = size(dfnc_corrs, 3);
        modelX = ones(size(dfnc_corrs, 1), 1);
        
        %% Initialize results
        t_u = cell(1, numClusters);
        p_u = t_u;
        stats_u = t_u;
        N = zeros(1, numClusters);
        mean_u = cell(1, numClusters);
        subject_indices = cell(1, numClusters);
        
        %% Compute and save
        for nC = 1:numClusters
            disp(['Computing one sample t-test on cluster state# ', num2str(nC), ' ...']);
            tmp = squeeze(dfnc_corrs(:, :, nC));
            chk = find(isfinite(tmp(:, 1)) == 1);
            if (length(chk) > 1)
                tmp = tmp(chk, :);
                N(nC) = size(tmp, 1);
                mean_u{nC} = mean(tmp);
                subject_indices{nC} = subInds(chk);
                [t_u{nC}, p_u{nC}, stats_u{nC}] = mT(tmp, modelX(chk), [], 0, {'verbose'});
            end
        end
        
        outFile = fullfile(cluster_stats_directory, [figData.prefix, '_one_sample_ttest_results.mat']);
        disp(['One sample t-test results are saved in ', outFile]);
        fprintf('\n\n');
        save(outFile, 't_u', 'p_u', 'stats_u', 'mean_u', 'N', 'groupNames', 'groupVals', 'subject_indices');
        
        drawnow;
        
        %% Plot results
        
        
    elseif (strcmpi(designCriteria, 'two sample t-test'))
        
        disp('Selected design criteria is two sample t-test');
        
        %% Two sample t-test
        if (size(dfnc_corrs, 1) > 1)
            dfnc_corrs = squeeze(mean(dfnc_corrs, 1));
        else
            dfnc_corrs = squeeze(dfnc_corrs);
        end
        
        numClusters = size(dfnc_corrs, 3);
        
        %% Initialize results
        t_u = cell(1, numClusters);
        p_u = t_u;
        stats_u = t_u;
        N = zeros(2, numClusters);
        mean_u = cell(2, numClusters);
        
        g1 = [figData.ttestOpts.t.val{1}];
        g2 = [figData.ttestOpts.t.val{2}];
        subject_indices = cell(2, numClusters);
        
        %% Compute and save
        for nC = 1:numClusters
            disp(['Computing two sample t-test on cluster state# ', num2str(nC), ' ...']);
            %tmp = squeeze(dfnc_corrs(:, :, nC));
            tmp1 =  squeeze(dfnc_corrs(g1, :, nC));
            tmp2 =  squeeze(dfnc_corrs(g2, :, nC));
            
            chk1 = find(isfinite(tmp1(:, 1)) == 1);
            chk2 = find(isfinite(tmp2(:, 1)) == 1);
            
            if (~isempty(chk1) && ~isempty(chk2))
                tmp1 = tmp1(chk1, :);
                tmp2 = tmp2(chk2, :);
                N1 = length(chk1);
                N2 = length(chk2);
                modelX = ones(N1 + N2, 1);
                modelX(N1 + 1:end) = 0;
                N(1, nC) = N1;
                N(2, nC) = N2;
                subject_indices{1, nC} = g1(chk1);
                subject_indices{2, nC} = g2(chk2);
                mean_u{1, nC} = mean(tmp1);
                mean_u{2, nC} = mean(tmp2);
                [t_u{nC}, p_u{nC}, stats_u{nC}] = mT([tmp1;tmp2], modelX, [], 1, {'verbose'});
            end
            
        end
        
        outFile = fullfile(cluster_stats_directory, [figData.prefix, '_two_sample_ttest_results.mat']);
        disp(['Two sample t-test results are saved in ', outFile]);
        fprintf('\n\n');
        save(outFile, 't_u', 'p_u', 'stats_u', 'mean_u', 'N', 'groupNames', 'groupVals', 'subject_indices');
        
        drawnow;
        
        %% Plot results
        
    elseif (strcmpi(designCriteria, 'paired t-test'))
        %% Paired sample t-test
        
        disp('Selected design criteria is paired sample t-test');
        
        % reshape to (sessions x subjects) x component pairs x cluster
        % states
        dfnc_corrs = reshape(dfnc_corrs, size(dfnc_corrs, 1)*size(dfnc_corrs, 2), size(dfnc_corrs, 3), size(dfnc_corrs, 4));
        
        numClusters = size(dfnc_corrs, 3);
        
        %% Initialize results
        t_u = cell(1, numClusters);
        p_u = t_u;
        stats_u = t_u;
        N = zeros(2, numClusters);
        mean_u = cell(2, numClusters);
        
        g1 = [figData.ttestOpts.t.val{1}];
        g2 = [figData.ttestOpts.t.val{2}];
        subject_indices = cell(2, numClusters);
        
        %% Compute and save
        for nC = 1:numClusters
            disp(['Computing paired sample t-test on cluster state# ', num2str(nC), ' ...']);
            %tmp = squeeze(dfnc_corrs(:, :, nC));
            tmp1 =  squeeze(dfnc_corrs(g1, :, nC));
            tmp2 =  squeeze(dfnc_corrs(g2, :, nC));
            
            chk = find(isfinite(tmp1(:, 1).*tmp2(:, 1)) == 1);
            
            if (length(chk) > 1)
                tmp1 = tmp1(chk, :);
                tmp2 = tmp2(chk, :);
                N(1, nC) = length(chk);
                N(2, nC) = N(1, nC);
                modelX = ones(N(1, nC), 1);
                mean_u{1, nC} = mean(tmp1);
                mean_u{2, nC} = mean(tmp2);
                subject_indices{1, nC} = g1(chk);
                subject_indices{2, nC} = g2(chk);
                [t_u{nC}, p_u{nC}, stats_u{nC}] = mT(tmp1 - tmp2, modelX, [], 0, {'verbose'});
            end
            
        end
        
        outFile = fullfile(cluster_stats_directory, [figData.prefix, '_paired_ttest_results.mat']);
        disp(['Paired t-test results are saved in ', outFile]);
        fprintf('\n\n');
        save(outFile, 't_u', 'p_u', 'stats_u', 'mean_u', 'N', 'groupNames', 'groupVals', 'subject_indices');
        
        drawnow;
        
        
        %     elseif (strcmpi(designCriteria, 'anova'))
        %         %% Anova
        %
        %
        %         anovaR = icatb_anova(tmpDat, anovaP, 'model', modelType, 'var_names', covNames);
        %         tbl = icatb_anova_table(anovaR);
        %
        %         clear anovaP tmpDat;
        %
        %         icatb_print_table(tbl, fid, 'a+', 0);
        %         fprintf(fid, '\n');
        
    end
    
    disp('');
    disp(['Stats completed. Please see results file ', outFile]);
    disp('');
    disp('');
    
    diary('off');
    
    set(handles, 'pointer', 'arrow');
    
catch
    
    icatb_errorDialog(lasterr, 'dFNC Stats');
    msg = lasterror;
    set(handles, 'pointer', 'arrow');
    rethrow(msg);
    diary('off');
    
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

if (~strcmpi(designCriteria, 'anova'))
    
    set(covH, 'enable', 'off');
    set(addCovH, 'enable', 'off');
    set(removeCovH, 'enable', 'off');
    
else
    
    set(covH, 'enable', 'on');
    set(addCovH, 'enable', 'on');
    set(removeCovH, 'enable', 'on');
    return;
    
end


figData = get(handles, 'userdata');

avgRuns = 1;
if (figData.numOfSess ~= 1)
    if (strcmpi(designCriteria, 'paired t-test'))
        avgRuns = 0;
    end
end

if (~avgRuns)
    subjectStr = cell( figData.numOfSub* figData.numOfSess, 1);
    count = 0;
    for nSub = 1: figData.numOfSub
        for nSess = 1: figData.numOfSess
            count = count + 1;
            subjectStr{count} = ['Subject ', num2str(nSub), ' Session ', num2str(nSess)];
        end
    end
else
    subjectStr = cell(figData.numOfSub, 1);
    count = 0;
    for nSub = 1:figData.numOfSub
        count = count + 1;
        subjectStr{count} = ['Subject ', num2str(nSub)];
    end
end

ttestOpts.des = designCriteria;

if (strcmpi(designCriteria, 'one sample t-test'))
    groupName = '';
    groupVal = [];
    try
        groupName = figData.ttestOpts.t.name{1};
    catch
    end
    try
        groupVal = figData.ttestOpts.t.val{1};
    catch
    end
    [groupName, groupVal] = icatb_select_groups_gui(subjectStr, 'Group', 'selGrp', groupName, groupVal);
    
    if (isempty(groupName))
        error('Data-sets are not selected');
    end
    
    ttestOpts.t.name = {groupName};
    ttestOpts.t.val = {groupVal};
    
elseif (strcmpi(designCriteria, 'two sample t-test'))
    
    group1Name = '';
    groupVal1 = [];
    try
        group1Name = figData.ttestOpts.t.name{1};
    catch
    end
    try
        groupVal1 = figData.ttestOpts.t.val{1};
    catch
    end
    
    group2Name = '';
    groupVal2 = [];
    try
        group2Name = figData.ttestOpts.t.name{2};
    catch
    end
    try
        groupVal2 = figData.ttestOpts.t.val{2};
    catch
    end
    
    [group1Name, groupVal1] = icatb_select_groups_gui(subjectStr, 'Group 1', 'selGrp1', group1Name, groupVal1);
    
    if (isempty(group1Name))
        error('Data-sets for group 1 are not selected');
    end
    
    
    [group2Name, groupVal2] = icatb_select_groups_gui(subjectStr, 'Group 2', 'selGrp2', group2Name, groupVal2);
    
    if (isempty(group2Name))
        error('Data-sets for group 2 are not selected');
    end
    
    ttestOpts.t.name = {group1Name, group2Name};
    ttestOpts.t.val = {groupVal1, groupVal2};
    
elseif (strcmpi(designCriteria, 'paired t-test'))
    
    group1Name = '';
    groupVal1 = [];
    try
        group1Name = figData.ttestOpts.t.name{1};
    catch
    end
    try
        groupVal1 = figData.ttestOpts.t.val{1};
    catch
    end
    
    group2Name = '';
    groupVal2 = [];
    try
        group2Name = figData.ttestOpts.t.name{2};
    catch
    end
    try
        groupVal2 = figData.ttestOpts.t.val{2};
    catch
    end
    
    [group1Name, groupVal1] = icatb_select_groups_gui(subjectStr, 'Condition 1', 'selGrp1', group1Name, groupVal1);
    
    if (isempty(group1Name))
        error('Data-sets for condition 1 are not selected');
    end
    
    [group2Name, groupVal2] = icatb_select_groups_gui(subjectStr, 'Condition 2', 'selGrp2', group2Name, groupVal2);
    
    if (isempty(group2Name))
        error('Data-sets for condition 2 are not selected');
    end
    
    if (length(groupVal1) ~= length(groupVal2))
        error('Please select same number of data-sets for both conditions');
    end
    
    ttestOpts.t.name = {group1Name, group2Name};
    ttestOpts.t.val = {groupVal1, groupVal2};
    
end

figData.ttestOpts = ttestOpts;

set(handles, 'userdata', figData);


function defaultsCallback(hObject, ed, handles)
%% Defaults menu callback

figData = get(handles, 'userdata');

threshdesc = figData.threshdesc;
p_thresh = figData.p_thresh;

opts = char('none', 'fdr');
val = strmatch(threshdesc, opts, 'exact');

numParameters = 1;

inputText(numParameters).promptString = 'Select threshold criteria';
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerString = opts;
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'threshdesc';
inputText(numParameters).value = val;
inputText(numParameters).enable = 'on';

numParameters = numParameters + 1;

inputText(numParameters).promptString = 'Enter p-threshold';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = num2str(p_thresh);
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'p_thresh';
inputText(numParameters).enable = 'on';


if ispc
    windowStyle = 'modal';
else
    windowStyle = 'normal';
end

answer = icatb_inputDialog('inputtext', inputText, 'Title', 'Stats defaults', 'handle_visibility',  'on', 'windowStyle', windowStyle);

if (isempty(answer))
    return;
end

threshdesc = answer{1};
p_thresh = answer{2};

figData.threshdesc = threshdesc;
figData.p_thresh = p_thresh;

set(handles, 'userdata', figData);


function displayCallback(hObject, ed, handles)
%% Display callback
%

set(handles, 'pointer', 'watch');

try
    %% Get info
    figData = get(handles, 'userdata');
    load(fullfile(figData.dfncOutputDir, [figData.prefix, '.mat'])); % dFNC file
    % Stats info
    statsInfo.threshdesc = figData.threshdesc;
    statsInfo.p_threshold = figData.p_thresh;
    cluster_stats_directory = figData.outputDir;
    statsInfo.outputDir = cluster_stats_directory;
    
    if (~exist(statsInfo.outputDir, 'dir') || isempty(icatb_listFiles_inDir(statsInfo.outputDir, '*.*')))
        error('No results to display. Use calculate button to compute stats');
    end
    
    %html_file = fullfile(cluster_stats_directory, 'html', [dfncInfo.prefix, '_stats.html']);
    
    %% Create figures
    helpStr = 'Creating HTML file. This will involve writing jpeg files to the disk. Please wait ...';
    helpH = helpdlg(helpStr);
    disp(helpStr);
    results = icatb_dfnc_results(dfncInfo, 'group comparisons', statsInfo);
    html_dir = fullfile( statsInfo.outputDir, 'html');
    html_file = fullfile(html_dir, [dfncInfo.prefix, '_stats.html']);
    
    %% HTML text
    start_string = '<html><head><title> DFNC Stats </title></head>';
    
    titleStr = 'Group Comparisons';
    
    results_string1 = ['<h1 align = "center">', titleStr, '</h2><p> </p><h2> Contents </h2><p><ul>'];
    
    
    for nR = 1:length(results)
        results_string1 = [results_string1, '<li><h3><a href="#results', num2str(nR), '">', results(nR).tag, '</a></h3></li>'];
    end
    
    results_string1 = [results_string1, '</ul></p><p> </p>'];
    
    for nR = 1:length(results);
        results_string1 = [results_string1, '<hr><h3 align = "left"><a name="results', num2str(nR), '">', results(nR).text, '</a></h3>'];
        results_string1 = [results_string1, '<p align = "center">'];
        for nF = 1:length(results(nR).file)
            if (~isempty(results(nR).file{nF}))
                [ppp, dddd, extn] = fileparts(results(nR).file{nF});
                if (~strcmpi(extn, '.txt'))
                    results_string1 = [results_string1, '<img src = "', results(nR).file{nF}, '" > </img>'];
                else
                    allLines = icatb_textscan(fullfile(html_dir, results(nR).file{nF}));
                    allLines = allLines(icatb_good_cells(allLines));
                    for nLine = 1:length(allLines)
                        
                       allCols = textscan(allLines{nLine}, '%s', 'whitespace', '\t');
                       allCols = allCols{1};
                        
                        if (nLine == 1)
                            tableStr = '<table align = "center" border = "1" width = "35%">';
                            tDStart = '<th align = "center">';
                            tDEnd = '</th>';
                        else
                            tDStart = '<td align = "center">';
                            tDEnd = '</td>';
                        end
                        
                        tableStr = [tableStr, '<tr>'];
                        
                        for nCols = 1:length(allCols)
                            tableStr = [tableStr, tDStart, allCols{nCols}, tDEnd];
                        end
                        
                        tableStr = [tableStr, '</tr>'];
                        
                    end
                    tableStr = [tableStr, '</table>'];
                    results_string1 = [results_string1, tableStr];
                end
            end
        end
        results_string1 = [results_string1, '<p> </p>'];
    end
    
    end_string =  '</html>';
    
    results_string = [start_string,  results_string1, end_string];
    
    dlmwrite(html_file, results_string, '');
    
    icatb_openHTMLHelpFile(html_file);
    
    try
        delete(helpH);
    catch
    end
    
    fprintf('Done\n');
    set(handles, 'pointer', 'arrow');
    
catch
    set(handles, 'pointer', 'arrow');
    icatb_errorDialog(lasterr, 'Display Error');
    try
        delete(helpH);
    catch
    end
    
    % rethrow(lasterror);
    
end