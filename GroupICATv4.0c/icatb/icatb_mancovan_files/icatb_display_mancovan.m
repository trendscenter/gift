function icatb_display_mancovan(paramFile)
%% Display Mancovan results
%

icatb_defaults;
global UI_FONTNAME;

if (~exist('paramFile', 'var'))
    paramFile = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Mancovan Parameter File', 'filter', '*mancovan.mat');
end

drawnow;

outputDir = fileparts(paramFile);
if (isempty(outputDir))
    outputDir = pwd;
end

cd(outputDir);
load(paramFile);

if (~exist('mancovanInfo', 'var'))
    error('Please select Mancovan parameter file');
end

if (~isfield(mancovanInfo, 'outputFiles'))
    error('Run setup features followed by mancova in order to visualize the results');
end

desCriteria = 'mancova';
try
    desCriteria = mancovanInfo.designCriteria;
catch
end

mancovanInfo.designCriteria = desCriteria;

%% Delete a previous figure of display GUI
checkDispGUI = findobj('tag', 'display_mancovan');

if ~isempty(checkDispGUI)
    for ii = 1:length(checkDispGUI)
        delete(checkDispGUI(ii));
    end
end

% display figure
graphicsHandle = icatb_getGraphics('Display Mancovan', 'displayGUI', 'display_mancovan', 'off');

set(graphicsHandle, 'CloseRequestFcn', @figCloseCallback);

% set graphics handle menu none
set(graphicsHandle, 'menubar', 'none');


% plot display defaults
resultsSummaryH = uimenu('parent', graphicsHandle, 'label', 'Results Summary (HTML or PDF)', 'callback', {@resultsSummaryCallback, graphicsHandle});

% offsets
xOffset = 0.025; yOffset = 0.075;

% title color
titleColor = [0 0.9 0.9];
% fonts
titleFont = 13;
axes('Parent', graphicsHandle, 'position', [0 0 1 1], 'visible', 'off');
xPos = 0.5; yPos = 0.97;
text(xPos, yPos, 'Display Mancovan Results', 'color', titleColor, 'FontAngle', 'italic', 'fontweight', 'bold', ...
    'fontsize', titleFont, 'HorizontalAlignment', 'center', 'FontName', UI_FONTNAME);

% plot display button
buttonWidth = 0.2; buttonHeight = 0.05;
displayButtonPos = [0.75 yOffset buttonWidth buttonHeight];


displayButtonH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', displayButtonPos, 'string', 'Display', 'tag', 'display_button', 'callback', ...
    {@displayButtonCallback, graphicsHandle});

% plot load anatomical
buttonWidth = 0.26; buttonHeight = 0.05;
loadAnatomicalPos = [0.05 yOffset buttonWidth buttonHeight];
loadAnatomicalH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', loadAnatomicalPos, 'string', 'Load Anatomical', 'tag', 'load_anatomical_button', ...
    'callback', {@loadAnatomicalCallback, graphicsHandle});


promptWidth = 0.52;
promptHeight = 0.05;

popupHeight = 0.05;
popupWidth = 0.2;

promptPos =  [xOffset, yPos - popupHeight - 0.01 - promptHeight, promptWidth, promptHeight];
textH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Select results to display', 'HorizontalAlignment', 'center');
% horizontal alignment - center, vertical alignment - middle
icatb_wrapStaticText(textH);
promptPos = get(textH, 'position');

popupPos = [promptPos(1) + promptPos(3) + xOffset, promptPos(2), popupWidth, popupHeight];
popupH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'popup', 'position', popupPos, 'string', char('Features (T-maps, Spectra, FNC)', 'Multivariate Results', 'Univariate Results'), 'tag', 'select_list', 'callback', ...
    {@selectCovariatesToPlot, graphicsHandle});

if ((icatb_get_matlab_version <= 2013) || strcmpi(version('-release'), '2014a'))
    
    extentPos = get(popupH(end), 'extent');
    extentPos(3) = extentPos(3) + 0.001;
    extentPos(4) = extentPos(4) + 0.001;
    popupPos(2) = popupPos(2) - 0.5*(extentPos(4) - popupPos(4));
    popupPos(3) = extentPos(3);
    popupPos(4) = extentPos(4);
    set(popupH(1), 'position', popupPos);
    
    if (length(popupH) == 2)
        set(popupH(2), 'position', [0, 0, 1, 1])
    end
    
end

promptPos =  [xOffset, promptPos(2) - yOffset - 0.5*promptHeight, promptWidth, promptHeight];
textH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'T-Threshold (Tmap)', 'HorizontalAlignment', 'center');
% horizontal alignment - center, vertical alignment - middle
icatb_wrapStaticText(textH);
promptPos = get(textH, 'position');

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = 0.12;
editH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', '1.0', 'tag', 't_threshold');

promptPos =  [xOffset, promptPos(2) - yOffset - 0.5*promptHeight, promptWidth, promptHeight];
textH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Select Image Values To Display (Tmap)', 'HorizontalAlignment', 'center');
% horizontal alignment - center, vertical alignment - middle
icatb_wrapStaticText(textH);
promptPos = get(textH, 'position');

popupPos = promptPos;
popupPos(1) = popupPos(1) + popupPos(3) + xOffset;
popupPos(3) = popupWidth;
popupH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'popup', 'position', popupPos, 'string', char('Positive and Negative', 'Positive', 'Absolute Value', 'Negative'), 'tag', 'image_values');

if ((icatb_get_matlab_version <= 2013) || strcmpi(version('-release'), '2014a'))
    
    extentPos = get(popupH(end), 'extent');
    extentPos(3) = extentPos(3) + 0.01;
    extentPos(4) = extentPos(4) + 0.01;
    popupPos(2) = popupPos(2) - 0.5*(extentPos(4) - popupPos(4));
    popupPos(3) = extentPos(3);
    popupPos(4) = extentPos(4);
    set(popupH(1), 'position', popupPos);
    
    if (length(popupH) == 2)
        set(popupH(2), 'position', [0, 0, 1, 1])
    end
    
end


promptPos =  [xOffset, promptPos(2) - yOffset - 0.005 - 0.5*promptHeight, promptWidth, promptHeight];
textH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Threshold Criteria (Univariate results)', 'HorizontalAlignment', 'center');
% horizontal alignment - center, vertical alignment - middle
icatb_wrapStaticText(textH);
promptPos = get(textH, 'position');

popupPos = promptPos;
popupPos(1) = popupPos(1) + popupPos(3) + xOffset;
popupPos(3) = 0.15;
popupH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'popup', 'position', popupPos, 'string', char('fdr', 'MAFDR', 'none'), 'tag', 'threshdesc', 'enable', 'on');



promptPos =  [xOffset, promptPos(2) - yOffset - 0.005 - 0.5*promptHeight, promptWidth, promptHeight];
textH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'P-Threshold (Univariate Results)', 'HorizontalAlignment', 'center');
% horizontal alignment - center, vertical alignment - middle
icatb_wrapStaticText(textH);
promptPos = get(textH, 'position');

popupPos = promptPos;
popupPos(1) = popupPos(1) + popupPos(3) + xOffset;
popupPos(3) = 0.15;
popupH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'edit', 'position', popupPos, 'string', num2str(mancovanInfo.userInput.p_threshold), 'tag', ...
    'thresh', 'enable', 'on');


promptPos =  [xOffset, promptPos(2) - yOffset - 0.5*promptHeight, promptWidth, promptHeight];
textH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Low and high frequency limits to compute fALFF', 'HorizontalAlignment', 'center');
% horizontal alignment - center, vertical alignment - middle
icatb_wrapStaticText(textH);
promptPos = get(textH, 'position');

popupPos = promptPos;
popupPos(1) = popupPos(1) + popupPos(3) + xOffset;
popupPos(3) = 0.16;
popupH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'edit', 'position', popupPos, 'string', '0.1, 0.15', 'tag', 'freq_limits', 'enable', 'on');

if ((icatb_get_matlab_version <= 2013) || strcmpi(version('-release'), '2014a'))
    
    extentPos = get(popupH(end), 'extent');
    extentPos(3) = extentPos(3) + 0.01;
    extentPos(4) = extentPos(4) + 0.01;
    popupPos(2) = popupPos(2) - 0.5*(extentPos(4) - popupPos(4));
    popupPos(3) = extentPos(3);
    popupPos(4) = extentPos(4);
    set(popupH(1), 'position', popupPos);
    
    if (length(popupH) == 2)
        set(popupH(2), 'position', [0, 0, 1, 1])
    end
    
end


promptPos =  [0.5 - 0.5*promptWidth, promptPos(2) - yOffset - 0.5*promptHeight, promptWidth, promptHeight];
checkH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'checkbox', 'position', promptPos, 'string', 'Display connectogram (univariate results)', ...
    'HorizontalAlignment', 'left', 'value', 1, 'tag', 'display_connectogram', 'callback', {@displayConnectCallback, graphicsHandle});


handles_data.mancovan = mancovanInfo;
handles_data.mancovan.outputDir = outputDir;
handles_data.structFile = fullfile(fileparts(which('gift.m')), 'icatb_templates', 'ch2bet.nii');
set(graphicsHandle, 'userdata', handles_data);

% promptPos =  [xOffset, popupPos(2) - yOffset - 0.5*promptHeight, promptWidth, promptHeight];
% textH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Select Covariates', 'HorizontalAlignment', 'center');
% % horizontal alignment - center, vertical alignment - middle
% icatb_wrapStaticText(textH);
% promptPos = get(textH, 'position');
%
% popupPos = promptPos;
% popupPos(1) = popupPos(1) + popupPos(3) + xOffset;
% popupPos(3) = 0.3;
% popupH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'popup', 'position', popupPos, 'string', 'Test');

set(graphicsHandle, 'visible', 'on');


function displayButtonCallback(hObject, event_data, handles)
%% Display button callback
%

icatb_defaults;
global UI_FONTNAME;
global UI_FS;

set(handles, 'pointer', 'watch');

handles_data = get(handles, 'userdata');

listH = findobj(handles, 'tag', 'select_list');
listVal = get(listH, 'value');
listStr = cellstr(get(listH, 'string'));
selListVal = listStr{listVal};

imValH = findobj(handles, 'tag', 'image_values');
imStr = cellstr(get(imValH, 'string'));
imStr = lower(imStr{get(imValH, 'value')});


threshDescH = findobj(handles, 'tag', 'threshdesc');
threshdesc = cellstr(get(threshDescH, 'string'));
threshdesc = lower(threshdesc{get(threshDescH, 'value')});
handles_data.mancovan.display.threshdesc = threshdesc;

t_threshold = str2num(get(findobj(handles, 'tag', 't_threshold'), 'string'));

freq_limits = str2num(get(findobj(handles, 'tag', 'freq_limits'), 'string'));

p_val_thresh = str2num(get(findobj(handles, 'tag', 'thresh'), 'string'));

structFile = handles_data.structFile;

display_connectogram = get(findobj(handles, 'tag', 'display_connectogram'),'value');
handles_data.mancovan.display.display_connectogram = display_connectogram;

handles_data.mancovan.display.image_values = imStr;
handles_data.mancovan.display.freq_limits = freq_limits;
handles_data.mancovan.display.structFile = structFile;
handles_data.mancovan.display.t_threshold = t_threshold;
handles_data.mancovan.display.p_threshold = p_val_thresh;

%load(fullfile(handles_data.mancovan.outputDir, [handles_data.mancovan.resultsFile]));

features = handles_data.mancovan.features;

sm_inds = strmatch('spatial maps', lower(features), 'exact');
spectra_inds = strmatch('timecourses spectra', lower(features), 'exact');
fnc_inds = strmatch('fnc correlations', lower(features), 'exact');
comps = handles_data.mancovan.comps;

load icatb_colors coldhot;
graphicsH = [];

try
    
    if (~isempty(findstr(lower(selListVal), 'features')))
        graphicsH = icatb_plot_mancova_features(handles_data.mancovan);
        
        %         plotSM = 0;
        %         plotTC = 0;
        %         if (~isempty(sm_inds))
        %             plotSM = 1;
        %             %sm_results = results{sm_inds};
        %             %tmap_files = cellstr(char(sm_results.tmap_file));
        %         end
        %
        %         if (~isempty(spectra_inds))
        %             plotTC = 1;
        %             %spectra_results = results{spectra_inds};
        %         end
        %
        %         if (plotSM || plotTC)
        %             %graphicsH = repmat(struct('H', []), 1, length(comps));
        %             %% Display ortho slices and spectra
        %             countF = 0;
        %             for nFiles = 1:length(comps)
        %                 countF = countF + 1;
        %
        %                 if (plotTC)
        %                     load(fullfile(handles_data.mancovan.outputDir, handles_data.mancovan.outputFiles(spectra_inds).filesInfo.result_files{nFiles}), 'spectra_tc', 'freq');
        %                     if (~exist('spectra_tc', 'var') || isempty(spectra_tc))
        %                         error('Please run setup features in order to view them.');
        %                     end
        %                     tc.data = spectra_tc;
        %                     tc.xAxis = freq;
        %                     tc.isSpectra = 1;
        %                     tc.xlabelStr = 'Frequency (Hz)';
        %                     tc.ylabelStr = 'Power';
        %                     dynamicrange = zeros(1, size(tc.data, 1));
        %                     fALFF = dynamicrange;
        %                     for nS = 1:size(tc.data, 1)
        %                         [dynamicrange(nS), fALFF(nS)] = icatb_get_spec_stats(tc.data(nS, :), tc.xAxis, freq_limits);
        %                     end
        %                     dynamicrange = mean(dynamicrange);
        %                     fALFF = mean(fALFF);
        %                     tc.titleStr = sprintf('Dynamic range: %0.3f, Power_L_F/Power_H_F: %0.3f', mean(dynamicrange), mean(fALFF));
        %                 end
        %
        %                 if (plotTC && plotSM)
        %                     H = icatb_orth_views(fullfile(handles_data.mancovan.outputDir, handles_data.mancovan.outputFiles(sm_inds).filesInfo.tmap_files{nFiles}), 'structfile', structFile, 'image_values', imStr, 'convert_to_zscores', 'no', 'set_to_max_voxel', 1, 'tc', tc, 'labels', ...
        %                         ['T-map ', icatb_returnFileIndex(comps(nFiles)), ' (T >= ', num2str(t_threshold), ')'], 'fig_title', ['Features (Tmap ', icatb_returnFileIndex(comps(nFiles)), ' and Spectra)'], 'threshold', t_threshold);
        %                 else
        %                     if (plotSM)
        %                         H = icatb_orth_views(fullfile(handles_data.mancovan.outputDir, handles_data.mancovan.outputFiles(sm_inds).filesInfo.tmap_files{nFiles}), 'structfile', structFile, 'image_values', imStr, 'convert_to_zscores', 'no', 'set_to_max_voxel', 1, 'labels', ...
        %                             ['T-map ', icatb_returnFileIndex(comps(nFiles)), ' (T >= ', num2str(t_threshold), ')'], 'fig_title', ['Features (Tmap ', icatb_returnFileIndex(comps(nFiles)), ')'], 'threshold', t_threshold);
        %                     else
        %                         tc.titleStr = sprintf('Comp %s Dynamic range: %0.3f, Power_L_F/Power_H_F: %0.3f', icatb_returnFileIndex(comps(nFiles)), mean(dynamicrange), mean(fALFF));
        %                         tc.fig_title = ['Power Spectra ', icatb_returnFileIndex(comps(nFiles))];
        %                         H = icatb_plot_spectra(tc);
        %                         axisH = axes('Parent', H, 'position', [0 0 1 1], 'visible', 'off');
        %                         xPos = 0.5; yPos = 0.97;
        %                         titleColor = 'c';
        %                         text(xPos, yPos, tc.fig_title, 'color', titleColor, 'fontweight', 'bold', 'fontsize', UI_FS + 2, 'HorizontalAlignment', 'center', 'FontName', UI_FONTNAME, 'parent', axisH);
        %                     end
        %                 end
        %
        %                 graphicsH(length(graphicsH) + 1).H = H;
        %             end
        %
        %         end
        %
        %         clear spectra_tc
        %
        %         if (~isempty(fnc_inds))
        %             load(fullfile(handles_data.mancovan.outputDir, handles_data.mancovan.outputFiles(fnc_inds).filesInfo.result_files{1}), 'fnc_corrs');
        %             if (~exist('fnc_corrs', 'var') || isempty(fnc_corrs))
        %                 error('Please run setup features in order to view them.');
        %             end
        %             M = icatb_vec2mat(icatb_z_to_r(squeeze(mean(fnc_corrs))), 1);
        %             CLIM = max(abs(M(:)));
        %             fig_title = 'Features (FNC Correlations)';
        %
        %
        %             network_values = zeros(1, length(handles_data.mancovan.userInput.comp));
        %             for nV = 1:length(network_values)
        %                 network_values(nV) = length(handles_data.mancovan.userInput.comp(nV).value);
        %             end
        %             network_names =  cellstr(char(handles_data.mancovan.userInput.comp.name));
        %
        %
        %             if (length(network_names) == 1)
        %
        %                 gH = icatb_plot_matrix(M, cellstr(num2str(handles_data.mancovan.comps(:))), cellstr(num2str(handles_data.mancovan.comps(:))), 'fig_title', fig_title, 'tag', 'FNC Correlations', 'title', ...
        %                     'FNC Correlations (Averaged over subjects)', 'cmap', coldhot, 'clim', [-CLIM, CLIM], 'ylabel', 'Components', 'xlabel', 'Components');
        %
        %             else
        %
        %                 gH = icatb_getGraphics(fig_title, 'graphics',  'FNC Correlations', 'on');
        %                 set(gH, 'resize', 'on');
        %                 axesH = axes('parent', gH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
        %                 icatb_plot_FNC(M, [-CLIM, CLIM], cellstr(num2str(handles_data.mancovan.comps(:))), (1:length(handles_data.mancovan.comps)), gH, fig_title, axesH, ...
        %                     network_values, network_names);
        %                 colormap(coldhot);
        %
        %             end
        %
        %             axisH = axes('Parent', gH, 'position', [0 0 1 1], 'visible', 'off');
        %             xPos = 0.5; yPos = 0.97;
        %             titleColor = 'c';
        %             text(xPos, yPos, fig_title, 'color', titleColor, 'fontweight', 'bold', 'fontsize', UI_FS + 2, 'HorizontalAlignment', 'center', 'FontName', UI_FONTNAME, 'parent', axisH);
        %
        %             graphicsH(length(graphicsH) + 1).H = gH;
        %
        %         end
        %
        %         clear fnc_corrs;
        
    elseif (strcmpi(selListVal, 'multivariate results'))
        load(fullfile(handles_data.mancovan.outputDir, handles_data.mancovan.outputFiles(1).filesInfo.result_files{1}), 'MULT');
        if (~exist('MULT', 'var') || isempty(MULT))
            % if (strcmpi(handles_data.mancovan.designCriteria, 'mancova'))
            %     error('Please run mancova in order to view results');
            % else
            disp('No multivariate results to display for t-tests. Use univariate results instead.');
            set(handles, 'pointer', 'arrow');
            return;
            % end
        end
        %% Display multi-variate results
        %graphicsH(length(graphicsH) + 1).H = icatb_plot_mult_mancovan(handles_data.mancovan);
        GH = icatb_plot_mult_mancovan(handles_data.mancovan);
        for nGH = 1:length(GH)
            graphicsH(length(graphicsH) + 1).H = GH(nGH);
        end
        if (isfield(handles_data.mancovan, 'time'))
            % Time 1
            graphicsH(length(graphicsH) + 1).H = icatb_plot_mult_mancovan(handles_data.mancovan, 1);
            % Time 2
            graphicsH(length(graphicsH) + 1).H = icatb_plot_mult_mancovan(handles_data.mancovan, 2);
        end
        
    else
        load(fullfile(handles_data.mancovan.outputDir, handles_data.mancovan.outputFiles(1).filesInfo.result_files{1}), 'UNI');
        if (~exist('UNI', 'var') || isempty(UNI))
            error('Please run mancova in order to view results');
        end
        
        %% Display univariate results
        handles_data.mancovan.structFile = structFile;
        handles_data.mancovan.image_values = imStr;
        graphicsH = icatb_plot_univariate_results(handles_data.mancovan);
        
        if (isfield(handles_data.mancovan, 'time'))
            % Time 1
            gH1 = icatb_plot_univariate_results(handles_data.mancovan, 1);
            % Time 2
            gH2 = icatb_plot_univariate_results(handles_data.mancovan, 2);
            graphicsH = [graphicsH, gH1, gH2];
        end
        
    end
    
    icatb_plotNextPreviousExitButtons(graphicsH);
    set(handles, 'pointer', 'arrow');
catch
    set(handles, 'pointer', 'arrow');
    disp(lasterr);
    icatb_errorDialog(lasterr, 'Display Error', 'modal');
end

function loadAnatomicalCallback(hObject, event_data, handles)
%% Load Anatomical callback

oldDir = pwd;

handles_data = get(handles, 'userdata');
startPath = fileparts(which('gift.m'));
startPath = fullfile(startPath, 'icatb_templates');
if (~exist(startPath, 'dir'))
    startPath = pwd;
end
structFile = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Structural File', 'filter', '*.img;*.nii', 'fileType', 'image', 'fileNumbers', 1, 'startpath', startPath);
drawnow;
if (isempty(structFile))
    error('Structural file is not selected');
end
handles_data.structFile = structFile;
set(handles, 'userdata', handles_data);

cd(oldDir);

function figCloseCallback(hObject, event_data, handles)
% figure close callback

delete(hObject);

function selectCovariatesToPlot(hObject, event_data, handles)
%% Select covariates to plot
%

strs = cellstr(get(hObject, 'string'));
val = get(hObject, 'value');

threshDescH = findobj(handles, 'tag', 'threshdesc');
set(threshDescH, 'enable', 'on');

threshH = findobj(handles, 'tag', 't_threshold');
set(threshH, 'enable', 'on');

if (strcmpi(strs{val}, 'univariate results'))
    set(threshH, 'enable', 'inactive');
    set(threshDescH, 'enable', 'on');
    handles_data = get(handles, 'userdata');
    outDir = handles_data.mancovan.outputDir;
    outputFiles = handles_data.mancovan.outputFiles;
    
    start_terms = {};
    for nO = 1:length(outputFiles)
        for nR = 1:length(outputFiles(nO).filesInfo.result_files)
            load(fullfile(outDir, outputFiles(nO).filesInfo.result_files{nR}), 'UNI');
            for nT = 1:length(UNI.tests)
                for nC = 1:length(UNI.stats{nT}.Contrast)
                    if (~isempty(UNI.stats{nT}.Contrast{nC}))
                        start_terms{end + 1} = UNI.stats{nT}.Contrast{nC};
                    else
                        start_terms{end + 1} = UNI.tests{nT};
                    end
                end
            end
        end
    end
    
    if (isempty(start_terms))
        error('No siginificant covariates found in all of the univariate tests');
    end
    
    [dd, inds] = unique(start_terms);
    
    start_terms = start_terms(sort(inds));
    
    %     load(fullfile(outDir, outputFiles(1).filesInfo.result_files{1}), 'MULT');
    %     start_terms = MULT.start_terms;
    
    titleFig = 'Select Covariates To Plot';
    if (length(start_terms) > 1)
        sel = icatb_listdlg('PromptString', titleFig, 'title_fig', 'Select Covariates To Plot', 'SelectionMode', 'multiple', 'ListString', start_terms, 'movegui', 'center', 'windowStyle', 'modal', 'help', struct('title', 'Covariates To Plot', 'str', ...
            'Covariates of interest will be included in the plots.'));
    else
        sel = 1;
    end
    
    drawnow;
    
    covariatesToPlot = start_terms(sel);
    handles_data.mancovan.covariatesToPlot = covariatesToPlot;
    set(handles, 'userdata', handles_data);
    
end


function resultsSummaryCallback(hObject, event_data, handles)
%% Results summary
%

formatName = questdlg('Select an option to view results summary', 'Mancova results summary', 'PDF', 'HTML', 'HTML');

if (isempty(formatName))
    return;
end

handles_data = get(handles, 'userdata');

mancovanInfo = handles_data.mancovan;

imValH = findobj(handles, 'tag', 'image_values');
imStr = cellstr(get(imValH, 'string'));
imStr = lower(imStr{get(imValH, 'value')});


p_val_thresh = str2num(get(findobj(handles, 'tag', 'thresh'), 'string'));
mancovanInfo.display.p_threshold = p_val_thresh;

threshDescH = findobj(handles, 'tag', 'threshdesc');
threshdesc = cellstr(get(threshDescH, 'string'));
threshdesc = lower(threshdesc{get(threshDescH, 'value')});
mancovanInfo.display.threshdesc = threshdesc;

t_threshold = str2num(get(findobj(handles, 'tag', 't_threshold'), 'string'));

freq_limits = str2num(get(findobj(handles, 'tag', 'freq_limits'), 'string'));

mancovanInfo.display.structFile = handles_data.structFile;

mancovanInfo.display.image_values = imStr;
mancovanInfo.display.freq_limits = freq_limits;
%handles_data.mancovan.structFile = structFile;
mancovanInfo.display.t_threshold = t_threshold;

outDir = fullfile(mancovanInfo.outputDir, [mancovanInfo.prefix, '_results_summary']);
opts.outputDir = outDir;
opts.showCode = false;
opts.useNewFigure = false;
opts.format = lower(formatName);
opts.createThumbnail = true;
if (strcmpi(opts.format, 'pdf'))
    opts.useNewFigure = false;
end
assignin('base', 'mancovanInfo', mancovanInfo);
opts.codeToEvaluate = 'icatb_mancovan_results_summary(mancovanInfo);';
disp('Generating reults summary. Please wait ....');
drawnow;
publish('icatb_mancovan_results_summary', opts);
close all;

if (strcmpi(opts.format, 'html'))
    icatb_openHTMLHelpFile(fullfile(outDir, 'icatb_mancovan_results_summary.html'));
else
    open(fullfile(outDir, 'icatb_mancovan_results_summary.pdf'));
end

disp('Done');

%icatb_mancovan_results_summary(handles_data.mancovan);


function displayConnectCallback(hObject, event_data, handles)


handles_data = get(handles, 'userdata');

display_connectogram = get(findobj(handles, 'tag', 'display_connectogram'),'value');
handles_data.mancovan.display.display_connectogram = display_connectogram;

set(handles, 'userdata', handles_data);