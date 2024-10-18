function icatb_auto_tb_gui()
    UI_FS=10;
    icatb_defaults;
    global PARAMETER_INFO_MAT_FILE;
    global UI_FS;
    s_anatomic_atlas = 'aal';
    s_function_atlas = 'yeo_buckner';
    clear sesInfo;

    %% Draw graphics % First Main Input GUI 
    figureTag = 'setup_nbic_gui';
    figHandle = findobj('tag', figureTag);
    if (~isempty(figHandle))
        delete(figHandle);
    end

    InputHandle = icatb_getGraphics('Auto Labeler', 'normal', figureTag, 'off');
    set(InputHandle, 'menubar', 'none');
    set(InputHandle, 'userdata', 5);

    promptWidth = 0.52; controlWidth = 0.32; promptHeight = 0.05;
    xOffset = 0.02; yOffset = promptHeight; listboxHeight = 0.04; listboxWidth = controlWidth; yPos = 0.9;
    okWidth = 0.12; okHeight = promptHeight;

    promptPos = [xOffset, yPos - 0.5*yOffset, promptWidth, promptHeight];

    listboxXOrigin = promptPos(1) + promptPos(3) + xOffset;
    listboxYOrigin = promptPos(2) - 0.5*listboxHeight;

    %% Get output dir
    s_folder = icatb_selectEntry('typeEntity', 'directory', 'title', 'Select output directory for analysis');

    %% Info that additional MATLAB packs needs to be installed
    textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', [0.0000    0.910    1.000    0.0500], ...
        'string', 'This will not run without adding BCT in your MATLAB path first', 'tag', 'prompt_components', 'fontsize', UI_FS + 2);
    icatb_wrapStaticText(textH);

     %%  Param vs Spatial Map
     d_ypos_plus=-12;
    textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos + [0 -(d_ypos_plus+12.75)*okHeight 0 0], ...
        'string', 'Select Paraneter File vs Spatial Map:', 'tag', 'prompt_components', 'fontsize', UI_FS - 1);
    icatb_wrapStaticText(textH);

    listboxYOrigin = promptPos(2) - 0.5*listboxHeight - (d_ypos_plus+12.75)*okHeight;
    listboxPos = [listboxXOrigin, listboxYOrigin, listboxWidth, listboxHeight];
    h_param_vs_smap = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'popup', 'position', listboxPos, ...
        'string', {'Parameter File', 'Spatial Map'}, 'tag', 'tag_param_vs_smap', 'fontsize', UI_FS - 1, 'callback', {@fun_param_vs_smap, InputHandle});

    rowdButPos =  [listboxPos(1) + 1*listboxPos(3) + 1*xOffset, promptPos(2)-(d_ypos_plus+13)*okHeight, promptHeight + 0.01, promptHeight - 0.01];           
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', rowdButPos , 'string', '?', 'tag', 'tag_help_mask', 'fontsize', ...
        UI_FS - 1, 'callback', {@fun_help_param_or_sm});

    %% Output folder
    d_ypos_plus = 2.25;       
    textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos + [0 -d_ypos_plus*okHeight 0 0], ...
        'string', 'Output Folder', 'tag', 'prompt_components', 'fontsize', UI_FS - 1);
    icatb_wrapStaticText(textH);
    editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', ...
        'position', promptPos + [+.54 -d_ypos_plus*okHeight -.2 0], 'String', s_folder, 'fontsize', UI_FS - 1, 'tag', 'tag_edit_output_folder');

    rowdButPos = [listboxPos(1) + listboxPos(3) + xOffset, promptPos(2)-d_ypos_plus*okHeight, promptHeight + 0.01, promptHeight - 0.01];
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', rowdButPos, 'string', '?', 'tag', 'tag_help_function_atlas', 'fontsize',...
        UI_FS - 1, 'callback', {@fun_help_output_folder});   


    %% Param File or Spatial Map
    d_ypos_plus = 3.5;       
    textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos + [0 -d_ypos_plus*okHeight 0 0], ...
        'string', '', 'tag', 'tag_text_param_or_smap', 'fontsize', UI_FS - 1);
    icatb_wrapStaticText(textH);

    editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', 'position', promptPos + [+.54 -d_ypos_plus*okHeight -.2 0], ...
        'String', '', 'fontsize', UI_FS - 1, 'tag', 'tag_edit_param_or_smap');

    rowdButPos = [listboxPos(1) + listboxPos(3) + xOffset, promptPos(2)-d_ypos_plus*okHeight, promptHeight + 0.01, promptHeight - 0.01];
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', rowdButPos, 'string', '?', 'tag', 'tag_help_function_atlas', 'fontsize',...
        UI_FS - 1, 'callback', {@fun_help_param_or_smap});       


    %% Mask File
    d_ypos_plus = 4.75;       
    textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos + [0 -d_ypos_plus*okHeight 0 0], ...
        'string', 'Do not use', 'tag', 'tag_text_mask', 'fontsize', UI_FS - 1);
    icatb_wrapStaticText(textH);
    editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', 'position', promptPos + [+.54 -d_ypos_plus*okHeight -.2 0], ...
        'String', '', 'fontsize', UI_FS - 1, 'tag', 'tag_edit_mask');

    rowdButPos = [listboxPos(1) + listboxPos(3) + xOffset, promptPos(2)-d_ypos_plus*okHeight, promptHeight + 0.01, promptHeight - 0.01];
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', rowdButPos, 'string', '?', 'tag', 'tag_help_function_atlas', 'fontsize',...
        UI_FS - 1, 'callback', {@fun_help_mask});   

    %%  Skip Artifact detection
    %promptPos(2) = listboxYOrigin - yOffset - 0.5*listboxHeight;
    d_ypos_plus = 4.5;
    textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos + [0 -(d_ypos_plus+1.25)*okHeight 0 0], 'string', 'Skip Artifact Detection:', 'tag', ...
        'prompt_components', 'fontsize', UI_FS - 1);
    icatb_wrapStaticText(textH);

    listboxYOrigin = promptPos(2) - 0.5*listboxHeight - 1*okHeight;
    listboxPos = [listboxXOrigin, listboxYOrigin-d_ypos_plus*okHeight, listboxWidth, listboxHeight];
    h_yn = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'popup', 'position', listboxPos, 'string', {'No', 'Yes'}, 'tag', ...
    'tag_detect_artifact', 'fontsize', UI_FS - 1, 'min', 0, 'max', 1, 'value', 1);

    rowdButPos =  [listboxPos(1) + 1*listboxPos(3) + 1*xOffset, promptPos(2)-(d_ypos_plus+1.5)*okHeight, promptHeight + 0.01, promptHeight - 0.01];    
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', rowdButPos , 'string', '?', 'tag', 'tag_help_arti', 'fontsize',...
        UI_FS - 1, 'callback', {@fun_callback_detect_artifact});
 
    
    %%  skip_anatomical
    d_ypos_plus = 2.5;    
    textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos + [0 -(d_ypos_plus+4.25)*okHeight 0 0], 'string', 'Skip Anatomical Labeling', 'tag', ...
        'prompt_components', 'fontsize', UI_FS - 1);
    icatb_wrapStaticText(textH);
    listboxYOrigin = promptPos(2) - 0.5*listboxHeight - (d_ypos_plus+4.25)*okHeight;
    listboxPos = [listboxXOrigin, listboxYOrigin, listboxWidth, listboxHeight];
    listCompH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'popup', 'position', listboxPos, 'string', {'No', 'Yes'}, 'tag', ...
        'tag_skip_anat', 'fontsize', UI_FS - 1, 'min', 0, 'max', 1, 'value', 1);
    rowdButPos =  [listboxPos(1) + 1*listboxPos(3) + 1*xOffset, promptPos(2)-(d_ypos_plus+4.5)*okHeight, promptHeight + 0.01, promptHeight - 0.01];
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', rowdButPos , 'string', '?', 'tag', 'tagHelpComponents', 'fontsize',...
        UI_FS - 1, 'callback', {@funCallbackHelpComponents});

    %% Anatomic prefix
    d_ypos_plus = 8.25;       
    textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos + [0 -d_ypos_plus*okHeight 0 0], ...
        'string', 'Prefix for Anatomical Atlas', 'tag', 'prompt_components', 'fontsize', UI_FS - 1);
    icatb_wrapStaticText(textH);
    editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', ...
        'position', promptPos + [+.54 -d_ypos_plus*okHeight -.2 0], 'String', s_anatomic_atlas, 'fontsize', UI_FS - 1, 'tag', 'tag_anatomic_atlas');

    rowdButPos = [listboxPos(1) + listboxPos(3) + xOffset, promptPos(2)-d_ypos_plus*okHeight, promptHeight + 0.01, promptHeight - 0.01];
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', rowdButPos, 'string', '?', 'tag', 'tag_help_anatiomic_atlas', 'fontsize',...
        UI_FS - 1, 'callback', {@fun_anatomic_atlas});    
    

    %%  skip_functional 
    d_ypos_plus = 0.75;
    textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos + [0 -(d_ypos_plus+8.5)*okHeight 0 0], 'string', 'Skip Functional Labeling', 'tag', ...
        'prompt_components', 'fontsize', UI_FS - 1);
    icatb_wrapStaticText(textH);
    listboxYOrigin = promptPos(2) - 0.5*listboxHeight - (d_ypos_plus+8.5)*okHeight;
    listboxPos = [listboxXOrigin, listboxYOrigin, listboxWidth, listboxHeight];
    listCompH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'popup', 'position', listboxPos, 'string', {'No', 'Yes'}, 'tag', ...
        'tag_skip_functional', 'fontsize', UI_FS - 1, 'min', 0, 'max', 1, 'value', 1);
    
    rowdButPos =  [listboxPos(1) + 1*listboxPos(3) + 1*xOffset, promptPos(2)-(d_ypos_plus+9.0)*okHeight, promptHeight + 0.01, promptHeight - 0.01];    
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', rowdButPos , 'string', '?', 'tag', 'tagHelpComponents', 'fontsize',...
        UI_FS - 1, 'callback', {@funCallbackHelpScoreCsv});


    %% Functional prefix
    d_ypos_plus = 11;       
    textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos + [0 -d_ypos_plus*okHeight 0 0], ...
        'string', 'Prefix for Functional Atlas', 'tag', 'prompt_components', 'fontsize', UI_FS - 1);
    icatb_wrapStaticText(textH);
    editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', ...
        'position', promptPos + [+.54 -d_ypos_plus*okHeight -.2 0], 'String', s_function_atlas, 'fontsize', UI_FS - 1, 'tag', 'tag_function_atlas');

    rowdButPos = [listboxPos(1) + listboxPos(3) + xOffset, promptPos(2)-d_ypos_plus*okHeight, promptHeight + 0.01, promptHeight - 0.01];
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', rowdButPos, 'string', '?', 'tag', 'tag_help_function_atlas', 'fontsize',...
        UI_FS - 1, 'callback', {@fun_function_atlas});    

    
    %% IC Threshold
    d_ypos_plus = 0.75; 
    textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos + [0 -(d_ypos_plus+12)*okHeight 0 0], ...
        'string', 'Independent component threshold', 'tag', 'prompt_components', 'fontsize', UI_FS - 1);
    icatb_wrapStaticText(textH);
    editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', ...
        'position', promptPos + [+.54 -(d_ypos_plus+12)*okHeight -.2 0], 'String', '3', 'fontsize', UI_FS - 1, 'tag', 'tag_threshold_ic');

    rowdButPos = [listboxPos(1) + listboxPos(3) + xOffset, promptPos(2)-(d_ypos_plus+12)*okHeight, promptHeight + 0.01, promptHeight - 0.01];
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', rowdButPos, 'string', '?', 'tag', 'tagHelpXSize', 'fontsize',...
        UI_FS - 1, 'callback', {@fun_ic_threshold});

    %% n cspatial correlations    
    d_ypos_plus = 1; 
    textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos + [0 -(d_ypos_plus+13)*okHeight 0 0], ...
        'string', 'Number of spatial correlation mappings', 'tag', 'prompt_components', 'fontsize', UI_FS - 1);
    icatb_wrapStaticText(textH);
    editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', ...
        'position', promptPos + [+.54 -(d_ypos_plus+13)*okHeight -.2 0], 'String', '3', 'fontsize', UI_FS - 1, 'tag', 'tag_n_corrs');

    rowdButPos = [listboxPos(1) + listboxPos(3) + xOffset, promptPos(2)-(d_ypos_plus+13)*okHeight, promptHeight + 0.01, promptHeight - 0.01];
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', rowdButPos , 'string', '?', 'tag', 'tagHelpYSize', 'fontsize',...
        UI_FS - 1, 'callback', {@fun_n_corrs});
 
    %% label file used for noise_training dataset
    d_ypos_plus = 1.25; 
    textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos + [0 -(d_ypos_plus+14)*okHeight 0 0], ...
        'string', 'label file prefix (for noise_training dataset)', 'tag', 'prompt_components', 'fontsize', UI_FS - 1);
    icatb_wrapStaticText(textH);  
    editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', ...
        'position', promptPos + [+.54 -(d_ypos_plus+14)*okHeight -.2 0], 'String', 'pre_fbirn_sub', 'fontsize', UI_FS - 1, 'tag', 'tag_prefix_label_train');
    rowdButPos = [listboxPos(1) + listboxPos(3) + xOffset, promptPos(2)-(d_ypos_plus+14)*okHeight, promptHeight + 0.01, promptHeight - 0.01];
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', rowdButPos , 'string', '?', 'tag', 'tagHelpTolerance', 'fontsize',...
        UI_FS - 1, 'callback', {@fun_prefix_label_train});
 





    



    
    
    %% Add cancel, save and run buttons
    cancelPos = [0.25 - 0.5*okWidth, promptPos(2) - yOffset - 15.5*okHeight, okWidth*2, okHeight];
    cancelPos(2) = cancelPos(2) - 0.5*cancelPos(4);
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', cancelPos, 'string', 'Cancel', 'tag', 'cancel_button', 'fontsize',...
        UI_FS - 1, 'callback', 'delete(gcbf);');

    okPos = [0.75 - 0.5*okWidth, promptPos(2) - yOffset - 15.5*okHeight, okWidth*2, okHeight];
    okPos(2) = okPos(2) - 0.5*okPos(4);
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', okPos, 'string', 'Run', 'tag', 'run_button', 'fontsize',...
        UI_FS - 1, 'callback', {@runCallback, InputHandle});

    try
        delete(msgH);
    catch
    end

    set(InputHandle, 'visible', 'on');
    drawnow;

    



function setCompCallback(hObject, event_data, compFigH, handles)
    %% Get fields from component

    icatbInfo = get(handles, 'userdata');
    networkNameH = findobj(compFigH, 'tag', 'comp_network_name');
    networkName = deblank(get(networkNameH, 'string'));

    try

        if (isempty(networkName))
            error('You must enter a component network name');
        end

        listH = findobj(compFigH, 'tag', 'components');
        comps = get(listH, 'value');

        if (isempty(comps))
            error('Components are not selected');
        end

        if (length(icatbInfo.userInput.comp) > 0)
            chk = strmatch(lower(networkName), lower(cellstr(char(icatbInfo.userInput.comp.name))), 'exact');
            if (~isempty(chk))
                ind = chk;
            end
        end

        if (~exist('ind', 'var'))
            ind = length(icatbInfo.userInput.comp) + 1;
        end

        %% Set user selected information in figure
        icatbInfo.userInput.comp(ind).name = networkName;
        icatbInfo.userInput.comp(ind).value =  comps(:)';
        set(handles, 'userdata', icatbInfo);
        compListH = findobj(handles, 'tag', 'comp');
        set(compListH, 'string', cellstr(char(icatbInfo.userInput.comp.name)));
        delete(compFigH);

    catch
        icatb_errorDialog(lasterr, 'Component Selection');
    end

function drawComp(hObject, event_data, figH, compFigHandle)
    %% Draw component

    icatb_defaults;
    global UI_FONTNAME;
    global FONT_COLOR;

    fontSizeText = 8;
    set(compFigHandle, 'pointer', 'watch');

    listH = findobj(figH, 'tag', 'comp');

    axesH = get(compFigHandle, 'currentaxes');


    sel_comp = get(findobj(compFigHandle, 'tag', 'components'), 'value');

    if (~isempty(sel_comp))
        DIM = [size(compData, 1), size(compData, 2), length(sel_comp)];
        [im, numImagesX, numImagesY, textToPlot] = icatb_returnMontage(compData(:, :, sel_comp), [], DIM, [1, 1, 1], sel_comp);
        image(im, 'parent', axesH, 'CDataMapping', 'scaled');
        set(axesH, 'clim', clim); % set the axis positions to the specified
        axis(axesH, 'off');
        axis(axesH, 'image');
        colormap(cmap);
        textCount = 0;
        dim = size(im);
        yPos = 1 + dim(1) / numImagesY;
        for nTextRows = 1:numImagesY
            xPos = 1;
            for nTextCols = 1:numImagesX
                textCount = textCount + 1;
                if textCount <= DIM(3)
                    text(xPos, yPos, num2str(round(textToPlot(textCount))), 'color', FONT_COLOR,  ...
                        'fontsize', fontSizeText, 'HorizontalAlignment', 'left', 'verticalalignment', 'bottom', ...
                        'FontName', UI_FONTNAME, 'parent', axesH);
                end
                xPos = xPos + (dim(2) / numImagesX);
            end
            % end for cols
            yPos = yPos + (dim(1) / numImagesY); % update the y position
        end
    else
        cla(axesH);
    end

    set(compFigHandle, 'pointer', 'arrow');


function runCallback(hObject, event_data, handles)
    % Run autolabeller
    
    % Get the variables
    h_tag_threshold_ic = findobj(handles, 'tag', 'tag_threshold_ic');
    d_tag_threshold_ic = str2num(h_tag_threshold_ic.String);
    if ~(d_tag_threshold_ic > -0.0001)  % rudimentory validation
        % failed to detect a positve number
        error('icatb_gui_al_toolbox: IC threshold was not validated as a number')
    end

    h_tag_n_corrs = findobj(handles, 'tag', 'tag_n_corrs');
    n_tag_n_corrs = str2num(h_tag_n_corrs.String);% 2;
    if mod(n_tag_n_corrs,1) ~= 0   % rudimentory validation
        % failed to detect integer
        error('icatb_gui_al_toolbox: Number of correlations was not validated to an integer')
    end

    h_tag_prefix_label_train = findobj(handles, 'tag', 'tag_prefix_label_train');
    s_tag_prefix_label_train = strtrim(h_tag_prefix_label_train.String);

    h_tag_detect_artifact = findobj(handles, 'tag', 'tag_detect_artifact');
    b_tag_detect_artifact = h_tag_detect_artifact.Value - 1;

    h_tag_skip_anat = findobj(handles, 'tag', 'tag_skip_anat');
    b_tag_skip_anat = h_tag_skip_anat.Value - 1;

    h_tag_skip_functional = findobj(handles, 'tag', 'tag_skip_functional');
    b_tag_skip_functional = h_tag_skip_functional.Value - 1;

    h_tag_anatomic_atlas = findobj(handles, 'tag', 'tag_anatomic_atlas');
    s_tag_anatomic_atlas = strtrim(h_tag_anatomic_atlas.String);    

    h_tag_function_atlas = findobj(handles, 'tag', 'tag_function_atlas');
    s_tag_function_atlas = strtrim(h_tag_function_atlas.String);  

    h_tag_edit_param_or_smap = findobj(handles, 'tag', 'tag_edit_param_or_smap');  
    s_tag_edit_param_or_smap = h_tag_edit_param_or_smap.String;
 
    h_tag_edit_output_folder = findobj(handles, 'tag', 'tag_edit_output_folder');
    s_tag_edit_output_folder = h_tag_edit_output_folder.String;

    h_tag_edit_mask = findobj(handles, 'tag', 'tag_edit_mask');
    s_tag_edit_mask = h_tag_edit_mask.String;

    h_param_vs_smap = findobj(handles, 'tag', 'tag_param_vs_smap');
    

    % GICA example with fbirn dataset
    clear params;
    if h_param_vs_smap.Value == 1
        % param file selected
        params.param_file = strtrim(s_tag_edit_param_or_smap);
    else
        % spat map sel        
        params.sm_path = strtrim(s_tag_edit_param_or_smap);
        params.mask_path = strtrim(s_tag_edit_mask);      
    end
    mkdir([s_tag_edit_output_folder filesep 'nc']); %folder needed for noisecloud
    params.outpath = strtrim(s_tag_edit_output_folder);
    params.n_corr = n_tag_n_corrs;
    params.skip_noise = b_tag_detect_artifact;
    params.skip_anatomical = b_tag_skip_anat;
    params.skip_functional = b_tag_skip_functional;
    params.noise_training_set = strtrim(s_tag_prefix_label_train);
    params.anatomical_atlas = strtrim(s_tag_anatomic_atlas); %aal
    params.threshold = d_tag_threshold_ic;
    params.functional_atlas = strtrim(s_tag_function_atlas); %'yeo_buckner'
    disp( 'Running the autolabeller on c dataset' )
    label_auto_main( params );

    disp('done'); 

function funCallbackHelpScoreCsv(hObject, event_data, handles)
    msg = sprintf(['Skip functional labels of your independent components']);
    disp(msg);
    msgH = helpdlg(msg, 'Auto Labeler Help');   
    
function fun_ic_threshold(hObject, event_data, handles)
    msg = 'Threshold level of independent components';
    disp(msg);
    msgH = helpdlg(msg, 'Auto Labeler Help');

function fun_prefix_label_train(hObject, event_data, handles)
    msg = 'Label file prefix to excel file for noise_training dataset and may be either pre_fbirn_sub or pre_aggregate';
    disp(msg);
    msgH = helpdlg(msg, 'Auto Labeler Help');

function fun_n_corrs(hObject, event_data, handles)
    msg = 'Number of spatial correlation mappings to match your components to';
    disp(msg);
    msgH = helpdlg(msg, 'Auto Labeler Help');

function fun_callback_detect_artifact(hObject, event_data, handles)
    msg = 'Option to ignore search for artifacts.';
    disp(msg);
    msgH = helpdlg(msg, 'Auto Labeler Help');

function funCallbackHelpComponents(hObject, event_data, handles)
    msg = 'Skip anatomical labels of your independent components';
    disp(msg);
    msgH = helpdlg(msg, 'Auto Labeler Help');

function fun_anatomic_atlas(hObject, event_data, handles)
    msg = 'Prefix of anatomic atlas picked from https://github.com/trendscenter/gift/tree/master/GroupICAT/icatb/toolbox/autolabeller/bin/autolabeller/data/Anatomical';
    disp(msg);
    msgH = helpdlg(msg, 'Auto Labeler Help');   

function fun_function_atlas(hObject, event_data, handles)
    msg = 'Prefix of functional atlas picked from https://github.com/trendscenter/gift/tree/master/GroupICAT/icatb/toolbox/autolabeller/bin/autolabeller/data/Functional';
    disp(msg);
    msgH = helpdlg(msg, 'Auto Labeler Help');       

function fun_help_param_or_sm(hObject, event_data, handles)
    msg = 'Decide if you will use a parameter file (*param_info.mat) after a group ICA or a spatial map (4D *.nii). For [Parameter File], popup window will let you select a *param_info.mat file resulting from a previous GIFT Group ICA run. For [Spatial Map] there are 2 steps, you first get a popup window that will let you select a spatial map (nii-file) and then another popup window will appear where you select a spatial mask file (img- or nii-file)';
    disp(msg);
    msgH = helpdlg(msg, 'Auto Labeler Help');     

function fun_param_vs_smap(hObject, event_data, handles)
    icatb_defaults;
    global PARAMETER_INFO_MAT_FILE;
    global UI_FS;
    h_param_vs_smap = findobj(handles, 'tag', 'tag_param_vs_smap');
    h_tag_text_param_or_smap = findobj(handles, 'tag', 'tag_text_param_or_smap');    
    h_tag_edit_param_or_smap = findobj(handles, 'tag', 'tag_edit_param_or_smap');        
    h_tag_text_mask = findobj(handles, 'tag', 'tag_text_mask');    
    h_tag_edit_mask = findobj(handles, 'tag', 'tag_edit_mask');  
    if h_param_vs_smap.Value == 1
        % parameter file
        h_tag_text_param_or_smap.String = 'Parameter File';
        h_tag_text_mask.String = 'Do Not Use';
        h_tag_edit_mask.String = '';
        try
            s_param = icatb_selectEntry('title', 'Select Parameter File', 'typeEntity', 'file', 'typeSelection', 'single', 'filter', ['*', PARAMETER_INFO_MAT_FILE, '.mat'], 'fileType', '*.mat');
            drawnow;
            [inDir, paramF, extn] = fileparts(s_param);
            if (isempty(inDir))
                inDir = pwd;
            end
            s_param = fullfile(inDir, [paramF, extn]);            
            if ~isempty(s_param)
                h_tag_edit_param_or_smap.String = s_param;
            else
                disp('Parameter file not selected. ');
                error('Parameter file not selected. ')
            end
        catch
            icatb_errorDialog(lasterr, 'Mask Error', 'modal');
        end
    else
        %spatial map
        h_tag_text_param_or_smap.String = 'Spatial Mask File';
        h_tag_text_mask.String = 'Mask File';        
        try
            [s_smap] = icatb_selectEntry('filter', '*.nii', 'title', 'Select a spatial map file in nifti format', ...
                'fileType', '*.nii', 'fileNumbers', 1);
            if ~isempty(s_smap)
                h_tag_edit_param_or_smap.String = s_smap;
            else
                disp('Spatial Map file not selected. ');
                error('Spatial Map file not selected. ')
            end
        catch
            icatb_errorDialog(lasterr, 'Mask Error', 'modal');
        end

        %mask
        try
            % select mask in 3D Analyze data
            [s_mask] = icatb_selectEntry('filter', '*.img;*.nii', 'title', 'Select a mask file in analyze or nifti format', ...
                'fileType', '*.img;*.nii', 'fileNumbers', 1);
            if ~isempty(s_mask)
                h_tag_edit_mask.String = s_mask;
            else
                disp('No mask file selected.');
                error('No mask file selected.')
            end
        catch
            icatb_errorDialog(lasterr, 'Mask Error', 'modal');
        end
    end

