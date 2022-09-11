function nbic_toolbox(param_file)
    icatb_defaults;
    global UI_FS;
    global PARAMETER_INFO_MAT_FILE;
    global DFNC_DEFAULTS;

    %% Select param file
    if (~exist('param_file', 'var'))
        param_file = icatb_selectEntry('title', 'Select Parameter File', 'typeEntity', 'file', 'typeSelection', 'single', 'filter', ['*', PARAMETER_INFO_MAT_FILE, '.mat;*dfnc.mat']);
        drawnow;
        if (isempty(param_file))
            error('Parameter file is not selected');
        end
    end
    
    drawnow;
    
    [inDir, paramF, extn] = fileparts(param_file);
    if (isempty(inDir))
        inDir = pwd;
    end
    
    param_file = fullfile(inDir, [paramF, extn]);
    
    load(param_file);

    if (exist('sesInfo', 'var'))
        sesInfo.userInput.pwd = inDir;
        sesInfo.outputDir = inDir;
        
        
        [modalityType, dataTitle, compSetFields] = icatb_get_modality;
        modalityType = 'smri';
        
        outputDir = icatb_selectEntry('typeEntity', 'directory', 'title', 'Select output directory to place NBIC results ...'); 
        drawnow; 
        if (isempty(outputDir))
            error('Output directory is not selected');
        end
        icatbInfo.userInput.ica_param_file = param_file;
        icatbInfo.userInput.outputDir = outputDir;
        icatbInfo.userInput.numOfSub = sesInfo.numOfScans; % number of subjects for sbm ;;  % number of subjects for sbm
        icatbInfo.userInput.numOfSess = 1; % only one subject for sbm ;; sesInfo.numOfSess
        icatbInfo.userInput.prefix = sesInfo.userInput.prefix; %
        compFiles = sesInfo.icaOutputFiles(1).ses(1).name;
        icatbInfo.userInput.compFiles = icatb_rename_4d_file(icatb_fullFile('directory', inDir, 'files', compFiles));
        icatbInfo.userInput.numICs = sesInfo.numComp;
        icatbInfo.userInput.HInfo = sesInfo.HInfo.V(1);
    elseif (exist('icatbInfo', 'var'))
        load (icatbInfo.userInput.ica_param_file);
        icatbInfo.userInput.numOfSess = sesInfo.numOfSess;
        clear sesInfo;
        outputDir = inDir;
        icatbInfo.userInput.outputDir = outputDir;
        inDir = fileparts(icatbInfo.userInput.ica_param_file);
    else
        error('Selected file is neither ICA parameter file nor dFNC parameter file');
    end

    drawnow;

    cd(outputDir);

    icatbInfo.userInput.outputDir = outputDir;
    load(icatbInfo.userInput.ica_param_file);

    % if (~exist(deblank(icatbInfo.userInput.compFiles(1, :)), 'file'))
    %     unzipFiles(sesInfo, sesInfo.icaOutputFiles(1).ses(1).name, fileparts(icatbInfo.userInput.ica_param_file));
    % end

    if (~exist(icatb_parseExtn(deblank(icatbInfo.userInput.compFiles(1, :))), 'file'))
        subjectICAFiles = icatb_parseOutputFiles('icaOutputFiles', sesInfo.icaOutputFiles, 'numOfSub', sesInfo.numOfSub, 'numOfSess', sesInfo.numOfSess);
        compFiles = subjectICAFiles(1).ses(1).name;
        icatbInfo.userInput.compFiles = icatb_rename_4d_file(icatb_fullFile('directory', inDir, 'files', compFiles));
        if (~exist(deblank(icatbInfo.userInput.compFiles(1, :)), 'file'))
            unzipFiles(sesInfo, compFiles, fileparts(icatbInfo.userInput.ica_param_file));
        end
    end


    msg = 'Opening N-BIC Setup GUI ...';

    disp(msg);

    msgH = helpdlg(msg, 'Setup NBIC');

    drawnow;

    % Loading components list init below
    matchCompNetworkNames = getCompData(icatbInfo.userInput.compFiles);
    covNames = ''; 
    compGroupNames = ''; 
    if (~isfield(icatbInfo.userInput, 'comp'))
        icatbInfo.userInput.comp = [];
    end
    try 
        compGroupNames = (cellstr(char(icatbInfo.userInput.comp.name))); 
    catch 
    end 
    
    % Loading Groups and subjects initiated below
    struSubjGrpData = 1:icatbInfo.userInput.numOfSub;
    csSubjGrpNames = ''; % Groups and subejcts list init
    if (~isfield(icatbInfo.userInput, 'group')) 
        icatbInfo.userInput.group = []; 
    end
    try 
        csSubjGrpNames = cellstr(char(icatbinfo.userInput.group.name));
    catch 
    end    
    
    % Import scores from CSV file initiated below
    csScoreNames = ''; % Scores (non ica loadings list init
    if (~isfield(icatbInfo.userInput, 'scores'))
        icatbInfo.userInput.scores = [];
    end
    try 
        csScoreNames = (cellstr(char(icatbInfo.userInput.scores.name))); 
    catch 
    end       
    
    clear sesInfo;
    %% Draw graphics % First Main Input GUI 
    figureTag = 'setup_nbic_gui';
    figHandle = findobj('tag', figureTag);
    if (~isempty(figHandle))
        delete(figHandle);
    end

    InputHandle = icatb_getGraphics('N-BIC Setup', 'normal', figureTag, 'off');
    set(InputHandle, 'menubar', 'none');
    set(InputHandle, 'userdata', icatbInfo);

    promptWidth = 0.52; controlWidth = 0.32; promptHeight = 0.05;
    xOffset = 0.02; yOffset = promptHeight; listboxHeight = 0.18; listboxWidth = controlWidth; yPos = 0.9;
    okWidth = 0.12; okHeight = promptHeight;

    promptPos = [xOffset, yPos - 0.5*yOffset, promptWidth, promptHeight];

    listboxXOrigin = promptPos(1) + promptPos(3) + xOffset;
    listboxYOrigin = promptPos(2) - 0.5*listboxHeight;

    %%  select subjects and groups listbox
    %promptPos(2) = listboxYOrigin - yOffset - 0.5*listboxHeight;
    textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Subjects & Groups', 'tag', ...
        'prompt_components', 'fontsize', UI_FS - 1);
    icatb_wrapStaticText(textH);
    listboxYOrigin = promptPos(2) - 0.5*listboxHeight;
    listboxPos = [listboxXOrigin, listboxYOrigin, listboxWidth, listboxHeight];

    listboxPos = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', listboxPos, 'string', csSubjGrpNames, 'tag', ...
    'group', 'fontsize', UI_FS - 1, 'min', 0, 'max', 1, 'value', 1, 'callback', {@fAddGroups, InputHandle});
    rowdButPos = [promptPos(1)+promptPos(3)/2, promptPos(2)-1*okHeight, promptHeight + 0.01, promptHeight - 0.01];
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', rowdButPos , 'string', '?', 'tag', 'tagHelpSubjects', 'fontsize',...
        UI_FS - 1, 'callback', {@funCallbackHelpSubjects});

    addButtonPos = [listboxXOrigin + listboxWidth + xOffset, listboxYOrigin + 0.5*listboxHeight + 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
    removeButtonPos = [listboxXOrigin + listboxWidth + xOffset, listboxYOrigin + 0.5*listboxHeight - 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', addButtonPos, 'string', '+', 'tag', 'tagAdd_group_button', 'fontsize',...
        UI_FS - 1, 'callback', {@fAddGroups, InputHandle});
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', removeButtonPos, 'string', '-', 'tag', 'tagRemove_group_button', 'fontsize',...
        UI_FS - 1, 'callback', {@fRemoveGroups, InputHandle});
    
    
    %%  Components listbox (Group components by name)
    textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos + [0 -4.25*okHeight 0 0], 'string', 'Add Components', 'tag', ...
        'prompt_components', 'fontsize', UI_FS - 1);
    icatb_wrapStaticText(textH);
    listboxYOrigin = promptPos(2) - 0.5*listboxHeight - 4.25*okHeight;
    listboxPos = [listboxXOrigin, listboxYOrigin, listboxWidth, listboxHeight];
    listCompH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', listboxPos, 'string', compGroupNames, 'tag', ...
        'comp', 'fontsize', UI_FS - 1, 'min', 0, 'max', 1, 'value', 1, 'callback', {@addCompNetwork, InputHandle}, 'userdata', matchCompNetworkNames);
    rowdButPos = [promptPos(1)+promptPos(3)/2, promptPos(2)-5.25*okHeight, promptHeight + 0.01, promptHeight - 0.01];
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', rowdButPos , 'string', '?', 'tag', 'tagHelpComponents', 'fontsize',...
        UI_FS - 1, 'callback', {@funCallbackHelpComponents});
    
    addButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) + 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
    removeButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) - 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', addButtonPos, 'string', '+', 'tag', 'add_comp_button', 'fontsize',...
        UI_FS - 1, 'callback', {@addCompNetwork, InputHandle});
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', removeButtonPos, 'string', '-', 'tag', 'remove_comp_button', 'fontsize',...
        UI_FS - 1, 'callback', {@removeCompNetwork, InputHandle});


    
    
    
     %%  CSV Scores listbox 
    textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos + [0 -8.5*okHeight 0 0], 'string', 'Import Scores from CSV-File', 'tag', ...
        'prompt_components', 'fontsize', UI_FS - 1);
    icatb_wrapStaticText(textH);
    listboxYOrigin = promptPos(2) - 0.5*listboxHeight - 8.5*okHeight;
    listboxPos = [listboxXOrigin, listboxYOrigin, listboxWidth, listboxHeight];
    listCompH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', listboxPos, 'string', csScoreNames, 'tag', ...
        'tagScore', 'fontsize', UI_FS - 1, 'min', 0, 'max', 1, 'value', 1, 'callback', {@fLoadCsvScore, InputHandle});
    
    rowdButPos = [promptPos(1)+promptPos(3)/2, promptPos(2)-9.5*okHeight, promptHeight + 0.01, promptHeight - 0.01];
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', rowdButPos , 'string', '?', 'tag', 'tagHelpComponents', 'fontsize',...
        UI_FS - 1, 'callback', {@funCallbackHelpScoreCsv});
    
    addButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) + 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
    removeButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) - 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', addButtonPos, 'string', 'L', 'tag', 'import_score_button', 'fontsize',...
        UI_FS - 1, 'callback', {@fLoadCsvScore, InputHandle});
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', removeButtonPos, 'string', '-', 'tag', 'remove_score_button', 'fontsize',...
        UI_FS - 1, 'callback', {@fRemoveCsvScore, InputHandle});
    

    textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos + [0 -12*okHeight 0 0], ...
        'string', 'X Size', 'tag', 'prompt_components', 'fontsize', UI_FS - 1);
    icatb_wrapStaticText(textH);
    editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', ...
        'position', promptPos + [+.54 -12*okHeight -.2 0], 'String', '50', 'fontsize', UI_FS - 1, 'tag', 'tagEditXSize');
    %inputText(numParameters).help = struct('title', 'Output Files', 'string', 'All the output files will be preprended with this prefix.');
    %promptPos = [xOffset, yPos - 0.5*yOffset, promptWidth, promptHeight];
    %addButtonPos = [listboxXOrigin + listboxWidth + xOffset, listboxYOrigin + 0.5*listboxHeight + 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
    rowdButPos = [listboxPos(1) + listboxPos(3) + xOffset, promptPos(2)-12*okHeight, promptHeight + 0.01, promptHeight - 0.01];
    % rowdButPos = promptPos + [+promptWidth -12*okHeight -promptWidth+.05 0];
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', rowdButPos, 'string', '?', 'tag', 'tagHelpXSize', 'fontsize',...
        UI_FS - 1, 'callback', {@funCallbackHelpXSize});

    textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos + [0 -13*okHeight 0 0], ...
        'string', 'Y Size', 'tag', 'prompt_components', 'fontsize', UI_FS - 1);
    icatb_wrapStaticText(textH);
    editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', ...
        'position', promptPos + [+.54 -13*okHeight -.2 0], 'String', '2', 'fontsize', UI_FS - 1, 'tag', 'tagEditYSize');
    rowdButPos = [listboxPos(1) + listboxPos(3) + xOffset, promptPos(2)-13*okHeight, promptHeight + 0.01, promptHeight - 0.01];
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', rowdButPos , 'string', '?', 'tag', 'tagHelpYSize', 'fontsize',...
        UI_FS - 1, 'callback', {@funCallbackHelpYSize});
 
    textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos + [0 -14*okHeight 0 0], ...
        'string', 'Tolerance', 'tag', 'prompt_components', 'fontsize', UI_FS - 1);
    icatb_wrapStaticText(textH);  
    editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', ...
        'position', promptPos + [+.54 -14*okHeight -.2 0], 'String', '20', 'fontsize', UI_FS - 1, 'tag', 'tagEditTolerance');
    rowdButPos = [listboxPos(1) + listboxPos(3) + xOffset, promptPos(2)-14*okHeight, promptHeight + 0.01, promptHeight - 0.01];
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', rowdButPos , 'string', '?', 'tag', 'tagHelpTolerance', 'fontsize',...
        UI_FS - 1, 'callback', {@funCallbackHelpTolerance});
 
    textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos + [0 -15*okHeight 0 0], ...
        'string', 'Repetitions', 'tag', 'prompt_components', 'fontsize', UI_FS - 1);
    icatb_wrapStaticText(textH);
    editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', ...
        'position', promptPos + [+.54 -15*okHeight -.2 0], 'String', '5', 'fontsize', UI_FS - 1, 'tag', 'tagEditReps');
    rowdButPos = [listboxPos(1) + listboxPos(3) + xOffset, promptPos(2)-15*okHeight, promptHeight + 0.01, promptHeight - 0.01];
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', rowdButPos , 'string', '?', 'tag', 'tagHelpReps', 'fontsize',...
        UI_FS - 1, 'callback', {@funCallbackHelpReps});
    


    
    
    %% Add cancel, save and run buttons
    cancelPos = [0.25 - 0.5*okWidth, promptPos(2) - yOffset - 15.5*okHeight, okWidth*2, okHeight];
    cancelPos(2) = cancelPos(2) - 0.5*cancelPos(4);
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', cancelPos, 'string', 'Cancel', 'tag', 'cancel_button', 'fontsize',...
        UI_FS - 1, 'callback', 'delete(gcbf);');

    okPos = [0.75 - 0.5*okWidth, promptPos(2) - yOffset - 15.5*okHeight, okWidth*2, okHeight];
    okPos(2) = okPos(2) - 0.5*okPos(4);
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', okPos, 'string', 'Run', 'tag', 'run_button', 'fontsize',...
        UI_FS - 1, 'callback', {@runCallback, InputHandle});

%     savePos = [0.5 - 0.5*okWidth, promptPos(2) - yOffset - 15.5*okHeight, okWidth, okHeight];
%     savePos(2) = savePos(2) - 0.5*savePos(4);
%     icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', savePos, 'string', 'Save', 'tag', 'save_button', 'fontsize',...
%         UI_FS - 1, 'callback', {@saveCallback, InputHandle});

    try
        delete(msgH);
    catch
    end

    set(InputHandle, 'visible', 'on');
    drawnow;

    
function matchCompNetworkNames = getCompData(file_names)

    structFile = deblank(file_names(1, :));

    %% Get colormap associated with the image values
    structData2 =  icatb_spm_read_vols(icatb_spm_vol(structFile));
    structData2(isfinite(structData2) == 0) = 0;
    structDIM = [size(structData2, 1), size(structData2, 2), 1];

    for nC = 1:size(file_names, 1)

        tmp = icatb_spm_read_vols(icatb_spm_vol(file_names(nC, :)));
        tmp(isfinite(tmp)==0) = 0;

        tmp(tmp ~= 0) = detrend(tmp(tmp ~= 0), 0) ./ std(tmp(tmp ~= 0));
        tmp(abs(tmp) < 1.0) = 0;

        if (nC == 1)
            compData = zeros(size(tmp, 1), size(tmp, 2), size(file_names, 1));
        end

        [dd, inds] = max(tmp(:));

        [x, y, z] = ind2sub(size(tmp), inds);

        [tmp, maxICAIM, minICAIM, minInterval, maxInterval] = icatb_overlayImages(reshape(tmp(:, :, z), [1, size(tmp, 1), size(tmp, 2), 1]), structData2(:, :, z), structDIM, structDIM, 1);

        compData(:, :, nC) = reshape(tmp, structDIM);

    end


    clim = [minInterval, 2*maxInterval];
    cmap = icatb_getColormap(1, 1, 1);

    matchCompNetworkNames.clim = clim;
    matchCompNetworkNames.cmap = cmap;
    matchCompNetworkNames.compData = compData;

function fLoadCsvScore(hObject, event_data, figH)
    % Imports Scores from CSV file

    %Initiations
    icatbInfo = get(figH, 'userdata');
    listH = findobj(figH, 'tag', 'tagScore');
    val = get(listH, 'value');

    icatb_defaults;
    global UI_FS;

    figureTag = 'tagImportCsvScore';
    figNbic = findobj('tag', figureTag);
    if (~isempty(figNbic))
        delete(figNbic);
    end
    
    ScoreVal = [];
    ScoreName = '';
    if (listH == hObject)
            return;
    end

    %Select CSV file
    chCsvFile = icatb_selectEntry('title', 'Select CSV file, with first row indicating score names and remaining rows loadings per subject', 'typeEntity', 'file', 'typeSelection', 'single', 'filter', ['*.csv;*.txt']);
    drawnow;
    if (isempty(chCsvFile))
        error('CSV file is not selected');
    end
    
    tabScoresCsv = readtable(chCsvFile); 

    [nCsvSubjs nScores] = size(tabScoresCsv);
    
    if ~(nCsvSubjs == icatbInfo.userInput.numOfSub)
        error(['CSV file list ' num2str(nCsvSubjs) ' subjects, but ICA parameter file indicates ' num2str(icatbInfo.userInput.numOfSub) ' subjects']);
    end
    
    try
        csScoreNames = tabScoresCsv.Properties.VariableNames;
        icatbInfo.userInput.scores = [];
        for nScore = 1:nScores
            ind = length(icatbInfo.userInput.scores) + 1;
            chScoreName = char(csScoreNames(nScore));
            coldScoreVals = table2array(tabScoresCsv(:,nScore));

            % Set user selected information in figure
            icatbInfo.userInput.scores(ind).name = chScoreName;
            icatbInfo.userInput.scores(ind).val =  coldScoreVals;
            set(figH, 'userdata', icatbInfo);
            groupListH = findobj(figH, 'tag', 'tagScore');
            set(groupListH, 'string', cellstr(char(icatbInfo.userInput.scores.name)));
        end
    catch
        icatb_errorDialog(lasterr, 'Problem reading in scores from CSV file');
    end
    
    %Save subjectID's to match if column is provided
    for n = 1:length(csScoreNames)
        if strcmpi(lower(char(csScoreNames(n))),'subjectidmatch')
            coldSubjectIdMatch = table2array(tabScoresCsv(:,n));
            save([icatbInfo.userInput.outputDir '/SubjectIdMatch.mat'],'coldSubjectIdMatch'); %save matchin info to sesInfo.userInput.files.name in *parameter_info.mat file
        end
    end

function fRemoveCsvScore(hObject, event_data, figH)
    % Remove single CSV score
    %

    % Initiations
    icatbInfo = get(figH, 'userdata');
    listH = findobj(figH, 'tag', 'tagScore');
    val = get(listH, 'value');
    strs = cellstr(get(listH, 'string'));

    if (~isempty(strs))
        check = icatb_questionDialog('title', 'Remove Score', 'textbody', ['Are you sure you want to remove the score ', strs{val}, ' from the list?']);
        if (~check)
            return;
        end
    end

    try
        strs = cellstr(char(icatbInfo.userInput.scores.name));
        icatbInfo.userInput.scores(val) = [];
        strs(val) = [];
        set(listH, 'value', 1);
        set(listH, 'string', strs);
        set(figH, 'userdata', icatbInfo);
    catch
    end
    
function fAddGroups(hObject, event_data, figH)
    % Add groups

    icatbInfo = get(figH, 'userdata');
    listH = findobj(figH, 'tag', 'group');
    val = get(listH, 'value');

    icatb_defaults;
    global UI_FS;

    figureTag = 'tagAddGroupsNbic';
    figNbic = findobj('tag', figureTag);
    if (~isempty(figNbic))
        delete(figNbic);
    end

    groupVal = [];
    groupName = '';
    if (listH == hObject)
        if (~strcmpi(get(figH, 'selectionType'), 'open'))
            return;
        end
        val = get(listH, 'value');
        try
            groupName = icatbInfo.userInput.group(val).name;
            groupVal = icatbInfo.userInput.group(val).val;
        catch
        end
    end

    subjectString = cellstr([repmat('Subject ', icatbInfo.userInput.numOfSub, 1), num2str((1:icatbInfo.userInput.numOfSub)')]);

    %[groupName, groupVal] = icatb_select_groups_gui(subjectString, groupName, 'select_subjects', groupVal);
    [groupName, groupVal] = icatb_select_groups_gui(subjectString, 'Group', 'select_subjects', groupName, groupVal);

    try

        if (isempty(groupName))
            error('Group name is not selected');
        end

        if (isempty(groupVal))
            error('Subjects are not selected');
        end

        if (length(icatbInfo.userInput.group) > 0)
            chk = strmatch(lower(groupName), lower(cellstr(char(icatbInfo.userInput.group.name))), 'exact');
            if (~isempty(chk))
                ind = chk;
            end
        end

        if (~exist('ind', 'var'))
            ind = length(icatbInfo.userInput.group) + 1;
        end

        %% Set user selected information in figure
        icatbInfo.userInput.group(ind).name = groupName;
        icatbInfo.userInput.group(ind).val =  groupVal;
        set(figH, 'userdata', icatbInfo);
        groupListH = findobj(figH, 'tag', 'group');
        set(groupListH, 'string', cellstr(char(icatbInfo.userInput.group.name)));

    catch
        icatb_errorDialog(lasterr, 'Group Selection');
    end



function fRemoveGroups(hObject, event_data, figH)
    % Remove groups
    %

    icatbInfo = get(figH, 'userdata');
    listH = findobj(figH, 'tag', 'group');
    val = get(listH, 'value');
    strs = cellstr(get(listH, 'string'));

    if (~isempty(strs))
        check = icatb_questionDialog('title', 'Remove groups', 'textbody', ['Do you want to remove the group ', strs{val}, ' from the list?']);
        if (~check)
            return;
        end
    end

    try
        strs = cellstr(char(icatbInfo.userInput.group.name));
        icatbInfo.userInput.group(val) = [];
        strs(val) = [];
        set(listH, 'value', 1);
        set(listH, 'string', strs);
        set(figH, 'userdata', icatbInfo);
    catch
    end
    
function addCompNetwork(hObject, event_data, figH)
    %% Add Component network
    %

    icatb_defaults;
    global UI_FS;

    figureTag = 'add_comp_dfnc';
    compFigHandle = findobj('tag', figureTag);
    if (~isempty(compFigHandle))
        delete(compFigHandle);
    end

    icatbInfo = get(figH, 'userdata');

    listH = findobj(figH, 'tag', 'comp');

    compVals = [];
    networkName = '';
    if (listH == hObject)
        if (~strcmpi(get(figH, 'selectionType'), 'open'))
            return;
        end
        val = get(listH, 'value');
        try
            networkName = icatbInfo.userInput.comp(val).name;
            compVals = icatbInfo.userInput.comp(val).value;
        catch
        end
    end

    compStr = num2str((1:icatbInfo.userInput.numICs)');

    compFigHandle = icatb_getGraphics('Select Component Networks', 'normal', figureTag);
    set(compFigHandle, 'menubar', 'none');

    promptWidth = 0.52; controlWidth = 0.32; promptHeight = 0.05;
    xOffset = 0.02; yOffset = promptHeight; listboxHeight = 0.6; yPos = 0.9;
    okWidth = 0.12; okHeight = promptHeight;

    %% Features text and listbox
    promptPos = [xOffset, yPos - 0.5*yOffset, promptWidth, promptHeight];
    textH = icatb_uicontrol('parent', compFigHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter Network Name', 'fontsize', UI_FS - 1);
    icatb_wrapStaticText(textH);

    promptPos = get(textH, 'position');

    editPos = promptPos;
    editPos(1) = editPos(1) + editPos(3) + xOffset;
    editPos(3) = controlWidth;
    icatb_uicontrol('parent', compFigHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', networkName, 'tag', 'comp_network_name', 'fontsize', UI_FS - 1);

    %% Right Listbox
    listbox2Wdith = 0.1;
    promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
    promptPos(1) = (1 - xOffset - 2*listbox2Wdith);
    promptPos(3) = 2*listbox2Wdith;
    textH = icatb_uicontrol('parent', compFigHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Components', 'fontsize', UI_FS - 1);
    icatb_wrapStaticText(textH);
    promptPos = get(textH, 'position');
    listboxYOrigin = promptPos(2) - 0.5*yOffset - listboxHeight;
    listboxXOrigin = promptPos(1) + 0.5*listbox2Wdith;
    listboxPos = [listboxXOrigin, listboxYOrigin, listbox2Wdith, listboxHeight];
    compListH = icatb_uicontrol('parent', compFigHandle, 'units', 'normalized', 'style', 'listbox', 'position', listboxPos, 'string', compStr, 'tag', 'components', 'fontsize', UI_FS - 1, 'min', 0, 'max', 2, 'value', compVals);

    %% Show components
    showWidth = 0.08; showHeight = 0.04;
    showButtonPos = [listboxXOrigin + 0.5*listbox2Wdith - 0.5*showWidth, listboxYOrigin - yOffset - 0.5*showHeight, showWidth, showHeight];
    showH = icatb_uicontrol('parent', compFigHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', showButtonPos, 'string', 'Show', 'fontsize', UI_FS - 1, 'callback', {@drawComp, figH, compFigHandle});

    %% Plot image on the left hand side
    axesPos = [xOffset, listboxYOrigin, listboxHeight, listboxHeight];
    axesH = axes('parent', compFigHandle, 'units', 'normalized', 'position', axesPos, 'tag', 'axes_display_comp');

    promptPos = axesPos;

    %% Add cancel and run buttons
    okPos = [0.5 - 0.5*okWidth, promptPos(2) - yOffset - 0.5*okHeight, okWidth, okHeight];
    okPos(2) = okPos(2) - 0.5*okPos(4);
    icatb_uicontrol('parent', compFigHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', okPos, 'string', 'Done', 'tag', 'done_button', 'fontsize', UI_FS - 1, 'callback', {@setCompCallback, ...
        compFigHandle, figH});

    %% Draw components on the left hand side
    drawComp(compListH, [], figH, compFigHandle);

function removeCompNetwork(hObject, event_data, figH)
    %% Remove Component network
    %

    icatbInfo = get(figH, 'userdata');
    listH = findobj(figH, 'tag', 'comp');
    val = get(listH, 'value');

    if (~isempty(val))
        check = icatb_questionDialog('title', 'Remove Component Network', 'textbody', 'Do you want to remove the component network from the list?');
        if (~check)
            return;
        end
    end

    try
        strs = cellstr(char(icatbInfo.userInput.comp.name));
        icatbInfo.userInput.comp(val) = [];
        strs(val) = [];
        set(listH, 'value', 1);
        set(listH, 'string', strs);
        set(figH, 'userdata', icatbInfo);
    catch
    end

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
    matchCompNetworkNames = get(listH, 'userdata');

    axesH = get(compFigHandle, 'currentaxes');

    clim = matchCompNetworkNames.clim;
    cmap = matchCompNetworkNames.cmap;
    compData = matchCompNetworkNames.compData;

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


function unzipFiles(sesInfo, compFiles, inDir)
%% Unzip files

zipContents.zipFiles = {};
zipContents.files_in_zip(1).name = {};
if isfield(sesInfo, 'zipContents')
    zipContents = sesInfo.zipContents;
end

zipFile = icatb_getViewingSet_zip(compFiles, [], 'real', zipContents);
if (~isempty(zipFile))
    icatb_unzip(fullfile(inDir, zipFile), inDir);
end

function saveCallback(hObject, event_data, handles)
% Save the icatbInfo
    disp('Not Implemented'); %ce071122

function runCallback(hObject, event_data, handles)
    % Run NBIC
    icatbInfo = get(handles, 'userdata');
    chLoadingsFile = [fileparts(icatbInfo.userInput.ica_param_file), '/' , icatbInfo.userInput.prefix, '_ica_c1-1.mat'];
    load(chLoadingsFile);
    
    % Get the NBIC init variables
    htagEditXSize = findobj(handles, 'tag', 'tagEditXSize');
    x_size = str2num(htagEditXSize.String); %50
    htagEditYSize = findobj(handles, 'tag', 'tagEditYSize');
    y_size = str2num(htagEditYSize.String);% 2;
    htagEditTolerance = findobj(handles, 'tag', 'tagEditTolerance');
    tolerance = str2num(htagEditTolerance.String) ; %2
    htagEditReps = findobj(handles, 'tag', 'tagEditReps');
    reps = str2num(htagEditReps.String); % repititions
    
    if size(int16(char(icatbInfo.userInput.comp.value)),2) > 1
        icatb_dialogBox('title', 'Incorrect Components', 'textBody', 'Detected that single component entry has more than one components in it which is not supported by NBIC. Please make sure each component entry only consists of a single component entry.', 'textType', 'large');
        disp('Detected that single component entry has more than one components in it which is not supported by NBIC. Please make sure each component entry only consists of a single component entry.');
        return;
    end
    
    matRawData = tc(:,int16(char(icatbInfo.userInput.comp.value))); %selects the components chosen
    csComNames = cellstr(char(icatbInfo.userInput.comp.name)); %store the chosen component names
    csSubjGrpNames = (cellstr(char(icatbInfo.userInput.group.name)));
      
    nScores = size(icatbInfo.userInput.scores,2);
    
    [nSubjects nIcaLoadings] = size(matRawData);
    for nScore = 1:nScores
        % Set user selected information in figure
        matRawData(:,nIcaLoadings+nScore) = icatbInfo.userInput.scores(nScore).val;
        csComNames(nIcaLoadings+nScore) = {icatbInfo.userInput.scores(nScore).name};
    end
    
    variable_Ids = 1:size(matRawData,2);
    
    biclusters = NBiC(x_size, y_size, tolerance, variable_Ids, reps, matRawData);
    
    % Sort NBIC numbers for presentation
    iGrps = length(csSubjGrpNames);
    iClusts = length(biclusters);
    % Get all Frequencies for sorting
    coliFreq = zeros(iClusts,1);
    for iClust = 1:iClusts % number of bi-clusters
        coliFreq(iClust) = sum(biclusters(iClust).freq);
    end
    [dummy, coliFreqSortIx] = sort(coliFreq, 'descend');
    coliFreqSortIx = coliFreqSortIx';
    matdFigLoadsClustGrpMean = zeros(iClusts, iGrps);
    matdFigLoadsClustGrpSd = zeros(iClusts, iGrps);
    matdFigLoadsClustGrpMeanFig = NaN(iClusts*(iGrps+1)-1,iGrps);
    matdFigLoadsClustGrpSdFig = NaN(iClusts*(iGrps+1)-1,iGrps);
    csXLabel = {};
    for iClust = coliFreqSortIx % Loop bi-clusters in freq order
        coliComClust = biclusters(iClust).comps;
        colixSubjsClust = biclusters(iClust).subs;
        colbSubjsAllClust = zeros(size(matRawData,1),1);
        colbSubjsAllClust(colixSubjsClust) = 1;
        for iGrp = 1:iGrps % number of bi-clusters
            chGroupName = csSubjGrpNames(iGrp);
            colixSubjsClustGrp = find(colbSubjsAllClust(icatbInfo.userInput.group(iGrp).val') == 1);
            if (min(size(matRawData(colixSubjsClustGrp, coliComClust))) > 1) %Only reshape if 2dim
                coldTmpFigLoadsClustGrp = reshape(matRawData(colixSubjsClustGrp, coliComClust),length(colixSubjsClustGrp)*length(coliComClust),1);
            else
                coldTmpFigLoadsClustGrp = matRawData(colixSubjsClustGrp, coliComClust);
            end      
            if ~isempty(colixSubjsClustGrp)
                matdFigLoadsClustGrpMean(iClust,iGrp) = mean(coldTmpFigLoadsClustGrp);
                matdFigLoadsClustGrpMeanFig(iGrp+((iClust-1)*(iGrps + 1)),iGrp) = mean(coldTmpFigLoadsClustGrp);
                matdFigLoadsClustGrpSd(iClust,iGrp) = std(coldTmpFigLoadsClustGrp);
                matdFigLoadsClustGrpSdFig(iGrp+((iClust-1)*(iGrps + 1)),iGrp) = mean(coldTmpFigLoadsClustGrp);
            else
                matdFigLoadsClustGrpMean(iClust,iGrp) = nan; 
                matdFigLoadsClustGrpSd(iClust,iGrp) = nan; 
            end
        end
        csXLabel{end + 1} = [num2str(coliFreq(iClust)) '; ' strjoin({csComNames(coliComClust')},',')];
    end  
    
    % Generate NBIC graph
    figNbic = figure;
    rowchColors = 'brymgcbrymgcbrymgcbrymgcbrymgcbrymgcbrymgcbrymgcbrymgcbrymgcbrymgcbrymgc'; %Only 6 colors for groups then repeat
    for iGrp = 1:iGrps
        errorbar(matdFigLoadsClustGrpMeanFig(:,iGrp),matdFigLoadsClustGrpSdFig(:,iGrp),['x' rowchColors(iGrp)]);
        hold on;
    end
    title('N-BIC Averages and SD per Group');
    ylabel('ICA Loadings and/or Scores');
    xlabel('Frequency; Cluster components and/or Score Types');
    xticks((0:(iGrps+1):iClusts*(iGrps+1)+1)+(iGrps+1)/2)
    xlim([0,(iClusts*(iGrps+1))]);
    grid on
    xticklabels(csXLabel);
    legend(csSubjGrpNames); 
    % Save nbic graph
    striTime = datestr(now,'yyyymmddHHMMss');
    saveas(figNbic, [icatbInfo.userInput.outputDir '/nbic' striTime '.png']); %worked w extra c
    savefig([icatbInfo.userInput.outputDir '/nbic' striTime '.fig']);
    % Save nbic data
    struNbic = struct;
    struNbic.grpInfo = icatbInfo.userInput.group;
    struNbic.biclusters = biclusters;
    struNbic.matdFigLoadsClustGrpMean = matdFigLoadsClustGrpMean;
    struNbic.matdFigLoadsClustGrpSd = matdFigLoadsClustGrpSd;
    save([icatbInfo.userInput.outputDir '/nbic' striTime '.mat'],'struNbic');
    
    disp('done'); 

function funCallbackHelpScoreCsv(hObject, event_data, handles)
    msg = sprintf(['[L] button loads CSV file for scores to complement ICA loadings. \n' ...
        'Scores from CSV file may be pruned away scores using the [-] button.\n' ...
        'The CSV file needs to have a first row with name of each score and \n' ...
        'then the imported scores will fill the remaining rows and have to be \n' ... 
        ' in the same subject order as GIFT read the image files.\n' ...
        'It may help to add a column, having score name SUBJECTIDMATCH, which\n'...
        ' will save the order of your subjectids in the file SubjectIdMatch.mat\n' ...
        ' and may be compared with variable sesInfo.userInput.files.name in\n' ...
        ' your *parameter_info.mat file.']);
    disp(msg);
    msgH = helpdlg(msg, 'NBIC Help');   
    
function funCallbackHelpXSize(hObject, event_data, handles)
    msg = 'Minimum number of components to be accepted for bicluster formation';
    disp(msg);
    msgH = helpdlg(msg, 'NBIC Help');
    
function funCallbackHelpReps(hObject, event_data, handles)
    msg = 'Number of expected repetitions. The input component ids would be randomly permutated this many times';
    disp(msg);
    msgH = helpdlg(msg, 'NBIC Help');

function funCallbackHelpTolerance(hObject, event_data, handles)
    msg = 'Allowed percentage of overlap between two biclusters ';
    disp(msg);
    msgH = helpdlg(msg, 'NBIC Help');

function funCallbackHelpYSize(hObject, event_data, handles)
    msg = 'Minimum number of subjects to be accepted for bicluster formation';
    disp(msg);
    msgH = helpdlg(msg, 'NBIC Help');

function funCallbackHelpSubjects(hObject, event_data, handles)
    msg = 'Define groups to add subjects into.';
    disp(msg);
    msgH = helpdlg(msg, 'NBIC Help');

function funCallbackHelpComponents(hObject, event_data, handles)
    msg = 'When naming components to add in the NBIC model you may have one component per name. If a component name contains more than one components, NBIC may fail at a later stage.';
    disp(msg);
    msgH = helpdlg(msg, 'NBIC Help');