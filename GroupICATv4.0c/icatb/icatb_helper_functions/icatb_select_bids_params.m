function bids_info = icatb_select_bids_params(bidsResults)
%% Select BIDS params using GUI
%
%
data_setDir = bidsResults.root_dir;

modality_dir = 'func';
try
    modality_dir = bidsResults.modality_dir;
catch
end

listFiles = 1;
try
    listFiles = bidsResults.listFiles;
catch
end

try
    bidsResults.subjects(1).name;
catch
    [~, bidsResults] = icatb_parseBIDS(struct('root_dir', data_setDir));
end

subjectList = cellstr(char(bidsResults.subjects.name));
[~, ids] = unique(subjectList);
subjectList = subjectList(sort(ids));
sessions = {};
try
    sessions = cellstr(char(bidsResults.subjects.session));
    [~, ids] = unique(sessions);
    sessions = sessions(sort(ids));
catch
end

if (~isempty(sessions))
    sess_chk_ids = ((icatb_good_cells(cellfun(@isempty, sessions, 'UniformOutput', false)))==0);
    sessions = sessions(sess_chk_ids);
end

tasks = {};

try
    tasks = cellstr(char(bidsResults.subjects(1).func.task));
    [~, ids] = unique(tasks);
    tasks = tasks(sort(ids));
catch
end

if (~isempty(tasks))
    task_chk_ids = ((icatb_good_cells(cellfun(@isempty, tasks, 'UniformOutput', false)))==0);
    tasks = tasks(task_chk_ids);
end

InputHandle = icatb_getGraphics('BIDS Parameters', 'normal', 'bids_params', 'on');

listboxWidth = 0.3;
listboxHeight = 0.2;

offset = 0.1;

yPos = 0.96 - 0.5*listboxHeight;
promptHeight = 0.05;
promptWidth = 0.5;
promptPos = [0.05, yPos - 0.5*promptHeight, promptWidth, promptHeight];

% Set no menu bar for the figure
set(InputHandle, 'Menubar', 'none');

tempH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', ...
    promptPos, 'string', 'Select subject ids', 'tag', 'prompt_subject_ids', 'panel_tag', ...
    'panel_prompt_subject_id', 'fontsize', 11, 'visible', 'on');

listboxPos = promptPos;
listboxPos(1) = listboxPos(1) + listboxPos(3) + 0.05;
listboxPos(2) = listboxPos(2) - 0.5*listboxHeight;
listboxPos(3) = listboxWidth;
listboxPos(4) = listboxHeight;
tempH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', ...
    listboxPos, 'string', subjectList, 'tag', 'subject_ids', 'panel_tag', ...
    'panel_subject_id', 'fontsize', 11, 'visible', 'on',  'max', 2, 'value', 1);

if (~isempty(sessions))
    promptPos(2) = promptPos(2) - offset - 0.5*listboxHeight - 0.5*promptHeight;
    tempH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', ...
        promptPos, 'string', 'Select sessions', 'tag', 'prompt_session_ids', 'panel_tag', ...
        'panel_prompt_session_id', 'fontsize', 11, 'visible', 'on');
    
    listboxPos = promptPos;
    listboxPos(1) = listboxPos(1) + listboxPos(3) + 0.05;
    listboxPos(2) = listboxPos(2) - 0.5*listboxHeight;
    listboxPos(3) = listboxWidth;
    listboxPos(4) = listboxHeight;
    tempH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', ...
        listboxPos, 'string', sessions, 'tag', 'session_ids', 'panel_tag', ...
        'panel_session_id', 'fontsize', 11, 'visible', 'on', 'max', 2, 'value', 1);
end


if (~isempty(tasks))
    promptPos(2) = promptPos(2) - offset - 0.5*listboxHeight - 0.5*promptHeight;
    tempH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', ...
        promptPos, 'string', 'Select tasks', 'tag', 'prompt_task_ids', 'panel_tag', ...
        'panel_prompt_task_id', 'fontsize', 11, 'visible', 'on');
    
    listboxPos = promptPos;
    listboxPos(1) = listboxPos(1) + listboxPos(3) + 0.05;
    listboxPos(2) = listboxPos(2) - 0.5*listboxHeight;
    listboxPos(3) = listboxWidth;
    listboxPos(4) = listboxHeight;
    tempH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', ...
        listboxPos, 'string', tasks, 'tag', 'task_ids', 'panel_tag', ...
        'panel_task_id', 'fontsize', 11, 'visible', 'on',  'max', 2, 'value', 1);
end

okWidth = 0.12;
okHeight = 0.05;
okPos = [0.5 - 0.5*okWidth, 0.1 - 0.5*okHeight, okWidth, okHeight];


% define push buttons
okHandle = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', ...
    okPos, 'string', 'Ok', 'Tag', 'bids_ok', 'callback', {@bidsGetParams, InputHandle});


try
    waitfor(InputHandle);
catch
end

appName = 'gift_bids_params';
if (isappdata(0, appName))
    
    bids_info = getappdata(0, appName);
    bids_info.modality_dir = modality_dir;
    bids_info.root_dir = data_setDir;
    rmappdata(0, appName);
    
    if (listFiles)
        bids_info.input_data_file_patterns = icatb_parseBIDS(bids_info);
    end
    
end

function bidsGetParams(hObject, event_data, handles)
%% Get Bids Params



subjectIDH = findobj(handles, 'tag', 'subject_ids');
subjects = get(subjectIDH, 'string');
subject_val = get(subjectIDH, 'value');
subjects = subjects(subject_val);

sessionIDH = findobj(handles, 'tag', 'session_ids');
sessions = {};
if (~isempty(sessionIDH))
    sessions = get(sessionIDH, 'string');
    sessions_val = get(sessionIDH, 'value');
    sessions = sessions(sessions_val);
end

taskIDH = findobj(handles, 'tag', 'task_ids');
tasks = {};
if (~isempty(taskIDH))
    tasks = get(taskIDH, 'string');
    tasks_val = get(taskIDH, 'value');
    tasks = tasks(tasks_val);
end


bids_info.subjects = subjects;

bids_info.sessions = sessions;

bids_info.tasks = tasks;


setappdata(0, 'gift_bids_params', bids_info);

delete(handles);