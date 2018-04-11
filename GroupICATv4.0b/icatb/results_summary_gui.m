function varargout = results_summary_gui(varargin)
% RESULTS_SUMMARY_GUI MATLAB code for results_summary_gui.fig
%      RESULTS_SUMMARY_GUI, by itself, creates a new RESULTS_SUMMARY_GUI or raises the existing
%      singleton*.
%
%      H = RESULTS_SUMMARY_GUI returns the handle to a new RESULTS_SUMMARY_GUI or the handle to
%      the existing singleton*.
%
%      RESULTS_SUMMARY_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RESULTS_SUMMARY_GUI.M with the given input arguments.
%
%      RESULTS_SUMMARY_GUI('Property','Value',...) creates a new RESULTS_SUMMARY_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before results_summary_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to results_summary_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help results_summary_gui

% Last Modified by GUIDE v2.5 10-Feb-2016 14:41:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @results_summary_gui_OpeningFcn, ...
    'gui_OutputFcn',  @results_summary_gui_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before results_summary_gui is made visible.
function results_summary_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to results_summary_gui (see VARARGIN)

% Choose default command line output for results_summary_gui
handles.output = hObject;

setInputs(handles);

num_subjects = 0;

for n = 1:2:length(varargin)
    if (strcmpi(varargin{n}, 'num_subjects'))
        num_subjects = varargin{n + 1};
    end
end

handles.num_subjects = num_subjects;

% if (length(varargin) > 0)
%     handles.compFiles = varargin{1}{1};
%     structFile = fullfile (fileparts(which('gift.m')), 'icatb_templates', 'nsingle_subj_T1_2_2_5.nii');
%     for nComp  = 1:size(handles.compFiles, 1)
%         hD = icatb_orth_views(deblank(handles.compFiles(nComp, :)), 'image_values', 'positive', ...
%             'convert_to_zscores', 'yes', 'get_interp_data', 1, 'set_to_max_voxel', 1, 'threshold', 1, 'structfile', structFile);
%         %'structfile', which('ch2bet.nii')
%         tmp = (squeeze(hD.data(:, :, hD.maxVoxelPos(end))));
%         if (nComp == 1)
%             compData = zeros([size(tmp), size(handles.compFiles, 1)]);
%         end
%         compData(:, :, nComp) = tmp;
%     end
%     handles.compData = compData;
%
%
%     getInputs(varargin{1}{2}, handles);
%
% end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes results_summary_gui wait for user response (see UIRESUME)
% uiwait(handles.results_summary_gui);


% --- Outputs from this function are returned to the command line.
function varargout = results_summary_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

waitfor(hObject);

results = [];
if (isappdata(0, 'resultsAppData'))
    results = getappdata(0, 'resultsAppData');
    rmappdata(0, 'resultsAppData');
end

varargout{1} = results;


% --- Executes on button press in results_format.
function results_format_Callback(hObject, eventdata, handles)
% hObject    handle to results_format (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% compFiles = handles.compFiles;
% if (isfield(handles, 'comps'))
%     comps = get(handles.results_format, 'userdata');
% else
%     comps = (1:size(compFiles, 1));
% end
% comps = selectComps(handles.compData, comps);
% set(handles.results_format, 'userdata', comps);
% set(handles.results_format, 'foregroundcolor', [0, 1, 0]);
% %handles.results_format = results_format;

function threshold_Callback(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshold as text
%        str2double(get(hObject,'String')) returns contents of threshold as a double


% --- Executes during object creation, after setting all properties.
function threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function anatomical_file_Callback(hObject, eventdata, handles)
% hObject    handle to anatomical_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of anatomical_file as text
%        str2double(get(hObject,'String')) returns contents of anatomical_file as a double


% --- Executes during object creation, after setting all properties.
function anatomical_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to anatomical_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in image_values.
function image_values_Callback(hObject, eventdata, handles)
% hObject    handle to image_values (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns image_values contents as cell array
%        contents{get(hObject,'Value')} returns selected item from image_values


% --- Executes during object creation, after setting all properties.
function image_values_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_values (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_file.
function browse_file_Callback(hObject, eventdata, handles)
% hObject    handle to browse_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

startPath = fileparts(which('groupica.m'));
startPath = fullfile(startPath, 'icatb_templates');

oldDir = pwd;

if (~exist(startPath, 'dir'))
    startPath = pwd;
end

structFile = icatb_selectEntry('typeEntity', 'file', 'title', 'Select anatomical file to overlay components ...', 'filter', '*.img;*.nii', 'filetype', 'image', 'typeSelection', 'single', ...
    'filenumbers', 1, 'startpath', startPath);
drawnow;

if (isempty(structFile))
    cd(oldDir);
    error('Structural file is not selected');
end

cd(oldDir);

set(handles.anatomical_file, 'string', structFile);

opts = lower(get(handles.anatomical_plane, 'string'));
val = get(handles.anatomical_plane, 'value');
updateSlices(handles, structFile, opts{val});

% --- Executes on selection change in convert_to_zscores.
function convert_to_zscores_Callback(hObject, eventdata, handles)
% hObject    handle to convert_to_zscores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns convert_to_zscores contents as cell array
%        contents{get(hObject,'Value')} returns selected item from convert_to_zscores


% --- Executes during object creation, after setting all properties.
function convert_to_zscores_CreateFcn(hObject, eventdata, handles)
% hObject    handle to convert_to_zscores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function slices_in_mm_Callback(hObject, eventdata, handles)
% hObject    handle to slices_in_mm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of slices_in_mm as text
%        str2double(get(hObject,'String')) returns contents of slices_in_mm as a double


% --- Executes during object creation, after setting all properties.
function slices_in_mm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slices_in_mm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in group_info.
function group_info_Callback(hObject, eventdata, handles)
% hObject    handle to group_info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns group_info contents as cell array
%        contents{get(hObject,'Value')} returns selected item from group_info


% --- Executes during object creation, after setting all properties.
function group_info_CreateFcn(hObject, eventdata, handles)
% hObject    handle to group_info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in anatomical_plane.
function anatomical_plane_Callback(hObject, eventdata, handles)
% hObject    handle to anatomical_plane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns anatomical_plane contents as cell array
%        contents{get(hObject,'Value')} returns selected item from anatomical_plane


structFile = get(handles.anatomical_file, 'string');
opts = lower(get(handles.anatomical_plane, 'string'));
val = get(handles.anatomical_plane, 'value');
updateSlices(handles, structFile, opts{val});


% --- Executes during object creation, after setting all properties.
function anatomical_plane_CreateFcn(hObject, eventdata, handles)
% hObject    handle to anatomical_plane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in apply_button.
function apply_button_Callback(hObject, eventdata, handles)
% hObject    handle to apply_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% Output format
opts = lower(get(handles.results_format, 'string'));
val = get(handles.results_format, 'value');
results.format = opts{val};

% Image values
opts = lower(get(handles.image_values, 'string'));
val = get(handles.image_values, 'value');
results.image_values = opts{val};

% Template name
results.anatomical_file = get(handles.anatomical_file, 'string');

% Slice plane
opts = lower(get(handles.anatomical_plane, 'string'));
val = get(handles.anatomical_plane, 'value');
results.anatomical_plane = opts{val};

% Slices in mm
slices_in_mm = str2num(get(handles.slices_in_mm, 'string'));
results.slices_in_mm = slices_in_mm;

% Convert to z-scores
opts = lower(get(handles.convert_to_zscores, 'string'));
val = get(handles.convert_to_zscores, 'value');
results.convert_to_zscores = opts{val};

% Threshold
threshold = str2num(get(handles.threshold, 'string'));
results.threshold = threshold;

chkG = get(handles.group_info, 'string');
chkVal = get(handles.group_info, 'value');
if strcmpi(chkG{chkVal}, 'yes')     
    groupsInfo = getGroupsInfo(handles.num_subjects);
    if (isempty(groupsInfo))
        error('Groups are not selected');
    end
    results.groupsInfo = groupsInfo;
end

setappdata(0, 'resultsAppData', results);

delete(gcbf);

function setInputs(handles)
% Set defaults

icatb_defaults;
global IMAGE_VALUES;
global CONVERT_Z;
global THRESHOLD_VALUE;
global ANATOMICAL_PLANE;


chk = strmatch('html', lower(get(handles.results_format, 'string')), 'exact');
if (~isempty(chk))
    set(handles.results_format, 'value', chk);
end

chk = strmatch(lower(IMAGE_VALUES), lower(get(handles.image_values, 'string')), 'exact');
if (~isempty(chk))
    set(handles.image_values, 'value', chk);
end

structFile = fullfile(fileparts(which('groupica.m')), 'icatb_templates', 'ch2bet.nii');
set(handles.anatomical_file, 'string', structFile);

convert_to_zscores = CONVERT_Z;
chk = strmatch(lower(convert_to_zscores), lower(get(handles.convert_to_zscores, 'string')), 'exact');
if (~isempty(chk))
    set(handles.convert_to_zscores, 'value', chk);
end

tval = THRESHOLD_VALUE;
if (~ischar(tval))
    tval = num2str(tval);
end

set(handles.threshold, 'string', tval);

chk = strmatch(ANATOMICAL_PLANE, lower(get(handles.anatomical_plane, 'string')), 'exact');
if (~isempty(chk))
    set(handles.anatomical_plane, 'value', chk);
end

updateSlices(handles, structFile, ANATOMICAL_PLANE);


function updateSlices(handles, structFile, slicePlane)
% Update slices in mm

imagVol = icatb_get_vol_nifti(structFile);
% get the slices in mm for the corresponding plane
[sliceParameters] = icatb_get_slice_def(imagVol, slicePlane);
% get the slices in mm
slices_in_mm = sliceParameters.slices;
clear sliceParameters;
% construct string
slices_in_mm = icatb_constructString(slices_in_mm);

set(handles.slices_in_mm, 'string', slices_in_mm);


function resultsInfo = getGroupsInfo(numOfSub)
%% Add groups
%

icatb_defaults;
global UI_FS;

resultsInfo.numOfSub = numOfSub;
resultsInfo.group = [];
figureTag = 'groupsinfo';
InputHandle = icatb_getGraphics('Enter groups info', 'normal', figureTag, 'on');
set(InputHandle, 'menubar', 'none');
set(InputHandle, 'userdata', resultsInfo);

%% Prompt
yPos = 0.92; yOffset = 0.05; xOffset  = 0.02;
promptHeight = 0.05;
promptWidth = 0.52;
promptPos = [0.5 - 0.5*promptWidth, yPos, promptWidth, promptHeight];
okWidth = 0.12; okHeight = 0.05;
textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Add Groups', 'tag', ...
    'prompt_groups', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);

%% Listbox
listboxWidth = 0.4;
listboxHeight = 0.4;
listboxPos(1) = 0.5 - 0.5*listboxWidth;
listboxPos(2) = promptPos(2) - yOffset - listboxHeight;
listboxPos(3) = listboxWidth;
listboxPos(4) = listboxHeight;
listCompH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', listboxPos, 'string', '', 'tag', ...
    'group', 'fontsize', UI_FS - 1, 'min', 0, 'max', 1, 'value', 1, 'callback', {@addGroups, InputHandle});

addButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) + 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
removeButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) - 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', addButtonPos, 'string', '+', 'tag', 'add_group_button', 'fontsize',...
    UI_FS - 1, 'callback', {@addGroups, InputHandle});
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', removeButtonPos, 'string', '-', 'tag', 'remove_group_button', 'fontsize',...
    UI_FS - 1, 'callback', {@removeGroups, InputHandle});


%% Ok button
okPos = [0.5 - 0.5*okWidth, listboxPos(2) - okHeight - yOffset, okWidth, okHeight];
okPos(2) = okPos(2) - 0.5*okPos(4);
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', okPos, 'string', 'Done', 'tag', 'run_button', 'fontsize',...
    UI_FS - 1, 'callback', {@runCallback, InputHandle});

waitfor(InputHandle);

resultsInfo = [];
if (isappdata(0, 'get_groups_data'))
    resultsInfo = getappdata(0, 'get_groups_data');
    rmappdata(0, 'get_groups_data');
end

function addGroups(hObject, event_data, figH)
%% Add groups

resultsInfo = get(figH, 'userdata');
listH = findobj(figH, 'tag', 'group');
val = get(listH, 'value');

icatb_defaults;
global UI_FS;

figureTag = 'add_groups';
compFigHandle = findobj('tag', figureTag);
if (~isempty(compFigHandle))
    delete(compFigHandle);
end

groupVal = [];
groupName = '';
if (listH == hObject)
    if (~strcmpi(get(figH, 'selectionType'), 'open'))
        return;
    end
    val = get(listH, 'value');
    try
        groupName = resultsInfo.group(val).name;
        groupVal = resultsInfo.group(val).val;
    catch
    end
end

subjectString = cellstr([repmat('Subject ', resultsInfo.numOfSub, 1), num2str((1:resultsInfo.numOfSub)')]);

%[groupName, groupVal] = icatb_select_groups_gui(subjectString, groupName, 'select_subjects', groupVal);
[groupName, groupVal] = icatb_select_groups_gui(subjectString, 'Group', 'select_subjects', groupName, groupVal);

try
    
    if (isempty(groupName))
        error('Group name is not selected');
    end
    
    if (isempty(groupVal))
        error('Subjects are not selected');
    end
    
    if (length(resultsInfo.group) > 0)
        chk = strmatch(lower(groupName), lower(cellstr(char(resultsInfo.group.name))), 'exact');
        if (~isempty(chk))
            ind = chk;
        end
    end
    
    if (~exist('ind', 'var'))
        ind = length(resultsInfo.group) + 1;
    end
    
    %% Set user selected information in figure
    resultsInfo.group(ind).name = groupName;
    resultsInfo.group(ind).val =  groupVal;
    set(figH, 'userdata', resultsInfo);
    groupListH = findobj(figH, 'tag', 'group');
    set(groupListH, 'string', cellstr(char(resultsInfo.group.name)));
    
catch
    icatb_errorDialog(lasterr, 'Group Selection');
end



function removeGroups(hObject, event_data, figH)
%% Remove groups
%

resultsInfo = get(figH, 'userdata');
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
    strs = cellstr(char(resultsInfo.group.name));
    resultsInfo.group(val) = [];
    strs(val) = [];
    set(listH, 'value', 1);
    set(listH, 'string', strs);
    set(figH, 'userdata', resultsInfo);
catch
end


function runCallback(hObject, event_data, handles)
%% Run callback
%

resultsInfo = get(handles, 'userdata');

if (isempty(resultsInfo.group))
    error('Group/Groups are not selected.');
end

groups_info = resultsInfo.group;

setappdata(0, 'get_groups_data', groups_info);

delete(handles);
