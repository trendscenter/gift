function varargout = sdfnc_post_process(varargin)
% SDFNC_POST_PROCESS MATLAB code for sdfnc_post_process.fig
%      SDFNC_POST_PROCESS, by itself, creates a new SDFNC_POST_PROCESS or raises the existing
%      singleton*.
%
%      H = SDFNC_POST_PROCESS returns the handle to a new SDFNC_POST_PROCESS or the handle to
%      the existing singleton*.
%
%      SDFNC_POST_PROCESS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SDFNC_POST_PROCESS.M with the given input arguments.
%
%      SDFNC_POST_PROCESS('Property','Value',...) creates a new SDFNC_POST_PROCESS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sdfnc_post_process_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sdfnc_post_process_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sdfnc_post_process

% Last Modified by GUIDE v2.5 19-Mar-2015 12:11:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @sdfnc_post_process_OpeningFcn, ...
    'gui_OutputFcn',  @sdfnc_post_process_OutputFcn, ...
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


% --- Executes just before sdfnc_post_process is made visible.
function sdfnc_post_process_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sdfnc_post_process (see VARARGIN)

% Choose default command line output for sdfnc_post_process
handles.output = hObject;

if (length(varargin) > 0)
    handles.compFiles = varargin{1}{1};
    structFile = fullfile (fileparts(which('gift.m')), 'icatb_templates', 'nsingle_subj_T1_2_2_5.nii');
    for nComp  = 1:size(handles.compFiles, 1)
        hD = icatb_orth_views(deblank(handles.compFiles(nComp, :)), 'image_values', 'positive', ...
            'convert_to_zscores', 'yes', 'get_interp_data', 1, 'set_to_max_voxel', 1, 'threshold', 1, 'structfile', structFile);
        %'structfile', which('ch2bet.nii')
        tmp = (squeeze(hD.data(:, :, hD.maxVoxelPos(end))));
        if (nComp == 1)
            compData = zeros([size(tmp), size(handles.compFiles, 1)]);
        end
        compData(:, :, nComp) = tmp;
    end
    handles.compData = compData;
    
    
    getInputs(varargin{1}{2}, handles);
    
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes sdfnc_post_process wait for user response (see UIRESUME)
% uiwait(handles.params_post_process_sdfnc);


% --- Outputs from this function are returned to the command line.
function varargout = sdfnc_post_process_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in comps.
function comps_Callback(hObject, eventdata, handles)
% hObject    handle to comps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

compFiles = handles.compFiles;
if (isfield(handles, 'comps'))
    comps = get(handles.comps, 'userdata');
else
    comps = (1:size(compFiles, 1));
end
comps = selectComps(handles.compData, comps);
set(handles.comps, 'userdata', comps);
set(handles.comps, 'foregroundcolor', [0, 1, 0]);
%handles.comps = comps;

% --- Executes on button press in apply_button.
function apply_button_Callback(hObject, eventdata, handles)
% hObject    handle to apply_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% get details
comps = get(handles.comps, 'userdata');

if (isempty(comps))
    error('Components are not selected');
end

%handlesData = guidata(handles);

params.comps = get(handles.comps, 'userdata');
params.num_clusters = str2num(get(handles.num_clusters, 'string'));
params.max_iter = str2num(get(handles.max_iter, 'string'));
opts = cellstr(get(handles.dmethod, 'string'));
val = get(handles.dmethod, 'value');
params.dmethod = lower(opts{val});
params.num_permutations = str2num(get(handles.num_permutations, 'string'));
params.transition_threshold = str2num(get(handles.transition_threshold, 'string'));

opts = cellstr(get(handles.threshold_type, 'string'));
val = get(handles.threshold_type, 'value');
params.threshold_type = lower(opts{val});
params.threshold = str2num(get(handles.threshold, 'string'));

setappdata(0, 'getParamsAppData', params);

delete(handles.output);

% --- Executes on selection change in threshold_type.
function threshold_type_Callback(hObject, eventdata, handles)
% hObject    handle to threshold_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns threshold_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from threshold_type


% --- Executes during object creation, after setting all properties.
function threshold_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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



function num_clusters_Callback(hObject, eventdata, handles)
% hObject    handle to num_clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_clusters as text
%        str2double(get(hObject,'String')) returns contents of num_clusters as a double


% --- Executes during object creation, after setting all properties.
function num_clusters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max_iter_Callback(hObject, eventdata, handles)
% hObject    handle to max_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_iter as text
%        str2double(get(hObject,'String')) returns contents of max_iter as a double


% --- Executes during object creation, after setting all properties.
function max_iter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in dmethod.
function dmethod_Callback(hObject, eventdata, handles)
% hObject    handle to dmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns dmethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dmethod



% --- Executes during object creation, after setting all properties.
function dmethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function comps = selectComps(compData, comps)

icatb_defaults;
global UI_FS;

compStr = num2str(((1:size(compData, 3))'));
figureTag = 'select_comps';

if (~exist('comps', 'var'))
    comps = [];
end

try
    delete(findobj(0, 'tag', figureTag));
catch
end

compFigHandle = icatb_getGraphics('Select Components', 'normal', figureTag);
set(compFigHandle, 'userdata', compData);
set(compFigHandle, 'menubar', 'none');

promptHeight = 0.05;
xOffset = 0.02; yOffset = promptHeight; listboxHeight = 0.72; yPos = 1;
okWidth = 0.12; okHeight = promptHeight;

%promptPos = [xOffset, yPos - 0.5*yOffset, promptWidth, promptHeight];

%% Right Listbox
listbox2Wdith = 0.1;
promptPos(2) = yPos - yOffset - 0.5*promptHeight;
promptPos(1) = (1 - xOffset - 2*listbox2Wdith);
promptPos(3) = 2*listbox2Wdith;
promptPos(4) = promptHeight;
textH = icatb_uicontrol('parent', compFigHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Components', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);
promptPos = get(textH, 'position');
listboxYOrigin = promptPos(2) - 0.5*yOffset - listboxHeight;
listboxXOrigin = promptPos(1) + 0.5*listbox2Wdith;
listboxPos = [listboxXOrigin, listboxYOrigin, listbox2Wdith, listboxHeight];
compListH = icatb_uicontrol('parent', compFigHandle, 'units', 'normalized', 'style', 'listbox', 'position', listboxPos, 'string', compStr, 'tag', 'components', 'fontsize', UI_FS - 1, ...
    'min', 0, 'max', 2, 'value', comps);

%% Show components
showWidth = 0.08; showHeight = 0.04;
showButtonPos = [listboxXOrigin + 0.5*listbox2Wdith - 0.5*showWidth, listboxYOrigin - yOffset - 0.5*showHeight, showWidth, showHeight];
showH = icatb_uicontrol('parent', compFigHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', showButtonPos, 'string', 'Show', 'fontsize', UI_FS - 1, 'callback', ...
    {@drawComp, compFigHandle});

%% Plot image on the left hand side
axesPos = [xOffset, listboxYOrigin, listboxHeight, listboxHeight];
axesH = axes('parent', compFigHandle, 'units', 'normalized', 'position', axesPos, 'tag', 'axes_display_comp');

promptPos = axesPos;

%% Add cancel and run buttons
okPos = [0.5 - 0.5*okWidth, promptPos(2) - yOffset - 0.5*okHeight, okWidth, okHeight];
okPos(2) = okPos(2) - 0.5*okPos(4);
icatb_uicontrol('parent', compFigHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', okPos, 'string', 'Done', 'tag', 'done_button', 'fontsize', UI_FS - 1, 'callback', ...
    {@getCompCallback, compFigHandle});

%% Draw components on the left hand side
drawComp(compListH, [], compFigHandle);

waitfor(compFigHandle);

appName = 'getCompAppData';
if (isappdata(0, appName))
    comps = getappdata(0, appName);
    rmappdata(0, appName);
end


function drawComp(compListH, event_data, handles)

icatb_defaults;
global UI_FONTNAME;
global FONT_COLOR;

sel_comp = get(findobj(handles, 'tag', 'components'), 'value');
compData = get(handles, 'userdata');
%axesH = findobj(handles, 'tag', 'axes_display_comp');
axesH = get(handles, 'currentaxes');

fontSizeText = 8;
cmap = icatb_getColormap(1, 2, 1);

if (~isempty(sel_comp))
    DIM = [size(compData, 1), size(compData, 2), length(sel_comp)];
    [im, numImagesX, numImagesY, textToPlot] = icatb_returnMontage(compData(:, :, sel_comp), [], DIM, [1, 1, 1], sel_comp);
    image(im, 'parent', axesH, 'CDataMapping', 'scaled');
    set(axesH, 'clim', [1, 200]); % set the axis positions to the specified
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

function getCompCallback(hObject, event_data, handles)

sel_comp = get(findobj(handles, 'tag', 'components'), 'value');
setappdata(0, 'getCompAppData', sel_comp);
delete(handles);


function getInputs(defaults, handles)

% compFiles = varargin{1};
% defaults = varargin{2};

set(handles.comps, 'userdata', defaults.comps);
if (~isempty(defaults.comps))
    set(handles.comps, 'foregroundcolor', [0, 1, 0]);
end
set(handles.num_clusters, 'string', num2str(defaults.num_clusters));
set(handles.max_iter, 'string', num2str(defaults.max_iter));
opts = cellstr(get(handles.dmethod, 'string'));
val = strmatch(lower(defaults.dmethod), lower(opts), 'exact');
if (isempty(val))
    val = 1;
end
set(handles.dmethod, 'value', val);
set(handles.num_permutations, 'string', num2str(defaults.num_permutations));
set(handles.transition_threshold, 'string', num2str(defaults.transition_threshold));

opts = cellstr(get(handles.threshold_type, 'string'));
val = strmatch(lower(defaults.threshold_type), lower(opts), 'exact');
if (isempty(val))
    val = 1;
end
set(handles.threshold_type, 'value', val);

set(handles.threshold, 'string', num2str(defaults.threshold));



function num_permutations_Callback(hObject, eventdata, handles)
% hObject    handle to num_permutations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_permutations as text
%        str2double(get(hObject,'String')) returns contents of num_permutations as a double


% --- Executes during object creation, after setting all properties.
function num_permutations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_permutations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function transition_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to transition_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of transition_threshold as text
%        str2double(get(hObject,'String')) returns contents of transition_threshold as a double


% --- Executes during object creation, after setting all properties.
function transition_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to transition_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
