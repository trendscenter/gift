function varargout = component_labeller(varargin)
% COMPONENT_LABELLER M-file for component_labeller.fig
%      COMPONENT_LABELLER, by itself, creates a new COMPONENT_LABELLER or raises the existing
%      singleton*.
%
%      H = COMPONENT_LABELLER returns the handle to a new COMPONENT_LABELLER or the handle to
%      the existing singleton*.
%
%      COMPONENT_LABELLER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COMPONENT_LABELLER.M with the given input arguments.
%
%      COMPONENT_LABELLER('Property','Value',...) creates a new COMPONENT_LABELLER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before component_labeller_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to component_labeller_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help component_labeller

% Last Modified by GUIDE v2.5 26-Oct-2011 15:58:08

if (icatb_get_matlab_version < 14)
    error('Component labeller requires R14 and higher');
end


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @component_labeller_OpeningFcn, ...
    'gui_OutputFcn',  @component_labeller_OutputFcn, ...
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


% --- Executes just before component_labeller is made visible.
function component_labeller_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to component_labeller (see VARARGIN)

icatb_defaults;
global CONVERT_Z;
global THRESHOLD_VALUE;
global IMAGE_VALUES;

handles.outDir = pwd;
% Choose default command line output for component_labeller
handles.output = hObject;

set(handles.anatomical_file, 'string', fullfile(fileparts(which('gift.m')), 'icatb_templates', 'nsingle_subj_T1_2_2_5.nii'));
set(handles.image_values, 'value', strmatch(lower(IMAGE_VALUES), lower({'Positive and Negative', 'Positive', 'Absolute Value', 'Negative'}), 'exact'));
set(handles.z_scores, 'value', strmatch(lower(CONVERT_Z), lower(get(handles.z_scores, 'string')), 'exact'));
set(handles.threshold, 'string', num2str(THRESHOLD_VALUE));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes component_labeller wait for user response (see UIRESUME)
% uiwait(handles.comp_labeller_figure);


% --- Outputs from this function are returned to the command line.
function varargout = component_labeller_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in select_comp_files.
function select_comp_files_Callback(hObject, eventdata, handles)
% hObject    handle to select_comp_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

getStyle = get(hObject, 'style');

selectFiles = 1;
if (strcmpi(getStyle, 'popup') || strcmpi(getStyle, 'popupmenu'))
    getString = cellstr(get(hObject, 'string'));
    if (strcmpi(getString{get(hObject, 'value')}, 'yes'))
        selectFiles = 0;
    end
end

if (selectFiles)

    files = icatb_selectEntry('typeSelection', 'multiple', 'filter', '*.img; *.nii', 'fileType', 'image', 'title', 'Select image files ...');

    if (isempty(files))
        error('Image files are not selected');
    end

    handles.image_files = files;

    set(hObject, 'style', 'popup');
    set(hObject, 'string', char('Yes', 'No'));
    set(hObject, 'value', strmatch('yes', lower(get(hObject, 'string')), 'exact'));

    guidata(hObject, handles);

end

%data = guidata(hObject);

% --- Executes on button press in select_template_files.
function select_template_files_Callback(hObject, eventdata, handles)
% hObject    handle to select_template_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

getStyle = get(hObject, 'style');

selectFiles = 1;
if (strcmpi(getStyle, 'popup') || strcmpi(getStyle, 'popupmenu'))
    getString = cellstr(get(hObject, 'string'));
    if (strcmpi(getString{get(hObject, 'value')}, 'yes'))
        selectFiles = 0;
    end
end

if (selectFiles)

    startPath = fullfile(fileparts(which('gift.m')), 'icatb_templates');

    templateFiles = icatb_selectEntry('typeSelection', 'multiple', 'filter', '*.img; *.nii; *.zip', 'title', 'Select template files ...', 'startPath', startPath);

    if (isempty(templateFiles))
        error('Template files are not selected');
    end

    chkZip = find(icatb_good_cells(regexp(cellstr(templateFiles), '\.zip$')) ~= 0);

    if (~isempty(chkZip))
        if (size(templateFiles, 1) > 1)
            error('Select atmost one zip file when selecting templates');
        end
    end

    handles.templateFiles = templateFiles;

    set(hObject, 'style', 'popup');
    set(hObject, 'string', char('Yes', 'No'));
    set(hObject, 'value', strmatch('yes', lower(get(hObject, 'string')), 'exact'));

    guidata(hObject, handles);

end

fn = deblank(handles.templateFiles(1, :));
[p, fn] = fileparts(fn);

labelsFile = fullfile(p, [fn, '.txt']);

if (~exist(labelsFile, 'file'))
    warnH = warndlg(['Labels file with the name ', labelsFile, ' not found. Please select the file under template labels'], 'Labels File', 'modal');
    waitfor(warnH);
    set(handles.prompt_template_labels, 'visible', 'on');
    set(handles.select_template_labels, 'visible', 'on');
    set(handles.select_template_labels, 'enable', 'on');
    set(handles.help_template_labels, 'visible', 'on');
    set(handles.help_template_labels, 'enable', 'on');
else
    handles.labelsFile = labelsFile;
    guidata(hObject, handles);
end

cd(handles.outDir);

% --- Executes on button press in select_template_labels.
function select_template_labels_Callback(hObject, eventdata, handles)
% hObject    handle to select_template_labels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

getStyle = get(hObject, 'style');

selectFiles = 1;
if (strcmpi(getStyle, 'popup') || strcmpi(getStyle, 'popupmenu'))
    getString = cellstr(get(hObject, 'string'));
    if (strcmpi(getString{get(hObject, 'value')}, 'yes'))
        selectFiles = 0;
    end
end

if (selectFiles)

    startPath = fullfile(fileparts(which('gift.m')), 'icatb_templates');

    labelsFile = icatb_selectEntry('typeEntity', 'file', 'typeselection', 'single', ...
        'title', 'Select component labels text file', 'filter', '*.txt', 'startPath', startPath);

    if (isempty(labelsFile))
        error('Labels file is not selected');
    end

    handles.labelsFile = labelsFile;

    set(hObject, 'style', 'popup');
    set(hObject, 'string', char('Yes', 'No'));
    set(hObject, 'value', strmatch('yes', lower(get(hObject, 'string')), 'exact'));

    guidata(hObject, handles);

end

cd(handles.outDir);

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


% --- Executes on button press in anatomical_button.
function anatomical_button_Callback(hObject, eventdata, handles)
% hObject    handle to anatomical_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

startPath = fullfile(fileparts(which('gift.m')), 'icatb_templates');

anatFile = icatb_selectEntry('typeSelection', 'single', 'filter', '*.img; *.nii', 'fileType', 'image', 'title', 'Select template files ...', 'startPath', startPath);

if (isempty(anatFile))
    error('anatomical file is not selected');
end

set(handles.anatomical_file, 'string', anatFile);

cd(handles.outDir);


% --- Executes on selection change in z_scores.
function z_scores_Callback(hObject, eventdata, handles)
% hObject    handle to z_scores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns z_scores contents as cell array
%        contents{get(hObject,'Value')} returns selected item from z_scores


% --- Executes during object creation, after setting all properties.
function z_scores_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z_scores (see GCBO)
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


% --- Executes on selection change in image_values.
function image_values_Callback(hObject, eventdata, handles)
% hObject    handle to image_values (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns image_values contents as cell array
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


% --- Executes on button press in run_labeller.
function run_labeller_Callback(hObject, eventdata, handles)
% hObject    handle to run_labeller (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (~isfield(handles, 'image_files') || isempty(handles.image_files))
    error('Image files are not selected');
end

if (~isfield(handles, 'templateFiles') || isempty(handles.templateFiles))
    error('Template files are not selected');
end

if (~isfield(handles, 'labelsFile') || isempty(handles.labelsFile))
    error('Labels file is not selected');
end

anatFile = get(handles.anatomical_file, 'string');
if (isempty(anatFile))
    error('Anatomical file is not selected');
end

imageStr = cellstr(get(handles.image_values, 'string'));
zStr = cellstr(get(handles.z_scores, 'string'));
dispDefs = struct('structFile', anatFile, 'threshold', str2num(get(handles.threshold, 'string')), 'image_values', lower(imageStr{get(handles.image_values, 'value')}), ...
    'convert_to_zscores', lower(zStr{get(handles.z_scores, 'value')}));

set(handles.comp_labeller_figure, 'pointer', 'watch');
drawnow;
icatb_componentLabeller(handles.image_files, handles.templateFiles, handles.labelsFile, dispDefs);
set(handles.comp_labeller_figure, 'pointer', 'arrow');



% --- Executes on button press in help_template_files.
function help_template_files_Callback(hObject, eventdata, handles)
% hObject    handle to help_template_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = icatb_dialogBox('title', 'Template Files', 'textbody', ...
    ['Select template files containing regions of interest. Resting state networks templates are provided in icatb/icatb_templates/RSN.zip.', ...
    ' Labels associated with the files are specified in RSN.txt.'], 'texttype', 'large');
waitfor(h);


% --- Executes on button press in help_template_labels.
function help_template_labels_Callback(hObject, eventdata, handles)
% hObject    handle to help_template_labels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = icatb_dialogBox('title', 'Template Labels', 'textbody', ['Select labels file. Use comma delimited or new line delimited text when specifying labels.', ...
    ' Example is shown in file icatb/icatb_templates/RSN.txt'], 'texttype', 'large');
waitfor(h);


% --- Executes on button press in help_anatomical_file.
function help_anatomical_file_Callback(hObject, eventdata, handles)
% hObject    handle to help_anatomical_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = icatb_dialogBox('title', 'Anatomical file', 'textbody', 'Select anatomical file or enter the file name in the text box next to anatomical file label.', 'texttype', 'large');
waitfor(h);


% --- Executes on button press in help_image_files.
function help_image_files_Callback(hObject, eventdata, handles)
% hObject    handle to help_image_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = icatb_dialogBox('title', 'Component images', 'textbody', 'Select images of interest. After the images are selected, a dropdown box will appear with Yes and No. If you need to reselect the images, select no.' , ...
    'texttype', 'large');
waitfor(h);

