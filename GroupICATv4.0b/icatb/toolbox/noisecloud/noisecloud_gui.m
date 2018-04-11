function varargout = noisecloud_gui(varargin)
%NOISECLOUD_GUI M-file for noisecloud_gui.fig
%      NOISECLOUD_GUI, by itself, creates a new NOISECLOUD_GUI or raises the existing
%      singleton*.
%
%      H = NOISECLOUD_GUI returns the handle to a new NOISECLOUD_GUI or the handle to
%      the existing singleton*.
%
%      NOISECLOUD_GUI('Property','Value',...) creates a new NOISECLOUD_GUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to noisecloud_gui_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      NOISECLOUD_GUI('CALLBACK') and NOISECLOUD_GUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in NOISECLOUD_GUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help noisecloud_gui

% Last Modified by GUIDE v2.5 03-Oct-2016 11:55:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @noisecloud_gui_OpeningFcn, ...
    'gui_OutputFcn',  @noisecloud_gui_OutputFcn, ...
    'gui_LayoutFcn',  [], ...
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


% --- Executes just before noisecloud_gui is made visible.
function noisecloud_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for noisecloud_gui
handles.output = hObject;

if (isempty(which('spm.m')))
    error('SPM does not exist on path');
end


chk = uigetdir(pwd, 'Select output directory ...');

if (isnumeric(chk) && ~chk)
    error('Output directory is not selected');
end


handles.outputDir = chk;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes noisecloud_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = noisecloud_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in run_class.
function run_class_Callback(hObject, eventdata, handles)
% hObject    handle to run_class (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


outputDir = pwd;

try
    outputDir = handles.outputDir;
catch
end

convert_to_z = 'yes';
try
    convert_to_z = handles.defaults.convert_to_z;
catch
end

permutations = 1;
try
    permutations = handles.defaults.permutations;
catch
end

cross_validations = 10;
try
    cross_validations = handles.defaults.cross_validations;
catch
end


training_opts.class_labels = handles.labels;

training_opts.TR = str2num(get(handles.train_tr, 'string'));
testing_opts.TR = str2num(get(handles.test_tr, 'string'));

try
    training_opts.sm = handles.training_opts.sm;
catch
    training_opts.sm = [];
end

try
    training_opts.tc = handles.training_opts.tc;
catch
    training_opts.tc = [];
end


try
    testing_opts.sm = handles.testing_opts.sm;
catch
    testing_opts.sm = [];
end

try
    testing_opts.tc = handles.testing_opts.tc;
catch
    testing_opts.tc = [];
end

try
    testing_opts.regress_cov = handles.testing_opts.regress_cov;
catch
    testing_opts.regress_cov = [];
end

try
    training_opts.regress_cov = handles.training_opts.regress_cov;
catch
    training_opts.regress_cov = [];
end

% Run classifier
[class_labels, fit_mdl, result_nc_classifier] = noisecloud_run(training_opts, testing_opts, 'convert_to_z', convert_to_z, 'outDir', outputDir, 'coregister', 1, ...
    'iterations', permutations, 'cross_validation', cross_validations);

disp('Saving results ...');
save(fullfile(outputDir, 'noise_cloud_results.mat'), 'class_labels', 'fit_mdl', 'result_nc_classifier', 'training_opts', 'testing_opts');



% --- Executes on selection change in convert_to_z.
function convert_to_z_Callback(hObject, eventdata, handles)
% hObject    handle to convert_to_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns convert_to_z contents as cell array
%        contents{get(hObject,'Value')} returns selected item from convert_to_z


% --- Executes during object creation, after setting all properties.
function convert_to_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to convert_to_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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



function cross_validations_Callback(hObject, eventdata, handles)
% hObject    handle to cross_validations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cross_validations as text
%        str2double(get(hObject,'String')) returns contents of cross_validations as a double


% --- Executes during object creation, after setting all properties.
function cross_validations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cross_validations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in help_convert_to_z.
function help_convert_to_z_Callback(hObject, eventdata, handles)
% hObject    handle to help_convert_to_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in help_num_permutations.
function help_num_permutations_Callback(hObject, eventdata, handles)
% hObject    handle to help_num_permutations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in help_cross_validations.
function help_cross_validations_Callback(hObject, eventdata, handles)
% hObject    handle to help_cross_validations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in browse_test_sm.
function browse_test_sm_Callback(hObject, eventdata, handles)
% hObject    handle to browse_test_sm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


sm_files = icatb_selectEntry('title', 'Select testing spatial maps set ...', 'typeSelection', 'multiple', 'filter', '*.img;*.nii');
drawnow;
handles.testing_opts.sm = sm_files;
guidata(hObject, handles);
if (~isempty(sm_files))
    set(handles.browse_test_sm, 'foregroundcolor', [0, 1, 0]);
end


% --- Executes on button press in browse_test_tc.
function browse_test_tc_Callback(hObject, eventdata, handles)
% hObject    handle to browse_test_tc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tc_files = icatb_selectEntry('title', 'Select testing timecourses set ...', 'typeSelection', 'multiple', 'filter', '*.img;*.nii');
drawnow;
handles.testing_opts.tc = tc_files;
guidata(hObject, handles);
if (~isempty(tc_files))
    set(handles.browse_test_tc, 'foregroundcolor', [0, 1, 0]);
end

function test_tr_Callback(hObject, eventdata, handles)
% hObject    handle to test_tr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of test_tr as text
%        str2double(get(hObject,'String')) returns contents of test_tr as a double


% --- Executes during object creation, after setting all properties.
function test_tr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to test_tr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in help_test_sm.
function help_test_sm_Callback(hObject, eventdata, handles)
% hObject    handle to help_test_sm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

H = helpdlg('Select subject component spatial maps (*sub*comp*nii)', 'Testing set (Spatial maps)');
waitfor(H);


% --- Executes on button press in help_test_tc.
function help_test_tc_Callback(hObject, eventdata, handles)
% hObject    handle to help_test_tc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


H = helpdlg('Select subject component timecourses (*sub*time*nii)', 'Testing set (Timecourses)');
waitfor(H);

% --- Executes on button press in help_test_tr.
function help_test_tr_Callback(hObject, eventdata, handles)
% hObject    handle to help_test_tr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

H = helpdlg('Experimental TR in seconds', 'Testing set');
waitfor(H);


% --- Executes on button press in browse_train_sm.
function browse_train_sm_Callback(hObject, eventdata, handles)
% hObject    handle to browse_train_sm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


sm_files = icatb_selectEntry('title', 'Select training spatial maps set ...', 'typeSelection', 'multiple', 'filter', '*.img;*.nii');
drawnow;
handles.training_opts.sm = sm_files;
guidata(hObject, handles);
if (~isempty(sm_files))
    set(handles.browse_train_sm, 'foregroundcolor', [0, 1, 0]);
end


% --- Executes on button press in browse_train_tc.
function browse_train_tc_Callback(hObject, eventdata, handles)
% hObject    handle to browse_train_tc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


tc_files = icatb_selectEntry('title', 'Select training timecourses set ...', 'typeSelection', 'multiple', 'filter', '*.img;*.nii');
drawnow;
handles.training_opts.tc = tc_files;
guidata(hObject, handles);
if (~isempty(tc_files))
    set(handles.browse_train_tc, 'foregroundcolor', [0, 1, 0]);
end


function train_tr_Callback(hObject, eventdata, handles)
% hObject    handle to train_tr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of train_tr as text
%        str2double(get(hObject,'String')) returns contents of train_tr as a double


% --- Executes during object creation, after setting all properties.
function train_tr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to train_tr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_train_labels.
function browse_train_labels_Callback(hObject, eventdata, handles)
% hObject    handle to browse_train_labels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


labelsFile = icatb_selectEntry('title', 'Select labels text file (2 levels: 1 - noise and 0 - network)', 'typeSelection', 'single', 'filter', '*txt');
drawnow;
handles.labels = labelsFile;
guidata(hObject, handles);
if (~isempty(labelsFile))
    set(handles.browse_train_labels, 'foregroundcolor', [0, 1, 0]);
end


% --- Executes on button press in help_train_sm.
function help_train_sm_Callback(hObject, eventdata, handles)
% hObject    handle to help_train_sm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

H = helpdlg('Select subject component spatial maps (*sub*comp*nii)', 'Training set (Spatial maps)');
waitfor(H);


% --- Executes on button press in help_train_tc.
function help_train_tc_Callback(hObject, eventdata, handles)
% hObject    handle to help_train_tc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

H = helpdlg('Select subject component timecourses (*sub*time*nii)', 'Training set (Timecourses)');
waitfor(H);


% --- Executes on button press in help_train_tr.
function help_train_tr_Callback(hObject, eventdata, handles)
% hObject    handle to help_train_tr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

H = helpdlg('Experimental TR in seconds', 'Training set');
waitfor(H);


% --- Executes on button press in help_train_labels.
function help_train_labels_Callback(hObject, eventdata, handles)
% hObject    handle to help_train_labels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

H = helpdlg('Select labels ascii file. Labels ascii file contains numbers 1 (Noise) and 0 (Network). Length of vector in the file must be subjects x components.', 'Labels file');
waitfor(H);


% --------------------------------------------------------------------
function defaults_menu_Callback(hObject, eventdata, handles)
% hObject    handle to defaults_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

convert_to_z = 'yes';
permutations = 1;
cross_validations = 10;

numParameters = 1;

inputText(numParameters).promptString = 'Do you want to convert components to Z-scores?';
inputText(numParameters).answerString = char('Yes', 'No');
inputText(numParameters).uiType = 'popupmenu';
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'convert_to_z';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;
inputText(numParameters).flag = 'scalar';
inputText(numParameters).help = struct('title', 'Z-scores', 'string', 'Components are scaled to z-scores. Mean is not removed while scaling as mean is used as one of the features.');


numParameters = numParameters + 1;
inputText(numParameters).promptString = 'Enter number of permutations';
inputText(numParameters).answerString = '1';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'permutations';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;
inputText(numParameters).flag = 'scalar';
inputText(numParameters).help = struct('title', 'Permutations', 'string', 'Number of times labels are randomly shuffled');

numParameters = numParameters + 1;
inputText(numParameters).promptString = 'Enter number of cross validations';
inputText(numParameters).answerString = '10';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'cross_validations';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;
inputText(numParameters).flag = 'scalar';
inputText(numParameters).help = struct('title', 'Cross validations', 'string', 'Number of cross validations');


answer = icatb_inputDialog('inputtext', inputText, 'title', 'Defaults ...', 'handle_visibility', 'on', 'windowstyle', 'modal');

if (~isempty(answer))
    convert_to_z = lower(answer{1});
    permutations = answer{2};
    cross_validations = answer{3};
end

handles.defaults.convert_to_z = convert_to_z;
handles.defaults.permutations = permutations;
handles.defaults.cross_validations = cross_validations;

guidata(hObject, handles);


% --- Executes on button press in test_covariates.
function test_covariates_Callback(hObject, eventdata, handles)
% hObject    handle to test_covariates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

regress_cov = icatb_selectEntry('title', 'Select covariates of testing set (optional)', 'typeSelection', 'multiple', 'filter', '*.txt');
drawnow;
handles.testing_opts.regress_cov = regress_cov;
guidata(hObject, handles);
if (~isempty(tc_files))
    set(handles.browse_test_tc, 'foregroundcolor', [0, 1, 0]);
end


% --- Executes on button press in help_test_covariates.
function help_test_covariates_Callback(hObject, eventdata, handles)
% hObject    handle to help_test_covariates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


H = helpdlg('Specify the covariates to be regressed from the timecourses. You could enter realignment parameters. Number of rows in each covariate file must match the number of timepoints. Number of files must match the number of subjects.', 'Regress Covariates (Testing set)');
waitfor(H);



% --- Executes on button press in train_covariates.
function train_covariates_Callback(hObject, eventdata, handles)
% hObject    handle to train_covariates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

regress_cov = icatb_selectEntry('title', 'Select covariates of training set (optional)', 'typeSelection', 'multiple', 'filter', '*.txt');
drawnow;
handles.training_opts.regress_cov = regress_cov;
guidata(hObject, handles);
if (~isempty(tc_files))
    set(handles.browse_test_tc, 'foregroundcolor', [0, 1, 0]);
end



% --- Executes on button press in help_train_covariates.
function help_train_covariates_Callback(hObject, eventdata, handles)
% hObject    handle to help_train_covariates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


H = helpdlg('Specify the covariates to be regressed from the timecourses. You could enter realignment parameters. Number of rows in each covariate file must match the number of timepoints. Number of files must match the number of subjects.', 'Regress Covariates (Training set)');
waitfor(H);
