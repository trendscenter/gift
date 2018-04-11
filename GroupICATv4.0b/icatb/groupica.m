function varargout = groupica(varargin)
% GROUPICA M-file for groupica.fig
%      GROUPICA, by itself, creates a new GROUPICA or raises the existing
%      singleton*.
%
%      H = GROUPICA returns the handle to a new GROUPICA or the handle to
%      the existing singleton*.
%
%      GROUPICA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GROUPICA.M with the given input arguments.
%
%      GROUPICA('Property','Value',...) creates a new GROUPICA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before groupica_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to groupica_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help groupica

% Last Modified by GUIDE v2.5 29-Jul-2011 15:26:52

icatb_delete_gui({'gift', 'eegift', 'sbm'});

try
    % add this statement to fix the button color
    feature('JavaFigures', 0);
catch
end

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @groupica_OpeningFcn, ...
    'gui_OutputFcn',  @groupica_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before groupica is made visible.
function groupica_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to groupica (see VARARGIN)

% Choose default command line output for groupica
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Check gift path
icatb_check_path;

% move the gui at the center of the screen
movegui(hObject, 'center');

if length(varargin) == 1
    if ischar(varargin{1})
        if strcmpi(varargin{1}, 'fmri')
            % fMRI Callback
            disp('Opening GIFT ...');
            fmri_button_Callback(hObject, eventdata, handles);
        elseif strcmpi(varargin{1}, 'eeg')
            % EEG Callback
            disp('Opening EEGIFT ...');
            eeg_button_Callback(hObject, eventdata, handles);
        elseif strcmpi(varargin{1}, 'smri')
            % sMRI Callback
            disp('Opening SBM ...');
            smri_button_Callback(hObject, eventdata, handles);
        elseif strcmpi(varargin{1}, 'defaults')
            % Defaults Callback
            disp('Opening Defaults ...');
            gica_defaults_Callback(hObject, eventdata, handles);
        elseif strcmpi(varargin{1}, 'exit') || strcmpi(varargin{1}, 'quit')
            disp('Quitting group ICA');
            % Exit callback
            exit_button_Callback(hObject, eventdata, handles);
        else
            % batch callback
            icatb_eval_script(varargin{1});
            exit_button_Callback(hObject, eventdata, handles);
        end
    end
end


%% Turn off finite warning
warning('off', 'MATLAB:FINITE:obsoleteFunction');

%% Turn off nifti class warning when loading SPM mat files
warning('off', 'MATLAB:unknownElementsNowStruc');

warning('off', 'MATLAB:dispatcher:ShadowedMEXExtension');

warning('off', 'MATLAB:pfileOlderThanMfile');


% UIWAIT makes groupica wait for user response (see UIRESUME)
% uiwait(handles.groupica);


% --- Outputs from this function are returned to the command line.
function varargout = groupica_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles)
    % Get default command line output from handles structure
    varargout{1} = handles.output;
else
    varargout{1} = [];
end


% --- Executes on button press in eeg_button.
function eeg_button_Callback(hObject, eventdata, handles)
% hObject    handle to eeg_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

eegift;

% --- Executes on button press in fmri_button.
function fmri_button_Callback(hObject, eventdata, handles)
% hObject    handle to fmri_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

gift;

% --- Executes on button press in smri_button.
function smri_button_Callback(hObject, eventdata, handles)
% hObject    handle to fmri_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sbm;

% --- Executes on button press in about_button.
function about_button_Callback(hObject, eventdata, handles)
% hObject    handle to about_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_titleDialog;

% --- Executes on button press in help_button.
function help_button_Callback(hObject, eventdata, handles)
% hObject    handle to help_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Help Callback
icatb_openHelp;

% --- Executes on button press in exit_button.
function exit_button_Callback(hObject, eventdata, handles)
% hObject    handle to exit_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_exit;

% --- Executes on button press in gica_defaults.
function gica_defaults_Callback(hObject, eventdata, handles)
% hObject    handle to gica_defaults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_defaults_gui;
