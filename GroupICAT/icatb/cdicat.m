function varargout = cdicat(varargin)
% cdicat MATLAB code for cdicat.fig
%      cdicat, by itself, creates a new cdicat or raises the existing
%      singleton*.
%
%      H = cdicat returns the handle to a new cdicat or the handle to
%      the existing singleton*.
%
%      cdicat('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in cdicat.M with the given input arguments.
%
%      cdicat('Property','Value',...) creates a new cdicat or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cdicat_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cdicat_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cdicat

% Last Modified by GUIDE v2.5 19-May-2025 23:30:47

% Begin initialization code - DO NOT EDIT
icatb_delete_gui({'groupica', 'gift', 'eegift', 'sbm', 'fnc'});

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @cdicat_OpeningFcn, ...
    'gui_OutputFcn',  @cdicat_OutputFcn, ...
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


% --- Executes just before cdicat is made visible.
function cdicat_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cdicat (see VARARGIN)


group_ica_modality = 'CONN';

setappdata(0, 'group_ica_modality', group_ica_modality);

% Choose default command line output for cdicat
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

movegui(hObject, 'center');

% UIWAIT makes cdicat wait for user response (see UIRESUME)
% uiwait(handles.cdicat);


% --- Outputs from this function are returned to the command line.
function varargout = cdicat_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in setup_conn_ica.
function setup_conn_ica_Callback(hObject, eventdata, handles)
% hObject    handle to setup_conn_ica (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_setup_conn_ica;


% --- Executes on button press in display_fnc_ica.
function display_fnc_ica_Callback(hObject, eventdata, handles)
% hObject    handle to display_fnc_ica (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



icatb_defaults;
global PARAMETER_INFO_MAT_FILE;
global GICA_PARAM_FILE;

filterP = ['*', PARAMETER_INFO_MAT_FILE, '*.mat'];
param_file = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Parameter File', 'filter', filterP);

if (isempty(param_file))
    error('Parameter file is not selected for analysis');
end

load(param_file);

if (~exist('sesInfo', 'var'))
    error('Selected file is not a valid parameter file');
end

temporal_stats_betas = [];

results = results_summary_gui('num_subjects', sesInfo.numOfSub);

formatName = results.format;
results.formatName = formatName;

drawnow;


% Run second matlab instance (matlab is installed on system)
disp('Generating summary with new matlab process ...');
giftPath = fileparts(which('gift.m'));
resultsFile = fullfile(fileparts(param_file), [sesInfo.userInput.prefix, '_tmp_results_struct.mat']);
save(resultsFile, 'results');

if ispc
    matlabExe = fullfile(matlabroot, 'bin', 'matlab.exe');
else
    matlabExe = fullfile(matlabroot, 'bin', 'matlab');
end

commandStr = ['"', matlabExe, '" -nodisplay -nodesktop -nosplash -r "', ...
    'addpath(genpath(''', giftPath, ''')); ', ...
    'icatb_report_generator(''', param_file, ''',''', resultsFile, '''); ', ...
    'exit;"'];
[status, message] = system(commandStr);

if (status ~= 0)
    error(message);
end

disp('Done');

% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(get(0, 'children'));


% --- Executes on button press in run_fnc_ica.
function run_fnc_ica_Callback(hObject, eventdata, handles)
% hObject    handle to run_fnc_ica (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%icatb_run_conn_ica;
param_file_x = handles.param_file;
load(param_file_x)
icatb_runAnalysis(sesInfo, 'All');

if isfield(sesInfo.userInput, 'b_enl2lin_sidecar')
    if sesInfo.userInput.b_enl2lin_sidecar
        % Copy ENL to LIN project and run it
        clear sesInfo;
        [pathstr, file_name, extn] = fileparts(param_file_x);
        load([pathstr filesep file_name(1:length(file_name)-23) '-lin_ica_parameter_info.mat']);
        icatb_runAnalysis(sesInfo, 'All');
        
        % MAtch and compare components between ENL to LIN project
        com_file1 = [pathstr filesep file_name(1:length(file_name)-23) '-enl_mean_component_ica_s_all_.nii'];
        com_file2 = [pathstr filesep file_name(1:length(file_name)-23) '-lin_mean_component_ica_s_all_.nii'];

        oc_sort = icatb_cls_greedy_sort_components([]); %initiates class
        s_file_name_greed = oc_sort.m_greedy_simple(com_file1, com_file2); % engages greedy sort
        
        load(s_file_name_greed);
        n_coms_above_04 = 0;
        for i_corr=1:size(o.ari_ordered_pairs_table,1)
            o.ari_ordered_pairs_table(i_corr,3) = o.ard_corrs_table(o.ari_ordered_pairs_table(i_corr,1),o.ari_ordered_pairs_table(i_corr,2));
            if o.ari_ordered_pairs_table(i_corr,3) > 0.4
                n_coms_above_04 = n_coms_above_04 + 1;
            end
        end
        
        T = array2table(o.ari_ordered_pairs_table, ...
            'VariableNames', {[file_name(1:length(file_name)-23) '-enl_mean_component_ica_s_all_.nii'], [file_name(1:length(file_name)-23) '-lin_mean_component_ica_s_all_.nii'], 'Corr'});
        
        s_corr_file = [file_name(1:length(file_name)-23) 'ENL_vs_LIN_IC_correlations' datestr(datetime('now'),'yyyymmddHHMMSS') '.tsv'];
        
        writetable(T, s_corr_file, ...
            'FileType', 'text', ...
            'Delimiter', '\t');
        
        msgH = msgbox([num2str(n_coms_above_04) ' components (of ' num2str(size(o.ari_ordered_pairs_table,1)) ') matches between explicitly nonlinear components and linear components (higher correlation than 0.4). Component numbers and correlations were saved in ' s_corr_file], 'Matching Components', 'modal');
        waitfor(msgH);

    end
end

disp('Please copy the input_mancovan_2ttest.m for both non-linear or linear and then run them for comparison across groups.');
disp('Done running Connectivity Domain Analysis');


% --- Executes on button press in stats.
function stats_Callback(hObject, eventdata, handles)
% hObject    handle to stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%icatb_stats_loadings;


% --- Executes on button press in import_data.
function import_data_Callback(hObject, eventdata, handles)
% hObject    handle to import_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

param_file = icatb_import_data_conn_ica;

handles.param_file = param_file;
guidata(hObject, handles);

% --------------------------------------------------------------------
function display_tools_Callback(hObject, eventdata, handles)
% hObject    handle to display_tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function image_viewer_Callback(hObject, eventdata, handles)
% hObject    handle to image_viewer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_image_viewer;
