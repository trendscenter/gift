function varargout = post_process_spatial_chronnectome(varargin)
%POST_PROCESS_SPATIAL_CHRONNECTOME MATLAB code file for post_process_spatial_chronnectome.fig
%      POST_PROCESS_SPATIAL_CHRONNECTOME, by itself, creates a new POST_PROCESS_SPATIAL_CHRONNECTOME or raises the existing
%      singleton*.
%
%      H = POST_PROCESS_SPATIAL_CHRONNECTOME returns the handle to a new POST_PROCESS_SPATIAL_CHRONNECTOME or the handle to
%      the existing singleton*.
%
%      POST_PROCESS_SPATIAL_CHRONNECTOME('Property','Value',...) creates a new POST_PROCESS_SPATIAL_CHRONNECTOME using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to post_process_spatial_chronnectome_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      POST_PROCESS_SPATIAL_CHRONNECTOME('CALLBACK') and POST_PROCESS_SPATIAL_CHRONNECTOME('CALLBACK',hObject,...) call the
%      local function named CALLBACK in POST_PROCESS_SPATIAL_CHRONNECTOME.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help post_process_spatial_chronnectome

% Last Modified by GUIDE v2.5 30-Oct-2020 12:37:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @post_process_spatial_chronnectome_OpeningFcn, ...
    'gui_OutputFcn',  @post_process_spatial_chronnectome_OutputFcn, ...
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


% --- Executes just before post_process_spatial_chronnectome is made visible.
function post_process_spatial_chronnectome_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for post_process_spatial_chronnectome
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes post_process_spatial_chronnectome wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = post_process_spatial_chronnectome_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;
waitfor(hObject);

results = [];
appName = 'pSchronnAppData';
if (isappdata(0, appName))
    results = getappdata(0, appName);
    rmappdata(0, appName);
end

varargout{1} = results;


% --- Executes on button press in done.
function done_Callback(hObject, eventdata, handles)
% hObject    handle to done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


results.num_intervals = str2num(deblank(get(handles.num_intervals, 'string')));
results.num_bins = str2num(deblank(get(handles.num_bins, 'string')));
results.tmap_transition_matrix_threshold = str2num(deblank(get(handles.tmap_transition_matrix, 'string')));
opts = cellstr(get(handles.est_clusters, 'string'));
val = get(handles.est_clusters, 'value');
results.est_clusters = lower(deblank(opts{val}));
results.num_clusters = str2num(deblank(get(handles.num_clusters, 'string')));

results.max_iter = str2num(deblank(get(handles.max_iter, 'string')));

opts = cellstr(get(handles.distance, 'string'));
val = get(handles.distance, 'value');
results.kmeans_distance_method = opts{val};
results.kmeans_num_replicates = str2num(deblank(get(handles.num_repetitions, 'string')));

opts = cellstr(get(handles.kmeans_init, 'string'));
val = get(handles.kmeans_init, 'value');
results.kmeans_init = opts{val};

use_tall_array = 'no';
if (isfield(handles, 'use_tall_array'))
    use_tall_array = handles.use_tall_array;
end

results.use_tall_array = use_tall_array;

setappdata(0, 'pSchronnAppData', results);

delete(gcbf);

drawnow;


% --- Executes on selection change in est_clusters.
function est_clusters_Callback(hObject, eventdata, handles)
% hObject    handle to est_clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns est_clusters contents as cell array
%        contents{get(hObject,'Value')} returns selected item from est_clusters


% --- Executes during object creation, after setting all properties.
function est_clusters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to est_clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in help_num_clusters.
function help_num_clusters_Callback(hObject, eventdata, handles)
% hObject    handle to help_num_clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgH = msgbox('Number of k-means clusters used in the analysis. If you selected estimate clusters, this option is not used. You can provide different clusters for different networks in a row vector like 5, 6.', 'No of clusters', 'modal');
waitfor(msgH);

% --- Executes on button press in help_est_clusters.
function help_est_clusters_Callback(hObject, eventdata, handles)
% hObject    handle to help_est_clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgH = msgbox(['Number of clusters is estimated using gap, silhouette, dunns, AIC/BIC and elbow methods. Estimated clusters are used in the analysis. ', ...
    'Please note that this approach involves running k-means several times to determine optimal number of clusters and can be slower for large data'], 'Estimate Clusters', 'modal');
waitfor(msgH);


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


% --- Executes on button press in help_max_iter.
function help_max_iter_Callback(hObject, eventdata, handles)
% hObject    handle to help_max_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


msgH = msgbox('Enter maximum number of iterations for the algorithm to converge.', 'Max Iterations', 'modal');
waitfor(msgH);

% --- Executes on selection change in distance.
function distance_Callback(hObject, eventdata, handles)
% hObject    handle to distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns distance contents as cell array
%        contents{get(hObject,'Value')} returns selected item from distance


% --- Executes during object creation, after setting all properties.
function distance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in help_distance.
function help_distance_Callback(hObject, eventdata, handles)
% hObject    handle to help_distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgH = msgbox('Select distance method. There are several options like City, sqEuclidean, Hamming, Correlation, Cosine distance measures.', 'Distance', 'modal');
waitfor(msgH);


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



function num_repetitions_Callback(hObject, eventdata, handles)
% hObject    handle to num_repetitions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_repetitions as text
%        str2double(get(hObject,'String')) returns contents of num_repetitions as a double


% --- Executes during object creation, after setting all properties.
function num_repetitions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_repetitions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in help_repetitions.
function help_repetitions_Callback(hObject, eventdata, handles)
% hObject    handle to help_repetitions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgH = msgbox('Enter number of times you want the algorithm to be run.', 'Num Repetitions', 'modal');
waitfor(msgH);

% --- Executes on button press in help_num_intervals.
function help_num_intervals_Callback(hObject, eventdata, handles)
% hObject    handle to help_num_intervals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgH = msgbox('Enter number of intervals.', 'Num intervals', 'modal');
waitfor(msgH);

% --- Executes on button press in help_num_bins.
function help_num_bins_Callback(hObject, eventdata, handles)
% hObject    handle to help_num_bins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgH = msgbox('Enter number of histogram bins.', 'Number of bins', 'modal');
waitfor(msgH);


function num_bins_Callback(hObject, eventdata, handles)
% hObject    handle to num_bins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_bins as text
%        str2double(get(hObject,'String')) returns contents of num_bins as a double


% --- Executes during object creation, after setting all properties.
function num_bins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_bins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_intervals_Callback(hObject, eventdata, handles)
% hObject    handle to num_intervals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_intervals as text
%        str2double(get(hObject,'String')) returns contents of num_intervals as a double


% --- Executes during object creation, after setting all properties.
function num_intervals_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_intervals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in help_tmap_transition.
function help_tmap_transition_Callback(hObject, eventdata, handles)
% hObject    handle to help_tmap_transition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgH = msgbox('T-map threshold. Significant voxels in atleast one of the states is used to compute spatial transition matrix.', 'T-map threshold (spatial transition matrix)', 'modal');
waitfor(msgH);


function tmap_transition_matrix_Callback(hObject, eventdata, handles)
% hObject    handle to tmap_transition_matrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tmap_transition_matrix as text
%        str2double(get(hObject,'String')) returns contents of tmap_transition_matrix as a double


% --- Executes during object creation, after setting all properties.
function tmap_transition_matrix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tmap_transition_matrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in kmeans_init.
function kmeans_init_Callback(hObject, eventdata, handles)
% hObject    handle to kmeans_init (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns kmeans_init contents as cell array
%        contents{get(hObject,'Value')} returns selected item from kmeans_init


% --- Executes during object creation, after setting all properties.
function kmeans_init_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kmeans_init (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in help_kmeans_init.
function help_kmeans_init_Callback(hObject, eventdata, handles)
% hObject    handle to help_kmeans_init (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


msgH = msgbox('Option is provided to initialize k-means using subject exemplars which are obtained using the global extrema points for each subject.', 'Distance', 'modal');
waitfor(msgH);


% --------------------------------------------------------------------
function options_menu_Callback(hObject, eventdata, handles)
% hObject    handle to options_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

use_tall_array = 'no';

if (~isempty(which('tall.m')))
    
    try
        use_tall_array = 'yes';
    catch
    end
    
    use_tall_array = questdlg('Do you want to use tall array option for large data?', ...
        'Tall array?', ...
        'No', 'Yes', 'No');
    
    if (isempty(use_tall_array))
        use_tall_array = 'no';
    end
    
end

handles.use_tall_array = use_tall_array;

guidata(hObject, handles);