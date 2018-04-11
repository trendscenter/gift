function varargout = post_process_dfnc(varargin)
% POST_PROCESS_DFNC MATLAB code for post_process_dfnc.fig
%      POST_PROCESS_DFNC, by itself, creates a new POST_PROCESS_DFNC or raises the existing
%      singleton*.
%
%      H = POST_PROCESS_DFNC returns the handle to a new POST_PROCESS_DFNC or the handle to
%      the existing singleton*.
%
%      POST_PROCESS_DFNC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POST_PROCESS_DFNC.M with the given input arguments.
%
%      POST_PROCESS_DFNC('Property','Value',...) creates a new POST_PROCESS_DFNC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before post_process_dfnc_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to post_process_dfnc_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help post_process_dfnc

% Last Modified by GUIDE v2.5 20-Feb-2017 09:26:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @post_process_dfnc_OpeningFcn, ...
    'gui_OutputFcn',  @post_process_dfnc_OutputFcn, ...
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


% --- Executes just before post_process_dfnc is made visible.
function post_process_dfnc_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to post_process_dfnc (see VARARGIN)

% Choose default command line output for post_process_dfnc
handles.output = hObject;

dfncInfo = varargin{1};

handles = setFieldVals(handles, dfncInfo);

movegui(hObject, 'center');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes post_process_dfnc wait for user response (see UIRESUME)
% uiwait(handles.post_process_dfnc);


% --- Outputs from this function are returned to the command line.
function varargout = post_process_dfnc_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;

waitfor(hObject);

results = [];
appName = 'pdFNCAppData';
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


results.regressCovFile = get(handles.regress_covariates, 'userdata');

opts = cellstr(get(handles.estimate_clusters, 'string'));
val = get(handles.estimate_clusters, 'value');
results.estimate_clusters = lower(deblank(opts{val}));
results.num_clusters = str2num(deblank(get(handles.num_clusters, 'string')));
results.num_comps_ica = str2num(deblank(get(handles.num_comps_ica, 'string')));

opts = cellstr(lower(get(handles.ica_algorithm, 'string')));
val = get(handles.ica_algorithm, 'value');
results.ica_algorithm = opts{val};

results.num_ica_runs = str2num(deblank(get(handles.num_ica_runs, 'string')));
results.kmeans_max_iter = handles.kmeans_max_iter; %str2num(deblank(get(handles.kmeans_max_iter, 'string')));

% opts = cellstr(lower(get(handles.dmethod, 'string')));
% val = get(handles.dmethod, 'value');
results.dmethod = handles.dmethod; %opts{val};
results.kmeans_num_replicates = handles.kmeans_num_replicates;
results.num_tests_est_clusters = handles.num_tests_est_clusters;

opts = cellstr(lower(get(handles.meta_method, 'string')));
val = get(handles.meta_method, 'value');
results.meta_method = opts{val};

setappdata(0, 'pdFNCAppData', results);

delete(gcbf);

drawnow;


% --- Executes on selection change in regress_covariates.
function regress_covariates_Callback(hObject, eventdata, handles)
% hObject    handle to regress_covariates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns regress_covariates contents as cell array
%        contents{get(hObject,'Value')} returns selected item from regress_covariates

val = get(hObject, 'value');

if (val == 2)
    regressCovFile = icatb_selectEntry('title', 'Select covariates file (continuous covariate or reduced model from mancova fnc results) ...', 'typeEntity', 'file', ...
        'typeSelection', 'single', 'filter', '*mancovan*results*fnc.mat;*txt');
    
    if (isempty(regressCovFile))
        set(hObject, 'value', 1);
    else
        set(hObject, 'userdata', regressCovFile);
    end
else
    set(hObject, 'userdata', '');
end


% --- Executes during object creation, after setting all properties.
function regress_covariates_CreateFcn(hObject, eventdata, handles)
% hObject    handle to regress_covariates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_comps_ica_Callback(hObject, eventdata, handles)
% hObject    handle to num_comps_ica (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_comps_ica as text
%        str2double(get(hObject,'String')) returns contents of num_comps_ica as a double


% --- Executes during object creation, after setting all properties.
function num_comps_ica_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_comps_ica (see GCBO)
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


% --- Executes on selection change in estimate_clusters.
function estimate_clusters_Callback(hObject, eventdata, handles)
% hObject    handle to estimate_clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns estimate_clusters contents as cell array
%        contents{get(hObject,'Value')} returns selected item from estimate_clusters

val = get(hObject, 'value');

set(handles.num_clusters, 'enable', 'on');
if (val == 2)
    set(handles.num_clusters, 'enable', 'inactive');
end


% --- Executes during object creation, after setting all properties.
function estimate_clusters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to estimate_clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in help_est_clusters.
function help_est_clusters_Callback(hObject, eventdata, handles)
% hObject    handle to help_est_clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgH = msgbox(['Number of clusters is estimated using gap and silhouette statistics. Estimated clusters are used in the standard dFNC analysis. ', ...
    'Please note that this approach involves running k-means several times to determine optimal number of clusters and can be slower for large data'], 'Estimate Clusters', 'modal');
waitfor(msgH);



% --- Executes on button press in help_num_clusters.
function help_num_clusters_Callback(hObject, eventdata, handles)
% hObject    handle to help_num_clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgH = msgbox('Number of k-means clusters used in the standard dFNC analysis.', 'No of clusters', 'modal');
waitfor(msgH);

% --- Executes on selection change in ica_algorithm.
function ica_algorithm_Callback(hObject, eventdata, handles)
% hObject    handle to ica_algorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ica_algorithm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ica_algorithm


% --- Executes during object creation, after setting all properties.
function ica_algorithm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ica_algorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_ica_runs_Callback(hObject, eventdata, handles)
% hObject    handle to num_ica_runs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_ica_runs as text
%        str2double(get(hObject,'String')) returns contents of num_ica_runs as a double


% --- Executes during object creation, after setting all properties.
function num_ica_runs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_ica_runs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in help_num_comps.
function help_num_comps_Callback(hObject, eventdata, handles)
% hObject    handle to help_num_comps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgH = msgbox(['Enter number of components/clusters to be used in the meta state analysis approach. K-means, PCA, Temporal ICA and spatial ICA are performed on windowed dFNCs. ', ...
    'It is recommended not to use more than 8 components/clusters.'], 'Number of components', 'modal');
waitfor(msgH);



% --- Executes on button press in help_ica_algorithms.
function help_ica_algorithms_Callback(hObject, eventdata, handles)
% hObject    handle to help_ica_algorithms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgH = msgbox('Specify ICA algorithm to use during temporal ICA or spatial ICA on windowed dFNC correlations.', 'ICA Algorithm', 'modal');
waitfor(msgH);



% --- Executes on button press in help_mst.
function help_mst_Callback(hObject, eventdata, handles)
% hObject    handle to help_mst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgH = msgbox('This option is used for temporal and spatial ICA only. ICA is run multiple times and the best stable run is determined using Minimum Spanning Tree approach', 'MST', 'modal');
waitfor(msgH);



function kmeans_max_iter_Callback(hObject, eventdata, handles)
% hObject    handle to kmeans_max_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kmeans_max_iter as text
%        str2double(get(hObject,'String')) returns contents of kmeans_max_iter as a double


% --- Executes during object creation, after setting all properties.
function kmeans_max_iter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kmeans_max_iter (see GCBO)
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


% --- Executes on button press in help_regress_cov.
function help_regress_cov_Callback(hObject, eventdata, handles)
% hObject    handle to help_regress_cov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgH = msgbox(['Variance associated with the covariates will be regressed out from the windowed dFNC correlations. ', ...
    'You could select reduced model from mancova (fnc_stats/*results*fnc*mat) or continuous covariates in a single ascii file.'], 'Regress covariates', 'modal');
waitfor(msgH);

function handles = setFieldVals(handles, dfncInfo)

if (~isempty(dfncInfo.regressCovFile))
    set(handles.regress_covariates, 'value', 2);
end

set(handles.num_clusters, 'enable', 'on');
if (strcmpi(dfncInfo.estimate_clusters, 'yes'))
    set(handles.estimate_clusters, 'value', 2);
    set(handles.num_clusters, 'enable', 'inactive');
end

set(handles.num_clusters, 'string', num2str(dfncInfo.num_clusters));

set(handles.num_comps_ica, 'string', num2str(dfncInfo.num_comps_ica));

chk = strmatch(lower(dfncInfo.ica_algorithm), lower(get(handles.ica_algorithm, 'string')), 'exact');
if (~isempty(chk))
    set(handles.ica_algorithm, 'value', chk);
end

set(handles.num_ica_runs, 'string', num2str(dfncInfo.num_ica_runs));

handles.kmeans_max_iter = dfncInfo.kmeans_max_iter;
handles.kmeans_methods = {'City', 'sqEuclidean', 'Hamming', 'Correlation', 'Cosine'};
handles.dmethod = lower(deblank(dfncInfo.dmethod));
kmeans_num_replicates = 5;
try
    kmeans_num_replicates = dfncInfo.kmeans_num_replicates;
catch
end

handles.kmeans_num_replicates = kmeans_num_replicates;


try
    num_tests_est_clusters  = dfncInfo.num_tests_est_clusters;
catch
    num_tests_est_clusters = 10;
end
handles.num_tests_est_clusters = num_tests_est_clusters;

% chk = strmatch(lower(dfncInfo.dmethod), lower(kmeans_methods), 'exact');
% if (~isempty(chk))
%     set(handles.dmethod, 'value', chk);
% end

% set(handles.kmeans_max_iter, 'string', num2str(dfncInfo.kmeans_max_iter));
%
% chk = strmatch(lower(dfncInfo.dmethod), lower(get(handles.dmethod, 'string')), 'exact');
% if (~isempty(chk))
%     set(handles.dmethod, 'value', chk);
% end

chk = strmatch(lower(dfncInfo.meta_method), lower(get(handles.meta_method, 'string')), 'exact');
if (~isempty(chk))
    set(handles.meta_method, 'value', chk);
end

meta_method_Callback(handles.meta_method, [], handles);

% --- Executes on selection change in meta_method.
function meta_method_Callback(hObject, eventdata, handles)
% hObject    handle to meta_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns meta_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from meta_method


val = get(hObject, 'value');
set(handles.num_ica_runs, 'enable', 'on');
set(handles.ica_algorithm, 'enable', 'on');
strs = cellstr(get(hObject, 'string'));
if (strcmpi(strs{val}, 'pca') || strcmpi(strs{val}, 'k-means'))
    set(handles.num_ica_runs, 'enable', 'inactive');
    set(handles.ica_algorithm, 'enable', 'inactive');
end


% --- Executes during object creation, after setting all properties.
function meta_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to meta_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in help_meta_method.
function help_meta_method_Callback(hObject, eventdata, handles)
% hObject    handle to help_meta_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


msgH = msgbox('Select the meta state method. If you select all, PCA, k-means, temporal ICA and spatial ICA will be performed on the dFNC correlations', 'Meta State Method', 'modal');
waitfor(msgH);


% --------------------------------------------------------------------
function k_means_options_Callback(hObject, eventdata, handles)
% hObject    handle to k_means_options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dlg_title = 'Select K-means options';

numParameters = 1;

inputText(numParameters).promptString = 'Enter maximum number of iterations';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = num2str(handles.kmeans_max_iter);
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'kmeans_max_iter';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;


chk = strmatch(lower(handles.dmethod), lower(handles.kmeans_methods), 'exact');
if (isempty(chk))
    chk = 1;
end


numParameters = numParameters + 1;

inputText(numParameters).promptString = 'Select distance method';
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerString = handles.kmeans_methods;
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'dmethod';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = chk;

numParameters = numParameters + 1;

inputText(numParameters).promptString = 'Number of times to repeat the clustering';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = num2str(handles.kmeans_num_replicates);
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'num_replicates';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;

numParameters = numParameters + 1;

inputText(numParameters).promptString = 'Number of reference data-sets for computing gap';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = num2str(handles.num_tests_est_clusters);
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'num_tests_est_clusters';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;



answers = icatb_inputDialog('inputtext', inputText, 'Title', dlg_title, 'handle_visibility', 'on');

if (~isempty(answers))
    
    handles.kmeans_max_iter = answers{1};
    handles.dmethod = answers{2};
    handles.kmeans_num_replicates = answers{3};
    try
        num_tests_est_clusters  = answers{4};
    catch
        num_tests_est_clusters = 10;
    end
    handles.num_tests_est_clusters = num_tests_est_clusters;
    guidata(handles.post_process_dfnc, handles);
    
end
