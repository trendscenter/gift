function varargout = fslica_id(varargin)
% fslica_id: Allows the user to select a group of ICA folders, and decide
% if each component is "good" or "bad" based on visual inspection of the
% component spatial map, timeseries, and power spectrum.  Note that this
% script "hijacks" some SPM file selection functions, so you can change the
% selectFoldersButton callback if you do not have SPM.  The output is:
%
% OUTPUT
% compnames: a vector of component labels, in format SUBID_IC#
% decision: a vector of [1,0] values, 1="bad" (flagged), 0="good"
%
% The purpose of this GUI is to allow researchers to make custom "gold
% standards" for evaluation of algorithms that identify bad/noisy networks
%
% RUN ON COMMAND LINE
% [compnames,labels] = fslica_id;
%
% STEPS TO USE
% 1) Run with command above
% 2) Select "Select IC Folders" to choose IC folders
% 3) For each component, select if it is "good" or "bad"
%  (in future do we want confidence rating?)
%  Making a selection will automatically advance you to the next
%  You can only click "Done" when all components have labels
%  You may go navigate around to change a label with the << and >> button
%  You cannot navigate away from a lab that does not have a rating, to
%  ensure that when you get to the end, all labels have ratings.
% 4) When all labels are done, click "done" to return to Matlab

% TO DO
% 1) Add a confidence interval to rating?
% 2) Give user ability to skip, have GUI return to ones that were skipped
%
% FSLICA_ID MATLAB code for fslica_id.fig
%      FSLICA_ID, by itself, creates a new FSLICA_ID or raises the existing
%      singleton*.
%
%      H = FSLICA_ID returns the handle to a new FSLICA_ID or the handle to
%      the existing singleton*.
%
%      FSLICA_ID('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FSLICA_ID.M with the given input arguments.
%
%      FSLICA_ID('Property','Value',...) creates a new FSLICA_ID or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fslica_id_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fslica_id_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fslica_id

% Last Modified by GUIDE v2.5 30-Jul-2012 21:06:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fslica_id_OpeningFcn, ...
                   'gui_OutputFcn',  @fslica_id_OutputFcn, ...
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


% --- Executes just before fslica_id is made visible.
function fslica_id_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fslica_id (see VARARGIN)

% Choose default command line output for fslica_id
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fslica_id wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = fslica_id_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
varargout{1} = handles.output;

% --- Executes on button press in selectFolderButton.
% This function should allow the user to select one or more folders
function selectFolderButton_Callback(hObject, eventdata, handles)
% hObject    handle to selectFolderButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Ask user to select ica folders
selected = 0;
while selected == 0
    [folders,selected] = spm_select([1 Inf],'dir','Select FSL ICA directories:','','','.ica');
end

% After we've selected, we might as well make it impossible to do again
set(handles.selectFolderButton,'Enable','off');

% Create list to hold ica folder paths and names
spatial_image_paths = {};  % will be master list of component image paths
timecourse_img_names = {}; % will be master list of all timeseries image files
frequency_img_names = {};  % will be master list of all frequency image files
spatial_image_labels = {};  % will be a master list of labels, subid_IC#
decision = []; % will be a master list of all decisions

% Put ica folders in a list, cycle through each one, extract imagenames,etc
for i=1:size(folders,1)
   [path, fname] = fileparts(folders(i,1:regexp(folders(i,:),'.ica')-1));
   icadir = strcat(path,'\',fname,'.ica\report\');

   % Get list of all zstat images
   zstatholder = dir(strcat(icadir,'IC_*_thresh.png'));
   compnames = cell(length(zstatholder),1);
   zstatimg = cell(length(zstatholder),1);
   tsimg = cell(length(zstatholder),1);
   fqimg = cell(length(zstatholder),1);
   decision_temp = zeros(length(zstatholder),1);
   for l = 1:length(zstatholder)
       % Save just the component name
       compnames{l} = [ fname '_IC' num2str(l) ];
       % Save the full png image path
       zstatimg{l} = [ icadir 'IC_' num2str(l) '_thresh.png' ];
       % Save the full timeseries image path
       tsimg{l} = [ icadir 't' num2str(l) '.png' ];
       % Save the full frequency image path
       fqimg{l} = [ icadir 'f' num2str(l) '.png' ];
   end
   
   % Put temporary lists into master variables
   spatial_image_labels = [ spatial_image_labels; compnames ];
   spatial_image_paths = [ spatial_image_paths; zstatimg ];
   timecourse_img_names = [ timecourse_img_names; tsimg ];
   frequency_img_names = [ frequency_img_names; fqimg ];
   decision = [ decision; decision_temp ];
end
  
% Place all variables into one structure, for later access
DATA.labels = spatial_image_labels;
DATA.paths = spatial_image_paths;
DATA.ts = timecourse_img_names;
DATA.freq = frequency_img_names;
DATA.numcomps = length(spatial_image_labels);

% Save DATA to user data of select button
setappdata(handles.selectFolderButton,'UserData',DATA);

% Save decision vector to User data of "Done" button
setappdata(handles.doneButton,'UserData',decision);

indices.curr = 1;     % Currently displayed image
indices.farthest = 1; % Farthest advanced

% Set userdata of ICimage to be the record of image indices
setappdata(handles.ICimage,'UserData',indices);

% Display the first set of images for the user
displayIC(handles,hObject,1)

% Enable the Good/Bad buttons
set(handles.goodButton,'Enable','on');
set(handles.BadButton,'Enable','on');

% Add total number of images
set(handles.last_Index,'String',num2str(length(spatial_image_labels)));


% --- Executes on button press in BadButton.
function BadButton_Callback(hObject, eventdata, handles)
% hObject    handle to BadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the current index that is being viewed
indices = getappdata(handles.ICimage,'UserData');
curr_index = indices.curr;
farthest_index = indices.farthest;
% Save a value of 1 (flagged as BAD) in this index in the output array
decision = getappdata(handles.doneButton,'UserData');
decision(curr_index) = 1;
setappdata(handles.doneButton,'UserData',decision);
% Update GUI to advance to next image, check button status, etc.
next_index = curr_index + 1;
advanceGUI(hObject,eventdata,handles,curr_index,next_index,farthest_index);

% --- Executes on button press in goodButton.
function goodButton_Callback(hObject, eventdata, handles)
% hObject    handle to goodButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get the current index that is being viewed
indices = getappdata(handles.ICimage,'UserData');
curr_index = indices.curr;
farthest_index = indices.farthest;
% Save a value of 0 (flagged as GOOD) in this index in the output array
decision = getappdata(handles.doneButton,'UserData');
decision(curr_index) = 0;
setappdata(handles.doneButton,'UserData',decision);
% Update GUI to advance to next image, check button status, etc.
next_index = curr_index + 1;
advanceGUI(hObject,eventdata,handles,curr_index,next_index,farthest_index);

% --- Executes on button press in nextButton.
function nextButton_Callback(hObject, eventdata, handles)
% hObject    handle to nextButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
indices = getappdata(handles.ICimage,'UserData');
curr_index = indices.curr;
farthest_index = indices.farthest;
% Update GUI to advance to next image, check button status, etc.
next_index = curr_index + 1;
advanceGUI(hObject,eventdata,handles,curr_index,next_index,farthest_index);

% --- Executes on button press in backButton.
function backButton_Callback(hObject, eventdata, handles)
% hObject    handle to backButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
indices = getappdata(handles.ICimage,'UserData');
curr_index = indices.curr;
farthest_index = indices.farthest;
% Update GUI to advance to next image, check button status, etc.
next_index = curr_index -1;
advanceGUI(hObject,eventdata,handles,curr_index,next_index,farthest_index);

% Display IC displays the desired index, and updates the ICimage UserData
% to reflect this index
function displayIC(handles,hObject,index_to_show)

% First get lists of files
DATA = getappdata(handles.selectFolderButton,'UserData');

% Load new ICimage, new Freq image, new TS image
ICimg = imread(DATA.paths{index_to_show});
[Freqimg,map,~] = imread(DATA.freq{index_to_show});
Freqimg = ind2rgb(Freqimg,map);
[TSimg,map,~] = imread(DATA.ts{index_to_show});
TSimg = ind2rgb(TSimg,map);

% Display selected ICimage in ICimage box
axes(handles.ICimage);
setappdata(handles.ICimage,'imageData',ICimg);
imshow(ICimg);

% Display selected TSimage in TSimage box
axes(handles.TSimage);
setappdata(handles.TSimage,'imageData',TSimg);
imshow(TSimg);

% Display selected Freqimage in FQimage box
axes(handles.FQimage);
setappdata(handles.FQimage,'imageData',Freqimg);
imshow(Freqimg);

% Update index that is currently being shown
indices = getappdata(handles.ICimage,'UserData');
indices.curr = index_to_show;
if index_to_show > indices.farthest
    indices.farthest = index_to_show;
end
setappdata(handles.ICimage,'UserData',indices);
set(handles.first_Index,'String',num2str(index_to_show));

% --- Executes on button press in doneButton.
function doneButton_Callback(hObject, eventdata, handles)
% hObject    handle to doneButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
decision = getappdata(handles.doneButton,'UserData');
DATA = getappdata(handles.selectFolderButton,'UserData');
output.decision = decision;
output.labels = DATA.labels;

% Get output name from user:
outname = inputdlg('Select output .mat name','Output File Name',1,{'fslICA_GS'});
save([ outname{1} '.mat' ],'output');
delete(handles.figure1);

fslica_id_OutputFcn(hObject, eventdata, handles)

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over selectFolderButton.
function selectFolderButton_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to selectFolderButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on mouse press over axes background.
function ICimage_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to ICimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Moves index forward / backward, checks to ensure OK button statuses
function advanceGUI(hObject,eventdata,handles,curr_index,next_index,farthest_index)

% Load DATA from selectFolderButton
DATA = getappdata(handles.selectFolderButton,'UserData');
% if the desired index is not outside the range, advance
if next_index <= DATA.numcomps
    displayIC(handles,hObject,next_index);
end

% If the current index is the last image, then allow the user to save
if curr_index == DATA.numcomps
    set(handles.doneButton,'Enable','on');
end

% If we are going to the first index (1), the back button should not be enabled
if next_index == 1
    set(handles.backButton,'Enable','off');
else
    set(handles.backButton,'Enable','on');
end

% If we are going to the last index, the forward button should not be enabled
if next_index == DATA.numcomps
    set(handles.nextButton,'Enable','off');
else
    set(handles.nextButton,'Enable','on');
end

% If the user is navigating around and curr_index < farthest_index, enable
% forward button
if curr_index < farthest_index
    set(handles.nextButton,'Enable','on');
else
    set(handles.nextButton,'Enable','off');
end
