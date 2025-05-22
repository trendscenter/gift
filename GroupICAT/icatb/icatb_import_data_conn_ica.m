function icatb_import_data_conn_ica(param_file)
%% Import data for running connectivity ICA
%

icatb_defaults;
global UI_FS;

if (~exist('param_file', 'var'))
    outputDir = icatb_selectEntry('typeEntity', 'directory', 'title', 'Select Analysis Output Directory');
    drawnow;
    sesInfo.userInput.outputDir = outputDir;
    sesInfo.userInput.pwd = outputDir;
    sesInfo.userInput.dataInfo.filesInfo = [];
    sesInfo.userInput.dataInfo.subsampling_depth = 1;
    sesInfo.userInput.dataInfo.maskFile = 'Default Mask';
else
    
    if (ischar(param_file))
        load(param_file);
        drawnow;
        if (~exist('sesInfo', 'var'))
            error('Selected file is not a valid ICA Parameter file');
        end
    else
        sesInfo = param_file;
    end
    outputDir = fileparts(param_file);
    if (isempty(outputDir))
        outputDir = pwd;
    end
    sesInfo.userInput.outputDir = outputDir;
    sesInfo.userInput.pwd = outputDir;
end

handles_data.sesInfo = sesInfo;

drawnow;
figTag = 'import_data_conn_ica';
graphicsHandle = icatb_getGraphics('Import data for Connectivity ICA', 'displaygui', figTag, 'on');
set(graphicsHandle, 'menubar', 'none');
set(graphicsHandle, 'CloseRequestFcn', @figCloseCallback);
set(graphicsHandle, 'userdata', handles_data);

% Offsets
xOffset = 0.05; yOffset = 0.05; yPos = 0.92;
buttonHeight = 0.052; promptHeight = 0.052;
promptWidth = 0.6;
editTextWidth = 0.2;
editTextHeight = 0.05;
promptTextPos = [xOffset, yPos - 0.5*yOffset - 0.5*promptHeight, promptWidth, promptHeight];


% Plot Text
promptH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', ...
    'position', promptTextPos, 'String', 'Enter output prefix to save analysis files', 'fontsize', UI_FS - 1, ...
    'horizontalalignment', 'center');

promptH = icatb_wrapStaticText(promptH);

editTextPos = get(promptH, 'position');
editTextPos(1) = editTextPos(1) + editTextPos(3) + xOffset;
editTextPos(3) = editTextWidth;
% Plot edit
editTextH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'edit', ...
    'position', editTextPos, 'String', '', 'fontsize', UI_FS - 1, ...
    'horizontalalignment', 'center',  'tag', 'output_prefix', 'callback', {@prefixCallback, graphicsHandle});


% Prompt
promptTextPos(2) = promptTextPos(2) - 1.2*yOffset - promptHeight;

% Plot Text
promptH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', ...
    'position', promptTextPos, 'String', 'Select fMRI data files for all subjects and sessions', 'fontsize', UI_FS - 1, ...
    'horizontalalignment', 'center');

promptH = icatb_wrapStaticText(promptH);

editTextPos = get(promptH, 'position');
editTextPos(1) = editTextPos(1) + editTextPos(3) + xOffset;
editTextPos(3) = editTextWidth;
buttonH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', editTextPos, 'String', 'Select', 'fontsize', UI_FS - 1, ...
    'horizontalalignment', 'center',  'tag', 'files', 'callback', {@dataCallback, graphicsHandle});


% Prompt
promptTextPos(2) = promptTextPos(2) - 1.2*yOffset - promptHeight;

maskOptions =  char('Default Mask', 'Average Mask', 'Select Mask', 'Default&ICV');

% Plot Text
promptH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', ...
    'position', promptTextPos, 'String', 'What Mask Do You Want To Use?', 'fontsize', UI_FS - 1, ...
    'horizontalalignment', 'center');

promptH = icatb_wrapStaticText(promptH);

popupTextPos = promptTextPos;
popupTextPos(1) = popupTextPos(1) + popupTextPos(3) + xOffset;
popupTextPos(3) = editTextWidth;

popupH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'popup', ...
    'position', popupTextPos, 'String', maskOptions, 'fontsize', UI_FS - 1, ...
    'horizontalalignment', 'left',  'tag', 'maskFile', 'value', 1, 'enable', 'inactive', 'callback', {@selectMaskCallback, graphicsHandle});

% Prompt
promptTextPos(2) = promptTextPos(2) - 1.2*yOffset - promptHeight;
promptH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', ...
    'position', promptTextPos, 'String', 'Enter sub-sampling depth to be applied on the mask', 'fontsize', UI_FS - 1, ...
    'horizontalalignment', 'center');

editTextPos = get(promptH, 'position');
editTextPos(1) = editTextPos(1) + editTextPos(3) + xOffset;
editTextPos(3) = editTextWidth;
editTextH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'edit', ...
    'position', editTextPos, 'String', '1', 'fontsize', UI_FS - 1, ...
    'horizontalalignment', 'center',  'tag', 'subsampling_depth');


% Prompt
promptTextPos(2) = promptTextPos(2) - 1.2*yOffset - promptHeight;
promptH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', ...
    'position', promptTextPos, 'String', 'Select connectivity type', 'fontsize', UI_FS - 1, ...
    'horizontalalignment', 'center');

popupTextPos = get(promptH, 'position');
popupTextPos(1) = popupTextPos(1) + popupTextPos(3) + xOffset;
popupTextPos(3) = editTextWidth;
popupH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'popup', ...
    'position', popupTextPos, 'String', 'ENLwFC', 'fontsize', UI_FS - 1, ...
    'horizontalalignment', 'left',  'tag', 'conn_type', 'value', 1, 'enable', 'inactive');


promptTextPos(2) = promptTextPos(2) - 1.5*yOffset - promptHeight;

% Plot done
buttonWidth = 0.12;
editTextPos(4) = buttonHeight;
editTextPos(1) = 0.5 - 0.5*buttonWidth;
editTextPos(2) = promptTextPos(2);
editTextPos(3) = buttonWidth;

% Plot done
icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', editTextPos, 'String', 'Done', 'fontsize', UI_FS - 1, 'horizontalalignment', 'center', 'callback', ...
    {@doneCallback, graphicsHandle});


function  prefixCallback (hObject, event_data, handles)
%% Prefix callback
%

icatb_defaults

handles_data = get(handles, 'userdata');

outputDir = handles_data.sesInfo.userInput.pwd;

prefix = get(hObject, 'string');
if (isempty(prefix))
    error('Enter output prefix');
end

param_file = [prefix, '_ica_parameter_info.mat'];
chkParamFile = fullfile(outputDir, param_file);
setC = 0;
if (exist(chkParamFile, 'file'))
    load(chkParamFile);
    handles_data.sesInfo = sesInfo;
    setC = 1;
end

drawnow;

handles_data.sesInfo.userInput.param_file = chkParamFile;
handles_data.sesInfo.userInput.prefix = prefix;
handles_data.sesInfo.userInput.pwd = outputDir;
handles_data.sesInfo.userInput.outputDir = outputDir;


set(handles, 'userdata', handles_data);

if (setC == 1)
    show_params(handles);
end



function dataCallback(hObject, event_data, handles)
%% Data callback

handles_data = get(handles, 'userdata');

try
    handles_data.sesInfo.userInput.prefix;
catch
    error('Output prefix is not entered');
end

set(handles, 'pointer', 'watch');


filesInfo = handles_data.sesInfo.userInput.dataInfo.filesInfo;

filesInfo = icatb_fileSelector(filesInfo, '*nii');

try
    filesInfo.filesList(1, :);
catch
    error ('Files are not selected');
end


handles_data.sesInfo.userInput.dataInfo.filesInfo = filesInfo;

set(handles, 'userdata', handles_data);

show_params(handles);

set(handles, 'pointer', 'arrow');

set(hObject, 'foregroundcolor', [0, 1, 0]);



function show_params(handles)

handles_data = get(handles, 'userdata');
sesInfo = handles_data.sesInfo;

%% Show selected params
filesH = findobj(handles, 'tag', 'files');
set(filesH, 'foregroundcolor', [0, 1, 0]);

% mask file
maskFile = 'Default Mask';
try
    maskFile = sesInfo.userInput.dataInfo.maskFile;
catch
end

maskH = findobj(handles, 'tag', 'maskFile');
maskOptions = cellstr(get(maskH, 'string'));

maskVal = strmatch(lower(maskFile), lower(maskOptions), 'exact');
if (isempty(maskVal))
    maskVal = 1;
end


try
    sesInfo.userInput.dataInfo.filesInfo.filesList(1, :);
catch
    error ('Files are not selected');
end

set(maskH, 'enable', 'on');
set(maskH, 'value', maskVal(1));

buttonH = findobj(handles, 'tag', 'files');
set(buttonH, 'foregroundcolor', [0, 1, 0]);

subsampling_depth = 1;
try
    subsampling_depth = sesInfo.userInput.dataInfo.subsampling_depth;
catch
end

subsampH = findobj(handles, 'tag', 'subsampling_depth');
set(subsampH, 'string', num2str(subsampling_depth));


maskH = findobj(handles, 'tag', 'maskFile');
maskOptions = cellstr(get(maskH, 'string'));


conn_type = 'ENLwFC';
try
    conn_type = sesInfo.userInput.dataInfo.conn_type;
catch

end

connH = findobj(handles, 'tag', 'conn_type');
connOptions = cellstr(get(connH, 'string'));
connVal = strmatch(lower(conn_type), lower(connOptions), 'exact');
if (isempty(connVal))
    connVal = 1;
end
set(connH, 'enable', 'on');
set(connH, 'value', connVal(1));



function selectMaskCallback(hObject, event_data, handles)
%% Mask callback
icatb_defaults;


handles_data = get(handles, 'userdata');

set(handles, 'pointer', 'watch');

filesList = [];

try
    filesList = handles_data.sesInfo.userInput.dataInfo.filesInfo.filesList;
catch
    
end

if (~isempty(filesList))
    [pathstr, fileN, extn] = fileparts(filesList(1, :));
    firstFile = icatb_listFiles_inDir(pathstr, [fileN, extn]);
    firstFile = icatb_fullFile('files', firstFile, 'directory', pathstr);
    firstFile = deblank(firstFile(1, :));
end

sesInfo = handles_data.sesInfo;
getString = get(hObject, 'string');
maskVal = get(hObject, 'value');
getMask = deblank(getString(maskVal, :));

maskFilter = '*.img;*.nii';
maskFigTitle = 'Select a mask file in analyze or nifti format';
fileType = 'image';

try
    
    maskFile = [];
    
    if (strcmpi(getMask, 'average mask'))
        maskFile = 'average mask';
    end
    
    % select mask in 3D Analyze data
    if (strcmpi(getMask, 'select mask'))
        [P] = icatb_selectEntry('filter', maskFilter, 'title', maskFigTitle, 'fileType', fileType, ...
            'fileNumbers', 1);
        if ~isempty(P)
            maskFile = P;
        else
            disp('No mask file selected. Setting to default mask ...');
            set(hObject, 'value', 1);
        end
    end
    
    
    getString = get(hObject, 'string');
    maskVal = get(hObject, 'value');
    getMask = deblank(getString(maskVal, :));
    
    %getMask = get(hObject, 'value');
    
    if (strcmpi(getMask, 'select mask'))
        if isempty(filesList)
            error('Data must be selected first before selecting mask');
        end
        [dT, extns, maskDim] = icatb_get_countTimePoints(icatb_parseExtn(maskFile));
        [dT, extns, dims] = icatb_get_countTimePoints(icatb_parseExtn(deblank(firstFile)));
        
        % If mask resolution do not match the mask will be resliced
        if length(find(maskDim == dims)) ~= length(maskDim)
            fprintf('Mask dimensions ([%s]) are not equal to image dimensions ([%s]). Resizing mask image/images to match functional image\n\n', ...
                num2str(maskDim), num2str(dims));
            
            % Handle gz files
            firstFileTmp = deblank(icatb_parseExtn(firstFile));
            if (strcmpi(firstFileTmp(end-2:end), '.gz'))
                gzfn = gunzip (firstFileTmp, tempdir);
                gzfn = char(gzfn);
                firstFile = icatb_rename_4d_file(gzfn);
                firstFile = deblank(firstFile(1, :));
            end
            
            sTemp = noisecloud_spm_coregister(firstFile, deblank(maskFile(1, :)), maskFile, sesInfo.userInput.pwd);
            [sFPath,sFName,sFExt] = fileparts(sTemp);
            % Rename tmp output to the mask name
            maskFile = [sesInfo.userInput.pwd filesep sesInfo.userInput.prefix 'Mask' sFExt];
            movefile(sTemp, maskFile, 'f')
        end
    end
    
    % create Mask
    sesInfo.userInput.dataInfo.maskFile = maskFile;
    handles_data.sesInfo = sesInfo;
    set(handles, 'userdata', handles_data);
    
    set(hObject, 'value', maskVal);
    
catch
    icatb_errorDialog(lasterr, 'Mask Error', 'modal');
end

set(handles, 'pointer', 'arrow');

function doneCallback(hObject, event_data, handles)
%% Done callback
%

handles_info = get(handles, 'userdata');

outputDir = handles_info.sesInfo.userInput.pwd;

if isempty(handles_info.sesInfo.userInput.dataInfo)
    error('Data is not selected for the analysis')
end

try
    handles_info.sesInfo.userInput.dataInfo.filesInfo.filesList(1,:);
catch
    error('Data is not selected for the analysis')
end

subsampH = findobj(handles, 'tag', 'subsampling_depth');
subsampling_depth = str2num(get(subsampH, 'string'));

handles_info.sesInfo.userInput.dataInfo.subsampling_depth = subsampling_depth;

conn_type = 'ENLwFC';
try
    conn_type = sesInfo.userInput.dataInfo.conn_type;
catch

end

connH = findobj(handles, 'tag', 'conn_type');
connOptions = cellstr(get(connH, 'string'));
connVal = strmatch(lower(conn_type), lower(connOptions), 'exact');
if (isempty(connVal))
    connVal = 1;
end

conn_type = deblank(connOptions{connVal});

handles_info.sesInfo.userInput.dataInfo.conn_type = conn_type;

sesInfo = handles_info.sesInfo;

disp('Saving parameters ...');
param_file = sesInfo.userInput.param_file;
%param_file = fullfile(outputDir, param_file); 
save(param_file, 'sesInfo');
disp('Done');
fprintf('\n');

try
    delete(handles);
catch
    
end

drawnow;

%% Generate connectivity matrices
icatb_gen_data_conn_ica(sesInfo);


function figCloseCallback(hObject, event_data, handles)
%% Fig close callback

delete(hObject);