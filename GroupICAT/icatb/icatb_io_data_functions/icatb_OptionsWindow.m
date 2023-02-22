function [output_parameters, valueButton] = icatb_OptionsWindow(inputParameters, defaults, figVisibility, varargin)
% Plot listbox on the left and frame on the right with the options that
% will be plotted on the frame. Set the fields to the parameters window
% using the tag information.
%
% OK: set the fields using the tag information of the objects that are
% plotted on the frame.
% Help: Show how to use the ICA options window using the text selected in
% the listbox.
% Defaults: Set defaults
% Cancel: Close the figure window and uses the parameters that were
% previously selected

if ~exist('figVisibility', 'var')
    figVisibility = 'off';
end

icatb_defaults;

% Screen Color Defaults
global BG_COLOR;
global BG2_COLOR;
global BUTTON_COLOR;
global BUTTON_FONT_COLOR;
global HELP_FONT_COLOR;
global FONT_COLOR;

% Fonts
global UI_FONTNAME;
global UI_FONTUNITS;
global UI_FS;

% title for the figure
titleFig = 'Options Window';


for n = 1:2:length(varargin)
    if (strcmpi(varargin{n}, 'help'))
        helpFile = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'title'))
        titleFig = varargin{n + 1};
    end
end


% plot parameters in a listbox
% display their options on the right frame with suitable information
% regarding uitype, strings, etc.

% Setup figure for GUI
[InputHandle] = icatb_getGraphics(titleFig, 'displaygui', titleFig, figVisibility);

% Set no menu bar for the figure
set(InputHandle, 'Menubar', 'None', 'userdata', inputParameters);

if (exist('helpFile', 'var'))
    % set the menu
    menu1H = uimenu('parent', InputHandle, 'label', 'GIFT-Help');
    menu2H = uimenu(menu1H, 'label', 'Options Window', 'callback', 'icatb_openHTMLHelpFile(helpFile);');
end

% title color
titleColor = [0 0.9 0.9];

% fonts
titleFont = 13;

axesH = axes('parent', InputHandle, 'position', [0 0 1 1], 'visible', 'off');

xPos = 0.5; yPos = 0.96;

% plot the title here
text(xPos, yPos, titleFig, 'color', titleColor, 'FontAngle', 'italic', ...
    'fontweight', 'bold', 'fontsize', titleFont, 'HorizontalAlignment', ...
    'center', 'FontName', UI_FONTNAME, 'parent', axesH);

xOffset = 0.02; yOffset = 0.032; % offsets

% Draw OK push button
okHeight = 0.05;
okWidth = 0.12;

% Draw help button
helpPos = [1 - xOffset - okWidth, 0.25*okHeight + yOffset okWidth okHeight];

% Draw default push button
defaultPos = [helpPos(1) - helpPos(3) - xOffset, 0.25*okHeight + yOffset okWidth okHeight];

% Draw cancel push button
cancelPos = [defaultPos(1) - defaultPos(3) - xOffset, 0.25*okHeight + yOffset okWidth okHeight];

% Draw Ok button
okPos = [cancelPos(1) - cancelPos(3) - xOffset, 0.25*okHeight + yOffset okWidth okHeight];

% Listbox positions
listboxXOrigin = xOffset; listboxYOrigin = okPos(2) + okPos(4) + yOffset;
listboxWidth = 0.32; listboxHeight = (yPos - 2*yOffset - listboxYOrigin);
listboxPos = [listboxXOrigin listboxYOrigin listboxWidth listboxHeight];

% determine frame Positions & Slider Positions
horzSliderHeight = 0.05; verSliderWidth = 0.04;
frameXOrigin = listboxPos(1) + listboxPos(3) + xOffset;
frameYOrigin = listboxPos(2) + horzSliderHeight + yOffset;
frameWidth = (1 - frameXOrigin - 2*xOffset - verSliderWidth);
frameHeight = listboxPos(2) + listboxPos(4) - frameYOrigin;

framePos = [frameXOrigin frameYOrigin frameWidth frameHeight];
horzSliderPos = [framePos(1) listboxPos(2) framePos(3) horzSliderHeight];
verSliderPos = [framePos(1) + framePos(3) + xOffset framePos(2) verSliderWidth frameHeight];

% Draw listbox
listH = uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', ...
    listboxPos, 'BackgroundColor', BG2_COLOR, 'ForegroundColor', FONT_COLOR, 'string', ...
    str2mat(inputParameters.listString), 'fontsize', UI_FS, 'fontname', UI_FONTNAME, 'fontunits', UI_FONTUNITS, 'value', 1, 'callback', ...
    {@optionslistboxCallback, InputHandle});

% Draw frame
frameH = uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'frame', 'position', framePos, ...
    'BackgroundColor', BG_COLOR, 'ForegroundColor', FONT_COLOR);

uiHeight = 0.05; % height of the uicontrol
uiWidth = (framePos(3) - 4*xOffset) / 2; % including prompt position and uicontrol position)

% Draw uicontrols on frame
intialXPosition = framePos(1) + xOffset;
intialYPosition = framePos(2) + framePos(4) - yOffset - uiHeight;

listboxValue = get(listH, 'value'); % get the current selection

count = 0;
checkPos = [intialXPosition intialYPosition uiWidth uiHeight];
promptPos = [intialXPosition intialYPosition uiWidth uiHeight];
controlPosition = checkPos;
% loop over number of parameters
for numParameters = 1:length(inputParameters)
    controlPosition(2) = promptPos(2);
    % check if the parameter value is equal to the listbox value
    if listboxValue == numParameters
        uiEnable = 'on';
        uiVisible = 'on';
    else
        uiEnable = 'on';
        uiVisible = 'off';
    end
    
    % loop over number of UICONTROLS
    for numUIControls = 1:length(inputParameters(numParameters).options)
        % xorigin and width
        controlPosition(1) = promptPos(1);
        controlPosition(3) = promptPos(3);
        % check if the control goes below the frame
        %         if controlPosition(2) < framePos(2) - yOffset
        %             uiEnable = 'off';
        %             uiVisible = 'off';
        %             plot_vertical_slider = 'on';
        %         end
        % check the prompt string
        if ~isempty(inputParameters(numParameters).options(numUIControls).promptString)
            count = count + 1;
            %% store the controls using tags
            tagPrefix(count).text = 'prefix';
            tag{count} = [tagPrefix(count).text, inputParameters(numParameters).options(numUIControls).tag];
            % Draw prompt string
            promptH = uicontrol('parent', InputHandle, 'units', 'normalized', 'style', ...
                'text', 'position', controlPosition, 'BackgroundColor', BG2_COLOR, ...
                'ForegroundColor', FONT_COLOR, 'string', ...
                inputParameters(numParameters).options(numUIControls).promptString, 'fontsize', UI_FS, 'fontname', UI_FONTNAME, 'fontunits', UI_FONTUNITS, ...
                'enable', uiEnable, 'visible', uiVisible);
            if ~iscell(inputParameters(numParameters).options(numUIControls).promptString)
                oldString = {inputParameters(numParameters).options(numUIControls).promptString};
            end
            % Adjust string in prompt string
            [newString, newPos] = textwrap(promptH, oldString);
            dd = controlPosition(4) - newPos(4);
            controlPosition(4) = newPos(4);
            controlPosition(2) = controlPosition(2) + dd;
            if (controlPosition(2) > intialYPosition)
                controlPosition(2) = intialYPosition;
            end
            set(promptH, 'string', newString, 'position', controlPosition, 'tag', tag{count});
            storePositions(count) = controlPosition(2);
            % update the control position here
            controlPosition(1) = controlPosition(1) + controlPosition(3) + xOffset;
        else
            % for checkboxes
            controlPosition(3) = 2*promptPos(3);
        end
        
        %controlPosition(4) = promptPos(4); % use the original uiHeight
        uiPos = inputParameters(numParameters).options(numUIControls).uiPos;
        answerPos = controlPosition;
        if ~isempty(uiPos)
            answerPos(3) = uiPos(1);
            answerPos(4) = uiPos(2);
        else
            answerPos(4) = promptPos(4);
        end
        
        tmpPos = controlPosition(2) + 0.5*controlPosition(4) - 0.5*answerPos(4);
        
        if (tmpPos > intialYPosition)
            tmpPos = intialYPosition;
        end
        
        answerPos(2) = tmpPos;
        
        count = count + 1;
        %% store the controls using tags
        tagPrefix(count).text = 'answer';
        tag{count} = [tagPrefix(count).text, inputParameters(numParameters).options(numUIControls).tag];
        storePositions(count) = controlPosition(2);
        
        % Draw UICONTROL
        callbackStr = '';
        if (isfield(inputParameters(numParameters).options(numUIControls), 'callback') && ~isempty(inputParameters(numParameters).options(numUIControls).callback))
            callbackStr = inputParameters(numParameters).options(numUIControls).callback;
        end                       
        
        uih = uicontrol('parent', InputHandle, 'units', 'normalized', 'style', ...
            inputParameters(numParameters).options(numUIControls).uiType, 'position', answerPos, ...
            'BackgroundColor', BG2_COLOR, 'ForegroundColor', FONT_COLOR, 'tag', ...
            tag{count}, 'string', inputParameters(numParameters).options(numUIControls).answerString, 'value', inputParameters(numParameters).options(numUIControls).value, ...
            'fontsize', UI_FS, 'fontname', UI_FONTNAME, 'fontunits', UI_FONTUNITS, 'callback', callbackStr, 'enable', 'on', 'visible', 'on');
        
        answerH = findobj(InputHandle, 'tag', tag{count});
        
        try
            if (~isempty(inputParameters(numParameters).options(numUIControls).userdata))
                set(answerH, 'userdata', inputParameters(numParameters).options(numUIControls).userdata);
            end
        catch
        end
        
        if (isfield(inputParameters(numParameters).options(numUIControls), 'contextmenu') && ~isempty(inputParameters(numParameters).options(numUIControls).contextmenu))
            set(answerH, 'tooltipString', 'Right click for more options ...');
            hcmenu = uicontextmenu;
            uimenu(hcmenu, 'Label', inputParameters(numParameters).options(numUIControls).contextmenu.label, 'Callback', ...
                inputParameters(numParameters).options(numUIControls).contextmenu.callback, 'tag', [tag{count}, 'ContextMenu']);
            set(answerH, 'uicontextmenu', hcmenu);
        end
        
        set(answerH, 'enable', uiEnable);
        set(answerH, 'visible', uiVisible);
              
        if (isfield(inputParameters(numParameters).options(numUIControls), 'enable') && ~isempty(inputParameters(numParameters).options(numUIControls).enable))
            set(answerH, 'enable', inputParameters(numParameters).options(numUIControls).enable);
        end
        
        % update the control position
        controlPosition(2) = controlPosition(2) - controlPosition(4) - yOffset;
    end
end

% Plot the pushbuttons here
okH = uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', okPos, ...
    'BackgroundColor', BUTTON_COLOR, 'ForegroundColor', BUTTON_FONT_COLOR, 'string', 'OK', 'fontsize', UI_FS, 'fontname', UI_FONTNAME, 'fontunits', UI_FONTUNITS, ...
    'userdata', inputParameters, 'callback', {@okCallback, InputHandle}, 'tooltipString', 'Uses the options and closes the figure ...');
cancelH = uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', cancelPos, ...
    'BackgroundColor', BUTTON_COLOR, 'ForegroundColor', BUTTON_FONT_COLOR, 'string', 'Cancel', 'fontsize', UI_FS, 'fontname', UI_FONTNAME, 'fontunits', UI_FONTUNITS, ...
    'callback', {@cancelCallback, InputHandle}, 'tooltipString', 'Deletes current Figure ...');
helpH = uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', helpPos, ...
    'BackgroundColor', BUTTON_COLOR, 'ForegroundColor', HELP_FONT_COLOR, 'string', '?', 'fontsize', UI_FS, 'fontname', UI_FONTNAME, 'fontunits', UI_FONTUNITS, 'fontweight', 'bold', ...
    'callback', {@helpCallback, InputHandle}, 'tooltipString', 'Help ...');
defaultH = uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', defaultPos, ...
    'BackgroundColor', BUTTON_COLOR, 'ForegroundColor', BUTTON_FONT_COLOR, 'string', 'Default', 'fontsize', UI_FS, 'fontname', UI_FONTNAME, 'fontunits', UI_FONTUNITS, ...
    'userdata', defaults, 'callback', {@defaultCallback, InputHandle}, 'tooltipstring', 'Show defaults ...');
% get the results
[results, inputParameters] = results_from_controls(inputParameters, InputHandle);
% check the figure visibility
if strcmpi(figVisibility, 'on')
    
    % set the figure visibility on
    try
        set(InputHandle, 'visible', 'on');
        waitfor(InputHandle);
    catch
        if ishandle(InputHandle)
            delete(InputHandle);
        end
    end
    
else
    delete(InputHandle);
end

% get the application data
if isappdata(0, 'optionsAppData')
    output_parameters = getappdata(0, 'optionsAppData');
    valueButton = 1;
    rmappdata(0, 'optionsAppData');
else
    valueButton = 0;
    output_parameters.results = results;
    output_parameters.inputParameters = inputParameters;
    clear inputParameters;
    clear results;
end

% listbox callback
function optionslistboxCallback(hObject, event_data, handles)
% enables and makes visible only the controls that are selected
listboxValue = get(hObject, 'value');
inputParameters = get(handles, 'userdata');

% loop over number of parameters
for ii = 1:length(inputParameters)
    % loop over number of options
    for jj = 1:length(inputParameters(ii).options)
        tag = inputParameters(ii).options(jj).tag;
        promptHandle = findobj(handles, 'tag', ['prefix', tag]);
        answerHandle = findobj(handles, 'tag', ['answer', tag]);
        % get prefix tag
        if listboxValue == ii
            if ~isempty(promptHandle)
                set(promptHandle, 'visible', 'on', 'enable', 'on'); % prompt handle
            end
            set(answerHandle, 'visible', 'on');
            %set(answerHandle, 'visible', 'on', 'enable', 'on'); % answer handle
        else
            if ~isempty(promptHandle)
                set(promptHandle, 'visible', 'off', 'enable', 'off'); % prompt handle
            end
            set(answerHandle, 'visible', 'off');
            %set(answerHandle, 'visible', 'off', 'enable', 'inactive'); % answer handle
        end
    end
end

% function callbacks
function okCallback(hObject, event_data, handles)

% Purpose: sets all the user seleceted parameters to the application data
inputParameters = get(hObject, 'userdata');

[results, inputParameters] = results_from_controls(inputParameters, handles); % get the results field
output_parameters.results = results;
% plot the user controls with the selected defaults
% loop over the number of input parameters
for ii = 1:length(inputParameters)
    % loop over number of controls for that list string
    for jj = 1:length(inputParameters(ii).options)
        objectTag = ['answer', inputParameters(ii).options(jj).tag]; % get the object tag
        if isempty(objectTag)
            error('Set the tags for the objects on the frame');
        end
        objectHandle = findobj(handles, 'tag', objectTag); % get the object handle on the frame
        inputParameters(ii).options(jj).answerString = get(objectHandle, 'string'); % get the corresponding string
        inputParameters(ii).options(jj).value = get(objectHandle, 'value'); % get the value
    end
end
output_parameters.inputParameters = inputParameters;
setappdata(0, 'optionsAppData', output_parameters); % application data
delete(handles);

drawnow;

% function callbacks
function cancelCallback(hObject, event_data, handles)

% delete the current figure
delete(handles);

% function callbacks
function defaultCallback(hObject, event_data, handles)

% Purpose: sets the default parameters
inputParameters = get(hObject, 'userdata'); % get the defaults
% plot the user controls with the selected defaults
% loop over the number of input parameters
for ii = 1:length(inputParameters)
    % loop over number of controls for that list string
    for jj = 1:length(inputParameters(ii).options)
        objectTag = ['answer', inputParameters(ii).options(jj).tag]; % get the object tag
        if isempty(objectTag)
            error('Set the tags for the objects on the frame');
        end
        objectHandle = findobj(handles, 'tag', objectTag); % get the object handle on the frame
        % Fix for R2009a (Enable controls and set the string)
        enableOn = get(objectHandle, 'enable');
        set(objectHandle, 'enable', 'on');
        set(objectHandle, 'string', inputParameters(ii).options(jj).answerString); % set the corresponding string
        set(objectHandle, 'value', inputParameters(ii).options(jj).value); % set the value
        set(objectHandle, 'enable', enableOn);
    end
end

function returnString = uiOptionsCallback(hObject, handles)
% callbacks of the objects plotted on the frame

getStyle = get(hObject, 'style'); % get the style

getStyle = lower(getStyle); % style of the uicontrol type
getString = get(hObject, 'string'); % get the object string
getValue = get(hObject, 'value'); % get the object value
getTag = get(hObject, 'tag'); % get the object tag value

% Check the uicontrols
if strcmpi(getStyle, 'edit')
    returnString = deblank(getString); % get the string
elseif strcmpi(getStyle, 'listbox') | findstr(lower(getStyle), 'popup')
    returnString = deblank(getString(getValue, :)); % get the selected strings
elseif strcmpi(getStyle, 'checkbox') | strcmpi(getStyle, 'radiobutton')
    returnString = num2str(getValue); % return value as a string
else
    returnString = deblank(getString);
end

% set the object userdata as the selected parameters
%set(hObject, 'userdata', returnString);

function [results, inputParameters] = results_from_controls(inputParameters, handles)
% get the results depending upon the numeric type
% set the field to the results structure

results = struct;
% loop over the number of input parameters
for ii = 1:length(inputParameters)
    % loop over number of controls for that list string
    for jj = 1:length(inputParameters(ii).options)
        objectTag = inputParameters(ii).options(jj).tag; % get the object tag
        if isempty(objectTag)
            error('Set the tags for the objects on the frame');
        end
        answerType = inputParameters(ii).options(jj).answerType; % get the answer type
        objectHandle = findobj(handles, 'tag', ['answer', objectTag]); % get the object handle on the frame
        %uiOptionsCallback(objectHandle, [], handles); % use the function callback
        %returnedString = get(objectHandle, 'userdata'); % get the selected parameters
        returnedString = uiOptionsCallback(objectHandle, handles);
        getValue = get(objectHandle, 'value');
        inputParameters(ii).options(jj).value = getValue;
        if strcmpi(inputParameters(ii).options(jj).uiType, 'edit')
            inputParameters(ii).options(jj).answerString = returnedString;
        end
        if strcmpi(answerType, 'numeric')
            returnedString = str2num(returnedString);
        end
        results = setfield(results, objectTag, returnedString); % set the field using the tag information
        ud = get(objectHandle, 'userdata');
        if (~isempty(ud))
            results = setfield(results, [objectTag, '_userdata'], ud);
            inputParameters(ii).options(jj).userdata = ud;
        end
    end
end


function helpCallback(hObject, event_data, handles)
% help button callback

D(1).string = 'Click on the left listbox and the options for that parameter will be displayed on the right in a frame. The functions for each button are listed below:';
D(size(D, 2) + 1).string = '';
D(size(D, 2) + 1).string = 'OK - The selected options are stored and the figure is closed. If the time course window is closed then the default data will be shown.';
D(size(D, 2) + 1).string = '';
D(size(D, 2) + 1).string = 'Cancel - The selected options will not be stored and the options figure is closed.';
D(size(D, 2) + 1).string = '';
D(size(D, 2) + 1).string = 'Default - Defaults will be shown. In order to use the defaults first default button should be pressed followed by OK.';
msgStr = str2mat(D.string);
dialogH = icatb_dialogBox('textBody', msgStr, 'title', 'Options Window', 'textType', 'large');
waitfor(dialogH);