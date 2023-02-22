function [answer, pos] = icatb_promptUI(uiType,promptString,choiceString,returnType,H)

icatb_defaults;
global UI_FONTNAME;
global UI_FONTUNITS;
global UI_FS;
if(~exist('returnType','var'))
    returnType = 'string';
end

if(~exist('H','var'))
    H=gcf;
end

% Which version
whichVersion = str2num(version('-release'));
if isempty(whichVersion)
    whichVersion = 14;
end

objHandle = findobj(H, 'tag', 'prompt');

if ~isempty(objHandle)
    existingUIPos = get(objHandle, 'position');
else
    existingUIPos = [];
end

% Get the panel position for popup
if whichVersion > 13
    if ~isempty(existingUIPos)
        if icatb_findstr(uiType, 'popup')
            panelHandle = get(objHandle, 'parent');
            existingUIPos = get(panelHandle, 'position');
        end
    end
end

%existingUIPos = get(findobj('tag','prompt'),'position');
if(~isempty(existingUIPos))
    if(iscell(existingUIPos))
        latestUI = existingUIPos{1};
        y0 = latestUI(2)-latestUI(4)-.02;
    else
        y0 = existingUIPos(2)-existingUIPos(4)-.02;
    end
else
    % changing offset
    y0 = .92;
end


if(strcmp(uiType,'popup'))

    string = str2mat(promptString,choiceString);
    pos = [.05 y0 .9 .05];
    [popupHandle]=icatb_getUIPopUp(H,string,pos,'','on','prompt');
    set(popupHandle,'FontSize',UI_FS);
    waitfor(popupHandle,'value');
    set(popupHandle,'enable','off');
    value = get(popupHandle,'value');
    if(strcmp(returnType,'string'))
        answer = deblank(choiceString(value-1,:));
    else
        answer = value-1;
    end
end


if(strcmp(uiType,'edit'))

    pos1 = [.05 y0 .6 .05];
    pos2 = [.75 y0 .2 .05];
    pos = pos1;


    keypressCallback = 'set(findobj(''tag'',''prompt''), ''UserData'', ''OK'');';

    [labelHandle]=icatb_getUILabel(H,promptString,pos1,'on','prompt');
    [editHandle]=icatb_getUIEdit(H,choiceString,pos2, keypressCallback, 'on','prompt');
    set(labelHandle,'FontSize',UI_FS);
    set(editHandle,'FontSize',UI_FS);
    waitfor(editHandle, 'UserData', 'OK');
    set(editHandle,'enable','off');
    value = get(editHandle,'string');
    if(strcmp(returnType,'string'))
        answer = value;
    else
        answer = str2num(value);
    end
end

%% Select Multiple entries using a multiple toolbox
if(strcmp(uiType,'listbox'))
    string = str2mat(promptString,choiceString);
    pos = [.05 y0 .9 .05];
    [listHandle]=icatb_getUIlistbox(H,string,pos,'','on','prompt');
    set(listHandle, 'min', 1, 'max', size(string, 1)-1);
    set(listHandle,'FontSize',UI_FS);
    waitfor(listHandle,'value');
    set(listHandle,'enable','off');
    value = get(listHandle,'value');
    if(strcmp(returnType,'string'))
        answer = deblank(choiceString(value-1,:));
    else
        answer = value-1;
    end

end