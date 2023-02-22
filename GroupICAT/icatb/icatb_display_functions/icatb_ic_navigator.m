function icatb_ic_navigator
% Navigator used for component explorer
% Borrowed code from  Bettyann Chodkowski at Kennedy Kreiger Institute

icatb_defaults;
global COMPONENT_NAMING;
global BG_COLOR;


figList = findall(0, 'Type', 'figure'); % figure list on command prompt

if ~isempty(figList)

    figNames = get(figList, 'name'); % figure names

    idx = find(icatb_regexpm(figNames, COMPONENT_NAMING));

    if isempty(idx)
        error('Cannot find component figures plotted using Component Explorer visualization method');
    end


    % Updated figure list
    figNames = str2mat(figNames{idx(:)});

    % Figure list
    figList = figList(idx);

    % Sort the figures based on their handle number
    [figList, newIndex] = sort(figList);

    figNames = figNames(newIndex(:), :);

    %%%%%% Draw IC Navigator %%%%%%%%
    figureTitle = 'IC Navigator';

    % delete a previous IC Navigator
    checkICGUI = findobj('tag', figureTitle);

    if ~isempty(checkICGUI)
        for ii = 1:length(checkICGUI)
            delete(checkICGUI(ii));
        end
    end

    [graphicsHandle] = icatb_getGraphics(figureTitle, 'normal', figureTitle, 'off');
    figPos = get(graphicsHandle, 'position');
    figPos = [figPos(1), figPos(2), 0.9*figPos(3), 0.85*figPos(4)];
    set(graphicsHandle, 'position', figPos);
    set(graphicsHandle, 'menubar', 'none');



    textPos = [0.05 0.9 0.9 0.05];
    textHandle = icatb_uicontrol('parent', graphicsHandle, 'style', 'text', 'position', textPos, 'units', ...
        'normalized', 'string', 'Component Figure Name: ', 'horizontalalignment', ...
        'center', 'tag', 'ic_navigator_text');

    extentPos = get(textHandle, 'extent');
    textPos(4) = extentPos(4) + 0.001;
    set(textHandle, 'position', textPos);

    yOffset = 0.03;

    xOrigin = textPos(1);
    yOrigin = 0.1;

    listWidth = 1 - 2*xOrigin;
    listHeight = textPos(2) - 0.5*textPos(4) - yOffset - yOrigin;

    % Listbox position
    listboxPos = [xOrigin, yOrigin, listWidth, listHeight];

    listData.figNames = figNames;
    listData.figHandles = figList;

    listHandle =  icatb_uicontrol('parent', graphicsHandle, 'style', 'listbox', 'units', ...
        'normalized', 'position', listboxPos, 'string', figNames, 'horizontalalignment', ...
        'center', 'min', 0, 'max', 1, 'callback', {@listboxCallback, graphicsHandle}, 'tag', ...
        'ic_navigator_listbox', 'userdata', listData);

    set(graphicsHandle, 'visible', 'on');

    % Make sure that when GUI is plotted axes is not plotted
    axis off;
    figure(graphicsHandle);
    %figure(findobj('tag', 'IC Navigator'));

    %%%% End for drawing IC Navigator %%%%%%

end



%%%%%% Listbox callback %%%%%%%%%
function listboxCallback(handleObj, event_data, handles)

listValue = get(handleObj, 'value');
listString = get(handleObj, 'string');
listData = get(handleObj, 'userdata');

figureName = deblank(listString(listValue, :));

% Figure names
figureNames = listString;
matchIndex = strmatch(lower(figureName), lower(figureNames), 'exact');

figHandles = listData.figHandles;

currentHandle = figHandles(matchIndex);

if ishandle(currentHandle)
    figure(currentHandle);
else
    error(['Cannot find figure with the name: ', figureName]);
end