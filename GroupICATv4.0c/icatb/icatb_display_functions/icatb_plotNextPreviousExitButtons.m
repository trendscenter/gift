function icatb_plotNextPreviousExitButtons(graphicsH)

%     % number of figures
xOffset = 0.02;
yOffset = 0.005;
buttonWidth = 0.04;
buttonHeight = 0.025;
xOrigin = 1 - 3*xOffset - 3*buttonWidth;

numFigures = length(graphicsH);

for numFig = numFigures:-1:1
    figureData = get(graphicsH(numFig).H, 'userdata');
    figure(graphicsH(numFig).H);
    figureData.index = numFig;
    figureData.numFigures = numFigures;
    figureData.GraphicsHandle = graphicsH;

    % --Buttons
    % back button and position
    BackButtonPos = [xOrigin yOffset buttonWidth buttonHeight];
    BackButton = icatb_uicontrol('parent', graphicsH(numFig).H, 'style', 'pushbutton', 'units', ...
        'normalized', 'position', BackButtonPos, 'string', '<-', 'callback', ...
        {@previousFigureCallback, graphicsH(numFig).H});
    % exit button and position
    ExitButtonPos = [BackButtonPos(1) + BackButtonPos(3) + xOffset BackButtonPos(2) buttonWidth buttonHeight];
    ExitButton = icatb_uicontrol('parent', graphicsH(numFig).H, 'style', 'pushbutton', 'units', 'normalized',...
        'position', ExitButtonPos, 'string', 'Exit', 'callback', {@exitCallback, graphicsH(numFig).H});
    % next button and position
    NextButtonPos =[ExitButtonPos(1) + ExitButtonPos(3) + xOffset BackButtonPos(2) buttonWidth buttonHeight];
    NextButton = icatb_uicontrol('parent', graphicsH(numFig).H, 'style', 'pushbutton', 'units', 'normalized', ...
        'position', NextButtonPos, 'string', '->', 'callback', {@nextFigureCallback, graphicsH(numFig).H});

    set(graphicsH(numFig).H, 'userdata', figureData);

    set(graphicsH(numFig).H, 'menubar', 'figure');

    try
        set(graphicsH(numFig).H, 'toolbar', 'figure');
    catch
    end
    clear figureData;
end


% previous figure callback
function previousFigureCallback(handleObj, event_data, handles)

% get the associated figure data
figureData = get(handles, 'userdata');
if(figureData.index == 1)
    disp('Can''t Display Previous Figure, This is the First Figure');
else
    figure(figureData.GraphicsHandle(figureData.index-1).H);
end

% next figure callback
function nextFigureCallback(handleObj, event_data, handles)

% get the associated figure data
figureData = get(handles, 'userdata');
if(figureData.index == figureData.numFigures)
    disp('Can''t Display Next Figure, This is the Last Figure');
else
    figure(figureData.GraphicsHandle(figureData.index+1).H);
end

% exit callback
function exitCallback(handleObj, event_data, handles)

% get the associated figure data
figureData = get(handles, 'userdata');
numOfFigs = size(figureData.GraphicsHandle,2);
for nFigs = 1:numOfFigs
    if(ishandle(figureData.GraphicsHandle(nFigs).H))
        close(figureData.GraphicsHandle(nFigs).H);
    end
end