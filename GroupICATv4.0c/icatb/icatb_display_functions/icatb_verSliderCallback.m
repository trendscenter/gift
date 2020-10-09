function icatb_verSliderCallback(handleObj, handles)
% vertical slider callback - Sets some of the controls invisible when
% scrolled up or down.
% 
% Input:
% 1. handleObj - Handle for object. Set user data for the object as a structure with the
% following fields 'maxHeight', 'minHeight', 'initialYPositions', 'tag'.
% 2. handles - Handle for the figure
%
% Output: none


% Get the user data
sliderData = get(handleObj, 'userdata');

% get the data information
initialYPos = sliderData.initialYPositions;
maxHeight = sliderData.maxHeight;
minHeight = sliderData.minHeight;

% Tags for the controls
tagFig = sliderData.tag;

for numSub = 1:length(tagFig)
    controlTag = findobj(handles, 'tag', tagFig{numSub});
    figSessPos = get(controlTag, 'position');
    % Get the slider value
    scrollValue = get(handleObj, 'value');
    % Determine the y position
    figSessPos(2) = initialYPos(numSub) - scrollValue;
    set(controlTag, 'position', figSessPos);
    
    if figSessPos(2) + figSessPos(4) > maxHeight
        set(controlTag, 'visible', 'off');
    else
        set(controlTag, 'visible', 'on');
    end
end