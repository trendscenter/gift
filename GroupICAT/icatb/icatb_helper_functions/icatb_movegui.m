function icatb_movegui(figHandle, option)
% movegui 

% get the screen size
screenSize = get(0, 'Screensize');

figUnits = get(figHandle, 'units');

if ~strcmpi(figUnits, 'pixels')
    error('Figure units should be in pixels');
end

% pixels
offset = 50;
% get the figure size
figPos = get(figHandle, 'position');

if ~exist('option', 'var')
    option = 'center';
end

xOrigin = (screenSize(3)/2) - (figPos(3)/2);
yOrigin = (screenSize(4)/2) - (figPos(4)/2);

% do the following things
switch lower(option)
    case 'center'
        figPos = [xOrigin, yOrigin, figPos(3), figPos(4)];
    case 'east'
        xOrigin = screenSize(3) - figPos(3) - offset;
        figPos = [xOrigin, yOrigin, figPos(3), figPos(4)];
    case 'west'
        xOrigin = screenSize(1) + offset;
        figPos = [xOrigin, yOrigin, figPos(3), figPos(4)];
    case 'north'
        yOrigin = screenSize(4) - figPos(4) - offset;
        figPos = [xOrigin, yOrigin, figPos(3), figPos(4)];
    case 'south'
        yOrigin = screenSize(2) + offset;
        figPos = [xOrigin, yOrigin, figPos(3), figPos(4)];
    otherwise
end

% set the new figure position
set(figHandle, 'position', figPos);