function [popupHandle]=icatb_getUIPopUp(Handle,labels,pos,callback_function,visible,Tag)
% [popupHandle]=icatb_getUIPopUp(Handle,labels,pos,callback_function,visible,Tag);
%--------------------------------------------------------------------------
% CREATED BY: Eric Egolf 
% LAST MODIFIED: 12-29-03
% ABOUT: Sets up a popup gui object
% 
% #########################################################################
% 
% USAGE: (* means optional)
% 
% - [popupHandle]=icatb_getUIPopUp(Handle,labels,pos,callback_function,visible,Tag);
%   INFO: Using parameters setup a popup gui object then returns the handle
%   PARAMETERS:
%         Handle = handle of figure you want to put object in
%         labels = labels for each of the choices
%         pos = position on the figure, normalized
%         callback_function = string containing callback function
%             (see matlab callback functions for help)
%         *visible = 'on' or 'off', sets visibility of object
%              default = 'on'
%         *Tag = string to set Tag property
%               default = 'popupMenu'
% #############################################################
% 
% LICENSING:
% 
%------------------------------------------------------

icatb_defaults;
global BG_COLOR;
global BG2_COLOR;
global FG_COLOR;
global AXES_COLOR;
global UI_FONTNAME;
global UI_FONTUNITS;
global UI_FS;
global FONT_COLOR;

if(~exist('Tag','var'));
    Tag = 'popupMenu';
end

if(~exist('visible','var'))
    visible='on'
end

whichVersion = str2num(version('-release'));
if isempty(whichVersion)
    whichVersion = 14;
end

% check the version of Matlab
if whichVersion <= 13    

    popupHandle = uicontrol(Handle , 'BackgroundColor', BG2_COLOR, 'ForegroundColor', FONT_COLOR, 'string', labels, 'units', ...
        'normalized', 'position', pos, 'fontunits', UI_FONTUNITS, 'fontname', UI_FONTNAME, 'FontSize', UI_FS, 'style', 'popupmenu', ...
        'callback', callback_function, 'Visible', visible, 'tag', Tag);
else

    panelHandle = uipanel('Parent', Handle, 'units', 'normalized', 'position', pos, 'BackgroundColor', BG2_COLOR, 'bordertype', 'none');

    popupHandle = uicontrol('Parent', panelHandle , 'BackgroundColor', BG2_COLOR, 'ForegroundColor', FONT_COLOR, 'string', labels, 'units', ...
        'normalized', 'position', [0 0 1 1], 'fontunits', UI_FONTUNITS, 'fontname', UI_FONTNAME, 'FontSize', UI_FS, 'style', 'popupmenu', ...
        'callback', callback_function, 'Visible', visible, 'tag', Tag);

end

