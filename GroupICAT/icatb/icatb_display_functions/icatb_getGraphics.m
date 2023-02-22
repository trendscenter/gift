function [GraphicsHandle]=icatb_getGraphics(TitleFig, GraphicsType, Tag, visible)
% [GraphicsHandle]=icatb_getGraphics(Title,Rect,Tag)
%--------------------------------------------------------------------------
% CREATED BY: Eric Egolf 
% LAST MODIFIED: 12-29-03
% ABOUT: Sets up a figure to be used as a graphics display figure. Returns
%        the Handle to the figure
% 
% #########################################################################
% 
% USAGE: (* means optional)
% 
% - [GraphicsHandle]=icatb_getGraphics(Title,Rect,Tag)
%   INFO: Using parameters setsup a figure
%   PARAMETERS:
%       *Title = Title for figure
%           default = 'Graphics'
%       *Rect = [x y height width]
%           default = [500 100 600 600];       
%       *Tag = Tag for figure
%           default = 'Graphics'
% #############################################################
% 
% LICENSING:
% 
%------------------------------------------------------

%get defaults
icatb_defaults;
global BG_COLOR;
global FG_COLOR;
global AXES_COLOR;
global FONT_COLOR;
global ASPECT_RATIO_FOR_SQUARE_FIGURE;
global S0;
global WS;
global MIN_SCREEN_DIM_IN_PIXELS;


if(~exist('TitleFig','var'))
    TitleFig = 'Graphics';
end
if(~exist('GraphicsType','var'))
    GraphicsType = 'normal';
end
if(~exist('Tag','var'))
    Tag = 'Graphics';
end

if ~exist('visible', 'var')
    visible = 'on';
end



%---------Get Coordinates For Figure
S=MIN_SCREEN_DIM_IN_PIXELS;
R=ASPECT_RATIO_FOR_SQUARE_FIGURE;
paperOrient = 'portrait';
switch lower(GraphicsType)
   
    case 'normal'
        extendLeft = round(S*R(1)*.55);
        extendUp = round(S*R(2)*.55);
        x0= (S0(3)/2)-(extendLeft/2);
        y0 = (S0(4)/2)-(extendUp/2);
        Rect = [x0 y0 extendLeft extendUp];
    case 'graphics'
        extendLeft = round(S*R(1)*.85);
        extendUp = round(S*R(2)*.85);
        x0= (S0(3)/2)-(extendLeft/2);
        y0 = (S0(4)/2)-(extendUp/2);
        Rect = [x0 y0 extendLeft extendUp];
    case 'displaygui'
        extendLeft = round(S*R(1)*.6);
        extendUp = round(S*R(2)*.6);
        x0= (S0(3)/2)-(extendLeft/2);
        y0 = (S0(4)/2)-(extendUp/2);
        Rect = [x0 y0 extendLeft extendUp];
    case 'statusbar'
        extendLeft = round(S*R(1)*.2);
        extendUp = round(S*R(2)*.3);
        x0= S0(3)*.1;
        y0 = S0(4)*.1;
        Rect = [x0 y0 extendLeft extendUp];
    case 'timecourse'
        extendLeft = round(S0(3)*R(1)*.9);
        extendUp = round(S*R(2)*.3);
        x0= S0(3)*.05;
        y0 = S0(4)*.5;
        Rect = [x0 y0 extendLeft extendUp];  
        paperOrient = 'landscape';
    otherwise        
        extendLeft = round(S*R(1)*.5);
        extendUp = round(S*R(2)*.5);
        x0= (S0(3)/2)-(extendLeft/2);
        y0 = (S0(4)/2)-(extendUp/2);
        Rect = [x0 y0 extendLeft extendUp];        
end


% %-----------Set Graphics Figure--------------------------     
GraphicsHandle = figure('Tag', Tag, 'units', 'pixels', 'Position', Rect, 'MenuBar', 'figure', ...
            'color', BG_COLOR, 'DefaultTextColor', FONT_COLOR, 'DefaultAxesColor', AXES_COLOR, ...
            'DefaultTextInterpreter', 'none', 'DefaultAxesYColor', 'k', ...
            'DefaultAxesZColor', 'k', 'DefaultPatchFaceColor', 'k', 'DefaultPatchEdgeColor', 'k',...
            'DefaultSurfaceEdgeColor', 'k', 'DefaultLineColor', 'k', 'DefaultUicontrolInterruptible', 'on', ...                                                     
            'Renderer', 'zbuffer', 'RendererMode' , 'manual', 'Visible', visible, 'Name', sprintf(TitleFig), ...
            'paperorientation', paperOrient, 'InvertHardcopy', 'off', 'resize', 'off', 'PaperPositionMode', 'auto');       