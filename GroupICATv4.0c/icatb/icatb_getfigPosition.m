function [Rect] = icatb_getfigPosition(GraphicsType)
% Purpose: Get the figure position in pixel units
%
% Input:
% GraphicsType - Options are 'normal', 'graphics', 'displaygui',
% 'statusbar', 'timecourse'
 
icatb_defaults;
global BG_COLOR;
global FG_COLOR;
global AXES_COLOR;
global FONT_COLOR;
global ASPECT_RATIO_FOR_SQUARE_FIGURE;
global S0;
global WS;
global MIN_SCREEN_DIM_IN_PIXELS;


%---------Get Coordinates For Figure
S=MIN_SCREEN_DIM_IN_PIXELS;
R=ASPECT_RATIO_FOR_SQUARE_FIGURE;
switch lower(GraphicsType)
   
    case 'normal'
        extendLeft = round(S*R(1)*.5);
        extendUp = round(S*R(2)*.5);
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
    otherwise
        
        extendLeft = round(S*R(1)*.5);
        extendUp = round(S*R(2)*.5);
        x0= (S0(3)/2)-(extendLeft/2);
        y0 = (S0(4)/2)-(extendUp/2);
        Rect = [x0 y0 extendLeft extendUp];
        
end