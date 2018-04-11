function icatb_titleDialog
% Title dialog displays the authors, organization and the collaborators

icatb_defaults;
global BG_COLOR;
global BG2_COLOR;
global BUTTON_COLOR;
global BUTTON_FONT_COLOR;
global FG_COLOR;
global AXES_COLOR;
global FONT_COLOR;

global UI_FONTNAME;
global UI_FONTUNITS;
global UI_FS;

keystr1 = sprintf(['keyVal = get(gcbf, ''currentcharacter''); \n']);

keystr2 = sprintf(['if keyVal == 13 \n delete(gcbf); clear keyVal; \n end\n']);

keypressCallback = [keystr1 keystr2];

% set the defaults for the figure window
figProperty = {'Resize','off', 'windowstyle', 'modal', ...
    'MenuBar', 'none',...
    'DefaultTextColor', FONT_COLOR,...
    'DefaultTextInterpreter', 'none',...
    'DefaultAxesColor', AXES_COLOR,...
    'DefaultAxesXColor', 'k',...
    'DefaultAxesYColor', 'k',...
    'DefaultAxesZColor', 'k',...
    'DefaultPatchFaceColor', 'k',...
    'DefaultPatchEdgeColor', 'k',...
    'DefaultSurfaceEdgeColor', 'k',...
    'DefaultLineColor', 'k',...
    'DefaultUicontrolInterruptible', 'on',...
    'PaperType', 'usletter', ...
    'PaperUnits', 'normalized', ...
    'PaperPositionmode', 'auto', ...
    'InvertHardcopy', 'off',...
    'Renderer', 'zbuffer',...
    'color', BG_COLOR, 'resize', 'off', 'keypressfcn', keypressCallback, 'numbertitle', 'off'};

h = figure('name', 'About GroupICAT', figProperty{1:length(figProperty)});

pos = get(h, 'position');

screenSize = get(0, 'screensize');
figurePos(1) = 0.5*(screenSize(1) + screenSize(3)) - 0.5*pos(3);
figurePos(2) = 0.5*(screenSize(2) + screenSize(4)) - 0.5*pos(4);
figurePos(3) = 0.95*pos(3);
figurePos(4) =  0.7*pos(4);

set(h, 'position', figurePos);

% positions of variables
axisPos = [0 0 1 1];
titlePos = [0.5 0.9];
normalPos = [0.04 0.6];
okPos(3) = 0.16; okPos(4) = 0.08;
okPos(1) = 0.75 - 0.5*okPos(3); okPos(2) = 0.1;
eegiftPos(3) = 0.16; eegiftPos(4) = 0.08;
eegiftPos(1) = 0.25 - 0.5*eegiftPos(3); eegiftPos(2) = okPos(2);

% title color
titleColor = [0 0.9 0.9];

% fonts
titleFont = 16;
subtitleFont = 13;
normalFont = 11;

% set axis handle to off
axisHandle = axes('Parent', h, 'Position', axisPos, 'Visible', 'off');


% Name of the toolbox
text('units', 'normalized', 'string', 'Group ICA/IVA Toolbox', 'position', titlePos, 'fontsize', titleFont, 'HorizontalAlignment', 'center', ...
    'fontweight', 'bold', 'FontName', UI_FONTNAME, 'color', titleColor, 'FontAngle', 'italic', 'parent', axisHandle);

% change the position of the text
titlePos(2) = titlePos(2) - 0.1;

% Caption of the toolbox
text('units', 'normalized', 'string', 'GroupICAT v4.0b', 'position', titlePos, 'fontsize', subtitleFont, 'HorizontalAlignment', 'center', ...
    'fontweight', 'normal', 'FontName', UI_FONTNAME, 'parent', axisHandle);

% Display the remaining things like the version number, organization,
% authors
giftInfo = {'', '', ['\bfRelease Date: \rm', num2str('20-Feb-2017')], '', str2mat('\bfAuthors: \rm\bfThe GIFT Team'), ...
    '', '\bfOrganization: \rmThe MIND Research Network', '', '\bfWebsite: \rmhttp://mialab.mrn.org'};

% change the position of the text
normalPos(2) = titlePos(2) - 0.25;


text('units', 'normalized', 'string', giftInfo, 'position', normalPos, 'fontsize', normalFont, 'HorizontalAlignment', 'left', ...
    'fontweight', 'normal', 'FontName', UI_FONTNAME, 'interpreter', 'tex', 'parent', axisHandle);

% Callback for the ok button
old_pushString = 'GIFT';
gift_pushbutton = icatb_uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', old_pushString, 'callback', {@giftMoreInfo, h, figProperty, 'GIFT'}, ...
    'fontunits', UI_FONTUNITS, 'fontname', UI_FONTNAME, 'FontSize', UI_FS, 'position', okPos, ...
    'backgroundcolor', BUTTON_COLOR, 'ForegroundColor', BUTTON_FONT_COLOR);

clear old_pushString;

old_pushString = 'EEGIFT';

% More information give the info about the collaborators
eegift_button= icatb_uicontrol('parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', ...
    old_pushString, 'position', eegiftPos, 'fontunits', UI_FONTUNITS, 'fontname', UI_FONTNAME, 'FontSize', UI_FS, 'ForegroundColor', ...
    BUTTON_FONT_COLOR, 'backgroundcolor', BUTTON_COLOR, 'callback', {@giftMoreInfo, h, figProperty, 'EEGIFT'});


% click on push button to display more information about GIFT
function giftMoreInfo(hObject, eventdata, handles, figProperty, toolbox_name)

%
icatb_defaults;
global BG_COLOR;
global BG2_COLOR;
global BUTTON_COLOR;
global BUTTON_FONT_COLOR;
global FG_COLOR;
global AXES_COLOR;
global UI_FONTNAME;
global UI_FONTUNITS;
global FONT_COLOR;
global UI_FS;

% set the defaults
figHandle = figure('name', ['More Info about ', toolbox_name], figProperty{1:length(figProperty)});

pos = get(figHandle, 'position');

screenSize = get(0, 'screensize');
figurePos(1) = 0.5*(screenSize(1) + screenSize(3)) - 0.5*pos(3);
figurePos(2) = 0.5*(screenSize(2) + screenSize(4)) - 0.5*pos(4);
figurePos(3) = 0.95*pos(3);
figurePos(4) =  0.7*pos(4);

set(figHandle, 'position', figurePos);

% set axis off
axisHandle = axes('Parent', figHandle, 'Position', [0 0 1 1], 'Visible', 'off');

normalFont = 10;

% Display the credits

if strcmpi(toolbox_name, 'gift')
    % Add changes here
    D(1).string = 'Additional contributions were provided by:';
    D(size(D,2)+1).string = '1. Tulay Adali at the University of Maryland Baltimore County.';
    D(size(D,2)+1).string = '2. Jim Pekar at the FM Kirby Center for Functional Brain Imaging.';
    D(size(D,2)+1).string = '3. Baoming Hong and Kent Kiehl at the ONRC.';
    D(size(D,2)+1).string = '4. Betty Ann Chodkowski at the Kennedy Krieger Institute.';
    D(size(D,2)+1).string = '5. Y. Li, Nicolle Correa, Sai Ma and Matthew Anderson at the University of Maryland Baltimore County.';
    
    %% end for changes
    D(size(D,2)+1).string = '';
    D(size(D,2)+1).string = 'Many additional ICA algorithms were generously contributed by Andrzej Cichocki.';
    D(size(D,2)+1).string = '';
    D(size(D,2)+1).string = 'ICA algorithms also available in the ICALab toolbox at http://www.bsp.brain.riken.jp/ICALAB/ICALABImageProc/.';
    D(size(D,2)+1).string = '';
    D(size(D,2)+1).string = 'More information on ICA can be found in the book "Adaptive Blind Signal and Image Processing" by Andrzej Cichocki and Shun-ichi Amari.';
    D(size(D,2)+1).string = '';
    D(size(D,2)+1).string = 'GIFT is released under the GNU General Public License. GIFT uses some functions from the SPM library mainly the spm_vol family of functions. Please visit SPM at http://www.fil.ion.ucl.ac.uk/spm/.';
    D(size(D,2)+1).string = '';
    D(size(D,2)+1).string = 'For more information please contact vcalhoun@unm.edu or visit GIFT project website http://icatb.sourceforge.net/.';
    
else
    %% end for changes
    D(1).string = 'Many additional ICA algorithms were generously contributed by Andrzej Cichocki.';
    D(size(D,2)+1).string = '';
    D(size(D,2)+1).string = 'ICA algorithms also available in the ICALab toolbox at http://www.bsp.brain.riken.jp/ICALAB/ICALABImageProc/.';
    D(size(D,2)+1).string = '';
    D(size(D,2)+1).string = 'More information on ICA can be found in the book "Adaptive Blind Signal and Image Processing" by Andrzej Cichocki and Shun-ichi Amari.';
    D(size(D,2)+1).string = '';
    D(size(D,2)+1).string = 'EEGIFT is released under the GNU General Public License. EEGIFT uses some functions from the EEGLAB library. Please visit EEGLAB at http://sccn.ucsd.edu/eeglab/.';
    D(size(D,2)+1).string = '';
    D(size(D,2)+1).string = 'For more information please contact tom.eichele@psybp.uib or visit EEGIFT project website http://icatb.sourceforge.net/.';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

str = {D.string};

position = get(figHandle, 'position');

xOffSet = 0.02; yOffSet = 0.04;

okPos(3) = 0.12; okPos(4) = 0.08;
okPos(1) = 0.5 - 0.5*okPos(3); okPos(2) = yOffSet;


% define the position of the listbox
position(1) = xOffSet; position(2) = okPos(2) + okPos(4) + yOffSet; position(3) = 1 - 2*xOffSet; position(4) = 1 - position(2) - yOffSet;


handle_scroll = icatb_uicontrol('parent', figHandle, 'units', 'normalized', 'style','listbox', ...
    'position', position, 'string', str, 'foregroundcolor', FONT_COLOR, ...
    'horizontalalignment','left', 'backgroundcolor', BG2_COLOR, 'FontSize', normalFont, 'fontunits', ...
    UI_FONTUNITS, 'fontname', UI_FONTNAME);

% Apply conditions for dialog box differently for different platforms
if ispc
    set(handle_scroll, 'enable', 'inactive');
else
    set(handle_scroll, 'enable', 'on');
end


%[newString, newPos] = textwrap(handle_scroll, str);
maxChars = 75;

% wrap the string inside the uicontrol
[newString, newPos] = textwrap(handle_scroll, str, maxChars);

set(handle_scroll, 'String', newString);

set(handle_scroll, 'min', 0, 'max', 2);

% make no selection
set(handle_scroll, 'value', []);

cancelCallback = 'delete(gcbf)';

old_pushString = 'Return';

ok_pushbutton = icatb_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'pushbutton', 'string', ...
    old_pushString, 'callback', cancelCallback, 'position', okPos, 'fontunits', UI_FONTUNITS, ...
    'fontname', UI_FONTNAME, 'FontSize', UI_FS, 'backgroundcolor', BUTTON_COLOR, 'ForegroundColor', BUTTON_FONT_COLOR);

