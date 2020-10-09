function icatb_legend(varargin)
% By default in Matlab 7 legend color is black
% This function adjusts the legend color such that it is white
% in color for GIFT toolbox

icatb_defaults;

global LEGENDCOLOR;

whichVersion = str2num(version('-release'));
if isempty(whichVersion)
    whichVersion = 14;
end

if whichVersion >= 14
    legendH = legend(varargin{:});
    set(legendH, 'textColor', LEGENDCOLOR);
    set(legendH, 'EdgeColor', LEGENDCOLOR);
    set(legendH, 'Location', 'best');
else
    legend(varargin{:});
end