function varargout = icatb_plot_spectra(axesH, tc)
%% Plot Timecourse or spectra

icatb_defaults;
global FONT_COLOR;

fig_title = 'Power Spectra';
if (nargin == 1)
    tc = axesH;
    clear axesH;
    figH = icatb_getGraphics(fig_title, 'graphics', 'spectral_power', 'on');
    axesH = axes('parent', figH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
end

if (isfield(tc, 'fig_title'))
    set(get(axesH, 'parent'), 'name', fig_title);
end

fnames = fieldnames(tc);

isSpectra = 0;
titleStr = '';
xlabelStr = '';
ylabelStr = '';
convertToSpectra = 0;
TR = 1;

for nF = 1:length(fnames)
    commandToEval = [fnames{nF}, '=tc.', fnames{nF}, ';'];
    eval(commandToEval);
end

if (length(data) == numel(data))
    data = data(:)';
end

if (isSpectra && convertToSpectra)
    [data, xAxis] = icatb_get_spectra(data, TR);
end

if (~isSpectra && (size(data, 1) > 1))
    data = mean(data);
end

if (~exist('xAxis', 'var'))
    xAxis = (1:size(data, 2));
end

if (isSpectra)
    icatb_plot_with_ste_area(axesH, xAxis, data, [], 'm', [.5 .5 1]);
else
    plot(xAxis, data, 'm', 'parent', axesH);
end

xlabel(xlabelStr, 'parent', axesH, 'Interpreter', 'tex');
ylabel(ylabelStr, 'parent', axesH, 'Interpreter', 'tex');
title(titleStr, 'parent', axesH, 'Interpreter', 'tex');

set(axesH, 'XColor', FONT_COLOR, 'YColor', FONT_COLOR);

axis(axesH, 'tight');

varargout{1} = get(axesH, 'parent');
