function gH = icatb_plot_composite_orth_views(fname, comp_network_names, varargin)
%% Orthogonal views of components are plotted as composite maps.
%
% Inputs:
% 1. fname - Nifti file name
% 2. comp_network_names - Cell array containing network names and
% components
% 3. varargin - Variable arguments passed as pairs.
%

icatb_defaults;
global FONT_COLOR;


%% Initialise vars
image_values = 'positive';
thresholds = 1;
structFile = fullfile(fileparts(which('groupica.m')), 'icatb_templates', 'ch2bet.nii');
convert_to_z = 'yes';

%% Parse arguments
for i = 1:2:length(varargin)
    if (strcmpi(varargin{i}, 'image_values'))
        image_values = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'convert_to_zscores') || strcmpi(varargin{i}, 'convert_to_z'))
        convert_to_z = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'threshold'))
        thresholds = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'anatomical_file'))
        structFile = varargin{i + 1};
    end
end

if (length(thresholds) == 1)
    thresholds = thresholds(1)*ones(1, size(comp_network_names, 1));
end

fname = icatb_rename_4d_file(fname);

lengths = cellfun(@length, comp_network_names(:,2), 'UniformOutput', false);
lengths = cellfun(@length, comp_network_names(:,2), 'UniformOutput', false);
[dd, inds]=max([lengths{:}]);


gH = icatb_getGraphics('Composite Maps', 'graphics', 'orthoviews', 'on');
sz = get(0, 'screensize');
pos = [50, 50, sz(3) - 100, sz(4) - 100];
set(gH, 'position', pos);
set(gH, 'resize', 'on');

nrows = ceil(size(comp_network_names, 1)/sqrt(size(comp_network_names, 1)));
ncols = ceil(size(comp_network_names, 1)/nrows);

maxSize = [0, 0];
ud = cell(1, size(comp_network_names, 1));
data = ud;
coordStrs = ud;

for nC = 1:size(comp_network_names, 1)
    comps = fname(comp_network_names{nC, 2}, :);
    comps = cellstr(comps);
    compVals = comp_network_names{nC, 2};
    valsIn = {comps, 'threshold', thresholds(nC), 'convert_to_zscores', convert_to_z, 'anatomical_file', structFile, ...
        'image_values', image_values, 'title', comp_network_names{nC, 1}, 'colorbar_label', cellstr(num2str(compVals(:)))};
    
    try
        tmpCoords = comp_network_names{nC, 3};
        valsIn{end + 1} = 'coords';
        valsIn{end + 1} = tmpCoords;
    catch
    end
    
    S1 = icatb_display_composite(valsIn{:});
    
    ud{nC} = valsIn;
    data{nC} = S1.slices;
    coordStrs{nC} = S1.xlabels;
    
    maxSize = [max([size(data{nC}, 1), maxSize(1)]), max([size(data{nC}, 2), maxSize(2)])];
    
    if (nC == inds(1))
        cmap = S1.cmap;
    end
    
end

clear valsIn;

colormap(cmap);

for nC = 1:size(comp_network_names, 1)
    
    tmp2 = data{nC};
    tmp = zeros([maxSize, 3]);
    midSize = ceil([size(tmp, 1), size(tmp, 2)]/2);
    midTmp = ceil([size(tmp2, 1), size(tmp2, 2)]/2);
    org = midSize - midTmp;
    tmp(org(1) + 1: org(1) + size(tmp2, 1), org(2) + 1: org(2) + size(tmp2, 2), :) = tmp2;
    
    
    valsIn = ud{nC};
    coordStr = coordStrs{nC};
    sh = subplot(nrows, ncols, nC);
    imagesc(tmp);
    title(char(comp_network_names{nC, 1}, coordStr), 'fontsize', 10, 'parent', sh);
    axis(sh, 'image');
    set(sh, 'XTick', []);
    set(sh, 'YTick', []);
    set(sh, 'box', 'off');
    set(sh, 'color', 'none');
    set(sh, 'userdata', valsIn);
    
    compNum = comp_network_names{nC, 2};
    labels = cellstr(num2str(compNum(:)));
    ch = colorbar('peer', sh);
    set(sh, 'clim', [1, size(cmap, 1)]);
    pos = get(ch, 'position');
    pos(1) = pos(1) + pos(3) + 0.03;
    colorLen = 64;
    colorMid = ceil(colorLen/2) + 1;
    xticks = colorMid:colorLen:(lengths{nC}*colorLen);
    set(ch, 'ylim', [1, 64*lengths{nC}]);
    set(ch, 'position', pos);
    set(ch, 'ytick', xticks);
    set(ch, 'YtickLabel', labels);
    set(ch, 'color', FONT_COLOR);
    
end

set(gH, 'visible', 'on');
set(gH, 'WindowButtonDownFcn', @plotInteractive);

function plotInteractive(hObject, event_data, handles)

if (strcmpi(get(hObject, 'selectionType'), 'open'))
    if (strcmpi(get(gco, 'type'), 'image'))
        valsIn = get(get(gco, 'parent'), 'userdata');
        if (iscell(valsIn))
            h = helpdlg('Opening interactive viewer. Please wait ...', 'Composite viewer');
            icatb_display_composite(valsIn{:});
            try
                delete(h);
            catch
            end
        end
    end
end