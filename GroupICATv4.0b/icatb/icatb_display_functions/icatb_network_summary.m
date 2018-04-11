function icatb_network_summary(network_opts)


icatb_defaults;
global UI_FS;
global IMAGE_VALUES;
global CONVERT_Z;
global THRESHOLD_VALUE;
global FONT_COLOR;


image_values = 'positive';
thresholds = 1;
structFile = fullfile(fileparts(which('groupica.m')), 'icatb_templates', 'ch2bet.nii');
convert_to_z = 'yes';
prefix = 'IC';
fnc_colorbar_label = 'Corr';
fnc_matrix_file = [];

outputDir = pwd;
try
    outputDir = network_opts.outputDir;
catch
end

try
    file_names = network_opts.file_names;
catch
end

try
    fnc_matrix_file = network_opts.fnc_matrix_file;
catch
end

try
    image_values = network_opts.image_values;
catch
end

try
    thresholds = network_opts.threshold;
catch
end

try
    structFile = network_opts.structFile;
catch
end

try
    convert_to_z = network_opts.convert_to_z;
catch
end

try
    prefix = network_opts.prefix;
catch
end

try
    fnc_colorbar_label = network_opts.fnc_colorbar_label;
catch
end

conn_threshold = [];
try
    conn_threshold = network_opts.conn_threshold;
catch
end

if (ischar(conn_threshold))
    conn_threshold = str2num(conn_threshold);
end

display_type = 'slices';
try
    display_type = network_opts.display_type;
catch
end

slice_plane = 'sagittal';
try
    slice_plane = network_opts.slice_plane;
catch
end

imWidth = [];
try
    imWidth = network_opts.imWidth;
catch
end

if (ischar(imWidth))
    imWidth = str2num(imWidth);
end

load icatb_colors coldhot;
cmap_corr = coldhot(1:4:end, :);
try
    cmap_corr = network_opts.cmap;
catch
end

CLIM_corr = [];
try
    CLIM_corr = network_opts.CLIM;
catch
end

if (ischar(CLIM_corr))
    CLIM_corr = str2num(CLIM_corr);
end

comp_network_names = network_opts.comp_network_names;

% for i = 1:2:length(varargin)
%     if (strcmpi(varargin{i}, 'image_values'))
%         image_values = varargin{i + 1};
%     elseif (strcmpi(varargin{i}, 'convert_to_z'))
%         convert_to_z = varargin{i + 1};
%     elseif (strcmpi(varargin{i}, 'threshold'))
%         thresholds = varargin{i + 1};
%     elseif (strcmpi(varargin{i}, 'anatomical_file'))
%         structFile = varargin{i + 1};
%     elseif (strcmpi(varargin{i}, 'prefix'))
%         prefix = varargin{i + 1};
%         elseif (strcmpi(varargin{i}, 'fnc_colorbar_label'))
%             fnc_colorbar_label = varargin{i + 1};
%     end
%     end
% end
%
% if (~exist('file_names', 'var') || isempty(file_names))
%     file_names = icatb_selectEntry('typeSelection', 'multiple', 'typeEntity', 'file', 'filter', '*.img;*.nii', 'title', 'Select component nifti images ...', 'fileType', 'image');
% end
%
%
% useGUI = 0;
% if (~exist('comp_network_names', 'var'))
%     useGUI = 1;
% end
%
% if (isempty(file_names))
%     error('Files are not selected');
% end
%
% if (~exist('fnc_matrix_file', 'var'))
%     fnc_matrix_file = icatb_selectEntry('typeSelection', 'single', 'typeEntity', 'file', 'filter', '*.mat;*.txt;*.dat', 'title', 'Select FNC correlations file');
% end

drawnow;

file_names = icatb_rename_4d_file(file_names);

% Load fnc matrix
FNCM = loadFNC(fnc_matrix_file);


drawnow;


fncComps = comp_network_names(:, 2);
fncComps = [fncComps{:}];
fncComps = fncComps(:);

if (~isempty(FNCM))
    if ((size(FNCM, 1) ~= length(fncComps)) && (size(FNCM, 1) ~= size(file_names, 1)))
        error('FNC matrix dimensions does not match the number of components entered');
    end
    
    if (size(FNCM, 1) == size(file_names, 1))
        % Truncate
        FNCM = FNCM(fncComps, fncComps);
    end
end

if (length(thresholds) == 1)
    thresholds = thresholds.*ones(1, size(comp_network_names, 1));
end

lengths = cellfun(@length, comp_network_names(:,2), 'UniformOutput', false);
lengths = cellfun(@length, comp_network_names(:,2), 'UniformOutput', false);
[dd, maxRInds] = max([lengths{:}]);

%% Rendering: Multiple components are rendered on surface of a 'standard' brain
nrows = ceil(size(comp_network_names, 1)/sqrt(size(comp_network_names, 1)));
ncols = ceil(size(comp_network_names, 1)/nrows);
for nF = 1:size(comp_network_names,1)
    fH = icatb_getGraphics('Render', 'normal', 'render_view', 'on');
    set(fH, 'resize', 'on');
    set(fH,'menubar','none');
    
    rData = icatb_display_composite(file_names(comp_network_names{nF, 2}, :), 'convert_to_zscores', convert_to_z, 'threshold', thresholds(1), 'image_values', image_values, ...
        'display_type', 'render');
    
    nLabels = comp_network_names{nF, 2};
    nLabels = cellstr(num2str(nLabels(:)));
    sh = axes('parent', fH, 'units', 'normalized', 'position', [0.16, 0.16, 0.72, 0.72]);
    image(rData.slices, 'parent', sh);
    colormap(rData.cmap);
    title(comp_network_names{nF, 1}, 'parent', sh, 'color', FONT_COLOR);
    try
        aspRatio = size(rData.slices, 1)/size(rData.slices, 2);
        set(sh, 'PlotBoxAspectRatio', [1, aspRatio, 1]);
    catch
        axis(sh, 'image;');
    end
    axis(sh, 'off');
    ch = colorbar('peer', sh);
    set(sh, 'clim', [1, size(rData.cmap, 1)]);
    pos = get(ch, 'position');
    origPos = pos(4);
    pos(4) = 0.7*pos(4);
    pos(2) = pos(2) + 0.5*(origPos - pos(4));
    pos(1) = pos(1) + pos(3) + 0.03;
    colorLen = rData.colorLen;
    colorMid = ceil(colorLen/2) + 1;
    xticks = colorMid:colorLen:(lengths{nF}*colorLen);
    set(ch, 'ylim', [1, colorLen*lengths{nF}]);
    set(ch, 'position', pos);
    set(ch, 'ytick', xticks);
    set(ch, 'YtickLabel', nLabels);
    set(ch, 'color', FONT_COLOR);
    
    outFNames = fullfile(outputDir, [prefix, '_render_', num2str(nF), '.fig']);
    saveas(fH, outFNames);
    
end


%% FNC correlations: Correlations are visualized in a matrix plot
if (~isempty(FNCM))
    %     load icatb_colors coldhot;
    %     cmap = coldhot(1:4:end, :);
    FNCM(isfinite(FNCM) == 0) = 0;
    CLIM = max(abs(FNCM(:)));
    network_values = cellfun(@length, comp_network_names(:,2), 'UniformOutput', false);
    network_values = [network_values{:}];
    %fH = figure('color', 'w');
    gH = icatb_getGraphics('FNC Correlations', 'graphics', 'fnc_corrs', 'on');
    pos = get(gH, 'position');
    %set(gH, 'position', pos);
    icatb_plot_FNC(FNCM, [-CLIM, CLIM], cellstr(num2str(fncComps)), (1:length(fncComps)), gH, fnc_colorbar_label, gca, ...
        network_values, comp_network_names(:,1));
    %     set(findobj(gH, 'type', 'line'), 'linewidth', 1.5);
    %     set(findobj(gH, 'type', 'line'), 'color', [0, 0, 0]);
    colormap(jet);
    outFNames = fullfile(outputDir, [prefix, '_FNC.fig']);
    saveas(gH, outFNames);
    
    drawnow;
    
    %% Connectogram view - FNC correlations are shown using bezier curves and thumbnails of spatial maps are shown in a circle. Components within the same network are shown in the same color.
    fH = icatb_plot_connectogram([], comp_network_names, 'C', FNCM, 'threshold', thresholds(1), 'image_file_names', file_names, 'colorbar_label', fnc_colorbar_label, 'cmap', cmap_corr, ...
        'slice_plane', slice_plane, 'conn_threshold', conn_threshold, 'imwidth', imWidth, 'display_type', display_type, 'CLIM', CLIM_corr);
    outFNames = fullfile(outputDir, [prefix, '_connectogram.fig']);
    saveas(fH, outFNames);
    
end

pause(1);

drawnow;

%% Multiple components are displayed in a composite plot. Orthogonal slices are shown.
gH = icatb_plot_composite_orth_views(file_names, comp_network_names, 'convert_to_zscores', convert_to_z, 'threshold', thresholds(1), 'image_values', image_values, ...
    'anatomical_file', structFile);
outFNames = fullfile(outputDir, [prefix, '_composite_orth_views.fig']);
saveas(gH, outFNames);

pause(1);

drawnow;

%% Stacked ortho slices are shown for each component in the network. Title shows component numbers plotted from top to bottom
icatb_groupNetworks(file_names, comp_network_names(:, 1), comp_network_names(:, 2), 'convert_to_z', convert_to_z, 'threshold', thresholds(1), 'image_values', image_values, ...
    'structfile', structFile, 'interMediatePlot', 0);
outFNames = fullfile(outputDir, [prefix, '_orth_views.fig']);
saveas(gcf, outFNames);


function FNCM = loadFNC(file_name)

if (isempty(file_name))
    FNCM = [];
end

if (ischar(file_name))
    [pp, fn, extn] = fileparts(file_name);
    
    if (strcmpi(extn, '.nii') || strcmpi(extn, '.img'))
        FNCM = icatb_loadData(file_name);
    else
        FNCM = icatb_load_ascii_or_mat(file_name);
    end
    
else
    FNCM = file_name;
end

if (numel(FNCM) == length(FNCM))
    FNCM = icatb_vec2mat(FNCM);
end