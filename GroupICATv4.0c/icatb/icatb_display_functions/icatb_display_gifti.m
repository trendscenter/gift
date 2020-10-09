function icatb_display_gifti(fname, varargin)
%% Display GIFTI images (code from CAT12 toolbox)
%
% 


warning off MATLAB:subscripting:noSubscriptsSpecified;

colorbar_title = '';
titleStr = '';

for n = 1:2:length(varargin)
    if (strcmpi(varargin{n}, 'colorbar_title'))
        colorbar_title = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'title'))
        titleStr = varargin{n + 1};
    end
end

global GIFTIH

GIFTIH.view = 1;
GIFTIH.show_transp = 1;
GIFTIH.clip         = [];
GIFTIH.clim         = [];
GIFTIH.XTick        = [];
GIFTIH.bkg_col      = [0 0 0];
GIFTIH.show_inv     = 0;
GIFTIH.no_neg       = 0;
GIFTIH.transp       = 1;
GIFTIH.col          = [.8 .8 .8; 1 .5 .5];
GIFTIH.FS           = (1:35);
GIFTIH.n_surf       = 1;
GIFTIH.thresh_value = 0;
GIFTIH.cursor_mode  = 1;
GIFTIH.text_mode    = 1;
GIFTIH.border_mode  = 0;
GIFTIH.is32k        = 0;
GIFTIH.str32k       = '';
GIFTIH.SPM_found    = 1;
GIFTIH.surf_sel     = 1;
GIFTIH.colorbar_title = colorbar_title;
GIFTIH.title = titleStr;


% positions
% ws = spm('Winsize', 'Graphics');
gH = icatb_getGraphics('Results', 'Graphics', 'graphics', 'on');
ws = get(gH, 'position');
delete(gH);
 ss = get(0, 'Screensize');
if 2.6 * ws(3) > ss(3)
    ws(3) = ws(3) / (2.6 * ws(3) / ss(3));
end



% result window with 5 surface views and alternative positions without top view and  only with lateral views
GIFTIH.viewpos = {[0.025 0.450 0.375 0.375;  0.025 0.450 0.375 0.375;  0.025 2.000 0.375 0.375],... % lh medial
    [0.025 0.025 0.375 0.375;  0.025 0.025 0.375 0.375;  0.175 0.350 0.175 0.350],... % lh lateral
    [0.600 0.450 0.375 0.375;  0.600 0.450 0.375 0.375;  0.600 2.000 0.375 0.375],... % rh medial
    [0.600 0.025 0.375 0.375;  0.600 0.025 0.375 0.375;  0.675 0.350 0.175 0.350],... % rh lateral
    [0.300 0.150 0.400 0.500;  0.300 2.000 0.400 0.500;  0.300 2.000 0.400 0.500],... % lh+rh top
    [0.400 0.750 0.200 0.225;  0.400 0.300 0.200 0.225;  0.400 0.750 0.200 0.225]};   % data plot

% change size and position of flatmaps for >= R20014b
if icatb_spm_check_version('matlab', '8.4') >= 0
    GIFTIH.viewpos{2}(3, :) = [-0.075 0.150 0.650 0.650]; % lh lateral
    GIFTIH.viewpos{4}(3, :) = [0.425 0.150 0.650 0.650]; % rh lateral
end

% figure 1 with result window
GIFTIH.pos{1} = struct( ...
    'fig', [10 10 round(2.6*ws(3)) ws(3)], ... % figure
    'cbar', [0.400 -0.150 0.200 0.300; 0.440 0.025 0.120 0.120]);% colorbar

% figure 2 with GUI
GIFTIH.pos{2} = struct(...
    'fig',   [2*ws(3)+10 10 0.6*ws(3) ws(3)],...
    'sel',   [0.290 0.930 0.425 0.050],...
    'nam',   [0.050 0.875 0.900 0.050],...
    'surf',  [0.050 0.800 0.425 0.050],'mview',   [0.525 0.800 0.425 0.050],...
    'text',  [0.050 0.725 0.425 0.050],'thresh',  [0.525 0.725 0.425 0.050],...
    'cmap',  [0.050 0.650 0.425 0.050],'atlas',   [0.525 0.650 0.425 0.050],...
    'cursor',[0.050 0.575 0.425 0.050],'border',  [0.525 0.575 0.425 0.050],...
    'info',  [0.050 0.500 0.425 0.050],'bkg',     [0.525 0.500 0.425 0.050],...
    'nocbar',[0.050 0.425 0.425 0.050],'transp',  [0.525 0.425 0.425 0.050],...
    'inv',   [0.050 0.350 0.425 0.050],'hide_neg',[0.525 0.350 0.425 0.050],...
    'ovmin', [0.050 0.175 0.425 0.070],'ovmax',   [0.525 0.175 0.425 0.070],...
    'save',  [0.050 0.050 0.425 0.050],'close',   [0.525 0.050 0.425 0.050]);

GIFTIH.figure = figure;
clf(GIFTIH.figure);

set(GIFTIH.figure, 'MenuBar', 'none', 'Position', GIFTIH.pos{1}.fig, ...
    'Name', GIFTIH.title, 'NumberTitle', 'off', 'Renderer', 'OpenGL');

GIFTIH.panel(1) = uipanel('Position',[0 0 1 1],'units','normalized','BackgroundColor',...
    GIFTIH.bkg_col,'BorderType','none');
%GIFTIH.panel(2) = uipanel('Position',[2/2.6 0 0.6/2.6 1],'units','normalized','BorderType','none','BackgroundColor',GIFTIH.col(1,:));

% define S structure that contains information for lh and rh
GIFTIH.S{1}.name = ''; GIFTIH.S{1}.side = 'lh';
GIFTIH.S{2}.name = ''; GIFTIH.S{2}.side = 'rh';







GIFTIH.S{1}.name = fname;
GIFTIH.S{2}.name = fname;

GIFTIH.logP = 0;
GIFTIH.title = titleStr;
meshes_merged = 0;

for ind = 1:2
    
    % read meshes
    GIFTIH.S{ind}.info = cat_surf_info(GIFTIH.S{ind}.name, 1);
    
    if GIFTIH.S{ind}.info(1).nvertices == 64984
        GIFTIH.str32k = '_32k';
        GIFTIH.is32k = 1;
    else
        GIFTIH.str32k = '';
        GIFTIH.is32k = 0;
    end
    
    if strcmp(GIFTIH.S{ind}.info(1).side, 'mesh')
        meshes_merged = 1;
        if ind == 1
            GIFTIH.S{ind}.info(1).side = 'lh';
            GIFTIH.S{ind}.info(1).Pmesh = fullfile(fileparts(which('gift')), 'icatb_templates', ...
                ['gifti_templates_surfaces' GIFTIH.str32k], 'lh.central.freesurfer.gii');
        else
            GIFTIH.S{ind}.info(1).side = 'rh';
            GIFTIH.S{ind}.info(1).Pmesh = fullfile(fileparts(which('gift')), 'icatb_templates',  ....
                ['gifti_templates_surfaces' GIFTIH.str32k], 'rh.central.freesurfer.gii');
        end
    end
    GIFTIH.S{ind}.M = gifti(GIFTIH.S{ind}.info(1).Pmesh);
    
    % get adjacency information
    GIFTIH.S{ind}.A = icatb_spm_mesh_adjacency(GIFTIH.S{ind}.M);
    
    % cdata found?
    if meshes_merged
        if ind == 1
            try
                Y = icatb_spm_data_read(icatb_spm_data_hdr_read(GIFTIH.S{ind}.name));
            catch
                error('No data in surfaces found.');
            end
            if GIFTIH.is32k
                GIFTIH.S{1}.Y = Y(1:32492, :);
                GIFTIH.S{2}.Y = Y(32493:end, :);
            else
                GIFTIH.S{1}.Y = Y(1:163842, :);
                GIFTIH.S{2}.Y = Y(163843:end, :);
            end
        end
    else
        gind = gifti(GIFTIH.S{ind}.name);
        if isfield(gind,'cdata')
            GIFTIH.S{ind}.Y = icatb_spm_data_read(icatb_spm_data_hdr_read(GIFTIH.S{ind}.name));
        else
            if ind == 1
                gind = gifti(GIFTIH.S{2}.name);
                if isfield(gind,'cdata')
                    if isnumeric(gind.cdata)
                        GIFTIH.S{1}.Y = zeros(size(gind.cdata));
                    else
                        GIFTIH.S{1}.Y = zeros(gind.cdata.dim);
                    end
                    GIFTIH.S{1}.name = '';
                else
                    error('No data in surfaces found.');
                end
            else
                gind = gifti(GIFTIH.S{1}.name);
                if isfield(gind,'cdata')
                    if isnumeric(gind.cdata)
                        GIFTIH.S{2}.Y = zeros(size(gind.cdata));
                    else
                        GIFTIH.S{2}.Y = zeros(gind.cdata.dim);
                    end
                    GIFTIH.S{2}.name = '';
                else
                    error('No data in surfaces found.');
                end
            end
        end
    end
end

display_results_all;


function display_results_all(obj, event_obj)
%-----------------------------------------------------------------------
global GIFTIH

if (size(GIFTIH.S{1}.Y) > 1 | size(GIFTIH.S{2}.Y) > 1) & min(min(GIFTIH.S{1}.Y(:)), min(GIFTIH.S{2}.Y(:))) < 0
    disp('Warning: Only results with positive values are displayed!');
end

% clear larger area and set background color to update labels and title
GIFTIH.Ha = axes('Parent', GIFTIH.panel(1), 'Position', [-.1 -.1 1.1 1.1], 'Color', GIFTIH.bkg_col);
cla(GIFTIH.Ha);

GIFTIH.renderer = get(GIFTIH.figure, 'Renderer');
set(GIFTIH.figure, 'Renderer', 'OpenGL');

%-Get mesh curvature and sulcal depth
%------------------------------------------------------------------
for i = 1:2
    g1 = gifti(fullfile(fileparts(which('gift')), 'icatb_templates', ['gifti_templates_surfaces' GIFTIH.str32k], [GIFTIH.S{i}.info(1).side '.mc.freesurfer.gii']));
    g2 = gifti(fullfile(fileparts(which('gift')), 'icatb_templates', ['gifti_templates_surfaces' GIFTIH.str32k], [GIFTIH.S{i}.info(1).side '.sqrtsulc.freesurfer.gii']));
    GIFTIH.S{i}.curv = cell(2, 1);
    GIFTIH.S{i}.curv{1} = g1.cdata;
    GIFTIH.S{i}.curv{2} = g2.cdata;
end

if GIFTIH.view == 1 % top view
    vv = [90 0; -90 0; -90 0; 90 0; 0 90];
else % bottom view
    vv = [90 0; -90 0; -90 0; 90 0; 0 -90];
end

for ind = 1:5
    display_results(ind, GIFTIH.viewpos{ind}(abs(GIFTIH.view), :), vv(ind, :));
end

% prepare dataplot axes
GIFTIH.dataplot = axes('Position', GIFTIH.viewpos{6}(abs(GIFTIH.view), :), 'Parent', GIFTIH.panel(1), 'Color', GIFTIH.bkg_col);
GIFTIH.figure = ancestor(GIFTIH.dataplot, 'figure');
try axes(GIFTIH.dataplot); end
axis off

% check whether data for left or right hemipshere are all non-zero
ind1 = find(GIFTIH.S{1}.Y(:) ~= 0);
ind2 = find(GIFTIH.S{2}.Y(:) ~= 0);

% estimate min value > 0 and min/max values
if ~isempty(ind1) & ~isempty(ind2)
    GIFTIH.S{1}.thresh = min(GIFTIH.S{1}.Y(GIFTIH.S{1}.Y(:) > 0));
    tmp = min(GIFTIH.S{2}.Y(GIFTIH.S{2}.Y(:) > 0));
    if ~isempty(tmp)
        GIFTIH.S{1}.thresh = min(GIFTIH.S{1}.thresh, tmp);
    end
    GIFTIH.S{1}.min = min(min(GIFTIH.S{1}.Y(~isinf(GIFTIH.S{1}.Y))), min(GIFTIH.S{2}.Y(~isinf(GIFTIH.S{2}.Y))));
    GIFTIH.S{1}.max = max(max(GIFTIH.S{1}.Y(~isinf(GIFTIH.S{1}.Y))), max(GIFTIH.S{2}.Y(~isinf(GIFTIH.S{2}.Y))));
elseif isempty(ind1)
    GIFTIH.S{1}.thresh = min(GIFTIH.S{2}.Y(GIFTIH.S{2}.Y(:) > 0));
    GIFTIH.S{1}.min = min(GIFTIH.S{2}.Y(~isinf(GIFTIH.S{2}.Y)));
    GIFTIH.S{1}.max = max(GIFTIH.S{2}.Y(~isinf(GIFTIH.S{2}.Y)));
elseif isempty(ind2)
    GIFTIH.S{1}.thresh = min(GIFTIH.S{1}.Y(GIFTIH.S{1}.Y(:) > 0));
    GIFTIH.S{1}.min = min(GIFTIH.S{1}.Y(~isinf(GIFTIH.S{1}.Y)));
    GIFTIH.S{1}.max = max(GIFTIH.S{1}.Y(~isinf(GIFTIH.S{1}.Y)));
end

% deal with neg. values
if GIFTIH.S{1}.min < 0
    mnx = max(abs([GIFTIH.S{1}.min, GIFTIH.S{1}.max]));
    GIFTIH.S{1}.min = - mnx;
    GIFTIH.S{1}.max = mnx;
end

% add 10% to min/max values
GIFTIH.S{1}.max = round(1100 * GIFTIH.S{1}.max) / 1000;
if GIFTIH.S{1}.min < 0
    GIFTIH.S{1}.min = round(1100 * GIFTIH.S{1}.min) / 1000;
else
    GIFTIH.S{1}.min = round(900 * GIFTIH.S{1}.min) / 1000;
end

GIFTIH.clim = [true GIFTIH.S{1}.min GIFTIH.S{1}.max];

% only apply thresholds that are slightly larger than zero
if GIFTIH.S{1}.thresh > 0.00015
    GIFTIH.clip = [true -GIFTIH.S{1}.thresh GIFTIH.S{1}.thresh];
end

for ind = 1:5
    if GIFTIH.S{1}.thresh > 0.00015
        setappdata(GIFTIH.patch(ind), 'clip', GIFTIH.clip);
    end
    setappdata(GIFTIH.patch(ind), 'clim', [true GIFTIH.S{1}.min GIFTIH.S{1}.max]);
    col = getappdata(GIFTIH.patch(ind), 'col');
    d = getappdata(GIFTIH.patch(ind), 'data');
    GIFTIH = updateTexture(GIFTIH, ind, d, col, GIFTIH.show_transp);
end



th = title(GIFTIH.title, 'parent', get(GIFTIH.patch(end),'parent'), 'color', [1, 1,1]);
set(th,'units', 'normalized');
thpos = get(th,'position');
thpos(2) = thpos(2) + 0.3;
set(th, 'position', thpos);

% % only show threshold popup if log-name was found and minimal value > 0 is < 1
% if H.logP & (H.S{1}.thresh < 1)
%     set(H.thresh, 'Enable', 'on');
% end

% if H.n_surf == 1
%     % get sure that image is thresholded and there are at least 20% zero/NaN areas
%     if (sum(d ~= 0) / numel(d) < 0.8)
%         set(H.atlas, 'Enable', 'on');
%     end
% end

%if ~H.disable_cbar
GIFTIH = show_colorbar(GIFTIH);
%end

% % show slider for range of results
% if H.n_surf == 1
%
%     % allow slider a more extended range
%     mnx = ceil(2 * max(abs([H.S{1}.min H.S{1}.max])));
%
%     H.slider_min = sliderPanel( ...
%         'Parent', H.panel(2), ...
%         'Title', 'Overlay min', ...
%         'Position', H.pos{2}.ovmin, ...
%         'Backgroundcolor', H.col(1,:), ...
%         'Min', - mnx, ...
%         'Max', mnx, ...
%         'Value', H.S{1}.min, ...
%         'FontName', 'Verdana', ...
%         'FontSize', 8, ...
%         'NumFormat', '%f', ...
%         'Callback', @slider_clim_min);
%
%     H.slider_max = sliderPanel( ...
%         'Parent', H.panel(2), ...
%         'Title', 'Overlay max', ...
%         'Position', H.pos{2}.ovmax, ...
%         'Backgroundcolor', H.col(1,:), ...
%         'Min', - mnx, ...
%         'Max', mnx, ...
%         'Value', H.S{1}.max, ...
%         'FontName', 'Verdana', ...
%         'FontSize', 8, ...
%         'NumFormat', '%f', ...
%         'Callback', @slider_clim_max);
% end

%-----------------------------------------------------------------------
function GIFTIH = show_colorbar(GIFTIH)
%-----------------------------------------------------------------------

% show colorbar
if GIFTIH.n_surf == 1
    
    if isfield(GIFTIH, 'cbar')
        try delete(GIFTIH.cbar); end
        GIFTIH = rmfield(GIFTIH, 'cbar');
    end
    
    GIFTIH.cbar = axes('Parent', GIFTIH.panel(1), 'Position', GIFTIH.pos{1}.cbar(1, :), 'Color', GIFTIH.bkg_col, 'Visible', 'off');
    GIFTIH.colourbar = colorbar('peer', GIFTIH.cbar, 'Northoutside');
    
    %if GIFTIH.logP, title(GIFTIH.cbar, 'p-value', 'Color', 1 - GIFTIH.bkg_col); end
    
    title(GIFTIH.cbar, GIFTIH.colorbar_title, 'Color', 1 - GIFTIH.bkg_col);
    clim = getappdata(GIFTIH.patch(1), 'clim');
    axis(GIFTIH.cbar, 'off');
    
    if clim(3) > clim(2)
        caxis([clim(2) clim(3)]);
    end
    
    col = getappdata(GIFTIH.patch(1), 'col');
    colormap(col);
    
    % Update colorbar colors if clipping is used
    clip = getappdata(GIFTIH.patch(1), 'clip');
    if ~isempty(clip)
        if ~isnan(clip(2)) & ~isnan(clip(3))
            ncol = length(col);
            col_step = (clim(3) - clim(2)) / ncol;
            cmin = max([1, ceil((clip(2) - clim(2)) / col_step)]);
            cmax = min([ncol, floor((clip(3) - clim(2)) / col_step)]);
            col(cmin:cmax, :) = repmat([0.5 0.5 0.5], (cmax - cmin + 1), 1);
            colormap(col);
        end
    end
    
    if GIFTIH.logP
        
        XTick = get(GIFTIH.colourbar, 'XTick');
        
        % save original XTick values
        if isempty(GIFTIH.XTick), GIFTIH.XTick = XTick; end
        
        % if threshold is between 1.3..1.4 (p<0.05) change XTick accordingly and correct by 0.3
        if ~isempty(clip)
            if clip(3) >= 1.3 & clip(3) <= 1.4
                XTick_step = ceil((clim(3) - clim(2)) / 5);
                if clip(2) <= - 1.3 & clip(2) >= - 1.4
                    XTick = [(round(clim(2)) - 0.3):XTick_step: - 1.3 0 1.3:XTick_step:(round(clim(3)) + 0.3)];
                else
                    XTick = [0 1.3:XTick_step:(round(clim(3)) + 0.3)];
                end
            else
                mn = floor(min(XTick));
                mx = ceil(max(XTick));
                
                % only allow integer values
                XTick = floor(mn:mx);
                %                if ~isempty(GIFTIH.XTick), XTick = GIFTIH.XTick; end
            end
        else
            % rescue original XThick values if clipping is changed
            if ~isempty(GIFTIH.XTick), XTick = GIFTIH.XTick; end
        end
        
        % change XTickLabel
        XTickLabel = [];
        for i = 1:length(XTick)
            if XTick(i) > 0
                XTickLabel = char(XTickLabel, remove_zeros(sprintf('%.g', 10^(-XTick(i)))));
            elseif XTick(i) < 0
                XTickLabel = char(XTickLabel, remove_zeros(sprintf('-%.g', 10^(XTick(i)))));
            else
                XTickLabel = char(XTickLabel, '');
            end
        end
        
        set(GIFTIH.colourbar, 'XTickLabel', XTickLabel(2:end, :), 'XTick', XTick);
    end % end GIFTIH.logP
    
    set(GIFTIH.colourbar, 'XColor', 1-GIFTIH.bkg_col, 'YColor', 1-GIFTIH.bkg_col, 'TickLength',0);
else
    
    if ~isfield(GIFTIH, 'cbar') || ~ishandle(GIFTIH.cbar)
        GIFTIH.cbar = axes('Parent', GIFTIH.panel(1), 'Position', GIFTIH.pos{1}.cbar(2, :), 'Color', GIFTIH.bkg_col, 'Enable', 'off');
    end
    
    % RGB colorbar
    if GIFTIH.n_surf == 3
        cb = [8 1 1 4 2 2 8; ...
            8 1 6 7 5 2 8; ...
            8 8 3 3 3 8 8];
    else %RG colorbar
        cb = [8 1 1 4 2 2 8; ...
            8 1 1 4 2 2 8];
    end
    imagesc(cb);
    colormap([1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 1 1 1; GIFTIH.bkg_col]);
    axis(GIFTIH.cbar, 'off'); axis('image');
end

%-----------------------------------------------------------------------
function display_results(ind, win, vw)
%-----------------------------------------------------------------------
global GIFTIH

% rescue old color before a new GIFTIH.patch is created
try
    col = getappdata(GIFTIH.patch(ind), 'col');
catch
    col = [];
end

if ind < 5 % single hemisphere views
    M = GIFTIH.S{round(ind / 2)}.M;
    Mc.cdata = GIFTIH.S{round(ind / 2)}.Y;
else
    Ml = GIFTIH.S{1}.M;
    Mr = GIFTIH.S{2}.M;
    Mcl.cdata = GIFTIH.S{1}.Y;
    Mcr.cdata = GIFTIH.S{2}.Y;
    
    % check whether number of data for lh/rh differ and fill with zeros
    diff_size_Y = size(GIFTIH.S{1}.Y, 2) - size(GIFTIH.S{2}.Y, 2);
    if diff_size_Y > 0
        Mcr.cdata = [Mcr.cdata zeros(size(GIFTIH.S{2}.Y, 1), 1)];
    end
    if diff_size_Y < 0
        Mcl.cdata = [Mcl.cdata; zeros(size(GIFTIH.S{1}.Y, 1), 1)];
    end
    
    M.faces = [Ml.faces; Mr.faces + size(Ml.vertices, 1)];
    M.vertices = [Ml.vertices; Mr.vertices];
    M.mat = Ml.mat;
    Mc.cdata = [Mcl.cdata; Mcr.cdata];
end

if isfield(Mc, 'cdata')
    M.cdata = Mc.cdata;
else
    M.cdata = [];
end

GIFTIH.axis = axes('Position', win, 'Parent', GIFTIH.panel(1), 'Visible', 'off');
GIFTIH.figure = ancestor(GIFTIH.axis, 'figure');
axes(GIFTIH.axis);

if isfield(M, 'facevertexcdata')
    GIFTIH.cdata = M.facevertexcdata;
else
    GIFTIH.cdata = [];
end

if ~isfield(M, 'vertices') || ~isfield(M, 'faces')
    error('cat_surf_results:nomesh', 'ERROR:cat_surf_render: No input mesh.');
end

%% -Patch
%------------------------------------------------------------------
P = struct('vertices', M.vertices, 'faces', double(M.faces));
GIFTIH.patch(ind) = patch(P, ...
    'FaceColor', [0.6 0.6 0.6], ...
    'EdgeColor', 'none', ...
    'FaceLighting', 'gouraud', ...
    'SpecularStrength', 0.1, ...
    'AmbientStrength', 1.0, ...
    'DiffuseStrength', 0.6, ...
    'SpecularExponent', 15, ...
    'Clipping', 'off', ...
    'Visible', 'off', ...
    'Tag', 'CATSurfRender', ...
    'Parent', GIFTIH.axis);
setappdata(GIFTIH.patch(ind), 'patch', P);
setappdata(GIFTIH.patch(ind), 'axis', GIFTIH.axis);

%-Apply texture to mesh
%------------------------------------------------------------------
if isfield(M, 'facevertexcdata')
    T = M.facevertexcdata;
elseif isfield(M, 'cdata')
    T = M.cdata;
else
    T = [];
end

if isempty(col)
    GIFTIH = updateTexture(GIFTIH, ind, T);
else
    GIFTIH = updateTexture(GIFTIH, ind, T, col);
end

axis(GIFTIH.axis, 'image');
axis(GIFTIH.axis, 'off');
view(GIFTIH.axis, vw);
material(GIFTIH.figure, 'dull');

% default lighting
GIFTIH.light(1) = camlight('headlight'); set(GIFTIH.light(1), 'Parent', GIFTIH.axis);
setappdata(GIFTIH.axis, 'handles', GIFTIH);
set(GIFTIH.patch(ind), 'Visible', 'on');
camlight(GIFTIH.light(1),'headlight')

%==========================================================================
function [GIFTIH, C] = updateTexture(GIFTIH, ind, v, col, transp)

%-Project data onto surface mesh
%--------------------------------------------------------------------------
if size(v, 2) < size(v, 1)
    v = v';
end
v(isinf(v)) = NaN;

%-Get colourmap
%--------------------------------------------------------------------------
if ~exist('col', 'var')
    if size(v, 1) == 1
        col = jet(256);
    else
        % use RGB colormap
        col = zeros(256, 3, size(v, 1));
        for i = 1:3
            col(:, i, i) = 1;
        end
    end
end

setappdata(GIFTIH.patch(ind), 'data', v);
setappdata(GIFTIH.patch(ind), 'col', col);

if ~exist('FaceColor', 'var') || isempty(FaceColor), FaceColor = 'interp'; end

%-Get curvature
%--------------------------------------------------------------------------
if ind < 5 % single hemisphere views
    curv = GIFTIH.S{round(ind / 2)}.curv{GIFTIH.text_mode};
else
    curv = [GIFTIH.S{1}.curv{GIFTIH.text_mode}; GIFTIH.S{2}.curv{GIFTIH.text_mode}];
end

if size(curv, 2) == 1
    
    % emphasize mean curvature values by using sqrt
    %  if GIFTIH.text_mode==1
    if 1
        indneg = find(curv < 0);
        curv(indneg) = - ((-curv(indneg)).^0.5);
        indpos = find(curv > 0);
        curv(indpos) = (curv(indpos).^0.5);
        curv = curv - min(curv);
    end
    
    curv = 0.5 + repmat(curv, 1, 3);
    curv = curv / max(curv(:));
    
    % for sulcal depth (with no neg. values) use inverted values
    if GIFTIH.text_mode == 2
        curv = 1 - curv;
    end
end

%-Create RGB representation of data according to colourmap
%--------------------------------------------------------------------------
C = zeros(size(v, 2), 3);
clim = getappdata(GIFTIH.patch(ind), 'clim');
if isempty(clim), clim = [false NaN NaN]; end
mi = clim(2); ma = clim(3);

if any(v(:))
    if ~clim(1), mi = min(v(:)); ma = max(v(:)); end
    % don't allow negative values for multiple maps
    if size(v, 1) > 1 & mi < 0
        if ~isempty(GIFTIH.clip)
            GIFTIH.clip(2) = - Inf;
        else
            GIFTIH.clip = [true -Inf 0];
        end
    end
    for i = 1:size(v, 1)
        C = C + squeeze(ind2rgb(floor(((v(i, :) - mi) / (ma - mi)) * size(col, 1)), col(:, :, i)));
    end
end

if ~isempty(GIFTIH.clip)
    v(v > GIFTIH.clip(2) & v < GIFTIH.clip(3)) = NaN;
    setappdata(GIFTIH.patch(ind), 'clip', [true GIFTIH.clip(2) GIFTIH.clip(3)]);
end

setappdata(GIFTIH.patch(ind), 'clim', [true mi ma]);
GIFTIH.clim = [true mi ma];

%-Build texture by merging curvature and data
%--------------------------------------------------------------------------
if size(v, 1) > 1 % RGB
    for i = 1:size(v, 1)
        C(:, i) = any(v(i, :), 1)' .* C(:, i);
    end
else
    C = repmat(any(v, 1), 3, 1)' .* C;
end

% add curvature pattern if transparency is defined
if nargin > 4
    if transp & size(C, 1) == size(curv, 1)
        C = (0.5 + 0.5 * curv) .* C;
    end
end

% replace regions below threshold by curvature
ind0 = repmat(~any(v, 1), 3, 1)';
if size(C, 1) == size(curv, 1)
    C(ind0) = curv(ind0);
else
    C(ind0) = 0;
end

%-Add atlas border
%--------------------------------------------------------------------------
if GIFTIH.border_mode
    if ind < 5 % single hemisphere views
        A = GIFTIH.S{round(ind / 2)}.A;
        A = sparse(1:size(GIFTIH.S{round(ind / 2)}.M.vertices, 1), 1:size(GIFTIH.S{round(ind / 2)}.M.vertices, 1), 1 ./ sum(A, 2)) * A;
        rdata = GIFTIH.rdata{GIFTIH.border_mode}(:, round(ind / 2));
        C0 = (A - speye(size(A))) * double(rdata);
        C(round(C0) ~= 0, :) = 0;
    else
        C0 = [];
        for i = 1:2
            A = GIFTIH.S{i}.A;
            A = sparse(1:size(GIFTIH.S{i}.M.vertices, 1), 1:size(GIFTIH.S{i}.M.vertices, 1), 1 ./ sum(A, 2)) * A;
            rdata = GIFTIH.rdata{GIFTIH.border_mode}(:, i);
            C0 = [C0; (A - speye(size(A))) * double(rdata)];
        end
        C(round(C0) ~= 0, :) = 0;
    end
end

set(GIFTIH.patch(ind), 'FaceVertexCData', C, 'FaceColor', FaceColor);


function s = remove_zeros(s)

pos = length(s);
while pos > 1
    if strcmp(s(pos), '0')
        s(pos) = '';
        pos = pos - 1;
    else break
    end
end



function [varargout] = cat_surf_info(P,readsurf,gui,verb)
% ______________________________________________________________________
% Extact surface information from filename.
%
% sinfo = cat_surf_info(P,readsurf,gui,verb)
%
%   P         .. surface filename
%   readsurf  .. read gifti or Freesurfer file to get more information
%   gui       .. interactive hemisphere selection
%   verb      .. verbose
%
% sinfo(i).
%   fname     .. full filename
%   pp        .. filepath
%   ff        .. filename
%   ee        .. filetype
%   exist     .. exist file?
%   fdata     .. structure from dir command
%   ftype     .. filetype [0=no surface,1=gifti,2=freesurfer]
%   statready .. ready for statistic (^s#.*.gii) [0|1]
%   side      .. hemisphere [lh|rh|lc|rc|mesh]
%   name      .. subject/template name
%   datatype  .. [-1=unknown|0=nosurf|1=mesh|2=data|3=surf]
%                only defined for readsurf==1 and surf=mesh+sidata
%   dataname  .. datafieldname [central|thickness|intensity...]
%   texture   .. textureclass [central|sphere|thickness|...]
%   label     .. labelmap
%   resampled .. resampled data [0|1]
%   template  .. template or individual mesh [0|1]
%   name      .. name of the dataset
%   roi       .. roi data
%   nvertices .. number vertices
%   nfaces    .. number faces
%   Pmesh     .. underlying meshfile
%   Psphere   .. sphere mesh
%   Pspherereg.. registered sphere mesh
%   Pdefects  .. topology defects mesh
%   Pdata     .. datafile
%   preside   .. prefix before hemi info (i.e. after smoothing)
%   posside   .. string after hemi info
%   smoothed  .. smoothing size
%   Phull     .. convex hull mesh
% ______________________________________________________________________
% Robert Dahnke
% $Id: cat_surf_info.m 1311 2018-04-26 08:16:03Z dahnke $

%#ok<*RGXP1>

if ~exist('P','var'), P=''; end
if strcmp(P,'selftest')
    pps = {
        fullfile(fileparts(which('gift')),'icatb_templates','gifti_templates_surfaces')
        fullfile('User','08.15','T1 T2','subs','mri')
        };
    ffs = {
        'lh.central.freesurfer'
        'lh.mymask'
        'output'
        'lh.texture.sub1.sub2'
        'lh.tex1.tex2.resampled.sub1.sub2'
        '08.15'
        's15mm.lh.tex03.33.resampled.S01.mri'
        's5mm.lh.t1.t2-3_3.resampled.S01_-.kdk.mri'
        'rh.s33mmtexture.S01.native.mri'
        'rh'
        'rh.'
        'rh.sphere.reg.sub1'
        'rc.defects.038.37.477'
        'lc.s33mmtexture.S01.native.mri'
        'rc.texture.sub1.sub2'
        };
    ees = {
        ''
        ... '.gii'
        ... '.annot'
        };
    varargout = cell(numel(pps),numel(ffs),numel(ees));
    for ppsi = 1:numel(pps)
        for ffsi = 1:numel(ffs)
            for eesi = 1:numel(ees)
                varargout{1}(ppsi,ffsi,eesi) = cat_surf_info(fullfile(pps{ppsi},[ffs{ffsi} ees{eesi}]),0,0,1);
            end
        end
    end
    return;
end



if nargin<2, readsurf = 0; end
if nargin<3, gui  = 0; end
if nargin<4, verb = 0; end

P = cellstr(P);

sinfo = struct(...
    'fname','',...      % full filename
    'pp','',...         % filepath
    'ff','',...         % filename
    'ee','',...         % filetype
    'exist','',...      % exist
    'fdata','',...      % datainfo (filesize)
    'ftype','',...      % filetype [0=no surface,1=gifti,2=freesurfer]
    'statready',0,...   % ready for statistic (^s#.*.gii)
    'side','',...       % hemishphere
    'name','',...       % subject/template name
    'datatype','',...   % datatype [0=nosurf/file|1=mesh|2=data|3=surf] with surf=mesh+data
    'dataname','',...   % datafieldname [central|thickness|s3thickness...]
    'texture','',...    % textureclass [central|sphere|thickness|...]
    'label','',...      % labelmap
    'resampled','',...  % dataspace
    'template','',...   % individual surface or tempalte
    'roi','',...        % roi data
    'nvertices',[],...  % number vertices
    'nfaces',[],...     % number faces
    'Pmesh','',...      % meshfile
    'Psphere','',...    % meshfile
    'Pspherereg','',... % meshfile
    'Pdefects','',...   % meshfile
    'Pdata','',...      % datafile
    'preside','', ...
    'posside','' ...
    );

if isempty(P{1}), varargout{1}=sinfo; return; end

for i=1:numel(P)
    [pp,ff,ee] = icatb_spm_fileparts(P{i});
    sinfo(i).fdata = dir(P{i});
    
    sinfo(i).fname = P{i};
    sinfo(i).exist = exist(P{i},'file') > 0;
    sinfo(i).pp = pp;
    switch ee
        case {'.xml','.txt','.html','.csv'}
            sinfo(i).ff = ff;
            sinfo(i).ee = ee;
            sinfo(i).ftype = 0;
            continue
        case '.gii'
            sinfo(i).ff = ff;
            sinfo(i).ee = ee;
            sinfo(i).ftype = 1;
            if sinfo(i).exist && readsurf
                S = gifti(P{i});
            end
        case '.annot'
            sinfo(i).ff = ff;
            sinfo(i).ee = ee;
            sinfo(i).ftype = 1;
            sinfo(i).label = 1;
            if sinfo(i).exist && readsurf
                clear S;
                try
                    S = cat_io_FreeSurfer('read_annotation',P{1});
                end
            end
            if exist('S','var')
                sinfo(i).ftype = 2;
            end
        otherwise
            sinfo(i).ff = [ff ee];
            sinfo(i).ee = '';
            sinfo(i).ftype = 0;
            if sinfo(i).exist && readsurf
                clear S;
                try
                    S = cat_io_FreeSurfer('read_surf',P{1});
                    if size(S.faces,2)~=3 || size(S.faces,1)<10000
                        clear S;
                    end
                end
                try
                    S.cdata = cat_io_FreeSurfer('read_surf_data',P{1});
                    if size(S.face,2)==3 || size(S.face,1)<10000
                        S = rmfield(S,'cdata');
                    end
                end
            end
            if exist('S','var')
                sinfo(i).ftype = 2;
            end
    end
    
    
    noname = sinfo(i).ff;
    
    % smoothed data
    sinfo(i).statready = ~isempty(regexp(noname,'^s(?<smooth>\d+)\..*'));
    
    % side
    if     strfind(noname,'lh.'),   sinfo(i).side='lh';   sidei = strfind(noname,'lh.');
    elseif strfind(noname,'rh.'),   sinfo(i).side='rh';   sidei = strfind(noname,'rh.');
    elseif strfind(noname,'mesh.'), sinfo(i).side='mesh'; sidei = strfind(noname,'mesh.');
    elseif strfind(noname,'lc.'),   sinfo(i).side='lc';   sidei = strfind(noname,'lc.');
    elseif strfind(noname,'rc.'),   sinfo(i).side='rc';   sidei = strfind(noname,'rc.');
    else
        % if SPM.mat exist use that for side information
        if exist(fullfile(pp,'SPM.mat'),'file')
            load(fullfile(pp,'SPM.mat'));
            [pp2,ff2]   = icatb_spm_fileparts(SPM.xY.VY(1).fname);
            
            % find mesh string
            hemi_ind = strfind(ff2,'mesh.');
            if ~isempty(hemi_ind)
                sinfo(i).side = ff2(hemi_ind(1):hemi_ind(1)+3);
            else
                % find lh|rh string
                hemi_ind = [strfind(ff2,'lh.') strfind(ff2,'rh.') strfind(ff2,'lc.') strfind(ff2,'rc.')];
                sinfo(i).side = ff2(hemi_ind(1):hemi_ind(1)+1);
            end
            
            sidei=[];
        else
            if gui
                if cat_get_defaults('extopts.expertgui')
                    sinfo(i).side = icatb_spm_input('Hemisphere',1,'lh|rh|lc|rc|mesh');
                else
                    sinfo(i).side = icatb_spm_input('Hemisphere',1,'lh|rh|mesh');
                end
            else
                sinfo(i).side = '';
            end
            sidei = strfind(noname,[sinfo(i).side '.']);
        end
    end
    if isempty(sidei), sidei = strfind(noname,sinfo(i).side); end
    if sidei>0
        sinfo(i).preside = noname(1:sidei-1);
        sinfo(i).posside = noname(sidei+numel(sinfo(i).side)+1:end);
    else
        sinfo(i).posside = noname;
    end
    
    % smoothed
    if isempty(sinfo(i).preside)
        sinfo(i).smoothed = 0;
    else
        sinfo(i).smoothed = max([0,double(cell2mat(textscan(sinfo(i).preside,'s%dmm.')))]);
    end
    
    % datatype
    if sinfo(i).exist && readsurf
        switch num2str([isfield(S,'vertices'),isfield(S,'cdata')],'%d%d')
            case '00',  sinfo(i).datatype  = 0;
            case '01',  sinfo(i).datatype  = 1;
            case '10',  sinfo(i).datatype  = 2;
            case '11',  sinfo(i).datatype  = 3;
        end
    else
        sinfo(i).datatype = -1;
    end
    
    
    % resampled
    sinfo(i).resampled = ~isempty(strfind(sinfo(i).posside,'.resampled'));
    % template
    sinfo(i).template  = ~isempty(strfind(lower(sinfo(i).ff),'.template'));
    if sinfo(i).template,  sinfo(i).resampled = 1; end
    
    
    % name / texture
    %  -----------------------------------------------------------------
    % ... name extraction is a problem, because the name can include points
    % and also the dataname / texture can include points ...
    resi = [strfind(sinfo(i).posside,'template.'),...
        strfind(sinfo(i).posside,'resampled.'),...
        strfind(sinfo(i).posside,'sphere.reg.')];
    if ~isempty(resi)
        sinfo(i).name = cat_io_strrep(sinfo(i).posside(max(resi):end),...
            {'template.','resampled.','sphere.reg'},''); %sinfo(i).posside,
        if ~isempty(sinfo(i).name) && sinfo(i).name(1)=='.', sinfo(i).name(1)=[]; end
        sinfo(i).texture = sinfo(i).posside(1:min(resi)-2);
    else
        % without no template/resampled string
        doti = strfind(sinfo(i).posside,'.');
        if numel(doti)==0
            % if not points exist that the string is the name
            sinfo(i).name    = '';
            sinfo(i).texture = sinfo(i).posside;
        elseif numel(doti)==1
            % if one point exist that the first string is the dataname and the second the subject name
            sinfo(i).name    = sinfo(i).posside(doti+1:end);
            sinfo(i).texture = sinfo(i).posside(1:doti-1);
        else
            % this is bad
            sinfo(i).name    = sinfo(i).posside(min(doti)+1:end);
            sinfo(i).texture = sinfo(i).posside(1:min(doti)-1);
        end
    end
    if verb
        fprintf('%50s: s%04.1f %2s ',sinfo(i).ff,sinfo(i).smoothed,sinfo(i).side);
        cat_io_cprintf([0.2 0.2 0.8],'%15s',sinfo(i).texture);
        cat_io_cprintf([0.0 0.5 0.2],'%15s',sinfo(i).name);
        fprintf('%4s\n',sinfo(i).ee);
    end
    % dataname
    sinfo(i).dataname  = cat_io_strrep(sinfo(i).posside,{sinfo(i).name,'template.','resampled.'},'');
    if ~isempty(sinfo(i).dataname) && sinfo(i).dataname(end)=='.', sinfo(i).dataname(end)=[]; end
    
    % ROI
    sinfo(i).roi = ~isempty(strfind(sinfo(i).posside,'.ROI'));
    
    
    
    % find Mesh and Data Files
    %  -----------------------------------------------------------------
    sinfo(i).Pmesh = '';
    sinfo(i).Pdata = '';
    % here we know that the gifti is a surf
    if sinfo(i).statready
        sinfo(i).Pmesh = sinfo(i).fname;
        sinfo(i).Pdata = sinfo(i).fname;
    end
    % if we have read the gifti than we can check for the fields
    if isempty(sinfo(i).Pmesh) && sinfo(i).exist && readsurf && isfield(S,'vertices')
        sinfo(i).Pmesh = sinfo(i).fname;
    end
    if isempty(sinfo(i).Pdata) && sinfo(i).exist && readsurf && isfield(S,'cdata')
        sinfo(i).Pdata = sinfo(i).fname;
    end
    % if the dataname is central we got a mesh or surf datafile
    if isempty(sinfo(i).Pdata) || isempty(sinfo(i).Pmesh)
        switch sinfo(i).texture
            case {'defects'} % surf
                sinfo(i).Pmesh = sinfo(i).fname;
                sinfo(i).Pdata = sinfo(i).fname;
            case {'central','inner','outer','sphere','hull'} % only mesh
                sinfo(i).Pmesh = sinfo(i).fname;
                sinfo(i).Pdata = '';
            case {'thickness','gyrification','frac','logsulc','GWMdepth','WMdepth','CSFdepth',...
                    'depthWM','depthGWM','depthCSF','depthWMg',...
                    'gyruswidth','gyruswidthWM','sulcuswidth'} % only thickness
                sinfo(i).Pdata = sinfo(i).fname;
        end
    end
    % if we still dont know what kind of datafile, we can try to find a
    % mesh surface
    if isempty(sinfo(i).Pmesh)
        if strcmp(ee,'.gii') && isempty(sinfo(i).side)
            sinfo(i).Pmesh = sinfo(i).fname;
            sinfo(i).Pdata = sinfo(i).fname;
        else
            % template mesh handling !!!
            Pmesh = char(cat_surf_rename(sinfo(i),'dataname','central','ee','.gii'));
            if exist(Pmesh,'file')
                sinfo(i).Pmesh = Pmesh;
                sinfo(i).Pdata = sinfo(i).fname;
            end
        end
    end
    % if we got still no mesh than we can use SPM.mat information or average mesh
    % ...
    if isempty(sinfo(i).Pmesh) %&& sinfo(i).ftype==1
        try
            if ischar(SPM.xVol.G)
                % data or analysis moved or data are on a different computer?
                if ~exist(SPM.xVol.G,'file')
                    [pp2,ff2,xx2] = icatb_spm_fileparts(SPM.xVol.G);
                    if strfind(ff2,'.central.freesurfer')
                        disp('FS')
                        if strfind(pp2,'templates_surfaces_32k')
                            SPM.xVol.G = fullfile(fileparts(which('gift')), 'icatb_templates', 'gifti_templates_surfaces_32k',[ff2 xx2]);
                        else
                            SPM.xVol.G = fullfile(fileparts(which('gift')), 'icatb_templates', 'gifti_templates_surfaces',[ff2 xx2]);
                        end
                    end
                end
                
                sinfo(i).Pmesh = SPM.xVol.G;
            else
                % 32k mesh?
                if SPM.xY.VY(1).dim(1) == 32492 || SPM.xY.VY(1).dim(1) == 64984
                    sinfo(i).Pmesh = fullfile(fileparts(which('gift')), 'icatb_templates', 'gifti_templates_surfaces_32k',...
                        [sinfo(i).side '.central.freesurfer.gii']);
                else
                    sinfo(i).Pmesh = fullfile(fileparts(which('gift')), 'icatb_templates','gifti_templates_surfaces',...
                        [sinfo(i).side '.central.freesurfer.gii']);
                end
            end
        catch
            % 32k mesh?
            switch sinfo(i).ee
                case '.gii'
                    if sinfo(i).exist && ~readsurf
                        S = gifti(P{i});
                    end
                case '.annot'
                    if sinfo(i).exist && ~readsurf
                        clear S;
                        try
                            S = cat_io_FreeSurfer('read_annotation',P{1});
                        end
                    end
            end
            
            if exist('S','var') && isfield(S,'cdata') && (length(S.cdata) == 32492 || length(S.cdata) == 64984)
                sinfo(i).Pmesh = fullfile(fileparts(which('gift')), 'icatb_templates','gifti_templates_surfaces_32k',...
                    [sinfo(i).side '.central.freesurfer.gii']);
            elseif exist('S','var') && isfloat(S) && (length(S) == 32492 || length(S) == 64984)
                sinfo(i).Pmesh = fullfile(fileparts(which('gift')), 'icatb_templates','gifti_templates_surfaces_32k',...
                    [sinfo(i).side '.central.freesurfer.gii']);
            else
                sinfo(i).Pmesh = fullfile(fileparts(which('gift')), 'icatb_templates', 'gifti_templates_surfaces',...
                    [sinfo(i).side '.central.freesurfer.gii']);
            end
        end
        sinfo(i).Pdata = sinfo(i).fname;
    end
    
    [ppm,ffm,eem]        = fileparts(sinfo(i).Pmesh);
    sinfo(i).Phull       = fullfile(ppm,strrep(strrep([ffm eem],'.central.','.hull.'),'.gii',''));
    sinfo(i).Psphere     = fullfile(ppm,strrep([ffm eem],'.central.','.sphere.'));
    sinfo(i).Pspherereg  = fullfile(ppm,strrep([ffm eem],'.central.','.sphere.reg.'));
    sinfo(i).Pdefects    = fullfile(ppm,strrep([ffm eem],'.central.','.defects.'));
    if ~exist(sinfo(i).Psphere ,'file'), sinfo(i).Psphere  = ''; end
    if ~exist(sinfo(i).Pdefects,'file'), sinfo(i).Pdefects = ''; end
    
    
    if sinfo(i).exist && readsurf
        if isfield(S,'vertices'),
            sinfo(i).nvertices = size(S.vertices,1);
        else
            if ~isempty(sinfo(i).Pmesh) && exist(sinfo(i).Pmesh,'file')
                S2 = gifti(sinfo(i).Pmesh);
                if ~isstruct(S), clear S; end
                if isfield(S2,'vertices'), S.vertices = S2.vertices; else S.vertices = []; end
                if isfield(S2,'faces'),    S.faces    = S2.faces;    else S.faces = []; end
            end
            if isfield(S,'vertices'),
                sinfo(i).nvertices = size(S.vertices,1);
            elseif isfield(S,'cdata'),
                sinfo(i).nvertices = size(S.cdata,1);
            else
                sinfo(i).nvertices = nan;
            end
        end
        if isfield(S,'faces'),    sinfo(i).nfaces    = size(S.faces,1); end
        if isfield(S,'cdata'),    sinfo(i).ncdata    = size(S.cdata,1); end
    end
    
    sinfo(i).catxml = fullfile(pp,['cat_' sinfo(i).name '*.xml']);
    if ~exist(sinfo(i).catxml,'file'), sinfo(i).catxml = ''; end
    
    if nargout>1
        varargout{2}{i} = S;
    else
        clear S
    end
end
varargout{1} = sinfo;





function MODIFIEDSTR = cat_io_strrep(ORIGSTR,OLDSUBSTR,NEWSUBSTR)
% _______________________________________________________________________
% cat_io_strrep replace strings by other strings. It based on the strrep 
% and allows to use cellstrings that were replace by another string or a
% similar number of cellstrings depending on their input order.
%
% claim{1} = 'This is a good example';
% claim{2} = 'This is a bad example';
% new_claimA = cat_io_strrep(claim,{' good',' bad'},'n')
% new_claimB = cat_io_strrep(claim,{'good','bad'},{'great','acceptable'})
%
% See also strrep, strfind, regexprep.
% _______________________________________________________________________
% Robert Dahnke
% $Id: cat_io_strrep.m 1017 2016-09-23 05:55:23Z dahnke $

  if nargin==0, help cat_io_strrep; return; end

  if iscell(ORIGSTR)
    MODIFIEDSTR = ORIGSTR; 
    for i=1:numel(ORIGSTR)
      MODIFIEDSTR{i} = cat_io_strrep(ORIGSTR{i},OLDSUBSTR,NEWSUBSTR);
    end
  else
    if iscell(OLDSUBSTR)
      if iscell(NEWSUBSTR)
        if numel(OLDSUBSTR)==numel(NEWSUBSTR)
          MODIFIEDSTR = ORIGSTR; 
          for i=1:numel(OLDSUBSTR) 
            MODIFIEDSTR = strrep(MODIFIEDSTR,OLDSUBSTR{i},NEWSUBSTR{i});
          end
        else
          error('cat_io_strrep:input',...
            'If multiple new stringswere used, their number must be equal to the number of old strings.\n'); 
        end
      else
        MODIFIEDSTR = ORIGSTR; 
        for i=1:numel(OLDSUBSTR) 
          MODIFIEDSTR = strrep(MODIFIEDSTR,OLDSUBSTR{i},NEWSUBSTR);
        end
      end
    else
      MODIFIEDSTR = strrep(ORIGSTR,OLDSUBSTR,NEWSUBSTR);
    end
  end

