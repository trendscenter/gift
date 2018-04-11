function fH = icatb_plot_connectogram(param_file, comp_network_names, varargin)
%% Functional network connectivity correlations are visualized as a connectogram plot.
%
% Inputs:
%
% 1. param_file - ICA parameter file (*ica*param*mat). By default, one sample t-test results of the FNC correlations are visualized. Also mean components are used as the component maps.
% If you want to pass the user defined correlation matrix pass after parameter 'C' and nifti file name after 'image_file_names' parameter.
% 2. comp_network_names - Network names and values are defined in a cell array of size number of networks by 2. First column corresponds to network names and second one corresponds to
% component numbers in the nifti file. You could optionally provide color values (Say [1, 1, 1] for white color) for each network in column 3.
% 3. Variable arguments are passed in pairs:
%   1. 'C' - Correlation matrix of size N x N where N is the number of
%   components defined in comp_network_names.
%   2. 'image_file_names' - Nifti file name containing spatial maps.
%   3. 'convert_to_zscores' - Convert images to z-scores. Options are 'yes' and 'no'.
%   4. 'image_values' - Image values. Options are 'positive and negative',
%   'postive', 'absolute value' and 'negative'.
%   5. 'template_file' - Anatomical file name.
%   6. 'slice_plane' - Anatomical plane to view.
%   7. 'threshold' - Threshold value.
%   8. 'colorbar_label' - Colorbar label.
%   9. cmap - Colormap of length 64 like hot(64).
%   10. imwidth - Image Width. You could customize the image width like 0.05
%   11. title - Figure title
%   12. display_type - By default slices are plotted. Option is provided to
%   use surface plots using 'render' option.
%   13. image_data - By default, image data is converted to RGB for each
%   component. You could provide custom RGB image in a cell array of length
%   equal to the number of components.
%   14. CLIM - Specify scaling like [-0.5, 0.5].
%   15. comp_labels - Override the default component labels string. Provide
%   a cell array of length equal to the number of components.
%   16. comp_txt_color - Component labels text color
%   17. conn_threshold - Connectivity threshold. Only values surpassing the threshold are used.
%   18. plot_arrows - Options are 0 and 1. Useful when plotting lag.
%   19. exclude_zeros - By default, only non-zero values are reported. You
%   could override this default using a value of 0.
%   20. fig_position - Figure position in a vector of length equal to 4 [x,
%   y, width, height] in pixels.
%   21. line_width - Connecting lines with.
%   22. radius_sm - Radius of spatial maps. Default is set to 1.7.
%
%
% Example is below:
%
%
% comp_network_names = {'BG', 21;                    % Basal ganglia 21st component
%     'AUD', 17;                   % Auditory 17th component
%     'SM', [7 23 24 29 38 56];    % Sensorimotor comps
%     'VIS', [39 46 48 59 64 67];  % Visual comps
%     'DMN', [25 50 53 68];        % DMN comps
%     'ATTN', [34 52 55 60 71 72]; % ATTN Comps
%     'FRONT', [20 42 47 49]};     % Frontal comps
%
% C = icatb_corr(rand(1000, 28)); C = C - eye(size(C)); C(abs(C) < 0.01) = 0;
% fname = 'F:\Example Subjects\mancova_sample_data\ica_output\rest_hcp_mean_component_ica_s_all_.nii';
% a. icatb_plot_connectogram([], comp_network_names, 'C', C, 'threshold', 1.5, 'image_file_names', fname, 'colorbar_label', 'Corr');
% b. icatb_plot_connectogram([], comp_network_names, 'C', C, 'threshold', 1.5, 'image_file_names', fname, 'colorbar_label', 'Corr', 'display_type', 'render');
%
% For better resolution, you could print to postscript file followed by
% postscript to pdf using adobe distiller or ghostscript.
% print -r300 -dpsc test.ps -bestfit


icatb_defaults;
global FONT_COLOR;
global CONNECTOGRAM_OPTIONS;
global CONNECTOGRAM_SM_WIDTH;
global DETRENDNUMBER;
global UI_FONTNAME;
global UI_FS;

%% Parse params
threshold = 1.5;
convert_to_zscores = 'yes';
image_values = 'positive';
load icatb_colors coldhot;
cmap = coldhot(1:4:end, :);
template_file = fullfile(fileparts(which('groupica.m')), 'icatb_templates', 'ch2bet.nii');
slice_plane = 'sagittal';
colorbar_label = 'Corr';

imWidth = CONNECTOGRAM_SM_WIDTH;
try
    imWidth = CONNECTOGRAM_OPTIONS.sm_width;
catch
end

titleStr = 'Connectogram';
comp_txt_color = [0, 0, 0];
display_type = 'slices';
conn_threshold = -Inf;
order_mtx = 0;
plotAnnotation = 0;
exclude_zeros = 1;
fig_position = [];
line_width = 2;
radius_sm = [];
try
    radius_sm = CONNECTOGRAM_OPTIONS.radius_sm;
catch
end

for nF = 1:2:length(varargin)
    if (strcmpi(varargin{nF}, 'C'))
        C = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'image_file_names'))
        fileNames = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'convert_to_zscores') || strcmpi(varargin{nF}, 'convert_to_z'))
        convert_to_zscores = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'image_values'))
        image_values = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'cmap'))
        cmap = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'template_file'))
        template_file = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'slice_plane'))
        slice_plane = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'threshold'))
        threshold = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'colorbar_label'))
        colorbar_label = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'imwidth'))
        imWidth = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'title'))
        titleStr = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'image_data'))
        I = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'CLIM'))
        FNC_CLIM = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'comp_labels'))
        compStr = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'comp_txt_color'))
        comp_txt_color = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'display_type'))
        display_type = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'conn_threshold'))
        conn_threshold = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'order_matrix'))
        % Order correlation matrix
        order_mtx = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'plot_arrows'))
        % plot arrows (useful in showing lag)
        plotAnnotation = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'exclude_zeros'))
        % exclude zeros in the connectivity matrix
        exclude_zeros = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'fig_position'))
        % figure position
        fig_position = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'line_width'))
        % Line width
        line_width = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'radius_sm'))
        % Radius for plotting spatial maps
        radius_sm = varargin{nF + 1};
    end
end


if (nargin == 0)
    %% Open GUI
    connectogramGUI;
    return;
end


disp('Displaying connectogram ...');

if (~exist('param_file', 'var'))
    param_file = [];
end


if (~isempty(param_file))
    load(param_file);
    sesInfo.outputDir = fileparts(param_file);
    sesInfo.userInput.pwd = sesInfo.outputDir ;
end

if (isempty(conn_threshold))
    conn_threshold = -Inf;
end


comps = [comp_network_names{:, 2}];
comps = comps(:)';

%% Check total number of components
if (~exist('I', 'var'))
    if (~exist('fileNames', 'var'))
        fileNames = fullfile(sesInfo.outputDir, sesInfo.icaOutputFiles(1).ses(1).name);
    end
    fileNames = cellstr(icatb_rename_4d_file(char(fileNames)));
    totalComps = length(fileNames);
else
    totalComps = length(I);
end


%% Compute t-test if correlation pairs not present
if (~exist('C', 'var'))
    
    %comps = [comp_network_names{:, 2}];
    %comps = comps(:)';
    
    nDataSet = 0;
    for nSub = 1:sesInfo.numOfSub
        for nSess = 1:sesInfo.numOfSess
            nDataSet = nDataSet + 1;
            TC = icatb_loadComp(sesInfo, comps, 'vars_to_load', 'tc', 'detrend_no', DETRENDNUMBER, 'subjects', nSub, 'sessions', nSess);
            cvals = icatb_corr(squeeze(TC));
            cpairs = icatb_r_to_z(icatb_mat2vec(cvals));
            if (nDataSet == 1)
                Ct = zeros(sesInfo.numOfSess, sesInfo.numOfSub, length(cpairs));
            end
            Ct(nSess, nSub, :) = cpairs;
        end
    end
    
    if (sesInfo.numOfSess == 1)
        Ct = squeeze(Ct);
    else
        Ct = squeeze(mean(Ct));
    end
    
    if (sesInfo.numOfSub > 1)
        [pi, tstat] = icatb_ttest(Ct);
    else
        tstat = squeeze(Ct);
    end
    C = icatb_vec2mat(tstat);
    
    order_mtx = 0;
    
else
    if (numel(C) == length(C))
        % Convert vector to matrix
        C = icatb_vec2mat(C);
    end
    
end

%% Error check for dimensions of Connectivity matrix
if ((size(C, 1) ~= length(comps)) && (size(C, 1) ~= totalComps))
    error('FNC matrix dimensions does not match the number of components entered');
end

if (size(C, 1) == totalComps)
    if ((length(comps) ~= totalComps) || order_mtx)
        % order matrix: Handle special case when all the components are selected for
        % connectogram plot (avoid reordering twice)
        C = C(comps, comps);
    end
end


% Change diagonals to zeros
C = C - diag(diag(C));

C(abs(C) < conn_threshold) = 0;

C(C==0) = NaN;

% Default Colors
frame_colors = getFrameColor;

for nC = 1:size(comp_network_names, 1)
    tmp = [];
    try
        tmp = comp_network_names{nC, 3};
    catch
    end
    if (isempty(tmp))
        if (nC <= size(frame_colors, 1))
            %tmp = colordg(nC);
            tmp = frame_colors(nC,:);
        else
            tmp = rand(1, 3);
        end
    end
    comp_network_names{nC, 3} = tmp;
    clear tmp;
end

if (exclude_zeros)
    %% Remove zeros in the connectivity matrix
    [C, comp_network_names, inc_comps] = removeZeros(C, comp_network_names);
else
    % plot all
    inc_comps = (1:size(C, 1));
    %[C, comp_network_names, inc_comps] = removeZeros(C, comp_network_names);
end

%% Important: Update components and its labels after removing zeros in the connectivity matrix
comps = [comp_network_names{:, 2}];
comps = comps(:)';

colorlims = max(abs(C(:)));
colorlims = [-colorlims, colorlims];

newC = cell(1, length(comps));

e = 0;
for nC = 1:size(comp_network_names, 1)
    s = e + 1;
    e = e + length(comp_network_names{nC, 2});
    tmp = comp_network_names{nC, 3};
    newC(s:e) = {tmp}; %color_values(nC);
end

clear tmp;

if (~exist('compStr', 'var'))
    prefix = '   IC  ';
    compStr = cell(length(comps), 1);
    for nComps = 1:length(comps)
        compStr{nComps} = charC(prefix, num2str(comps(nComps)));
    end
else
    compStr = cellstr(compStr);
    % Update user defined labels
    compStr = compStr(inc_comps);
end


if (~exist('I', 'var'))
    
    I = cell(1, length(comps));
    
    structVol = icatb_spm_vol(template_file);
    structVol = structVol(1);
    structDIM = structVol.dim(1:3);
    
    returnValue = strmatch(lower(image_values), {'positive and negative', 'positive', 'absolute value', 'negative'}, 'exact');
    
    if (strcmpi(display_type, 'slices'))
        
        imagesCmap = icatb_getColormap(1, returnValue, 1);
        
        for nComp = 1:length(comps)
            
            fileName = fileNames{comps(nComp)};
            
            [tmpF, tmpN] = icatb_parseExtn(fileName);
            
            imageTmp = icatb_resizeImage(structVol, tmpF, 'axial', [], tmpN);
            
            structData = reshape(imageTmp(1, :, :, :), structVol.dim(1:3));
            imageTmp = reshape(imageTmp(end, :, :, :), structVol.dim(1:3));
            
            imageTmp = icatb_applyDispParameters_comp(imageTmp(:)', strcmpi(convert_to_zscores, 'yes'), returnValue, threshold);
            imageTmp = reshape(imageTmp, structDIM);
            
            [dd, tmp_inds] = max(abs(imageTmp(:)));
            [pixX, pixY, pixZ] = ind2sub(structDIM, tmp_inds);
            
            if (strcmpi(slice_plane, 'sagittal'))
                cdata = (rot90(squeeze(imageTmp(pixX,:,:))));
                structData = (rot90(squeeze(structData(pixX,:,:))));
            elseif (strcmpi(slice_plane, 'coronal'))
                cdata = (rot90(squeeze(imageTmp(:,pixY,:))));
                structData = (rot90(squeeze(structData(:,pixY,:))));
            else
                cdata = (rot90(squeeze(imageTmp(:,:,pixZ))));
                structData = (rot90(squeeze(structData(:,:,pixZ))));
            end
            
            cdata = makeCompositeMap(structData, cdata, returnValue);
            
            cdata = returnRGB(cdata, imagesCmap);
            
            I{nComp} = cdata;
        end
        
    else
        
        % render
        for nComp = 1:length(comps)
            fileName = fileNames{comps(nComp)};
            hOut = icatb_display_composite(fileName, 'image_values', image_values, 'convert_to_z', convert_to_zscores, 'threshold', threshold, 'display_type', 'render');
            [ddOut, indsOut] = max(hOut.countVoxelsRend);
            I{nComp} = hOut.slices_all{indsOut(1)};
        end
        
        
    end
    
else
    
    % Use only included components
    I = I(inc_comps);
    
    % If image data is provided in RGB
    if (length(I) ~= length(comps))
        error('Image data must be in a cell array of length equal to the number of components');
    end
    
end


if (exist('FNC_CLIM', 'var') && ~isempty(FNC_CLIM))
    [C, colorlims]  = checkCmap(C, FNC_CLIM);
end


if (isempty(radius_sm))
    if (size(C, 1) < 10)
        radius_sm = 2;
    else
        radius_sm = 1.7;
    end
end


[h, fH] = plot_connecto_gram(C, I, cmap, compStr, newC, comp_network_names, imWidth, colorlims, comp_txt_color, plotAnnotation, fig_position, line_width, radius_sm);
set(fH, 'name', titleStr);


titleColor = [0 0.9 0.9];

% fonts
titleFont = 11;

axisH = axes('Parent', fH, 'position', [0 0 1 1], 'visible', 'off');

text(0.5, 0.97, titleStr, 'color', titleColor, 'FontAngle', 'italic', 'fontweight', 'bold', 'fontsize', titleFont, 'HorizontalAlignment', ...
    'center', 'FontName', UI_FONTNAME, 'parent', axisH);

xlims = [1, size(cmap, 1)];

colormap(cmap);
set(gca, 'clim', [1, size(cmap, 1)]);

cbWidth = 0.12;
cbHeight = 0.025;
ch = colorbar('horiz');
cpos = get(ch, 'position');
cpos(1) = 0.75 - 0.5*cbWidth;
cpos(2) = 0.062;
cpos(3) = cbWidth;
cpos(4) = cbHeight;
set(ch, 'position', cpos);
%set(ch, 'xlim', xlims);
%xlims = get(ch,'xlim');
set(ch, 'xtick', [xlims(1), xlims(end)]);
set(ch, 'xticklabel', num2str(colorlims','%0.1f'));
xlabel(colorbar_label, 'parent', ch);
set(ch, 'XCOLOR', FONT_COLOR);
set(ch, 'YCOLOR', FONT_COLOR);
set(ch, 'fontsize', UI_FS - 2);

handles_data = get(fH, 'userdata');
%uistack(handles_data.PH, 'top');
cns = comp_network_names(:,1);
lH = legend(handles_data.PH, cns{:}, 'location', 'eastoutside');
lpos = get(lH, 'position');
lpos(1) = lpos(1)+ 0.045;
set(lH, 'position', lpos);
set(lH, 'fontsize', UI_FS - 1);


fprintf('Done\n\n');


function [h, fH] = plot_connecto_gram(C, I, cmap_corr, compNames, colorFrames, comp_network_names, imWidth, FNC_CLIM, comp_txt_color, plotAnnotation, fig_position, line_width, R)
% Plot connectogram

icatb_defaults;
global UI_FS;
global CONNECTOGRAM_FIG_POS;
global CONNECTOGRAM_OPTIONS;

if (~exist('fig_position', 'var'))
    fig_position = CONNECTOGRAM_FIG_POS;
    try
        fig_position = CONNECTOGRAM_OPTIONS.fig_pos;
    catch
    end
end


% C = C./(max(abs(C(:))) + eps);
C(C == 0) = NaN;


if (~exist('colorFrames', 'var'))
    colorFrames = repmat({[1, 1, 0]}, length(compNames), 1);
end

if (~exist('cmap_corr', 'var'))
    cmap_corr = hsv(64);
end

if (~exist('imWidth', 'var'))
    imWidth = [];
end

if (~exist('plotAnnotation', 'var'))
    plotAnnotation = 0;
end

numSlices = size(C, 2);
fH = icatb_getGraphics('Fig', 'graphics', 'gg', 'on');
fig_color = get(fH, 'color');
set(fH, 'resize', 'on');
%set(fH,'menubar','none');
sz = get(0, 'ScreenSize');
if (isempty(fig_position))
    figPos = [50, 50, sz(3) - 100, sz(4) - 100];
else
    figPos = fig_position;
end
try
    set(fH, 'position', figPos);
catch
end
if (numSlices < 10)
    axesWidth = 0.4;
else
    axesWidth = 0.62;
end
axesPos = [0.5 - 0.5*axesWidth, 0.5 - 0.5*axesWidth, axesWidth, axesWidth];
aH = axes('parent', fH, 'units', 'normalized', 'position', axesPos, 'NextPlot','add', 'tag', 'main_axes');
set(aH, 'color', fig_color);
pi_incr = 2*pi/numSlices;
rin = 1;
rout = 1.3;
% if (numSlices < 10)
%     rout = 1.3;
% else
%     rout = 1.5;
% end
startAngle = -pi/2;

mpX = zeros(1, numSlices);
mpY = zeros(1, numSlices);
th = zeros(1, numSlices);
for n = 1:numSlices
    tmpC = colorFrames{n};
    endAngle = startAngle + pi_incr;
    t2 = linspace(startAngle, endAngle, 256);
    th(n) = t2(ceil(length(t2)/2));
    startAngle = endAngle;
    xin = rin*cos(t2(end:-1:1));
    yin = rin*sin(t2(end:-1:1));
    xout = rout*cos(t2);
    yout = rout*sin(t2);
    
    PH(n) = patch([xout, xin],[yout, yin], tmpC, 'edgecolor', 'k', 'parent', aH);
    xx = (rin + rout)*cos(th(n))/2;
    yy = (rin + rout)*sin(th(n))/2;
    %trot(n) = 180/pi* (atan((newY(2) - newY(1)) / (newX(2) - newX(1))));
    %thA=text(xx,yy,compNames{n},'horizontalalignment','center','rotation',(180/pi)*th(n),'fontsize',11,'fontweight','bold', 'color', 'k', 'tag', ['text_', num2str(n)]);
    mpX(n) = xx;
    mpY(n) = yy;
    hold(aH, 'on');
end

mpX(end + 1) = mpX(1);
mpY(end + 1) = mpY(1);
axis (aH, 'equal');
axis(aH, 'off');



phInc = [];
e = 0;
for n = 1:size(comp_network_names, 1)
    nIn = length(comp_network_names{n, 2});
    s = e + ceil(nIn/2);
    e = e + nIn;
    %inc = 0;
    %textDeg = 180/pi* (atan((mpY(s + 1) - mpY(s)) / ( mpX(s + 1) - mpX(s) )));
    %text(xtext(s) + inc, ytext(s) + inc, comp_network_names{n, 1}, 'parent', aH, 'color', colorFrames{s}, 'fontsize', UI_FS + 2, 'fontweight', 'bold', 'rotation', ...
    %   textDeg, 'tag', ['Label_', num2str(n)]);
    phInc= [phInc, s];
end


PH = PH(phInc);


% Set figure position
try
    set(fH, 'position', figPos);
catch
end

% Make sure that axes x origin is visible
chkAxesPos = get(aH, 'position');
chkAxesPos(1) = 0.5 - 0.5*chkAxesPos(3);
set(aH, 'position', chkAxesPos);

textHandles = [];
th_deg = th/pi*180;
for n = 1:numSlices
    %textDeg = 180/pi* (atan((mpY(n + 1) - mpY(n)) / ( mpX(n + 1) - mpX(n) )));
    % Rotation tweak (Ashkan)
    textDeg = th_deg(n) + 90;
    textDeg = rem(textDeg + 3600, 360);
    if (textDeg > 80 && textDeg < 270)
        textDeg = textDeg + 180;
    end
    thA = text(mpX(n), mpY(n), compNames{n}, 'parent', aH, 'horizontalalignment','center','rotation',textDeg,'fontsize', UI_FS - 3, 'fontweight', 'bold', 'color', comp_txt_color, ...
        'tag', ['text_', num2str(n)], 'verticalalignment','middle');
    textHandles = [textHandles, thA];
end

handles_data.textHandles = textHandles;
handles_data.PH = PH;
set(fH, 'userdata', handles_data);

if (numSlices < 10)
    RText = 2.8;
else
    RText = 2;
end
[xtext,ytext]=pol2cart(th, RText);

if (~exist('R', 'var'))
    if (numSlices < 10)
        R = 2;
    else
        R = 1.7;
    end
end

[Xna, Yna]= pol2cart(th, R);

[x0, y0] = getAxesPos(Xna, Yna, aH);

rL = axesPos(4)/(2*rout);

rL = max(abs([rL*cos(2*pi/numSlices), rL*sin(2*pi/numSlices)]));

if (isempty(imWidth))
    imWidth = (2*pi*rL/numSlices);
    imWidth = min([imWidth, 0.16]);
    
    if (imWidth < 0.04)
        imWidth = 0.04;
    end
end

imHeight = imWidth;
offset = 0.5*imWidth;

ahs = zeros(1, length(Xna));
for n = 1:length(Xna)
    xr = -offset;
    yr = -offset;
    ah = axes('parent', fH, 'units', 'normalized', 'position',[x0(n) + xr, y0(n) + yr, imWidth, imHeight]);
    image(I{n}, 'parent', ah);
    axis(ah, 'image');
    axis(ah, 'off');
    ahs(n) = ah;
end

labels = cellstr(strcat('', num2str((1:size(C,1))')));
[h,X,Y,cmap_corr] = schemaball(C, cmap_corr, aH, th, FNC_CLIM, plotAnnotation, line_width);



function [mx, my] = midpoint(x, y)

newX = [x(end), x];
newY = [y(end), y];

my = zeros(1, length(x));
mx = my;

for n = 2:length(newX)
    
    mx(n - 1) = mean([newX(n - 1), newX(n)]);
    my(n - 1) = mean([newY(n - 1), newY(n)]);
end



function [h,x,y,ccolor] = schemaball(r, ccolor, gH, theta, FNC_CLIM, plotAnnotation, line_width)
% Adapted version of schemaball https://www.mathworks.com/matlabcentral/fileexchange/42279-okomarov-schemaball?requestedDomain=www.mathworks.com
%

% Number of color shades/buckets (large N simply creates many perceptually indifferent color shades)
N = size(ccolor, 1)/2;

sz = size(r);

% Points in [0, 1] for bezier curves: leave space at the extremes to detach a bit the nodes.
% Smaller step will use more points to plot the curves.
t      = (0.025: 0.05 :1)';

% Create figure
axes(gH);

% Index only low triangular matrix without main diag
tf        = tril(true(sz),-1);

% Index correlations into bucketed colormap to determine plotting order (darkest to brightest)
N2        = 2*N;
%[n, isrt] = histc(r(tf), linspace(-1,1 + eps(100),N2 + 1));
rtf = r(tf);
[n, isrt] = histc(rtf, linspace(min(FNC_CLIM), max(FNC_CLIM) + eps(100), N2 + 1));
plotorder = reshape([N:-1:1; N+1:N2],N2,1);

% Retrieve pairings of nodes
[row,col] = find(tf);

% Positions of nodes on the circle starting from (0,-1), useful later for label orientation
%step  = tau/sz(1);
%theta = -.25*tau : step : .75*tau - step;
% Get cartesian x-y coordinates of the nodes
x     = cos(theta);
y     = sin(theta);

% PLOT BEZIER CURVES
% Calculate Bx and By positions of quadratic Bezier curves with P1 at (0,0)
% B(t) = (1-t)^2*P0 + t^2*P2 where t is a vector of points in [0, 1] and determines, i.e.
% how many points are used for each curve, and P0-P2 is the node pair with (x,y) coordinates.
t2  = [1-t, t].^2;
s.l = NaN(N2,1);
%plotAnnnotation = 1;
% LOOP per color bucket
for c = 1:N2
    pos = plotorder(c);
    idx = isrt == pos;
    if nnz(idx)
        Bx     = [t2*[x(col(idx)); x(row(idx))]; NaN(1,n(pos))];
        By     = [t2*[y(col(idx)); y(row(idx))]; NaN(1,n(pos))];
        xx = Bx(:);
        yy = By(:);
        [x0, y0] = getAxesPos(xx, yy, gH);
        s.l(c) = plot(Bx(:), By(:), 'Color', ccolor(pos, :), 'LineWidth', line_width, 'parent', gH);
        
        if (plotAnnotation)
            
            [newX, newY]= getBezXY(x0, y0);
            
            %chk = find(isnan(xx)==1);
            tmpcvals = rtf(idx);
            for na = 1:length(newX)
                tmpx = newX{na};
                tmpy = newY{na};
                hold on;
                if (sign(tmpcvals(na)) == 1)
                    annotation('arrow', [tmpx(end-1), tmpx(end)], [tmpy(end-1), tmpy(end)], 'color', ccolor(pos, :));
                else
                    annotation('arrow', [tmpx(2), tmpx(1)], [tmpy(2), tmpy(1)], 'color', ccolor(pos, :));
                end
            end
        end
        
    end
end

h = s;


function scaledData = scaleIm(tmin, tmax, imageData, returnValue, data_range)
%% Scale images
%

if (~exist('returnValue', 'var'))
    returnValue = 2;
end

if (~exist('data_range', 'var'))
    minVal = min(imageData);
    maxVal = max(imageData);
else
    minVal = min(data_range);
    maxVal = max(data_range);
end

if (returnValue == 1)
    maxVal = max(abs([minVal, maxVal]));
    minVal = -maxVal;
end

rangeVal = (maxVal-minVal) + eps;
trange = tmax-tmin;
if (rangeVal == 0)
    rangeVal = eps;
end
scaledData = (((imageData-minVal)./rangeVal)./(1/trange))+tmin;


function compositeMap = makeCompositeMap(structData, imageTmp, returnValue)
%% Make composite map
%
mask_inds = (abs(imageTmp) > eps);

compositeMap = structData(:);
structInds = (abs(structData) > eps);
compositeMap(structInds) = scaleIm(65, 128, structData(structInds));
compositeMap(mask_inds) = scaleIm(1, 64, imageTmp(mask_inds), returnValue);

compositeMap = reshape(compositeMap, size(structData));

function RGB = returnRGB(compositeMap, cmap)
%% Get RGB values
%

dims = [size(compositeMap, 1), size(compositeMap, 2), 3];
compositeMap = ceil(compositeMap);
compositeMap(compositeMap == 0) = size(cmap, 1) + 1;
cmap(end + 1, :) = [0, 0, 0];
RGB = reshape(cmap(compositeMap(:), :), dims);



function [x0, y0] = getAxesPos(Xna, Yna, aH)
%% Get axes position in figure units
%

xlim = get(aH, 'xlim');
ylim = get(aH, 'ylim');
pos = get(aH, 'position');

mn = min(xlim);
rn = range(xlim);
x0 = pos(1) + (pos(3))*(Xna - mn)./rn;

mn = min(ylim);
rn = range(ylim);
y0 = pos(2) + (pos(4) )*(Yna - mn)./rn;


function rangeVal = range(imageData)
%% Get range

minVal = min(imageData(:));
maxVal = max(imageData(:));
rangeVal = maxVal-minVal;



function [C, comp_network_names, inc] = removeZeros(C, comp_network_names)
%% Cleanup the correlation matrix
%

C(eye(size(C)) == 1) = 0;
C(isfinite(C)==0) = 0;
rows = (1:size(C, 1));
inc = [];
for n = 1:size(C, 1)
    tmp = C(n, :);
    if length(find(abs(tmp) > eps)) > 0
        inc = [inc, n];
    end
end

C = C(inc, inc);

if (isempty(C))
    error('Correlation matrix has all zeros');
end

Nrows = size(comp_network_names, 1);
RowsToexclude = [];
endLen = 0;
for nC = 1:Nrows
    startLen = endLen + 1;
    tmpC = comp_network_names{nC, 2};
    endLen = endLen + length(tmpC);
    compN = (startLen:endLen); %comp_network_names{nC, 2};
    [chk, ia] = intersect(compN, inc);
    if (isempty(chk))
        comp_network_names{nC, 2} = [];
    else
        ia = sort(ia);
        comp_network_names{nC, 2} = tmpC(ia);
    end
    %if (isempty(chk))
    %RowsToexclude = [RowsToexclude, nC];
    %end
end

chk = cellfun('isempty', comp_network_names(:,2));
comp_network_names(chk, :) = [];

function [C, colorlims] = checkCmap(C, FNC_CLIM)
% Apply scaling

%maxAbsVal = max(abs(FNC_CLIM(:)));
maxVal = max(FNC_CLIM(:));
minVal = min(FNC_CLIM(:));

c_inds = (abs(C) > eps);
CN = C(c_inds);
CN(CN > maxVal) = maxVal;
CN(CN < minVal) = minVal;

C(c_inds) = CN;

colorlims = [minVal, maxVal];


function c = charC(a, b)
%% Character array. Strings are centered.
%

len = max([length(a), length(b)]);

c = repmat(' ', 2, len);

s = ceil(len/2 - length(a)/2) + 1;
c(1, s:s+length(a) - 1) = a;
s = ceil(len/2 - length(b)/2) + 1;
c(2, s:s + length(b) - 1) = b;


function connectogramGUI
%% Connectogram GUI
%

icatb_defaults;
global UI_FS;
global CONVERT_Z;
global THRESHOLD_VALUE;
global IMAGE_VALUES;
global PARAMETER_INFO_MAT_FILE;



%% Select component networks
filterP = ['*', PARAMETER_INFO_MAT_FILE, '*.mat'];
file_names = icatb_selectEntry('typeSelection', 'multiple', 'typeEntity', 'file', 'filter', [filterP, ';*comp*.img;*comp*.nii'], 'title', 'Select ICA parameter file/ component spatial maps', 'filetype', 'image');
drawnow;

if (isempty(file_names))
    error('Files are not selected');
end


[pathstr, fN, extn] = fileparts(deblank(file_names(1, :)));
isParamFile = 0;
if (~strcmpi(extn, '.mat'))
    file_names = icatb_rename_4d_file(file_names);
    %% Load fnc matrix file
    fnc_matrix_file = icatb_selectEntry('typeSelection', 'single', 'typeEntity', 'file', 'filter', '*.mat;*.txt;*.dat', 'title', 'Select FNC correlations file (2D matrix: comps x comps)');
    FNCM = loadFNC(fnc_matrix_file);
    numComp = size(FNCM, 2);
else
    param_file = deblank(file_names(1, :));
    load(param_file);
    if (~exist('sesInfo', 'var'))
        error('Selected file is not a valid ICA parameter file');
    end
    numComp = sesInfo.numComp;
    isParamFile = 1;
    file_names = fullfile(pathstr, sesInfo.icaOutputFiles(1).ses(1).name);
    file_names = icatb_rename_4d_file(file_names);
end

MenuOptions = {'image_values', IMAGE_VALUES, 'convert_to_z', CONVERT_Z, 'threshold', str2num(THRESHOLD_VALUE), 'slice_plane', 'sagittal', 'imWidth', [], 'cmap', 'coldhot', ...
    'CLIM', []};


compGroupNames = '';

figData.Options = MenuOptions;
figData.numComps = numComp;
figData.numICs = size(file_names, 1); %% Draw graphics
figData.comp = [];
figureTag = 'setup_comp_networks_gui';
figHandle = findobj('tag', figureTag);
if (~isempty(figHandle))
    delete(figHandle);
end

% Setup figure for GUI
InputHandle = icatb_getGraphics('Connectogram GUI', 'displaygui', figureTag, 'on');
set(InputHandle, 'menubar', 'none');
set(InputHandle, 'userdata', figData);
optionsMenuH = uimenu('parent', InputHandle, 'label', 'Defaults', 'callback', {@options_callback, InputHandle});

promptWidth = 0.52; controlWidth = 0.32; promptHeight = 0.05;
xOffset = 0.02; yOffset = promptHeight; listboxHeight = 0.18; listboxWidth = controlWidth; yPos = 0.9;
okWidth = 0.12; okHeight = promptHeight;

promptPos = [xOffset, yPos - 0.5*yOffset, promptWidth, promptHeight];

compNetworkNameData = getCompData(file_names);
drawnow;



listboxXOrigin = promptPos(1) + promptPos(3) + xOffset;
listboxYOrigin = promptPos(2) - 0.5*listboxHeight;

%%  Components listbox (Group components by name)
%promptPos(2) = listboxYOrigin - yOffset - 0.5*listboxHeight;
textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Add Components', 'tag', ...
    'prompt_components', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);
listboxYOrigin = promptPos(2) - 0.5*listboxHeight;
listboxPos = [listboxXOrigin, listboxYOrigin, listboxWidth, listboxHeight];
listCompH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', listboxPos, 'string', compGroupNames, 'tag', ...
    'comp', 'fontsize', UI_FS - 1, 'min', 0, 'max', 2, 'value', 1, 'callback', {@addCompNetwork, InputHandle}, 'userdata', compNetworkNameData);

addButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) + 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
removeButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) - 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', addButtonPos, 'string', '+', 'tag', 'add_comp_button', 'fontsize',...
    UI_FS - 1, 'callback', {@addCompNetwork, InputHandle});
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', removeButtonPos, 'string', '-', 'tag', 'remove_comp_button', 'fontsize',...
    UI_FS - 1, 'callback', {@removeCompNetwork, InputHandle});


promptPos(2) = listboxYOrigin - yOffset - 0.5*promptHeight;
promptPos(3) = promptWidth;
textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter FNC threshold', 'tag', 'prompt_fnc_thresh', ...
    'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', '', 'tag', 'conn_threshold', 'fontsize', UI_FS - 1, 'tooltipstring', ...
    'Values less than the connectivity thresholded are not included in the plot ...');


promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
promptPos(3) = promptWidth;
textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter FNC label (in Colorbar)', 'tag', 'prompt_fnc_label', ...
    'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', 'T-value', 'tag', 'fnc_label', 'fontsize', UI_FS - 1);


promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
promptPos(3) = promptWidth;
textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Select display type', 'tag', 'prompt_display_type', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
popupH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'popup', 'position', editPos, 'string', {'Slices', 'Render'}, 'tag', 'display_type', ...
    'fontsize', UI_FS - 1, 'value', 1, 'tooltipstring', 'You could select either surface or slices to be plotted in a circle...');

anatWidth = 0.2;
anatHeight = 0.05;
anatPos = [0.25 - 0.5*okWidth, promptPos(2) - yOffset - 0.5*okHeight, anatWidth, anatHeight];
anatPos(2) = anatPos(2) - 0.5*anatPos(4);
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', anatPos, 'string', 'Select Anatomical', 'tag', 'anat_button', 'fontsize',...
    UI_FS - 1, 'callback', {@selectAnatomical, InputHandle});


%% Add cancel, save and run buttons
okPos = [0.75 - 0.5*okWidth, promptPos(2) - yOffset - 0.5*okHeight, okWidth, okHeight];
okPos(2) = okPos(2) - 0.5*okPos(4);
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', okPos, 'string', 'Ok', 'tag', 'run_button', 'fontsize',...
    UI_FS - 1, 'callback', {@runCallback, InputHandle});

waitfor(InputHandle);


% Plot connectogram
appName = 'inputConnGramData';
if isappdata(0, appName)
    % get the application data
    answers = getappdata(0, appName);
    rmappdata(0, appName); % remove the application data
    comp_network_names = answers.comp_network_names;
    threshold = answers.threshold;
    image_values = answers.image_values;
    convert_to_z = answers.convert_to_z;
    slice_plane = answers.slice_plane;
    colorbar_label = answers.colorbar_label;
    imWidth = answers.imWidth;
    CLIM = answers.CLIM;
    display_type = answers.display_type;
    conn_threshold = answers.conn_threshold;
    cmap = answers.cmapN;
    
    structFile = fullfile(fileparts(which('groupica.m')), 'icatb_templates', 'ch2bet.nii');
    try
        structFile = answers.structFile;
    catch
    end
    
    drawnow;
    
    % params
    if (~isParamFile)
        connectogram_params = {[], comp_network_names, 'C', FNCM, 'image_file_names', file_names, 'convert_to_z', convert_to_z, 'image_values', image_values, ...
            'threshold', threshold, 'cmap', cmap, 'template_file', structFile, 'slice_plane', slice_plane, 'colorbar_label', colorbar_label, 'imwidth', imWidth, ...
            'CLIM', CLIM, 'display_type', display_type, 'conn_threshold', conn_threshold, 'order_matrix', 1};
    else
        connectogram_params = {param_file, comp_network_names, 'convert_to_z', convert_to_z, 'image_values', image_values, ...
            'threshold', threshold, 'cmap', cmap, 'template_file', structFile, 'slice_plane', slice_plane, 'colorbar_label', colorbar_label, 'imwidth', imWidth, ...
            'CLIM', CLIM, 'display_type', display_type, 'conn_threshold', conn_threshold, 'order_matrix', 0};
    end
    
    icatb_plot_connectogram(connectogram_params{:});
    
end



function compNetworkNameData = getCompData(file_names)

structFile = deblank(file_names(1, :));

%% Get colormap associated with the image values
structData2 =  icatb_spm_read_vols(icatb_spm_vol(structFile));
structData2(isfinite(structData2) == 0) = 0;
structDIM = [size(structData2, 1), size(structData2, 2), 1];

for nC = 1:size(file_names, 1)
    
    tmp = icatb_spm_read_vols(icatb_spm_vol(file_names(nC, :)));
    tmp(isfinite(tmp)==0) = 0;
    
    tmp(tmp ~= 0) = detrend(tmp(tmp ~= 0), 0) ./ std(tmp(tmp ~= 0));
    tmp(abs(tmp) < 1.0) = 0;
    
    if (nC == 1)
        compData = zeros(size(tmp, 1), size(tmp, 2), size(file_names, 1));
    end
    
    [dd, inds] = max(tmp(:));
    
    [x, y, z] = ind2sub(size(tmp), inds);
    
    [tmp, maxICAIM, minICAIM, minInterval, maxInterval] = icatb_overlayImages(reshape(tmp(:, :, z), [1, size(tmp, 1), size(tmp, 2), 1]), structData2(:, :, z), structDIM, structDIM, 1);
    
    compData(:, :, nC) = reshape(tmp, structDIM);
    
end


clim = [minInterval, 2*maxInterval];
cmap = icatb_getColormap(1, 1, 1);

compNetworkNameData.clim = clim;
compNetworkNameData.cmap = cmap;
compNetworkNameData.compData = compData;


function addCompNetwork(hObject, event_data, figH)
%% Add Component network
%

icatb_defaults;
global UI_FS;

figureTag = 'add_comp_cngm';
compFigHandle = findobj('tag', figureTag);
if (~isempty(compFigHandle))
    delete(compFigHandle);
end

figData = get(figH, 'userdata');

listH = findobj(figH, 'tag', 'comp');

compVals = [];
networkName = '';
network_color = [];
if (listH == hObject)
    if (~strcmpi(get(figH, 'selectionType'), 'open'))
        return;
    end
    val = get(listH, 'value');
    try
        networkName = figData.comp(val).name;
        compVals = figData.comp(val).value;
        network_color = figData.comp(val).network_color;
    catch
    end
end

compStr = num2str((1:figData.numICs)');

compFigHandle = icatb_getGraphics('Select Component Networks', 'normal', figureTag);
set(compFigHandle, 'menubar', 'none');

promptWidth = 0.52; controlWidth = 0.28; promptHeight = 0.05;
xOffset = 0.02; yOffset = promptHeight; listboxHeight = 0.6; yPos = 0.9;
okWidth = 0.12; okHeight = promptHeight;

%% Features text and listbox
promptPos = [xOffset, yPos - 0.5*yOffset, promptWidth, promptHeight];
textH = icatb_uicontrol('parent', compFigHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter Network Name', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);

promptPos = get(textH, 'position');

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
icatb_uicontrol('parent', compFigHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', networkName, 'tag', 'comp_network_name', 'fontsize', UI_FS - 1);


%editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = 0.1;
icatb_uicontrol('parent', compFigHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', editPos, 'string', 'color', 'tag', 'comp_network_color', 'fontsize', ...
    UI_FS - 2, 'userdata', network_color, 'callback', {@setColor, figH}, 'tooltipstring', 'You could select custom color for plotting networks. If nothing is specified, default colors are used. ...');

%% Right Listbox
listbox2Wdith = 0.1;
promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
promptPos(1) = (1 - xOffset - 2*listbox2Wdith);
promptPos(3) = 2*listbox2Wdith;
textH = icatb_uicontrol('parent', compFigHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Components', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);
promptPos = get(textH, 'position');
listboxYOrigin = promptPos(2) - 0.5*yOffset - listboxHeight;
listboxXOrigin = promptPos(1) + 0.5*listbox2Wdith;
listboxPos = [listboxXOrigin, listboxYOrigin, listbox2Wdith, listboxHeight];
compListH = icatb_uicontrol('parent', compFigHandle, 'units', 'normalized', 'style', 'listbox', 'position', listboxPos, 'string', compStr, 'tag', 'components', 'fontsize', UI_FS - 1, ...
    'min', 0, 'max', 2, 'value', compVals);

%% Show components
showWidth = 0.08; showHeight = 0.04;
showButtonPos = [listboxXOrigin + 0.5*listbox2Wdith - 0.5*showWidth, listboxYOrigin - yOffset - 0.5*showHeight, showWidth, showHeight];
showH = icatb_uicontrol('parent', compFigHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', showButtonPos, 'string', 'Show', 'fontsize', UI_FS - 1, 'callback', ...
    {@drawComp, figH, compFigHandle});

%% Plot image on the left hand side
axesPos = [xOffset, listboxYOrigin, listboxHeight, listboxHeight];
axesH = axes('parent', compFigHandle, 'units', 'normalized', 'position', axesPos, 'tag', 'axes_display_comp');

promptPos = axesPos;

%% Add cancel and run buttons
okPos = [0.5 - 0.5*okWidth, promptPos(2) - yOffset - 0.5*okHeight, okWidth, okHeight];
okPos(2) = okPos(2) - 0.5*okPos(4);
icatb_uicontrol('parent', compFigHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', okPos, 'string', 'Done', 'tag', 'done_button', 'fontsize', UI_FS - 1, 'callback', ...
    {@setCompCallback, compFigHandle, figH});

%% Draw components on the left hand side
drawComp(compListH, [], figH, compFigHandle);


function removeCompNetwork(hObject, event_data, figH)
%% Remove Component network
%

figData = get(figH, 'userdata');
listH = findobj(figH, 'tag', 'comp');
val = get(listH, 'value');

if (~isempty(val))
    check = icatb_questionDialog('title', 'Remove Component Network', 'textbody', 'Do you want to remove the component network from the list?');
    if (~check)
        return;
    end
end

try
    strs = cellstr(char(figData.comp.name));
    figData.comp(val) = [];
    strs(val) = [];
    set(listH, 'value', 1);
    set(listH, 'string', strs);
    set(figH, 'userdata', figData);
catch
end

function drawComp(hObject, event_data, figH, compFigHandle)
%% Draw component

icatb_defaults;
global UI_FONTNAME;
global FONT_COLOR;

fontSizeText = 8;
set(compFigHandle, 'pointer', 'watch');

listH = findobj(figH, 'tag', 'comp');
compNetworkNameData = get(listH, 'userdata');

axesH = get(compFigHandle, 'currentaxes');

clim = compNetworkNameData.clim;
cmap = compNetworkNameData.cmap;
compData = compNetworkNameData.compData;

sel_comp = get(findobj(compFigHandle, 'tag', 'components'), 'value');

if (~isempty(sel_comp))
    DIM = [size(compData, 1), size(compData, 2), length(sel_comp)];
    [im, numImagesX, numImagesY, textToPlot] = icatb_returnMontage(compData(:, :, sel_comp), [], DIM, [1, 1, 1], sel_comp);
    image(im, 'parent', axesH, 'CDataMapping', 'scaled');
    set(axesH, 'clim', clim); % set the axis positions to the specified
    axis(axesH, 'off');
    axis(axesH, 'image');
    colormap(cmap);
    textCount = 0;
    dim = size(im);
    yPos = 1 + dim(1) / numImagesY;
    for nTextRows = 1:numImagesY
        xPos = 1;
        for nTextCols = 1:numImagesX
            textCount = textCount + 1;
            if textCount <= DIM(3)
                text(xPos, yPos, num2str(round(textToPlot(textCount))), 'color', FONT_COLOR,  ...
                    'fontsize', fontSizeText, 'HorizontalAlignment', 'left', 'verticalalignment', 'bottom', ...
                    'FontName', UI_FONTNAME, 'parent', axesH);
            end
            xPos = xPos + (dim(2) / numImagesX);
        end
        % end for cols
        yPos = yPos + (dim(1) / numImagesY); % update the y position
    end
else
    cla(axesH);
end

set(compFigHandle, 'pointer', 'arrow');


function FNCM = loadFNC(file_name)
%% Load FNC matrix

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



function setCompCallback(hObject, event_data, compFigH, handles)
%% Get fields from component

figData = get(handles, 'userdata');
networkNameH = findobj(compFigH, 'tag', 'comp_network_name');
networkName = deblank(get(networkNameH, 'string'));

try
    
    if (isempty(networkName))
        error('You must enter a component network name');
    end
    
    listH = findobj(compFigH, 'tag', 'components');
    comps = get(listH, 'value');
    
    if (isempty(comps))
        error('Components are not selected');
    end
    
    if (length(figData.comp) > 0)
        chk = strmatch(lower(networkName), lower(cellstr(char(figData.comp.name))), 'exact');
        if (~isempty(chk))
            ind = chk;
        end
    end
    
    if (~exist('ind', 'var'))
        ind = length(figData.comp) + 1;
    end
    
    colorH = findobj(compFigH, 'tag', 'comp_network_color');
    colorValues = get(colorH, 'userdata');
    
    %% Set user selected information in figure
    figData.comp(ind).name = networkName;
    figData.comp(ind).value =  comps(:)';
    figData.comp(ind).network_color =  colorValues;
    set(handles, 'userdata', figData);
    compListH = findobj(handles, 'tag', 'comp');
    set(compListH, 'string', cellstr(char(figData.comp.name)));
    delete(compFigH);
    
catch
    icatb_errorDialog(lasterr, 'Component Selection');
end


function setColor(hObject, event_data, handles)
%% set color values
%

figData = get(handles, 'userdata');

networkNameH = findobj(gcbf, 'tag', 'comp_network_name');
networkName = deblank(get(networkNameH, 'string'));


if (length(figData.comp) > 0)
    chk = strmatch(lower(networkName), lower(cellstr(char(figData.comp.name))), 'exact');
    if (~isempty(chk))
        ind = chk;
    end
end

if (~exist('ind', 'var'))
    ind = length(figData.comp) + 1;
end

ud = get(hObject, 'userdata');

if (numel(ud) == 1)
    if (ud == 0)
        ud = [];
    end
end

if (isempty(ud))
    try
        ud = getFrameColor(ind);
    catch
    end
end

colorVal = uisetcolor(ud);
set(hObject, 'userdata', colorVal);



function runCallback(hObject, event_data, handles)
%% Select values
%

load('icatb_colors.mat', 'coldhot', 'coldhot_sensitive');

figData = get(handles, 'userdata');

for n = 1:2:length(figData.Options)
    figData.(figData.Options{n}) = figData.Options{n+1};
end


cmap = coldhot;
try
    cmap = eval(figData.cmap);
catch
end

if (mod(size(cmap, 1), 64) == 0)
    numToSkip = size(cmap, 1)/64;
    cmap = cmap(1:numToSkip:end, :);
end

figData.cmapN = cmap;

if (isempty(figData.comp))
    error('Please select components using + button');
end

network_names = cellstr(char(figData.comp.name));
network_values = cell(size(network_names));
network_colors = cell(size(network_names));

for n = 1:length(network_values)
    network_values{n} = figData.comp(n).value;
end

for n = 1:length(network_values)
    tmp = figData.comp(n).network_color;
    if (numel(tmp) == 1)
        if (tmp == 0)
            tmp = [];
        end
    end
    if (isempty(tmp))
        tmp = getFrameColor(n);
    end
    network_colors{n} = tmp;
end

%% Network names and values
figData.comp_network_names = [network_names, network_values, network_colors];

displayTypeH = findobj(handles, 'tag', 'display_type');
displayStr = cellstr(get(displayTypeH, 'string'));
displayVal = get(displayTypeH, 'value');
figData.display_type = lower(displayStr{displayVal});

connH = findobj(handles, 'tag', 'conn_threshold');
figData.conn_threshold = str2num(get(connH, 'string'));

fncLabelH = findobj(handles, 'tag',  'fnc_label');
figData.colorbar_label = get(fncLabelH, 'string');

setappdata(0, 'inputConnGramData', figData);

drawnow;

delete(handles);


function frame_colors = getFrameColor(N)
%% Default frame colors


% Default Colors used from
% (http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=12)
frame_colors = [166, 206, 227; 31,120,180; 178,223,138; 51,160,44; 251,154,153;
    227,26,28; 253,191,111; 255,127,0; 202,178,214; 106,61,154; 255,255,153; 177,89,40];

frame_colors = frame_colors/256;
frame_colors = [frame_colors(1:2:end,:);frame_colors(2:2:end,:)];

if (~exist('N', 'var'))
    return;
end

try
    frame_colors = frame_colors(N, :);
catch
    frame_colors = rand(1, 3);
end

function options_callback(hObject, event_data, handles)
%% Options callback
%

figData = get(handles, 'userdata');
Options = figData.Options;

imOptions = char('Positive and Negative', 'Positive', 'Absolute value', 'Negative');
numParameters = 1;

matchedInd = strmatch('image_values', lower(Options(1:2:end)), 'exact');
image_values = Options{2*matchedInd};
returnValue = strmatch(lower(image_values), lower(imOptions), 'exact');
inputText(numParameters).promptString = 'Select Image Values';
inputText(numParameters).answerString = imOptions;
inputText(numParameters).uiType = 'popup';
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'image_values';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = returnValue;
inputText(numParameters).help = struct('title', 'Image Values', 'string', 'Option is provided to display positive and negative, positive, negative and absolute image values.');

numParameters = numParameters + 1;
matchedInd = strmatch('convert_to_z', lower(Options(1:2:end)), 'exact');
convert_to_z = Options{2*matchedInd};
zOptions = char('Yes', 'No');
chkInd = strmatch(lower(convert_to_z), lower(zOptions), 'exact');
inputText(numParameters).promptString = 'Convert Spatial Maps To Z-scores?';
inputText(numParameters).answerString = zOptions;
inputText(numParameters).uiType = 'popup';
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'convert_to_z';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = chkInd;
inputText(numParameters).help = struct('title', 'Threshold', 'string', 'If you have selected option yes, images are converted to z-scores.');

numParameters = numParameters + 1;
matchedInd = strmatch('threshold', lower(Options(1:2:end)), 'exact');
threshold = Options{2*matchedInd};
inputText(numParameters).promptString = 'Enter Threshold Value';
inputText(numParameters).answerString = num2str(threshold);
inputText(numParameters).uiType = 'edit';
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'threshold';
inputText(numParameters).enable = 'on';
inputText(numParameters).help = struct('title', 'Threshold', 'string', 'Threshold is applied on the spatial maps.');

numParameters = numParameters + 1;
sliceOptions = {'Axial', 'Coronal', 'Sagittal'};
matchedInd = strmatch('slice_plane', lower(Options(1:2:end)), 'exact');
slice_plane = Options{2*matchedInd};
chkInd = strmatch(lower(slice_plane), lower(sliceOptions), 'exact');
inputText(numParameters).promptString = 'Select Anatomical Plane';
inputText(numParameters).answerString = sliceOptions;
inputText(numParameters).uiType = 'popup';
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'slice_plane';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = chkInd;
inputText(numParameters).help = struct('title', 'Slice Plane', 'string', ...
    'Thumbnails are plotted using the selected slice plane. If you have selected render option, rendered images are shown.');

numParameters = numParameters + 1;
matchedInd = strmatch('cmap', lower(Options(1:2:end)), 'exact');
inputText(numParameters).promptString = 'Enter colormap for correlations';
inputText(numParameters).answerString = Options{2*matchedInd};
inputText(numParameters).uiType = 'edit';
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'cmap';
inputText(numParameters).enable = 'on';
inputText(numParameters).help = struct('title', 'colormap', 'string', 'You could use any valid command to get the colormap like coldhot, jet, winter.');


numParameters = numParameters + 1;
matchedInd = strmatch('clim', lower(Options(1:2:end)), 'exact');
inputText(numParameters).promptString = 'Enter min and max of correlations.';
inputText(numParameters).answerString = num2str(Options{2*matchedInd});
inputText(numParameters).uiType = 'edit';
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'CLIM';
inputText(numParameters).enable = 'on';
inputText(numParameters).help = struct('title', 'Colorbar range', 'string', 'Enter min and max of correlations. You could leave it as empty if max and min are determined using data.');

numParameters = numParameters + 1;
matchedInd = strmatch('imwidth', lower(Options(1:2:end)), 'exact');
inputText(numParameters).promptString = 'Enter image width of the spatial maps';
inputText(numParameters).answerString = num2str(Options{2*matchedInd});
inputText(numParameters).uiType = 'edit';
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'imWidth';
inputText(numParameters).enable = 'on';
inputText(numParameters).help = struct('title', 'Image Width', 'string', ...
    'You could provide a specified value like 0.06. If nothing is specified, image width is automatically determined.');


answer = icatb_inputDialog('inputtext', inputText, 'Title', 'Select connectogram options', 'handle_visibility',  'on');

if (~isempty(answer))
    % ICA options with flags and the values corresponding to it
    Options = cell(1, 2*length(answer));
    
    for i = 1:length(answer)
        Options{2*i - 1} = inputText(i).tag;
        Options{2*i} = answer{i};
    end
    figData.Options = Options;
    set(handles, 'userdata', figData);
end



function selectAnatomical(hObject, event_data, figH)
%% Anatomical callback
%

figData = get(figH, 'userdata');

startPath = fileparts(which('groupica.m'));
startPath = fullfile(startPath, 'icatb_templates');

oldDir = pwd;

if (~exist(startPath, 'dir'))
    startPath = pwd;
end

% get the structural file
structFile = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Structural File', 'filter', ...
    '*.img;*.nii', 'fileType', 'image', 'fileNumbers', 1, 'startpath', startPath);

drawnow;

cd(oldDir);

if (~isempty(structFile))
    figData.structFile = structFile;
    set(figH, 'userdata', figData);
end


function [newX, newY]= getBezXY(x0, y0)

chk = find(isnan(x0) == 1);

if (isempty(chk))
    newX{1} = x0;
    newY{1} = y0;
else
    
    newX = cell(1, length(chk));
    newY = newX;
    
    s = 1;
    for n = 1:length(chk)
        if (n > 1)
            s = chk(n-1) + 1;
        end
        e = chk(n) - 1;
        newX{n} = x0(s:e);
        newY{n} = y0(s:e);
    end
    
end
