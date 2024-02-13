%--------------------------------------------------------------------------
% Sample Matlab code for creating images with hue and alpha color-mapping.
%
%  Notes for using this code:
%  You must have OpenGL available on your system to use transparency (alpha). 
%  When rendering transparency MATLAB automatically uses OpenGL if it is
%  available. If it is not available, transparency will not display. 
%  See the figure property RendererMode for more information.
%
% EA Allen August 30, 2011
% eallen@mrn.org
%--------------------------------------------------------------------------

%% 1. Load the AOD_data.mat file with sample data from the fMRI AOD experiment
%--------------------------------------------------------------------------
%load AOD_data.mat
% For a single axial slice (Z = 2 mm) of data, you should have:
% Bmap_N_S: 'Difference between Novel and Standard betas averaged over 28 subjects'
% Tmap_N_S: 'T-statistics for the paired t-test comparing Novel and Standard betas'
% Pmap_N_S: 'Binary map indicating significance at P<0.001 (fdr corrected)'
% Underlay: 'Structural image ch2bet from MRIcron, warped to functional data'
%--------------------------------------------------------------------------

function cdata = icatb_dualcodeImage(Bmap_N_S, Tmap_N_S, Pmap_N_S, Underlay, cmap)


%% 2. Set some defaults that will affect the appearance of the image
%--------------------------------------------------------------------------
% Set the Min/Max values for hue coding
absmax = max(abs(Bmap_N_S(:))); 
H_range = [-absmax absmax]; % The colormap is symmetric around zero

% Set the Min/Max T-values for alpha coding
absmaxT = max(abs(Tmap_N_S(:))); 
A_range = [-absmaxT, absmaxT];
% Voxels with t-stats of 0 will be completely transparent; 
% voxels with t-stat magnitudes greater or equal than 5 will be opaque.

% Set the labels for the colorbar
hue_label = 'Beta Difference (Novel - Standard)';
alpha_label = '|t|';

CM_under = cmap(size(cmap, 1)/2 + 1:end, :);
CM_over = cmap(1:size(cmap, 1)/2, :);
% Choose a colormap for the underlay
%CM_under = gray(256);

% Choose a colormap for the overlay
%CM_over = jet(256);
%--------------------------------------------------------------------------

%% 3. Do the actual plotting
%--------------------------------------------------------------------------
% Make a figure and set of axes
F = figure('Color', 'k', 'Units', 'Normalized', 'Position', [0.3, 0.4, 0.2, 0.35], 'visible', 'off'); 
axesH = axes('Position', [0 0 1 1], 'parent', F); 

% Transform the underlay and beta map to RGB values, based on specified colormaps
% See function convert_to_RGB() for more information
U_RGB = convert_to_RGB(Underlay, CM_under);
O_RGB = convert_to_RGB(Bmap_N_S, CM_over, H_range);

% Plot the underlay
layer1 = image(U_RGB); axis image
hold on;
% Now, add the Beta difference map as an overlay
layer2 = image(O_RGB); axis image

% Use the T-statistics to create an alpha map (which must be in [0,1])
alphamap = abs(Tmap_N_S);
alphamap(alphamap > A_range(2)) = A_range(2);
alphamap(alphamap < A_range(1)) = 0;
alphamap = alphamap/A_range(2);

% Adjust the alpha values of the overlay 
set(layer2, 'alphaData', alphamap);

% Add some (black) contours to annotate nominal significance
if (exist('Pmap_N_S', 'var') && ~isempty(Pmap_N_S))
    hold on;
    if (sum(Pmap_N_S(:)) ~= 0)
        [C, CH] = contour(Pmap_N_S, 1, 'k');
    end
end
%--------------------------------------------------------------------------

fd = getframe(axesH);
cdata = fd.cdata;
delete(F);
drawnow;

% %% 4. Create a 2D colorbar for the dual-coded overlay
% %--------------------------------------------------------------------------
% G = figure('color', 'k', 'Units', 'Normalized', 'Position', [0.5, 0.4, 0.06, 0.35]);
% x = linspace(A_range(1), A_range(2), 256); 
% % x represents the range in alpha (abs(t-stats))
% y = linspace(H_range(1), H_range(2), size(CM_over,1));
% % y represents the range in hue (beta weight difference)
% [X,Y] = meshgrid(x,y); % Transform into a 2D matrix
% imagesc(x,y,Y); axis xy; % Plot the colorbar
% set(gca, 'Xcolor', 'w', 'Ycolor', 'w')
% colormap(CM_over); 
% alpha(X);
% alpha('scaled');  
% xlabel(alpha_label)
% set(gca, 'YAxisLocation', 'right')
% ylabel(hue_label)
%--------------------------------------------------------------------------

function IMrgb = convert_to_RGB(IM, cm, cmLIM)
% convert_to_RGB - converts any image to truecolor RGB using a specified colormap  
%
% USAGE: IMrgb = convert_to_RGB(IM, cm, cmLIM)
%
% INPUTS: 
%    IM    = the image [m x n]
%    cm    = the colormap [p x 3] <optional; default = jet(256)>
%    cmLIM = the data limits [min max] to be used in the color-mapping 
%            <optional; default = [min(IM) max(IM)]>
% OUTPUTS: 
%    IMrgb = the truecolor RGB image [m x n x 3]
%
% Based on ind2rgb from the Image Processing Toolbox
% EA Allen August 30, 2011
% eallen@mrn.org
%--------------------------------------------------------------------------
if nargin < 2
    cm = jet(256);
end

if nargin < 3
    cmLIM = [min(IM(:)) max(IM(:))];
end

IM = IM-cmLIM(1);
IM = IM/(cmLIM(2)-cmLIM(1)+eps);
nIND = size(cm,1);
IM = round(IM*(nIND-1));

IM = double(IM)+1;
r = zeros(size(IM)); r(:) = cm(IM,1);
g = zeros(size(IM)); g(:) = cm(IM,2);
b = zeros(size(IM)); b(:) = cm(IM,3);

IMrgb = zeros([size(IM),3]);
% Fill in the r, g, and b channels
IMrgb(:,:,1) = r;
IMrgb(:,:,2) = g;
IMrgb(:,:,3) = b;
