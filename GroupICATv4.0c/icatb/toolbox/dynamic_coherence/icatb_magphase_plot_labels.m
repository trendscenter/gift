function [ output_args ] = icatb_magphase_plot_labels(Wxy, CM_over, coin)

% [F,A,C,I] = plot_FNC(FNC, CLIM, LABEL, RSN_I, F, T)
%if nargin < 2 || ~isvarname('MOD')
%[LABEL, RSN_I, aa] = comp_labels;
%end
% define the size of each "MODULE"
MOD = [5, 2, 11, 6, 13, 8, 2];
cMOD = cumsum(MOD);

Alpha_Image = abs(Wxy);
%Overlay_Image = exp(1i*angle(Wxy));
Overlay_Image = angle(Wxy);


if ~exist('coin', 'var') || isempty(coin)
    Underlay_Image = ones(size(Overlay_Image));
else
    Underlay_Image = coin;
end

if ~exist('CLIM', 'var') || isempty(CLIM)
     CLIM = [-pi, pi];
end

if ~exist('T', 'var') || isempty(T)
    T = 'correlation';
end

% Set the Min/Max T-values for alpha coding

A_range = [0 max(Alpha_Image(:))];
%fprintf('%1.2f\n', max(Alpha_Image(:)));
%A_range = [0 2.9];


% Choose a colormap for the underlay
CM_under = gray(256);

U_RGB = convert_to_RGB(Underlay_Image, CM_under, [0 1]);

layer1 = image(U_RGB); %axis image

hold on;

Overlay_Image = convert_to_RGB(Overlay_Image, CM_over, [-pi pi]);


layer2 = icatb_simtb_pcolor(1:size(Overlay_Image,2), 1:size(Overlay_Image,1), Overlay_Image);

axis(gca, 'square');

    
alphamap = abs(Alpha_Image);
alphamap(alphamap > A_range(2)) = A_range(2);
alphamap(alphamap < A_range(1)) = 0;
alphamap = alphamap/A_range(2);


% Adjust the alpha values of the overlay
set(layer2, 'alphaData', alphamap);



function IMrgb = convert_to_RGB(IM, cm, cmLIM)
% convert_to_RGB - converts any image to truecolor RGB using a specified colormap  
% USAGE: IMrgb = convert_to_RGB(IM, cm, cmLIM)
% INPUTS: 
%    IM    = the image [m x n]
%    cm    = the colormap [p x 3], optional; default = jet(256)
%    cmLIM = the data limits [min max] to be used in the color-mapping 
%            optional; default = [min(IM) max(IM)]
% OUTPUTS: 
%    IMrgb = the truecolor RGB image [m x n x 3]
% Based on ind2rgb from the Image Processing Toolbox
% EA Allen August 30, 2011
% eallen@mrn.org
%--------------------------------------------------------------------------
if nargin < 2, cm = jet(256); end
if nargin < 3, cmLIM = [min(IM(:)) max(IM(:))]; end

IM = IM-cmLIM(1);
IM = IM/(cmLIM(2)-cmLIM(1));
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
