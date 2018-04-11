function [v, offset] = icatb_recenter_image(v)
%% Center image distribution to zero (Based on code from Elena Allen)
%
% Inputs:
% v - Input image
%
% Outputs:
% v - Centered image
%

icatb_defaults;
global SMOOTHINGVALUE;

if (isempty(SMOOTHINGVALUE))
    SMOOTHINGVALUE = 1.1;
end

mask_ind = (v ~= 0);
v2 = v(mask_ind);

if (isempty(mask_ind))
    error('No non-zero voxels found');
end

v2 = v2(:);
nPoints = ceil(length(v2)/10);
hist_bins = linspace(min(v2), max(v2), nPoints);
%hist_bins = min(v2):.1:max(v2);
[dist,bins] = hist(v2, hist_bins);
dist = icatb_gauss_smooth1D(dist, SMOOTHINGVALUE);
[mv, mind] = max(dist);
offset = bins(mind);
v(mask_ind) = v2-offset;