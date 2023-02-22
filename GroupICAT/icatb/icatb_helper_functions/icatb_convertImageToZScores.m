function [im] = icatb_convertImageToZScores(icasig, center_image)
% converts images to z scores


if ~exist('center_image', 'var')
    center_image = 1;
end


% Number of images
numberOfImages = size(icasig,1);
% loop over images
for i=1:numberOfImages

    % get the current image
    v = icasig(i, :);
    mask_ind = (v ~= 0);
    v2 = v(mask_ind);

    if (~center_image)
        v2 = detrend(v2, 0);
    end

    vstd = norm(v2, 2) ./ sqrt(length(v2) - 1);

    if (vstd ~= 0)
        v=v./(eps + vstd);
    else
        disp('-- Not converting to z-scores as division by zero warning may occur.');
    end

    icasig(i, :) = v;
end

im = icasig;


