function [newData] = icatb_resizeData(templateImg, files, flipWarning)
% Resize files to the template image


if ~exist('flipWarning', 'var')
    flipWarning = 0;
end

if ~isstruct(templateImg)
    % Template volume
    VT = icatb_spm_vol(deblank(templateImg(1, :)));
    VT = VT(1);
else
    VT = templateImg;
    clear templateImg;
end

% Form new file name
files = icatb_rename_4d_file(files);

newData = zeros([size(files, 1), VT.dim(1:3)]);

% Loop over files
for nn = 1:size(files, 1)
    V = icatb_spm_vol(deblank(files(nn, :)));

    % Interpolate images
    check = V.dim(1:3) ./ VT.dim(1:3);
    if any(check ~=1 )
        data = zeros(VT.dim(1:3));
        for slice = 1:VT.dim(3)
            M  = inv(icatb_spm_matrix([0 0 -slice 0 0 0 1 1 1])*inv(VT.mat)*V.mat);
            i1 = icatb_spm_slice_vol(V, M, VT.dim(1:2), 1);
            if flipWarning
                if sign(VT.mat(1, 1)) ~= sign(V.mat(1, 1))
                    warning('warning:flip', 'File %s is flipped w.r.t \n %s', V.fname, VT.fname);
                end
            end
            data(:, :, slice) = i1;
        end
    else
        data = icatb_read_vols(V);
    end
    % End for interpolating

    newData(nn, :, :, :) = data;
    clear data;
    clear V
end
% End loop over files

newData(isnan(newData)) = 0;