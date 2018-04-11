function icatb_write_nifti_data(outputFileName, V, data, description)
% Use volume information to create a 4d nifti file.

if ~exist('description', 'var')
    description = '4D Nifti Data';
end

N    = cat(1, V.private);

% Initialise maximum and minimum
mx   = -Inf;
mn   = Inf;

% Compute scaling factor
for i = 1:size(data, 4)
    dat = data(:, :, :, i);
    mx = max(mx, max(dat(:)));
    mn = min(mn, min(dat(:)));
end

sf         = max(mx, -mn) / 32767;
% End for computing scaling factor

ni         = icatb_nifti;
ni.dat     = icatb_file_array(outputFileName, [V(1).dim size(data, 4)], 'INT16-BE', 0, sf, 0);
ni.mat     = N(1).mat;
ni.mat0    = N(1).mat;
ni.descrip = description;

create(ni);

% Write data
for i = 1:size(ni.dat, 4)
    ni.dat(:, :, :, i) = data(:, :, :, i);
    icatb_spm_get_space([ni.dat.fname ',' num2str(i)], V(i).mat);
end
% End for writing data
